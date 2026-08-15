"""
03_explore_nifecgdb.py — look before building. Answer one question first.

WHY THIS EXISTS
NIFECGDB is the only route to testing maternal QRS detection on a PREGNANT
ABDOMEN without ethics approval and a participant. But the published
literature contradicts itself about what its .qrs annotation files contain:

    "only annotations specifying the location of maternal R-waves"
    "only have QRS time samples for FECG"       (i.e. FETAL)
    "without reference annotations"

Three papers, three claims. Everything downstream depends on which is right,
so this script settles it EMPIRICALLY before a single line of pipeline code
gets written.

THE TEST IS TRIVIAL AND DECISIVE
    maternal heart rate in pregnancy : ~60-100 BPM (elevated 10-20 over
                                       non-pregnant baseline)
    fetal heart rate                 : ~110-160 BPM
Median RR from the annotations separates these with no overlap whatsoever.
One number, no ambiguity, no assumption.

    -> if MATERNAL: the dataset answers our question directly.
    -> if FETAL   : it does NOT serve the MHR question. Stop and rethink
                    rather than forcing it.

WHAT THE DATABASE IS
    55 recordings, ONE subject, gestational weeks 21-40
    EDF+ format, 1 kHz, 16-bit, bandpass filtered 0-100 Hz at acquisition
    2 thoracic (maternal ECG) + 3-4 abdominal (maternal + fetal) channels
    annotations in <record>.edf.qrs (WFDB annotation format)

THE PAIRED COMPARISON THIS ENABLES
Thoracic and abdominal channels are recorded SIMULTANEOUSLY in every record.
Thoracic is analogous to our AD8232 smoke-test placement; abdominal is the
belt scenario. Running the frozen detector on both isolates the PLACEMENT
effect with subject, session, gestational age and equipment all held
constant. That is a far cleaner experiment than device-vs-database, and it
is the actual question: what does moving electrodes from chest to abdomen
cost the maternal QRS detector?

TWO CONFOUNDS TO RECORD NOW, NOT LATER
  - Single subject. This tests abdominal placement in pregnancy, NOT
    population variation. It cannot support any claim about generalisation
    across women.
  - Electrode positions were NOT fixed across recordings and were sometimes
    moved to improve SNR. So placement varies alongside gestational age, and
    any trend across weeks cannot be attributed to gestation alone.

FORMAT NOTES (wfdb 4.3.1)
    signals : wfdb.io.convert.read_edf(path)   -- NOT wfdb.rdrecord
    beats   : wfdb.rdann('<...>/ecgca102.edf', 'qrs')
              the base name already ends in '.edf', which is unusual and
              trips the obvious guess
    rdedfann() reads EDF-EMBEDDED annotations. That is a different thing and
    is not what we want here.

STATUS: written without access to the actual files (PhysioNet is not
reachable from the authoring sandbox). The EDF load path and the annotation
filename convention are the parts most likely to need adjustment. Failures
are caught and reported per record rather than raising, so a partial run
still tells you what works.
"""
import os
import glob
import argparse
import collections

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import wfdb
from wfdb.io.convert import read_edf

DATASET = 'nifecgdb'
DATA_SUBPATH = os.path.join('data', DATASET)
PLOT_DIR = os.path.join('plots', DATASET)

SEARCH_DEPTH = 5          # parent directories to walk when locating data/


def resolve_data_dir(subpath=DATA_SUBPATH, depth=SEARCH_DEPTH):
    """Find data/nifecgdb/ in cwd or an ancestor.

    This script lives in hardware/ecg_bringup/ but the reference datasets sit
    at the repo root alongside mitdb, ctu-uhb-ctgdb and the rest. Walking up
    means the script runs correctly from either location, matching the same
    resolver 02_detect_device.py uses for the frozen pipeline files.
    """
    here = os.path.abspath(os.getcwd())
    for _ in range(depth + 1):
        cand = os.path.join(here, subpath)
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    return subpath          # fall through; main() reports the miss clearly


ANN_EXT = 'qrs'
FS_EXPECTED = 1000.0

# Physiological bands used ONLY to classify what the annotations are.
# Maternal upper bound is 110, not 105: maternal HR in late pregnancy commonly
# reaches 100-110, while fetal sits near 140. The first version used 105 and
# labelled two records (ecgca699 at 108.7, ecgca868 at 107.1) as FETAL, which
# was an artifact of the threshold rather than a finding.
MATERNAL_BPM = (50.0, 110.0)
FETAL_BPM = (110.0, 200.0)

# Annotation-quality thresholds. These reference sets are NOT exhaustive the
# way MIT-BIH's are, and that has to be measured before scoring.
GAP_FACTOR = 1.75        # RR above this x median => a beat is probably MISSING
REFRACTORY_S = 0.200     # RR below this => physiologically impossible, spurious
EDF_ANN_LABEL = 'EDF Annotations'   # the embedded EDF+ text channel

PLOT_SECONDS = 8.0


def list_records(data_dir):
    """Record base names, e.g. 'ecgca102.edf' (the '.edf' IS part of it)."""
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(data_dir, '*.edf')))


def load_annotations(data_dir, rec):
    path = os.path.join(data_dir, rec)
    if not os.path.exists(path + '.' + ANN_EXT):
        return None
    return wfdb.rdann(path, ANN_EXT)


def read_edf_quiet(path, header_only=False):
    """read_edf with its stdout captured, and the header_only dict normalised.

    TWO WFDB BEHAVIOURS THIS WRAPS.

    1. read_edf prints, once per matching channel:
           *** This may be an EDF+ Annotation file instead, please see the
               `rdedfann` function. ***
       It is a bare print(), not a warning or an exception, emitted from
       wfdb/io/convert/edf.py when a channel LABEL equals 'EDF Annotations'.
       It says nothing about the .qrs files — those are read by rdann, a
       different code path entirely. Over 55 records it buries the output.

       It does carry one real consequence: that annotation channel is counted
       as a SIGNAL. It occupies a column of p_signal holding text bytes
       reinterpreted as numbers. Any later loop over channels must exclude it
       by name or it will feed garbage to the detector.

    2. With header_only=True, read_edf returns a plain DICT, not a Record
       (see the `if header_only:` branch in edf.py). Attribute access such as
       r.sig_name raises AttributeError. Normalised here to a dict either way
       so callers have one shape to handle.
    """
    import io as _io
    import contextlib
    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf):
        r = read_edf(path, header_only=header_only)
    noisy = EDF_ANN_LABEL in buf.getvalue() or 'EDF+ Annotation' in buf.getvalue()
    if header_only:
        return dict(r), noisy
    return {'fs': float(r.fs), 'sig_len': int(r.sig_len), 'n_sig': int(r.n_sig),
            'sig_name': [str(n).strip() for n in r.sig_name],
            'record': r}, noisy


def annotation_quality(rr, fs):
    """Measure how complete and clean an annotation set is.

    MIT-BIH annotations are exhaustive and cardiologist-corrected. These are
    not, and scoring against an incomplete reference produces a meaningless
    PPV: every real beat the reference omits becomes a false positive that
    says nothing about the detector.

        suspected_missed : RR above GAP_FACTOR x median. A skipped beat
                           doubles the interval, so these cluster near 2x.
        suspected_spurious: RR below the 200 ms refractory floor. No heart
                           produces these; they are annotation noise.
    """
    med = float(np.median(rr))
    missed = int((rr > GAP_FACTOR * med).sum())
    spurious = int((rr < REFRACTORY_S).sum())
    near_double = int(((rr > 1.8 * med) & (rr < 2.2 * med)).sum())
    return {'median_rr_s': med,
            'n_intervals': int(rr.size),
            'suspected_missed': missed,
            'pct_missed': 100.0 * missed / rr.size if rr.size else 0.0,
            'near_double': near_double,
            'suspected_spurious': spurious,
            'pct_spurious': 100.0 * spurious / rr.size if rr.size else 0.0}


def classify(bpm):
    if MATERNAL_BPM[0] <= bpm < MATERNAL_BPM[1]:
        return 'MATERNAL'
    if FETAL_BPM[0] <= bpm < FETAL_BPM[1]:
        return 'FETAL'
    return 'UNCLEAR'


def summarise_record(data_dir, rec, want_signal=False):
    """Returns a dict of findings, or {'error': ...}. Never raises."""
    out = {'record': rec}
    path = os.path.join(data_dir, rec)

    # --- annotations: the load-bearing part ---
    try:
        ann = load_annotations(data_dir, rec)
    except Exception as e:
        return {'record': rec, 'error': f'rdann failed: {e}'}
    if ann is None:
        return {'record': rec, 'error': 'no .qrs file'}

    fs = float(ann.fs) if getattr(ann, 'fs', None) else FS_EXPECTED
    out['ann_fs'] = fs
    out['n_beats'] = int(ann.sample.size)
    if ann.sample.size < 3:
        out['error'] = f'only {ann.sample.size} annotations'
        return out

    rr = np.diff(ann.sample) / fs
    rr = rr[rr > 0]
    med = float(np.median(rr))
    out['median_rr_s'] = med
    out['median_bpm'] = 60.0 / med
    out['bpm_p5'] = 60.0 / float(np.percentile(rr, 95))
    out['bpm_p95'] = 60.0 / float(np.percentile(rr, 5))
    out['verdict'] = classify(out['median_bpm'])
    out['symbols'] = dict(collections.Counter(ann.symbol).most_common(5))
    out.update(annotation_quality(rr, fs))
    out['ann_span_s'] = float(ann.sample[-1] - ann.sample[0]) / fs

    # --- signals: channel inventory ---
    try:
        hdr, noisy = read_edf_quiet(path, header_only=True)
        out['sig_name'] = [str(n).strip() for n in hdr['sig_name']]
        out['n_sig'] = int(hdr['n_sig'])
        out['sig_fs'] = float(hdr['fs'])
        out['duration_s'] = (float(hdr['sig_len']) / float(hdr['fs'])
                             if hdr.get('sig_len') else None)
        out['has_edf_ann_channel'] = EDF_ANN_LABEL in out['sig_name']
        if out['duration_s']:
            out['ann_coverage_pct'] = 100.0 * out['ann_span_s'] / out['duration_s']
    except Exception as e:
        out['sig_error'] = f'read_edf(header_only) failed: {e}'

    if want_signal:
        try:
            full, _ = read_edf_quiet(path, header_only=False)
            out['record_obj'] = full['record']
        except Exception as e:
            out['sig_error'] = f'read_edf failed: {e}'
    return out


def plot_one(res, data_dir, rec):
    r = res.get('record_obj')
    ann = load_annotations(data_dir, rec)
    if r is None or ann is None:
        print("  (skipping plot: signal or annotations unavailable)")
        return
    os.makedirs(PLOT_DIR, exist_ok=True)

    fs = float(r.fs)
    n = min(int(PLOT_SECONDS * fs), r.p_signal.shape[0])
    t = np.arange(n) / fs
    sel = ann.sample[ann.sample < n]

    fig, axes = plt.subplots(r.n_sig, 1, figsize=(14, 2.0 * r.n_sig), sharex=True)
    axes = np.atleast_1d(axes)
    for i, ax in enumerate(axes):
        sig = r.p_signal[:n, i]
        ax.plot(t, sig, lw=0.7, color='0.35')
        if sel.size:
            ax.vlines(sel / fs, np.nanmin(sig), np.nanmax(sig),
                      colors='tab:red', alpha=0.45, lw=0.9)
        ax.set_ylabel(r.sig_name[i].strip()[:14], fontsize=8)
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel('seconds')
    fig.suptitle(f'{rec} — all channels, {sel.size} annotations in first '
                 f'{PLOT_SECONDS:.0f} s  ({res.get("verdict","?")})')
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f'03_explore_{rec.replace(".edf","")}.png')
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  plot -> {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', default=None)
    ap.add_argument('--plot-record', default=None,
                    help="record to plot, e.g. ecgca102.edf (default: first OK)")
    args = ap.parse_args()

    data_dir = args.data_dir or resolve_data_dir()
    recs = list_records(data_dir)
    if not recs:
        print(f"No .edf files in {os.path.abspath(data_dir)}.")
        print("Run download_nifecgdb.py from the repo root first —")
        print("download_data.py will NOT work here: dl_database() infers files")
        print("from WFDB headers, and this database has none (it is EDF+).")
        raise SystemExit(1)

    print(f"=== {DATASET}: {len(recs)} records in {data_dir} ===\n")

    rows = [summarise_record(data_dir, r) for r in recs]
    ok = [r for r in rows if 'error' not in r]
    bad = [r for r in rows if 'error' in r]

    print(f"{'record':>16} {'beats':>7} {'med BPM':>8} {'5-95% BPM':>14}  verdict")
    for r in rows:
        if 'error' in r:
            print(f"{r['record']:>16}  {r['error']}")
            continue
        print(f"{r['record']:>16} {r['n_beats']:>7} {r['median_bpm']:>8.1f} "
              f"{r['bpm_p5']:>6.1f}-{r['bpm_p95']:<7.1f} {r['verdict']}")

    if not ok:
        print("\nNo record produced usable annotations. Check the .qrs filename")
        print("convention: the base name should already end in '.edf'.")
        raise SystemExit(1)

    # ---------------------------------------------------------- THE VERDICT
    verdicts = collections.Counter(r['verdict'] for r in ok)
    all_bpm = np.array([r['median_bpm'] for r in ok])
    print("\n" + "=" * 62)
    print("  WHAT ARE THE ANNOTATIONS?")
    print("=" * 62)
    print(f"  records with usable annotations : {len(ok)}/{len(recs)}")
    print(f"  median BPM across records       : {np.median(all_bpm):.1f}")
    print(f"  range                           : {all_bpm.min():.1f}-{all_bpm.max():.1f}")
    print(f"  verdicts                        : {dict(verdicts)}")
    print(f"\n  reference bands: maternal {MATERNAL_BPM[0]:.0f}-{MATERNAL_BPM[1]:.0f}, "
          f"fetal {FETAL_BPM[0]:.0f}-{FETAL_BPM[1]:.0f} BPM")

    top = verdicts.most_common(1)[0]
    frac = top[1] / len(ok)
    print()
    if top[0] == 'MATERNAL' and frac >= 0.9:
        print("  => MATERNAL R-peak annotations. The dataset ANSWERS our question.")
        print("     Next: pre-register the thoracic-vs-abdominal comparison, then")
        print("     run the frozen detector on both channel groups per record.")
    elif top[0] == 'FETAL' and frac >= 0.9:
        print("  => FETAL annotations. This dataset does NOT serve the MHR")
        print("     question. Do not force it. Reconsider the source before")
        print("     writing any pipeline code.")
    else:
        print(f"  => MIXED or UNCLEAR ({top[0]} in only {100*frac:.0f}% of records).")
        print("     Do not proceed on an assumption. Inspect individual records")
        print("     and the database documentation before going further.")

    # ---------------------------------------------------------- quality
    print("\n" + "=" * 62)
    print("  ANNOTATION QUALITY — measure this BEFORE scoring anything")
    print("=" * 62)
    print("  MIT-BIH annotations are exhaustive and cardiologist-corrected.")
    print("  These are not. Every real beat the reference OMITS becomes a")
    print("  false positive that says nothing about the detector, so PPV is")
    print("  uninterpretable until this is characterised.\n")
    print(f"{'record':>16} {'ints':>6} {'missed%':>8} {'~2x':>5} {'spur%':>7} "
          f"{'cover%':>7}")
    q_missed, q_spur = [], []
    for r in ok:
        cov = r.get('ann_coverage_pct')
        print(f"{r['record']:>16} {r['n_intervals']:>6} "
              f"{r['pct_missed']:>8.2f} {r['near_double']:>5} "
              f"{r['pct_spurious']:>7.2f} "
              f"{(f'{cov:.1f}' if cov else '?'):>7}")
        q_missed.append(r['pct_missed'])
        q_spur.append(r['pct_spurious'])

    q_missed = np.array(q_missed)
    q_spur = np.array(q_spur)
    print(f"\n  suspected MISSED beats  : median {np.median(q_missed):.2f}% of "
          f"intervals, worst {q_missed.max():.2f}%")
    print(f"  suspected SPURIOUS      : median {np.median(q_spur):.2f}%, "
          f"worst {q_spur.max():.2f}%")
    print(f"  records with >5% missed : "
          f"{int((q_missed > 5).sum())}/{len(ok)}")
    print(f"  records with any spurious: "
          f"{int((q_spur > 0).sum())}/{len(ok)}")
    print("\n  A cluster near 2x the median RR is the signature of a SKIPPED")
    print("  beat. Intervals below the 200 ms refractory floor are impossible")
    print("  and are annotation noise.")
    print("\n  IMPLICATION: report SENSITIVITY as the headline (did the detector")
    print("  find the annotated beats?) and treat PPV as a lower bound, since")
    print("  unannotated real beats are counted against it. Decide and")
    print("  PRE-REGISTER this before running the detector, not after seeing")
    print("  a disappointing number.")

    # ---------------------------------------------------------- channels
    names = collections.Counter()
    counts = collections.Counter()
    n_edf_ann = 0
    for r in ok:
        if 'sig_name' in r:
            names.update(r['sig_name'])
            counts[r['n_sig']] += 1
            if r.get('has_edf_ann_channel'):
                n_edf_ann += 1
    if names:
        print("\n=== Channel inventory ===")
        print(f"  channel counts per record : {dict(counts)}")
        print("  channel names seen (name: records):")
        for n, c in names.most_common():
            flag = '   <-- EXCLUDE: EDF+ text channel, not a signal' \
                if n == EDF_ANN_LABEL else ''
            print(f"    {n:<24} {c}{flag}")
        if n_edf_ann:
            print(f"\n  {n_edf_ann}/{len(ok)} records carry an '{EDF_ANN_LABEL}'")
            print("  channel. This is the source of the repeated wfdb notice.")
            print("  It is TEXT stored in the signal stream, counted in n_sig")
            print("  and occupying a p_signal column as reinterpreted bytes.")
            print("  Any loop over channels MUST skip it by name, or the")
            print("  detector will be handed garbage.")
        print("\n  Identify which are THORACIC (maternal reference, analogous to")
        print("  our smoke-test placement) and which are ABDOMINAL (the belt")
        print("  scenario). That split defines the paired comparison.")
    else:
        errs = {r.get('sig_error') for r in ok if 'sig_error' in r}
        if errs:
            print("\n=== Channel inventory unavailable ===")
            for e in errs:
                print(f"  {e}")

    # ---------------------------------------------------------- one plot
    target = args.plot_record or ok[0]['record']
    print(f"\n=== Plotting {target} ===")
    res = summarise_record(data_dir, target, want_signal=True)
    if 'error' in res:
        print(f"  {res['error']}")
    else:
        if 'duration_s' in res and res['duration_s']:
            print(f"  duration {res['duration_s']:.1f} s, {res.get('n_sig','?')} channels, "
                  f"fs {res.get('sig_fs', '?')} Hz")
        plot_one(res, data_dir, target)

    if bad:
        print(f"\n  {len(bad)} record(s) unusable — list them as a "
              f"pre-registered exclusion before scoring anything.")

    print("\n  Nothing has been detected or scored. This step only establishes")
    print("  what the data IS. Pre-register the comparison before running the")
    print("  detector, same as the device phase.")


if __name__ == '__main__':
    main()
