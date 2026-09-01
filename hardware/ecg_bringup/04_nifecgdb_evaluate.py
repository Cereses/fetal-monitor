"""
04_nifecgdb_evaluate.py — what does abdominal placement cost the QRS detector?

THE QUESTION
Phase A ran the frozen detector on an AD8232 recording from a CHEST/TORSO
placement on a non-pregnant male: 31/31, Se 100%, PPV 100%. The device is for
a pregnant woman, and a belt puts the electrodes on the ABDOMEN, where maternal
QRS is uV rather than mV and fetal ECG plus uterine EMG sit inside the QRS band.

NIFECGDB records both placements SIMULTANEOUSLY in every recording, so subject,
session, gestational age, amplifier and electrodes are all held constant and
only placement varies. The deliverable is the GAP between the two groups, not
either number on its own.

Governed by notes/PREREGISTRATION_nifecgdb.md, committed before this file
existed. Read that first; this script only executes it.

WHAT IS IMPORTED, NEVER COPIED
    detect_qrs                <- 02_qrs_detect.py     (repo root)
    match_beats, TOLERANCE_S  <- 03_qrs_evaluate.py   (repo root)
    load_record, ...          <- nifecgdb_io.py       (this directory)
Those first two produced the MIT-BIH headline (Se 99.52%, PPV 99.35%) and the
Phase A device result. NO PARAMETER IN THEM MAY BE CHANGED IN RESPONSE TO
THESE NUMBERS. Editing them retroactively invalidates both prior results.
Their SHA-256 is recorded per run so that claim is checkable.

REPORTING CONVENTION, FIXED IN THE PRE-REGISTRATION
The .qrs reference is ~6-7% INCOMPLETE (two independent estimates: 7.33% from
RR gaps per record, 5.84% pooled from beats/duration). NIFECGDB exists to
support FETAL ECG research; the maternal annotations are supplied so the
maternal component can be located and cancelled. They are a means, not the
product, and are not the exhaustive cardiologist-corrected reference MIT-BIH
provides.

    Sensitivity  = HEADLINE.    Unaffected by omissions in the reference.
    PPV          = LOWER BOUND. Every real beat the reference omits is charged
                   to the detector as a false positive.
    Predicted PPV ceiling from a PERFECT detector: 92-94%.

THE TEST THAT MAKES THAT CONVENTION HONEST (pre-registered P4)
A lower bound is an excuse unless the unmatched detections are shown to land
where the reference is missing beats. So every false positive is classified:

    in_gap   : inside a reference interval longer than 1.75x that record's
               median RR — i.e. exactly where a beat was skipped
    outside  : before the first or after the last annotation, where the
               reference does not cover the signal at all
    genuine  : everywhere else. THESE are real false positives.

P4 passes if >= 70% of in-span false positives are in_gap. If it fails, the
lower-bound framing collapses and the false positives must be reported as real.

RUNTIME
313 channel-records over 9.16 h of 1000 Hz signal. The longest record is 2.78M
samples and detect_qrs runs sosfiltfilt plus a 150-sample moving-window
convolution per channel. Expect 10-60 minutes for the full run. Use --limit to
validate the plumbing on a few records first; a trial run is NOT a result and
is marked as such in the output.

USAGE
    python 04_nifecgdb_evaluate.py --limit 3     # plumbing check
    python 04_nifecgdb_evaluate.py               # full run
"""
import os
import csv
import json
import time
import hashlib
import argparse
import importlib.util
import collections

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import nifecgdb_io as nio

RESULTS_DIR = 'results'
PLOT_DIR = 'plots'
NOTES_DIR = 'notes'

DETECTOR_FILE = '02_qrs_detect.py'
EVALUATOR_FILE = '03_qrs_evaluate.py'
PREREG_FILE = os.path.join(NOTES_DIR, 'PREREGISTRATION_nifecgdb.md')
SEARCH_DEPTH = 5

GAP_FACTOR = 1.75          # matches 03_explore; a skipped beat doubles the RR
P4_THRESHOLD = 0.70        # pre-registered
P1_THRESHOLD = 0.98        # thoracic pooled sensitivity
P2_THRESHOLD = 0.75        # fraction of records where abdominal < thoracic
P3_RANGE = (0.92, 0.94)    # predicted PPV band
SUPPLEMENTARY = ['ecgca998.edf']    # 10.12% sub-refractory annotations


# ---------------------------------------------------------------- plumbing
def resolve_upward(filename, depth=SEARCH_DEPTH):
    here = os.path.abspath(os.getcwd())
    for _ in range(depth + 1):
        cand = os.path.join(here, filename)
        if os.path.exists(cand):
            return cand
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    raise FileNotFoundError(
        f"could not find {filename} in cwd or {depth} parents "
        f"(from {os.path.abspath(os.getcwd())}). This script must import the "
        f"ORIGINAL frozen file, never a copy.")


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- scoring
def reference_gaps(ref, gap_factor=GAP_FACTOR):
    """Intervals where the reference probably skipped a beat.

    Returns a list of (start, end) sample pairs. A skipped beat doubles the
    interval, so these cluster near 2x the record median — verified in
    03_explore, where the flagged count and the count between 1.8x and 2.2x
    matched almost exactly.
    """
    if ref.size < 3:
        return []
    rr = np.diff(ref)
    med = float(np.median(rr))
    idx = np.where(rr > gap_factor * med)[0]
    return [(int(ref[i]), int(ref[i + 1])) for i in idx]


def classify_fps(fps, ref, gaps):
    """Split false positives into in_gap / outside_span / genuine.

    This is the P4 test. Without it, 'PPV is a lower bound' is an unfalsifiable
    excuse rather than a measured claim.
    """
    if fps.size == 0:
        return {'in_gap': 0, 'outside_span': 0, 'genuine': 0,
                'in_span': 0, 'frac_in_gap': float('nan')}
    lo, hi = int(ref[0]), int(ref[-1])
    outside = int(((fps < lo) | (fps > hi)).sum())
    in_span = fps[(fps >= lo) & (fps <= hi)]
    in_gap = 0
    if gaps and in_span.size:
        mask = np.zeros(in_span.size, dtype=bool)
        for a, b in gaps:
            mask |= (in_span >= a) & (in_span <= b)
        in_gap = int(mask.sum())
    genuine = int(in_span.size) - in_gap
    return {'in_gap': in_gap, 'outside_span': outside, 'genuine': genuine,
            'in_span': int(in_span.size),
            'frac_in_gap': in_gap / in_span.size if in_span.size else float('nan')}


def score_channel(qrs, ev, sig, ref, fs, gaps):
    det, _ = qrs.detect_qrs(sig.astype(float), fs)
    tol = int(round(ev.TOLERANCE_S * fs))
    tp, fp, fn, errors, claimed = ev.match_beats(ref, det, tol)
    fps = det[~claimed]
    fpc = classify_fps(fps, ref, gaps)
    se = tp / (tp + fn) if (tp + fn) else 0.0
    ppv = tp / (tp + fp) if (tp + fp) else 0.0
    return {
        'detected': int(det.size), 'TP': int(tp), 'FP': int(fp), 'FN': int(fn),
        'sensitivity': se, 'PPV': ppv,
        'DER': (fp + fn) / ref.size if ref.size else 0.0,
        'timing_mae_ms': (float(np.mean(np.abs(errors)) / fs * 1000)
                          if errors.size else float('nan')),
        **{f'fp_{k}': v for k, v in fpc.items()},
    }


def group_of(name):
    if name.startswith(nio.THORACIC_PREFIX):
        return 'thoracic'
    if name.startswith(nio.ABDOMINAL_PREFIX):
        return 'abdominal'
    return 'other'


# ---------------------------------------------------------------- reporting
def pooled(rows):
    tp = sum(r['TP'] for r in rows)
    fp = sum(r['FP'] for r in rows)
    fn = sum(r['FN'] for r in rows)
    return {
        'n': len(rows), 'TP': tp, 'FP': fp, 'FN': fn,
        'sensitivity': tp / (tp + fn) if (tp + fn) else 0.0,
        'PPV': tp / (tp + fp) if (tp + fp) else 0.0,
    }


def plot_results(rows, out_path):
    thor = [r for r in rows if r['group'] == 'thoracic']
    abdo = [r for r in rows if r['group'] == 'abdominal']

    per_rec = collections.defaultdict(dict)
    for r in rows:
        per_rec[r['record']].setdefault(r['group'], []).append(r)
    paired = []
    for rec, g in per_rec.items():
        if 'thoracic' in g and 'abdominal' in g:
            paired.append((
                rec,
                float(np.median([x['sensitivity'] for x in g['thoracic']])),
                float(np.median([x['sensitivity'] for x in g['abdominal']])),
                g['thoracic'][0]['gestation_weeks']))

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    if paired:
        t = [100 * p[1] for p in paired]
        a = [100 * p[2] for p in paired]
        ax.scatter(t, a, s=28, alpha=0.75, color='tab:blue')
        lim = [min(t + a) - 2, 100.5]
        ax.plot(lim, lim, ls='--', lw=1, color='0.4', label='equal')
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.legend(fontsize=8)
    ax.set_xlabel('thoracic sensitivity (%)')
    ax.set_ylabel('abdominal sensitivity (%)')
    ax.set_title('Paired by record — below the line = abdominal is worse',
                 loc='left', fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    for grp, col, lbl in ((thor, 'tab:green', 'thoracic'),
                          (abdo, 'tab:orange', 'abdominal')):
        if grp:
            ax.scatter([r['gestation_weeks'] for r in grp],
                       [100 * r['sensitivity'] for r in grp],
                       s=16, alpha=0.6, color=col, label=lbl)
    ax.set_xlabel('gestational age (weeks)')
    ax.set_ylabel('sensitivity (%)')
    ax.set_title('EXPLORATORY — confounded with electrode repositioning',
                 loc='left', fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    data = [[100 * r['sensitivity'] for r in thor],
            [100 * r['sensitivity'] for r in abdo]]
    keep = [(l, d) for l, d in zip(['thoracic', 'abdominal'], data) if d]
    if keep:
        # 'labels' was renamed 'tick_labels' in matplotlib 3.9. Set the ticks
        # explicitly so this works on either side of that change.
        ax.boxplot([d for _, d in keep])
        ax.set_xticks(range(1, len(keep) + 1))
        ax.set_xticklabels([l for l, _ in keep])
    ax.set_ylabel('sensitivity (%)')
    ax.set_title('Per-channel sensitivity by placement', loc='left', fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    cats = ['in_gap', 'genuine', 'outside_span']
    x = np.arange(len(cats))
    w = 0.38
    for off, grp, col, lbl in ((-w / 2, thor, 'tab:green', 'thoracic'),
                               (w / 2, abdo, 'tab:orange', 'abdominal')):
        vals = [sum(r[f'fp_{c}'] for r in grp) for c in cats]
        ax.bar(x + off, vals, w, color=col, label=lbl)
    ax.set_xticks(x)
    ax.set_xticklabels(['in reference gap', 'genuine FP', 'outside span'])
    ax.set_ylabel('false positives')
    ax.set_title('P4 — where do unmatched detections fall?', loc='left', fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle('NIFECGDB — frozen Pan-Tompkins, thoracic vs abdominal placement')
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=0,
                    help='score only the first N records (0 = all). '
                         'A limited run is a PLUMBING CHECK, not a result.')
    ap.add_argument('--data-dir', default=None)
    args = ap.parse_args()

    t0 = time.time()
    print("=== 04_nifecgdb_evaluate.py ===\n")

    det_path = resolve_upward(DETECTOR_FILE)
    ev_path = resolve_upward(EVALUATOR_FILE)
    qrs = load_module(det_path, 'qrs_detect')
    ev = load_module(ev_path, 'qrs_evaluate')

    print("=== Frozen code (imported, not copied) ===")
    print(f"  detector : {det_path}")
    print(f"             sha256 {sha256_file(det_path)[:16]}...")
    print(f"  scorer   : {ev_path}")
    print(f"             sha256 {sha256_file(ev_path)[:16]}...")
    prereg_hash = None
    if os.path.exists(PREREG_FILE):
        prereg_hash = sha256_file(PREREG_FILE)
        print(f"  prereg   : {PREREG_FILE}")
        print(f"             sha256 {prereg_hash[:16]}...")
    else:
        print(f"  [warn] {PREREG_FILE} not found — the predictions this run")
        print("         resolves are supposed to be committed BEFORE it.")
    print(f"  tolerance: +/-{ev.TOLERANCE_S*1000:.0f} ms\n")

    paths = nio.list_records(args.data_dir)
    if args.limit:
        paths = paths[:args.limit]
        print(f"  *** LIMITED RUN: {len(paths)} record(s). Plumbing check "
              f"only, NOT a result. ***\n")

    rows = []
    for i, p in enumerate(paths, 1):
        name = os.path.basename(p)
        tr = time.time()
        try:
            rec = nio.load_record(p)
            ref, ann_fs, scale = nio.load_annotations(p, rec['fs'], rec['n_samples'])
            ok, msg = nio.check_alignment(ref, rec['n_samples'])
            if not ok:
                print(f"  [{i:2d}/{len(paths)}] {name}: SKIPPED — {msg}")
                continue
        except Exception as e:
            print(f"  [{i:2d}/{len(paths)}] {name}: LOAD FAILED — {e}")
            continue

        gaps = reference_gaps(ref)
        g = rec['gestation']
        gw = g[0] + g[1] / 7.0 if g else float('nan')

        for c, cname in enumerate(rec['sig_name']):
            grp = group_of(cname)
            if grp == 'other':
                continue
            r = score_channel(qrs, ev, rec['sig'][:, c], ref, rec['fs'], gaps)
            r.update({'record': name, 'channel': cname, 'group': grp,
                      'gestation_weeks': gw, 'ref_beats': int(ref.size),
                      'units': rec['units'][c], 'fs': rec['fs'],
                      'n_gaps': len(gaps),
                      'supplementary': name in SUPPLEMENTARY})
            rows.append(r)

        el = time.time() - tr
        done = time.time() - t0
        eta = done / i * (len(paths) - i)
        se_t = [x['sensitivity'] for x in rows
                if x['record'] == name and x['group'] == 'thoracic']
        se_a = [x['sensitivity'] for x in rows
                if x['record'] == name and x['group'] == 'abdominal']
        print(f"  [{i:2d}/{len(paths)}] {name}  ref {ref.size:>5}  gaps {len(gaps):>4}  "
              f"Se thor {100*np.median(se_t) if se_t else float('nan'):6.2f}%  "
              f"abd {100*np.median(se_a) if se_a else float('nan'):6.2f}%  "
              f"({el:.1f}s, eta {eta/60:.1f}m)")

    if not rows:
        print("\n  No channel-records scored. Nothing to report.")
        raise SystemExit(1)

    head = rows
    supp = [r for r in rows if r['supplementary']]
    thor = [r for r in head if r['group'] == 'thoracic']
    abdo = [r for r in head if r['group'] == 'abdominal']
    pt, pa = pooled(thor), pooled(abdo)

    print("\n" + "=" * 66)
    print("  HEADLINE — pooled by placement")
    print("=" * 66)
    print(f"  {'group':>10} {'chans':>6} {'ref beats':>10} {'TP':>8} {'FP':>7} "
          f"{'FN':>7} {'Se':>8} {'PPV':>8}")
    for lbl, s in (('thoracic', pt), ('abdominal', pa)):
        print(f"  {lbl:>10} {s['n']:>6} {s['TP']+s['FN']:>10,} {s['TP']:>8,} "
              f"{s['FP']:>7,} {s['FN']:>7,} {100*s['sensitivity']:>7.2f}% "
              f"{100*s['PPV']:>7.2f}%")
    gap_se = pt['sensitivity'] - pa['sensitivity']
    print(f"\n  PLACEMENT COST (the deliverable): "
          f"{100*gap_se:+.2f} percentage points of sensitivity")
    print(f"  thoracic {100*pt['sensitivity']:.2f}% -> abdominal "
          f"{100*pa['sensitivity']:.2f}%")

    # --- paired ---
    per_rec = collections.defaultdict(dict)
    for r in head:
        per_rec[r['record']].setdefault(r['group'], []).append(r['sensitivity'])
    paired = [(rec, float(np.median(g['thoracic'])), float(np.median(g['abdominal'])))
              for rec, g in per_rec.items()
              if 'thoracic' in g and 'abdominal' in g]
    worse = sum(1 for _, t, a in paired if a < t)
    print(f"\n  Paired by record: abdominal worse in {worse}/{len(paired)} "
          f"({100*worse/len(paired) if paired else 0:.0f}%)")
    if paired:
        d = [100 * (t - a) for _, t, a in paired]
        print(f"  per-record gap: median {np.median(d):+.2f} pp, "
              f"range {min(d):+.2f} to {max(d):+.2f}")

    # --- P4 ---
    print("\n" + "=" * 66)
    print("  P4 — where do the unmatched detections fall?")
    print("=" * 66)
    for lbl, grp in (('thoracic', thor), ('abdominal', abdo)):
        ig = sum(r['fp_in_gap'] for r in grp)
        gn = sum(r['fp_genuine'] for r in grp)
        os_ = sum(r['fp_outside_span'] for r in grp)
        insp = ig + gn
        frac = ig / insp if insp else float('nan')
        print(f"  {lbl:>10}: in reference gap {ig:>7,}   genuine {gn:>7,}   "
              f"outside span {os_:>6,}   -> {100*frac:5.1f}% in gap")
    tot_ig = sum(r['fp_in_gap'] for r in head)
    tot_gn = sum(r['fp_genuine'] for r in head)
    p4 = tot_ig / (tot_ig + tot_gn) if (tot_ig + tot_gn) else float('nan')
    print(f"\n  overall: {100*p4:.1f}% of in-span false positives fall inside a "
          f"reference gap")
    p4_pass = p4 >= P4_THRESHOLD
    print(f"  P4 threshold {100*P4_THRESHOLD:.0f}% -> "
          f"{'PASS' if p4_pass else 'FAIL'}")
    if p4_pass:
        print("  The detector is finding beats the reference omitted. Reporting")
        print("  PPV as a lower bound is justified BY MEASUREMENT.")
    else:
        print("  The lower-bound framing does NOT hold. These false positives")
        print("  are real and must be reported as such, not explained away.")

    # --- predictions ---
    print("\n" + "=" * 66)
    print("  PRE-REGISTERED PREDICTIONS")
    print("=" * 66)
    p1 = pt['sensitivity'] >= P1_THRESHOLD
    p2 = (worse / len(paired) >= P2_THRESHOLD) if paired else False
    in_band = [P3_RANGE[0] <= s['PPV'] <= P3_RANGE[1] for s in (pt, pa)]
    print(f"  P1 thoracic Se >= {100*P1_THRESHOLD:.0f}%            : "
          f"{100*pt['sensitivity']:.2f}%  {'HELD' if p1 else 'DID NOT HOLD'}")
    print(f"  P2 abdominal worse in >= {100*P2_THRESHOLD:.0f}% recs : "
          f"{100*worse/len(paired) if paired else 0:.0f}%  "
          f"{'HELD' if p2 else 'DID NOT HOLD'}")
    print(f"  P3 PPV in {100*P3_RANGE[0]:.0f}-{100*P3_RANGE[1]:.0f}% both groups : "
          f"thor {100*pt['PPV']:.2f}%, abd {100*pa['PPV']:.2f}%  "
          f"{'HELD' if all(in_band) else 'DID NOT HOLD'}")
    print(f"  P4 >= {100*P4_THRESHOLD:.0f}% FPs in gaps        : "
          f"{100*p4:.1f}%  {'HELD' if p4_pass else 'DID NOT HOLD'}")
    print("\n  A failed prediction is a RESULT. Report it and diagnose the")
    print("  mechanism. Do not revise the prediction after the fact, and do")
    print("  not touch 02_qrs_detect.py or 03_qrs_evaluate.py.")

    if supp:
        ps = pooled(supp)
        print(f"\n  SUPPLEMENTARY — {', '.join(sorted({r['record'] for r in supp}))} "
              f"(10.12% sub-refractory annotations)")
        print(f"    Se {100*ps['sensitivity']:.2f}%  PPV {100*ps['PPV']:.2f}%  "
              f"({ps['n']} channels). Retained in the headline; shown here so")
        print("    nothing is hidden and no record list was tuned.")

    # --- persist ---
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(PLOT_DIR, exist_ok=True)
    csv_path = os.path.join(RESULTS_DIR, 'nifecgdb_channel_results.csv')
    with open(csv_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    summary = {
        'limited_run': bool(args.limit),
        'records_scored': len({r['record'] for r in rows}),
        'channel_records': len(rows),
        'detector_sha256': sha256_file(det_path),
        'evaluator_sha256': sha256_file(ev_path),
        'prereg_sha256': prereg_hash,
        'tolerance_ms': ev.TOLERANCE_S * 1000,
        'tolerance_samples': int(round(ev.TOLERANCE_S * rows[0]['fs'])),
        'thoracic': pt, 'abdominal': pa,
        'placement_cost_se_pp': 100 * gap_se,
        'paired_records': len(paired),
        'abdominal_worse_records': worse,
        'p4_frac_in_gap': p4,
        'predictions': {'P1': bool(p1), 'P2': bool(p2),
                        'P3': bool(all(in_band)), 'P4': bool(p4_pass)},
        'runtime_s': time.time() - t0,
    }
    json_path = os.path.join(RESULTS_DIR, 'nifecgdb_summary.json')
    with open(json_path, 'w') as fh:
        json.dump(summary, fh, indent=2)

    plot_path = os.path.join(PLOT_DIR, 'nifecgdb_placement_comparison.png')
    plot_results(head, plot_path)

    print(f"\n  results -> {csv_path}")
    print(f"  summary -> {json_path}")
    print(f"  plot    -> {plot_path}")
    print(f"  runtime : {(time.time()-t0)/60:.1f} min")

    print("\n  SCOPE: ONE subject, 55 recordings. This measures the cost of")
    print("  abdominal placement in a pregnant woman with everything else held")
    print("  constant. It says NOTHING about variation across women, nothing")
    print("  about fetal heart rate, and nothing about this project's hardware.")


if __name__ == '__main__':
    main()
