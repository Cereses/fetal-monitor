"""
05_mhr_from_ecg.py — does the pipeline report the right maternal BPM?

THE GAP THIS CLOSES
Everything in the ECG work so far scores BEAT POSITIONS: sensitivity, PPV,
timing error. None of it computes the number the device would actually show a
pregnant woman. Those are different questions:

    Se/PPV : did we find each beat?
    BPM    : is the RATE right?

BPM is derived from RR INTERVALS, not from individual detections. A detector at
99.69% sensitivity still reports a bad rate if its misses CLUSTER, and a worse
detector can report a fine rate if its errors cancel. The detection metrics make
a good BPM LIKELY; they do not measure it. This does.

WHY MEDIAN RR, NOT MEAN
Within an 8 s window at ~92 BPM there are roughly 12 intervals. The NIFECGDB
reference omits ~6-7% of beats, so about one interval per window is DOUBLE
length. A mean is dragged down by that; a median is not — it would take more
than half the intervals being wrong to move it.

That matters more than it sounds: it means the reference incompleteness which
caps PPV at ~94% does NOT contaminate a median-based rate. Rate and PPV are
limited by different things, which is exactly why both get reported.

WHY 8-SECOND WINDOWS
Matches the frozen PPG pipeline (1.13 BPM MAE on BIDMC), so the two maternal
heart rate paths — AD8232 electrical and MAX30102 optical — become directly
comparable in the same unit on the same window length. Without that the two
results cannot be put in one sentence.

Windows are NON-OVERLAPPING. Overlapping windows share beats, so their errors
are correlated and any spread statistic computed from them is optimistic.

WHAT ARBITRATES
The same annotations the detector is scored against, put through the IDENTICAL
windowing function. Only the beat positions differ between the two sides, so
the comparison isolates detection and nothing else. A separate windowing path
for the reference would fold windowing differences into the reported error.

SOURCES
  nifecgdb thoracic : the target condition — pregnant woman, CHEST electrodes.
                      This is the placement the project has now committed to.
                      Abdominal channels are deliberately NOT included: the
                      placement study (04) measured a 12% catastrophic failure
                      rate there, and the design decision was taken because of
                      it. Rate accuracy on a placement we are not using would
                      be a number with no home.
  device capture    : the actual AD8232 hardware, against its frozen
                      hand-reviewed reference. Only ~31 beats, so a handful of
                      windows — reported for completeness, not for weight.

USAGE
    python 05_mhr_from_ecg.py
    python 05_mhr_from_ecg.py --limit 5
    python 05_mhr_from_ecg.py --window 10
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
PLOT_DIR = os.path.join('plots', 'nifecgdb')
CAPTURE_DIR = 'captures'

DETECTOR_FILE = '02_qrs_detect.py'
SEARCH_DEPTH = 5

WINDOW_S = 8.0             # matches the frozen PPG pipeline
MIN_INTERVALS = 4          # below this a median RR is not stable
BPM_LO, BPM_HI = 30.0, 220.0   # physiological sanity bound on a window

DEVICE_CAPTURE = 'raw_ecg_capture_20260729_180936.csv'
DEVICE_REF = os.path.join(
    RESULTS_DIR, 'raw_ecg_capture_20260729_180936_annotations_frozen.csv')

CLOSE_BPM = (2.0, 5.0)     # report fraction of windows within these


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
    raise FileNotFoundError(f"could not find {filename} from "
                            f"{os.path.abspath(os.getcwd())}")


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


# ---------------------------------------------------------------- rate
def windowed_bpm(samples, fs, n_total, window_s=WINDOW_S, apply_bounds=True):
    """Median-RR heart rate per non-overlapping window.

    Returns (window_index, bpm) for windows carrying enough intervals.
    Applied IDENTICALLY to detector output and reference annotations except
    for `apply_bounds`, so any difference between the two is detection.

    `apply_bounds` MUST BE FALSE FOR THE DETECTOR. The first version applied
    a 30-220 BPM physiological filter to BOTH sides:

        if not (BPM_LO <= r <= BPM_HI):
            continue

    On ecgca699 Thorax_2 the detector reports ~268 BPM — above the ceiling —
    so 24 of its 41 windows were silently discarded and the channel appeared
    to have only 17. The filter removed precisely the windows where the
    detector failed worst, and the reported MAE described performance only on
    windows where the detector was already roughly sane.

    The bound is legitimate on the REFERENCE: annotations implying 480 BPM
    are corrupt and cannot be scored against. It is not legitimate on the
    DETECTOR: a detector reporting 268 BPM has given a wrong answer, and the
    honest treatment is to score it as one.

    Windows dropped for too few intervals are still counted and reported as
    coverage, since a detector that reports nothing cannot be wrong.
    """
    samples = np.asarray(samples, dtype=np.int64)
    w = int(round(window_s * fs))
    n_win = int(n_total // w)
    idx, bpm, n_out = [], [], 0
    for k in range(n_win):
        a, b = k * w, (k + 1) * w
        s = samples[(samples >= a) & (samples < b)]
        if s.size < MIN_INTERVALS + 1:
            continue
        rr = np.diff(s) / fs
        rr = rr[rr > 0]
        if rr.size < MIN_INTERVALS:
            continue
        r = 60.0 / float(np.median(rr))
        if not (BPM_LO <= r <= BPM_HI):
            n_out += 1
            if apply_bounds:
                continue
        idx.append(k)
        bpm.append(r)
    return np.array(idx, dtype=int), np.array(bpm, dtype=float), n_win, n_out


def compare(det_samples, ref_samples, fs, n_total, window_s=WINDOW_S):
    """Paired per-window BPM. Only windows valid on BOTH sides are compared.

    Bounds are applied to the reference and NOT to the detector — see
    windowed_bpm. An out-of-range detector window is an error to be reported,
    not an inconvenience to be filtered.
    """
    di, db, n_win, d_out = windowed_bpm(det_samples, fs, n_total, window_s,
                                        apply_bounds=False)
    ri, rb, _, r_out = windowed_bpm(ref_samples, fs, n_total, window_s,
                                    apply_bounds=True)
    dmap = dict(zip(di.tolist(), db.tolist()))
    rmap = dict(zip(ri.tolist(), rb.tolist()))
    common = sorted(set(dmap) & set(rmap))
    if not common:
        return None
    d = np.array([dmap[k] for k in common])
    r = np.array([rmap[k] for k in common])
    return {
        'n_windows_total': int(n_win),
        'n_windows_ref': int(ri.size),
        'n_windows_det': int(di.size),
        'n_windows_paired': len(common),
        'n_det_out_of_range': int(d_out),
        'n_ref_out_of_range': int(r_out),
        'bpm_det': d, 'bpm_ref': r, 'error': d - r,
    }


def stats(err, ref):
    """MAE, bias, Bland-Altman limits, and closeness fractions."""
    out = {
        'n': int(err.size),
        'mae_bpm': float(np.mean(np.abs(err))),
        'bias_bpm': float(np.mean(err)),
        'sd_bpm': float(np.std(err, ddof=1)) if err.size > 1 else float('nan'),
        'median_abs_bpm': float(np.median(np.abs(err))),
        'max_abs_bpm': float(np.max(np.abs(err))),
        'ref_bpm_mean': float(np.mean(ref)),
        'ref_bpm_range': [float(ref.min()), float(ref.max())],
    }
    # Bland-Altman: agreement between two methods, not correlation. Two series
    # can correlate at 0.99 and still disagree by 10 BPM at every point.
    if err.size > 1:
        out['loa_lower'] = out['bias_bpm'] - 1.96 * out['sd_bpm']
        out['loa_upper'] = out['bias_bpm'] + 1.96 * out['sd_bpm']
    for c in CLOSE_BPM:
        out[f'within_{c:g}_bpm'] = float(np.mean(np.abs(err) <= c))
    return out


# ---------------------------------------------------------------- sources
def run_nifecgdb(qrs, data_dir, limit, window_s):
    paths = nio.list_records(data_dir)
    if limit:
        paths = paths[:limit]
    rows, errs, refs = [], [], []
    print(f"\n=== NIFECGDB — THORACIC channels only ({len(paths)} records) ===")
    print("  Abdominal channels excluded by design: 04 measured a 12% "
          "catastrophic")
    print("  failure rate there and the project committed to chest placement.\n")
    for i, p in enumerate(paths, 1):
        name = os.path.basename(p)
        try:
            rec = nio.load_record(p)
            ref, _, _ = nio.load_annotations(p, rec['fs'], rec['n_samples'])
            ok, msg = nio.check_alignment(ref, rec['n_samples'])
            if not ok:
                print(f"  [{i:2d}/{len(paths)}] {name}: SKIPPED — {msg}")
                continue
        except Exception as e:
            print(f"  [{i:2d}/{len(paths)}] {name}: LOAD FAILED — {e}")
            continue

        g = rec['gestation']
        gw = g[0] + g[1] / 7.0 if g else float('nan')
        per_chan = []
        for c, cname in enumerate(rec['sig_name']):
            if not cname.startswith(nio.THORACIC_PREFIX):
                continue
            det, _ = qrs.detect_qrs(rec['sig'][:, c].astype(float), rec['fs'])
            cmp = compare(det, ref, rec['fs'], rec['n_samples'], window_s)
            if cmp is None:
                continue
            st = stats(cmp['error'], cmp['bpm_ref'])
            st.update({'record': name, 'channel': cname,
                       'gestation_weeks': gw,
                       'windows_total': cmp['n_windows_total'],
                       'windows_paired': cmp['n_windows_paired'],
                       'coverage': cmp['n_windows_paired'] / cmp['n_windows_total']
                       if cmp['n_windows_total'] else 0.0})
            rows.append(st)
            per_chan.append(st['mae_bpm'])
            errs.append(cmp['error'])
            refs.append(cmp['bpm_ref'])
        if per_chan:
            print(f"  [{i:2d}/{len(paths)}] {name}  gest {gw:5.1f}w  "
                  f"MAE {np.mean(per_chan):5.2f} BPM  "
                  f"({rows[-1]['windows_paired']} windows)")
    if not errs:
        return rows, None, None
    return rows, np.concatenate(errs), np.concatenate(refs)


def run_device(qrs, window_s):
    cap = os.path.join(CAPTURE_DIR, DEVICE_CAPTURE)
    if not (os.path.exists(cap) and os.path.exists(DEVICE_REF)):
        print("\n=== Device capture: not found, skipping ===")
        return None
    ms, adc = [], []
    with open(cap, newline='') as fh:
        r = csv.reader(fh)
        next(r)
        for row in r:
            if row:
                ms.append(int(row[0]))
                adc.append(int(row[1]))
    ms, adc = np.array(ms), np.array(adc, dtype=float)
    fs = 1000.0 / float(np.median(np.diff(ms)))
    ref = []
    with open(DEVICE_REF, newline='') as fh:
        for row in csv.DictReader(fh):
            ref.append(int(row['sample']))
    ref = np.array(sorted(ref))
    det, _ = qrs.detect_qrs(adc, fs)
    cmp = compare(det, ref, fs, adc.size, window_s)
    print(f"\n=== Device capture — AD8232 hardware ===")
    if cmp is None:
        print("  no window valid on both sides (record is only 30 s)")
        return None
    st = stats(cmp['error'], cmp['bpm_ref'])
    print(f"  {DEVICE_CAPTURE}  fs {fs:.0f} Hz  {ref.size} reference beats")
    print(f"  windows paired : {cmp['n_windows_paired']} of "
          f"{cmp['n_windows_total']}")
    print(f"  MAE            : {st['mae_bpm']:.2f} BPM")
    print(f"  bias           : {st['bias_bpm']:+.2f} BPM")
    print(f"  max abs error  : {st['max_abs_bpm']:.2f} BPM")
    print(f"  reference rate : {st['ref_bpm_mean']:.1f} BPM "
          f"({st['ref_bpm_range'][0]:.1f}-{st['ref_bpm_range'][1]:.1f})")
    print(f"  Only {cmp['n_windows_paired']} windows from one 30 s capture. "
          f"Reported for")
    print("  completeness, not for weight.")
    return st


# ---------------------------------------------------------------- plots
def plot_agreement(err, ref, rows, out_path, window_s):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    det = ref + err
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    ax.scatter(ref, det, s=6, alpha=0.25, color='tab:blue')
    lim = [min(ref.min(), det.min()) - 3, max(ref.max(), det.max()) + 3]
    ax.plot(lim, lim, ls='--', lw=1, color='0.4')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('reference BPM')
    ax.set_ylabel('detector BPM')
    ax.set_title(f'Per-window rate ({window_s:.0f} s windows)',
                 loc='left', fontsize=9)
    ax.grid(alpha=0.3)

    # Bland-Altman: mean of the two methods against their difference.
    ax = axes[0, 1]
    mean = (ref + det) / 2
    bias, sd = float(np.mean(err)), float(np.std(err, ddof=1))
    ax.scatter(mean, err, s=6, alpha=0.25, color='tab:blue')
    ax.axhline(bias, color='tab:red', lw=1.2, label=f'bias {bias:+.2f}')
    ax.axhline(bias + 1.96 * sd, color='tab:orange', ls='--', lw=1,
               label=f'+1.96SD {bias+1.96*sd:+.2f}')
    ax.axhline(bias - 1.96 * sd, color='tab:orange', ls='--', lw=1,
               label=f'-1.96SD {bias-1.96*sd:+.2f}')
    ax.set_xlabel('mean of detector and reference (BPM)')
    ax.set_ylabel('detector - reference (BPM)')
    ax.set_title('Bland-Altman — agreement, not correlation',
                 loc='left', fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    lim = float(np.percentile(np.abs(err), 99.5)) + 0.5
    ax.hist(err, bins=60, range=(-lim, lim), color='tab:blue', alpha=0.8)
    ax.axvline(0, color='0.3', lw=1)
    ax.set_xlabel('error (BPM)')
    ax.set_ylabel('windows')
    ax.set_title('Error distribution', loc='left', fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    if rows:
        ax.scatter([r['gestation_weeks'] for r in rows],
                   [r['mae_bpm'] for r in rows], s=20, alpha=0.7,
                   color='tab:green')
    ax.set_xlabel('gestational age (weeks)')
    ax.set_ylabel('per-channel MAE (BPM)')
    ax.set_title('EXPLORATORY — confounded with electrode repositioning',
                 loc='left', fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle('Maternal BPM accuracy — NIFECGDB thoracic, frozen detector')
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--window', type=float, default=WINDOW_S)
    ap.add_argument('--data-dir', default=None)
    args = ap.parse_args()

    t0 = time.time()
    print("=== 05_mhr_from_ecg.py — maternal BPM accuracy ===\n")
    det_path = resolve_upward(DETECTOR_FILE)
    qrs = load_module(det_path, 'qrs_detect')
    print(f"  detector : {det_path}")
    print(f"             sha256 {sha256_file(det_path)[:16]}...")
    print(f"  window   : {args.window:.0f} s, non-overlapping, "
          f"median RR, >= {MIN_INTERVALS} intervals")

    data_dir = args.data_dir or nio.resolve_data_dir()
    rows, err, ref = run_nifecgdb(qrs, data_dir, args.limit, args.window)
    if err is None:
        print("\n  No windows produced. Nothing to report.")
        raise SystemExit(1)

    st = stats(err, ref)
    print("\n" + "=" * 62)
    print("  MATERNAL BPM ACCURACY — NIFECGDB thoracic, pooled")
    print("=" * 62)
    print(f"  channels        : {len(rows)}")
    print(f"  windows compared: {st['n']:,}")
    print(f"  reference rate  : {st['ref_bpm_mean']:.1f} BPM "
          f"({st['ref_bpm_range'][0]:.1f}-{st['ref_bpm_range'][1]:.1f})")
    print(f"\n  MAE             : {st['mae_bpm']:.3f} BPM")
    print(f"  median |error|  : {st['median_abs_bpm']:.3f} BPM")
    print(f"  bias            : {st['bias_bpm']:+.3f} BPM")
    print(f"  SD of error     : {st['sd_bpm']:.3f} BPM")
    print(f"  95% limits of agreement: [{st['loa_lower']:+.2f}, "
          f"{st['loa_upper']:+.2f}] BPM")
    print(f"  max |error|     : {st['max_abs_bpm']:.2f} BPM")
    for c in CLOSE_BPM:
        print(f"  within {c:g} BPM      : {100*st[f'within_{c:g}_bpm']:.2f}% "
              f"of windows")
    cov = float(np.mean([r['coverage'] for r in rows]))
    print(f"  window coverage : {100*cov:.1f}% "
          f"(windows usable on both sides)")

    print("\n  Limits of agreement, not MAE, is the number to quote for a")
    print("  device: it states the range within which the reported rate sits")
    print("  95% of the time. MAE alone hides the tail.")

    dev = run_device(qrs, args.window)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(os.path.join(RESULTS_DIR, 'mhr_ecg_channel_results.csv'),
              'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    summary = {
        'window_s': args.window,
        'min_intervals': MIN_INTERVALS,
        'rate_estimator': 'median RR within window',
        'detector_sha256': sha256_file(det_path),
        'nifecgdb_thoracic': {k: v for k, v in st.items()
                              if not isinstance(v, np.ndarray)},
        'device_capture': dev,
        'mean_window_coverage': cov,
        'runtime_s': time.time() - t0,
    }
    with open(os.path.join(RESULTS_DIR, 'mhr_ecg_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=2, default=float)

    plot_path = os.path.join(PLOT_DIR, 'mhr_ecg_agreement.png')
    plot_agreement(err, ref, rows, plot_path, args.window)

    print(f"\n  results -> {os.path.join(RESULTS_DIR, 'mhr_ecg_channel_results.csv')}")
    print(f"  summary -> {os.path.join(RESULTS_DIR, 'mhr_ecg_summary.json')}")
    print(f"  plot    -> {plot_path}")
    print(f"  runtime : {(time.time()-t0)/60:.1f} min")

    print("\n  COMPARABILITY: the frozen PPG pipeline reports 1.13 BPM MAE on")
    print("  BIDMC. Same window length, same estimator, so these two maternal")
    print("  heart rate paths can now be stated in one sentence. Different")
    print("  databases and subjects, so it is a comparison of paths, not a")
    print("  controlled head-to-head.")
    print("\n  SCOPE: ONE pregnant subject, chest electrodes, a research")
    print("  amplifier with a 50 Hz notch this project's hardware does not")
    print("  have. Rate accuracy of the ALGORITHM on pregnant chest ECG —")
    print("  not of the device on a pregnant woman.")


if __name__ == '__main__':
    main()
