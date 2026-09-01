"""
06_diagnose_window_bpm.py — why does one record report a 49 BPM error?

THE OBSERVATION
05 reported, pooled over 8,166 windows of NIFECGDB thoracic ECG:
    median |error| 0.145 BPM, 98.76% of windows within 2 BPM
but also
    SD 4.988 BPM, max |error| 110.58 BPM, reference max 213.9 BPM
and one record, ecgca699, at 49.19 BPM MAE — roughly 0.4% of the windows
carrying about 37% of the pooled MAE.

An SD 34x the median error is not a spread. It is a tight core plus a few
catastrophic outliers, and those outliers must be explained before any
agreement statistic is quoted.

TWO HYPOTHESES, AND THEY POINT OPPOSITE WAYS
  H1  THE REFERENCE IS DOUBLED. Each beat carries two annotations, halving
      the median RR and doubling the apparent rate. The detector is right.
      Supporting arithmetic: in the worst window the reference reads 213.9
      BPM while the detector reports ~103.3 — a ratio of 2.07.
  H2  THE DETECTOR IS HALVING. It misses alternate beats, so its rate is
      half the truth. The reference is right and this is a real failure.

These are not stylistic alternatives. Under H1 the pipeline is fine and the
database annotation has a defect; under H2 the pipeline has a serious flaw
that the pooled numbers are hiding. The record cannot be excluded, and the
result cannot be written up, until one of them is established.

IS 108.7 BPM ITSELF SUSPICIOUS? NO.
ecgca699's record-level median is 108.7 BPM at 32.4 weeks. Sinus tachycardia
begins above 100, pregnancy raises resting rate by 10-20, and anaemia,
dehydration or movement raise it further. An elevated median is plausible
physiology and is NOT what is being questioned here. What is questioned is
individual windows reading ~214 BPM, which is not sinus tachycardia in a
healthy pregnancy at rest and would not appear only in isolated 8 s windows.

Note also that the record's 95th-percentile rate is 117.4 BPM, not 214. If
most intervals were doubled the 95th percentile would sit near 214. So any
doubling is CONCENTRATED in short stretches rather than spread throughout —
which is exactly why a record-level median hides it and a windowed comparison
exposes it.

THE TESTS, IN ORDER OF DECISIVENESS

  1. ALTERNATION. Genuine doubling places a spurious mark between real beats
     (typically on the S or T deflection), so intervals alternate
     short-long-short-long. Consecutive RR intervals then correlate
     STRONGLY NEGATIVELY. Real tachycardia has no such structure; its
     consecutive intervals are near-independent or mildly positively
     correlated through respiratory sinus arrhythmia.
     A lag-1 correlation near -0.5 or below is close to conclusive.

  2. HALF-MEDIAN CLUSTER. Doubling creates a population of intervals near
     half the record median. Measured as the fraction of intervals within
     0.4-0.6x the median.

  3. THE PICTURE. The raw signal across the worst window with reference marks
     and detections overlaid. If two marks sit on one QRS, H1. If the trace
     genuinely beats twice as fast as the detector claims, H2. Every
     statistic above is a proxy for this plot.

WHY THIS RECORD IS NOT BEING EXCLUDED
Dropping a record after seeing that it hurts the headline is exactly what the
flutter audit refused to do and what the ecgca998 handling avoided. The
outcome here is a MEASURED FINDING about the reference, reported alongside
the result — not a quietly shorter record list.

USAGE
    python 06_diagnose_window_bpm.py                       # ecgca699
    python 06_diagnose_window_bpm.py --record ecgca416.edf
"""
import os
import json
import argparse
import importlib.util

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import nifecgdb_io as nio

RESULTS_DIR = 'results'
PLOT_DIR = os.path.join('plots', 'nifecgdb')
DETECTOR_FILE = '02_qrs_detect.py'
SEARCH_DEPTH = 5

DEFAULT_RECORD = 'ecgca699.edf'
WINDOW_S = 8.0
MIN_INTERVALS = 4
N_WORST = 5

ALTERNATION_THRESHOLD = -0.30   # lag-1 RR correlation below this = alternating
HALF_LO, HALF_HI = 0.40, 0.60   # window around half the median RR


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
    raise FileNotFoundError(f"could not find {filename}")


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def window_bpm(samples, fs, n_total, window_s=WINDOW_S):
    """Median-RR rate per non-overlapping window. Mirrors 05 exactly."""
    samples = np.asarray(samples, dtype=np.int64)
    w = int(round(window_s * fs))
    out = {}
    for k in range(int(n_total // w)):
        a, b = k * w, (k + 1) * w
        s = samples[(samples >= a) & (samples < b)]
        if s.size < MIN_INTERVALS + 1:
            continue
        rr = np.diff(s) / fs
        rr = rr[rr > 0]
        if rr.size < MIN_INTERVALS:
            continue
        out[k] = 60.0 / float(np.median(rr))
    return out


def alternation(rr):
    """Lag-1 correlation of consecutive RR intervals.

    Strongly negative means short-long-short-long: the signature of an extra
    mark inserted between real beats. Real rate changes are smooth, so genuine
    tachycardia does not alternate.
    """
    if rr.size < 3:
        return float('nan')
    a, b = rr[:-1], rr[1:]
    if np.std(a) == 0 or np.std(b) == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def half_cluster(rr, med):
    """Fraction of intervals sitting near HALF the record median."""
    if rr.size == 0:
        return 0.0
    return float(np.mean((rr > HALF_LO * med) & (rr < HALF_HI * med)))


def extra_detection_phase(det, ref, fs, tol_s=0.150):
    """Where in the cardiac cycle do the EXTRA detections sit?

    For every detection with no reference beat within +/-tol_s, measure its
    offset from the PRECEDING reference beat, expressed both in seconds and as
    a fraction of that beat's RR interval.

    This exists because the plot cannot settle it. refine_to_r_peak returns the
    index of maximum |bandpassed| signal, and the marker is drawn at the RAW
    amplitude there, so a detection genuinely on a QRS can render well below
    the visible R peak. Vertical marker positions are misleading and reading a
    mechanism off them is guesswork.

    The histogram is not:
        mode near 0.20-0.30 s, ~0.35-0.50 of RR  -> T WAVE
        mode near 0.04-0.09 s, ~0.10 of RR       -> S-wave rebound / QRS
                                                    counted twice
        broad, no mode                           -> neither; the mechanism is
                                                    something else and must be
                                                    found before anything is
                                                    claimed
    """
    det = np.asarray(det, dtype=np.int64)
    ref = np.asarray(ref, dtype=np.int64)
    if det.size == 0 or ref.size < 2:
        return np.array([]), np.array([])
    tol = int(round(tol_s * fs))
    nearest = np.searchsorted(ref, det)
    nearest = np.clip(nearest, 1, ref.size - 1)
    d_prev = det - ref[nearest - 1]
    d_next = ref[np.minimum(nearest, ref.size - 1)] - det
    unmatched = (np.minimum(d_prev, d_next) > tol)

    offs, frac = [], []
    for i in np.where(unmatched)[0]:
        j = int(np.searchsorted(ref, det[i])) - 1
        if j < 0 or j + 1 >= ref.size:
            continue
        rr = (ref[j + 1] - ref[j]) / fs
        if rr <= 0:
            continue
        off = (det[i] - ref[j]) / fs
        if 0 < off < rr:
            offs.append(off)
            frac.append(off / rr)
    return np.array(offs), np.array(frac)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--record', default=DEFAULT_RECORD)
    ap.add_argument('--data-dir', default=None)
    ap.add_argument('--window', type=float, default=WINDOW_S)
    ap.add_argument('--channel', default=None,
                    help='diagnose one named channel instead of the worst')
    args = ap.parse_args()

    data_dir = args.data_dir or nio.resolve_data_dir()
    path = os.path.join(data_dir, args.record)
    if not os.path.exists(path):
        print(f"Missing {path}")
        raise SystemExit(1)

    qrs = load_module(resolve_upward(DETECTOR_FILE), 'qrs_detect')
    rec = nio.load_record(path)
    ref, ann_fs, scale = nio.load_annotations(path, rec['fs'], rec['n_samples'])
    fs = rec['fs']
    g = rec['gestation']

    print(f"=== 06_diagnose_window_bpm.py — {args.record} ===")
    print(f"  gestation  : {g[0]}w{g[1]}d" if g else "  gestation  : unknown")
    print(f"  duration   : {rec['duration_s']:.1f} s   fs {fs:.0f} Hz")
    print(f"  channels   : {rec['sig_name']}")
    print(f"  reference  : {ref.size} beats\n")

    # ---------------------------------------------------------- record level
    rr = np.diff(ref) / fs
    rr = rr[rr > 0]
    med = float(np.median(rr))
    bpm = 60.0 / rr
    print("=== Record-level reference rate ===")
    print(f"  median RR   : {med:.3f} s  -> {60/med:.1f} BPM")
    print(f"  BPM percentiles  5th {np.percentile(bpm,5):.1f}  "
          f"25th {np.percentile(bpm,25):.1f}  50th {np.median(bpm):.1f}  "
          f"75th {np.percentile(bpm,75):.1f}  95th {np.percentile(bpm,95):.1f}")
    print(f"  BPM max     : {bpm.max():.1f}   (RR min {rr.min():.3f} s)")
    print(f"  BPM min     : {bpm.min():.1f}   (RR max {rr.max():.3f} s)")
    print(f"\n  An elevated MEDIAN is plausible physiology: sinus tachycardia")
    print("  starts above 100 BPM and pregnancy adds 10-20. It is the extreme")
    print("  windows, not the median, that this script is testing.")

    r_alt = alternation(rr)
    r_half = half_cluster(rr, med)
    print(f"\n=== Test 1: alternation (record-wide) ===")
    print(f"  lag-1 RR correlation : {r_alt:+.3f}")
    print(f"  threshold for 'alternating' : <= {ALTERNATION_THRESHOLD:+.2f}")
    print(f"\n=== Test 2: half-median cluster (record-wide) ===")
    print(f"  intervals in {HALF_LO:.1f}-{HALF_HI:.1f}x median : "
          f"{100*r_half:.2f}%")

    # ---------------------------------------------------------- windows
    # EVERY thoracic channel is evaluated, then the WORST is diagnosed.
    # The first version took thor[0] unconditionally. On ecgca699 that is
    # Thorax_1, which scores 0.22 BPM MAE — while the record's reported 49.19
    # BPM came from Thorax_2. The failing channel was never examined. A
    # diagnostic that can silently inspect the healthy half of a record is
    # worse than none, because it returns a confident wrong answer.
    thor = [i for i, n in enumerate(rec['sig_name'])
            if n.startswith(nio.THORACIC_PREFIX)]
    if args.channel:
        sel = [i for i, n in enumerate(rec['sig_name']) if n == args.channel]
        if not sel:
            print(f"  channel {args.channel!r} not in {rec['sig_name']}")
            raise SystemExit(1)
        thor = sel
    if not thor:
        thor = [0]

    print(f"\n=== Per-channel windowed comparison ===")
    print(f"  {'channel':>12} {'windows':>8} {'MAE':>9} {'median|e|':>10} "
          f"{'max|e|':>9} {'mean ratio':>11}")
    chan_stats = []
    for c in thor:
        d, _ = qrs.detect_qrs(rec['sig'][:, c].astype(float), fs)
        rw = window_bpm(ref, fs, rec['n_samples'], args.window)
        dw = window_bpm(d, fs, rec['n_samples'], args.window)
        com = sorted(set(rw) & set(dw))
        if not com:
            print(f"  {rec['sig_name'][c]:>12}  no paired windows")
            continue
        e = np.array([dw[k] - rw[k] for k in com])
        ratio = float(np.mean([rw[k] / dw[k] for k in com if dw[k]]))
        chan_stats.append({'idx': c, 'name': rec['sig_name'][c], 'det': d,
                           'ref_w': rw, 'det_w': dw, 'common': com,
                           'errs': e, 'mae': float(np.mean(np.abs(e))),
                           'ratio': ratio})
        print(f"  {rec['sig_name'][c]:>12} {len(com):>8} "
              f"{np.mean(np.abs(e)):>9.2f} {np.median(np.abs(e)):>10.2f} "
              f"{np.max(np.abs(e)):>9.2f} {ratio:>11.2f}")

    if not chan_stats:
        print("  No channel produced paired windows.")
        raise SystemExit(1)

    worst = max(chan_stats, key=lambda s: s['mae'])
    ch = worst['idx']
    det = worst['det']
    ref_w, det_w, common = worst['ref_w'], worst['det_w'], worst['common']
    errs = worst['errs']
    print(f"\n  diagnosing WORST channel: {worst['name']} "
          f"(MAE {worst['mae']:.2f} BPM)")
    if worst['ratio'] < 0.65:
        print(f"  mean reference/detector ratio {worst['ratio']:.2f} — the")
        print("  DETECTOR is reporting roughly DOUBLE the reference rate.")
        print("  That is over-detection (T waves counted as beats), not a")
        print("  reference defect. Neither H1 nor H2 as originally framed.")
    elif worst['ratio'] > 1.6:
        print(f"  mean reference/detector ratio {worst['ratio']:.2f} — the")
        print("  reference is roughly double the detector.")

    print(f"\n=== Windowed comparison on {rec['sig_name'][ch]} ===")
    print(f"  windows paired : {len(common)}")
    print(f"  MAE            : {np.mean(np.abs(errs)):.2f} BPM")
    print(f"  median |error| : {np.median(np.abs(errs)):.2f} BPM")

    order = np.argsort(-np.abs(errs))[:N_WORST]
    w = int(round(args.window * fs))
    print(f"\n  {N_WORST} worst windows:")
    print(f"  {'win':>4} {'t (s)':>8} {'ref BPM':>9} {'det BPM':>9} "
          f"{'ratio':>7} {'err':>9} {'n_ref':>6} {'n_det':>6} "
          f"{'lag1':>7} {'half%':>7}")
    worst_rows = []
    for j in order:
        k = common[j]
        a, b = k * w, (k + 1) * w
        rs = ref[(ref >= a) & (ref < b)]
        ds = det[(det >= a) & (det < b)]
        wrr = np.diff(rs) / fs
        wrr = wrr[wrr > 0]
        alt = alternation(wrr)
        hc = half_cluster(wrr, med)
        ratio = ref_w[k] / det_w[k] if det_w[k] else float('nan')
        print(f"  {k:>4} {a/fs:>8.1f} {ref_w[k]:>9.1f} {det_w[k]:>9.1f} "
              f"{ratio:>7.2f} {errs[j]:>+9.1f} {rs.size:>6} {ds.size:>6} "
              f"{alt:>+7.2f} {100*hc:>6.1f}%")
        worst_rows.append({'window': int(k), 't_start_s': a / fs,
                           'ref_bpm': ref_w[k], 'det_bpm': det_w[k],
                           'ratio': ratio, 'error_bpm': float(errs[j]),
                           'n_ref': int(rs.size), 'n_det': int(ds.size),
                           'lag1_corr': alt, 'half_cluster': hc})

    # ------------------------------------------- what ARE the extra beats?
    offs, frac = extra_detection_phase(det, ref, fs)
    print(f"\n=== Test 3: phase of UNMATCHED detections (objective) ===")
    print(f"  unmatched detections analysed : {offs.size}")
    phase_verdict = 'insufficient data'
    if offs.size >= 20:
        med_off = float(np.median(offs))
        med_frac = float(np.median(frac))
        iqr = float(np.percentile(offs, 75) - np.percentile(offs, 25))
        print(f"  offset from preceding reference R peak:")
        print(f"    median {med_off*1000:.0f} ms  "
              f"({100*med_frac:.0f}% of the RR interval)")
        print(f"    IQR    {iqr*1000:.0f} ms   "
              f"[{np.percentile(offs,25)*1000:.0f}-"
              f"{np.percentile(offs,75)*1000:.0f} ms]")
        tight = iqr < 0.080
        if tight and 0.15 <= med_off <= 0.35:
            phase_verdict = 'T wave'
            print("    => T WAVE. A tight mode at this offset is the T wave.")
        elif tight and med_off < 0.12:
            phase_verdict = 'S rebound / QRS twice'
            print("    => S-WAVE REBOUND. The QRS itself is being counted twice.")
        elif tight:
            phase_verdict = f'tight mode at {med_off*1000:.0f} ms'
            print(f"    => tight mode at {med_off*1000:.0f} ms, but not where a")
            print("       T wave or an S rebound sits. Investigate directly.")
        else:
            phase_verdict = 'no mode'
            print(f"    => NO TIGHT MODE (IQR {iqr*1000:.0f} ms). The extra")
            print("       detections are scattered through the cycle, so")
            print("       neither T-wave nor double-QRS explains them.")
            print("       Do not claim a mechanism that has not been measured.")
    else:
        print("  too few unmatched detections to identify a mode")

    # ---------------------------------------------------------- verdict
    bad = [r for r in worst_rows if r['ratio'] > 1.6]
    alts = [r['lag1_corr'] for r in worst_rows
            if not np.isnan(r['lag1_corr'])]
    mean_alt = float(np.mean(alts)) if alts else float('nan')
    mean_half = float(np.mean([r['half_cluster'] for r in worst_rows]))

    print("\n" + "=" * 62)
    print("  VERDICT")
    print("=" * 62)
    print(f"  worst windows with reference >1.6x detector : "
          f"{len(bad)}/{len(worst_rows)}")
    print(f"  mean lag-1 RR correlation in those windows  : {mean_alt:+.3f}")
    print(f"  mean half-median cluster                    : {100*mean_half:.1f}%")
    print()
    over = worst['ratio'] < 0.65
    if over:
        print("  => DETECTOR OVER-DETECTION. The detector reports roughly")
        print("     DOUBLE the reference rate on this channel, so it is")
        print("     counting non-QRS deflections — T waves being the usual")
        print("     culprit, since stage 5 only rejects a peak inside the")
        print("     360 ms window when its slope is under half the previous")
        print("     beat's. Neither original hypothesis: the reference is")
        print("     sound and the detector is adding beats, not missing them.")
        print("     Characterise it. Do NOT retune 02_qrs_detect.py — that")
        print("     file produced the MIT-BIH headline and the Phase A result.")
    elif len(bad) >= 3 and (mean_alt <= ALTERNATION_THRESHOLD or mean_half > 0.25):
        print("  => H1 SUPPORTED: the REFERENCE is doubled in these windows.")
        print("     Reference rate is ~2x the detector, intervals alternate")
        print("     short-long, and a population sits near half the record")
        print("     median. The detector is reporting the physiological rate.")
        print("     Report this as a measured property of the annotations.")
        print("     Do NOT drop the record and do NOT change the detector.")
    elif len(bad) >= 3:
        print("  => AMBIGUOUS: the reference reads ~2x the detector, but the")
        print("     alternation and cluster tests do not confirm doubling.")
        print("     The plot is now the arbiter — inspect it before writing")
        print("     anything up.")
    else:
        print("  => H1 NOT SUPPORTED. The discrepancy is not a 2x pattern, so")
        print("     doubling does not explain it. Treat H2 as live: the")
        print("     detector may be failing on this record, which is a")
        print("     serious finding the pooled numbers were masking.")
    print("\n  The plot settles it. Two marks on one QRS = doubled reference.")
    print("  A trace genuinely beating twice as fast as the detector claims")
    print("  = detector failure.")

    # ---------------------------------------------------------- plots
    os.makedirs(PLOT_DIR, exist_ok=True)
    k = common[int(order[0])]
    a, b = k * w, (k + 1) * w
    t = np.arange(a, b) / fs
    sig = rec['sig'][a:b, ch]
    rs = ref[(ref >= a) & (ref < b)]
    ds = det[(det >= a) & (det < b)]

    fig, axes = plt.subplots(4, 1, figsize=(14, 13))

    ax = axes[0]
    ax.plot(t, sig, lw=0.8, color='0.35')
    ax.plot(rs / fs, rec['sig'][rs, ch], 'o', ms=10, mfc='none',
            color='tab:green', label=f'reference ({rs.size})')
    ax.plot(ds / fs, rec['sig'][ds, ch], 'v', ms=7, color='tab:red',
            label=f'detected ({ds.size})')
    ax.set_title(f"WORST WINDOW {k} ({a/fs:.0f}-{b/fs:.0f} s) — "
                 f"ref {ref_w[k]:.0f} BPM vs det {det_w[k]:.0f} BPM. "
                 f"Two green marks on one QRS = doubled reference.",
                 loc='left', fontsize=9)
    ax.set_ylabel(f"{rec['sig_name'][ch]} ({rec['units'][ch]})")
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(alpha=0.3)

    ax = axes[1]
    zt = min(2.0, args.window)
    m = t < (a / fs + zt)
    rs2 = rs[rs < a + int(zt * fs)]
    ds2 = ds[ds < a + int(zt * fs)]
    ax.plot(t[m], sig[m], lw=1.2, color='0.35')
    ax.plot(rs2 / fs, rec['sig'][rs2, ch], 'o', ms=12, mfc='none', mew=2,
            color='tab:green')
    ax.plot(ds2 / fs, rec['sig'][ds2, ch], 'v', ms=9, color='tab:red')
    ax.set_title(f'Same window, first {zt:.0f} s — this is the decisive view',
                 loc='left', fontsize=9)
    ax.set_ylabel(rec['units'][ch])
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(ref[1:] / fs, rr, lw=0.6, color='tab:blue')
    ax.axhline(med, ls='--', lw=1, color='tab:green', label=f'median {med:.3f}s')
    ax.axhline(med / 2, ls='--', lw=1, color='tab:red',
               label=f'half median {med/2:.3f}s')
    ax.axvspan(a / fs, b / fs, color='tab:orange', alpha=0.25,
               label='worst window')
    ax.set_xlabel('seconds')
    ax.set_ylabel('RR interval (s)')
    ax.set_title('RR over the record — a band at half the median is doubling',
                 loc='left', fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[3]
    if offs.size:
        ax.hist(offs * 1000, bins=60, color='tab:purple', alpha=0.8)
        ax.axvspan(150, 350, color='tab:orange', alpha=0.15, label='T-wave band')
        ax.axvspan(40, 120, color='tab:blue', alpha=0.15, label='S-rebound band')
        ax.axvline(float(np.median(offs)) * 1000, color='tab:red', lw=1.2,
                   label=f'median {1000*np.median(offs):.0f} ms')
        ax.legend(fontsize=8)
    ax.set_xlabel('offset of unmatched detection from preceding reference R peak (ms)')
    ax.set_ylabel('count')
    ax.set_title('Test 3 — WHAT the extra detections are. This, not the '
                 'marker positions above, is the evidence.',
                 loc='left', fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle(f'{args.record} — reference doubling vs detector failure')
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f'06_diagnose_{args.record.replace(".edf","")}.png')
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n  plot -> {out}")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    res = os.path.join(RESULTS_DIR,
                       f'06_diagnose_{args.record.replace(".edf","")}.json')
    with open(res, 'w') as fh:
        json.dump({
            'record': args.record,
            'gestation': list(g) if g else None,
            'channel': rec['sig_name'][ch],
            'reference_beats': int(ref.size),
            'median_rr_s': med, 'median_bpm': 60 / med,
            'bpm_max': float(bpm.max()), 'bpm_min': float(bpm.min()),
            'bpm_p95': float(np.percentile(bpm, 95)),
            'record_lag1_corr': r_alt,
            'record_half_cluster': r_half,
            'windows_paired': len(common),
            'window_mae_bpm': float(np.mean(np.abs(errs))),
            'unmatched_detections': int(offs.size),
            'unmatched_offset_median_s': float(np.median(offs)) if offs.size else None,
            'unmatched_offset_iqr_s': float(np.percentile(offs, 75)
                                            - np.percentile(offs, 25)) if offs.size else None,
            'unmatched_frac_of_rr_median': float(np.median(frac)) if frac.size else None,
            'phase_verdict': phase_verdict,
            'worst_windows': worst_rows,
        }, fh, indent=2, default=float)
    print(f"  results -> {res}")


if __name__ == '__main__':
    main()
