"""
04_flutter_audit.py — test one hypothesis, numerically, and end the argument.

HYPOTHESIS
  Record 207's 202 false positives are largely an artifact of OUR annotation
  filtering, not a detector failure.

  Reasoning: '!' (ventricular flutter wave) is in NON_BEAT, so we DELETE those
  472 annotations from the reference. But the detector still fires on them,
  because a flutter wave is a real deflection with QRS-like slope. Every such
  firing becomes a guaranteed false positive. We would be penalising the
  detector for finding something genuinely present in the signal.

  MIT-BIH brackets these episodes explicitly: '[' opens, ']' closes.

THREE DEFENSIBLE CONVENTIONS
  (a) exclude the whole '['..']' region from scoring   <- ANSI/AAMI EC57 aligned
  (b) treat '!' as reference beats                     <- defensible
  (c) drop '!' from reference, keep detections         <- what 03 currently does

  EC57 excludes ventricular flutter/fibrillation episodes from beat-by-beat
  comparison, because "should a detector fire during VF?" has no agreed answer.
  VERIFY THIS CLAUSE AGAINST THE STANDARD before citing it in the writeup.

WHAT THIS DECIDES
  If most of 207's FPs fall inside the bracketed episodes, (a) is justified and
  the fix is principled. If they do not, the hypothesis is wrong, we keep (c),
  and we move on. Either outcome is a result. This script does not assume the
  answer.

Detector and matching are IMPORTED, never copied.
"""
import os
import importlib.util
import collections

import numpy as np
import pandas as pd

import wfdb

EVAL_FILE = '03_qrs_evaluate.py'
DATASET = 'mitdb'
DATA_DIR = os.path.join('data', DATASET)
RESULTS_DIR = 'results'

RECORD = '207'


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def list_records(data_dir):
    return sorted(f[:-4] for f in os.listdir(data_dir) if f.endswith('.hea'))


def scan_all_records(recs, ev):
    """Confirm which records carry episode brackets at all.

    Discovery reported 6 '[' and 6 ']' database-wide. If they are all in 207,
    this is a single-record methodology question, not a database-wide one.
    """
    print("=== Records containing episode brackets ===")
    found = {}
    for name in recs:
        try:
            ann = wfdb.rdann(os.path.join(DATA_DIR, name), 'atr')
        except Exception:
            continue
        c = collections.Counter(ann.symbol)
        n_open, n_close = c.get(ev.EPISODE_OPEN, 0), c.get(ev.EPISODE_CLOSE, 0)
        if n_open or n_close:
            found[name] = (n_open, n_close)
            print(f"  {name}: {n_open} x '{ev.EPISODE_OPEN}', "
                  f"{n_close} x '{ev.EPISODE_CLOSE}'")
    if not found:
        print("  none")
    print(f"  -> {len(found)} record(s) affected\n")
    return found


def main():
    ev = load_module(EVAL_FILE, 'qrs_evaluate')
    qrs = ev.load_detector()

    # 03 owns the convention: episode helpers and bracket symbols come from it.
    find_episodes, in_any = ev.find_episodes, ev.in_any

    recs = list_records(DATA_DIR)
    scan_all_records(recs, ev)

    path = os.path.join(DATA_DIR, RECORD)
    rec = wfdb.rdrecord(path)
    fs = rec.fs
    ch = qrs.find_channel(rec, 'MLII', 'II')
    sig = rec.p_signal[:, ch]
    ann = wfdb.rdann(path, 'atr')

    # --- episodes ---
    eps = find_episodes(ann, sig.size)
    print(f"=== Record {RECORD}: flutter episodes ===")
    total_dur = 0
    for i, (a, b) in enumerate(eps, 1):
        dur = (b - a) / fs
        total_dur += dur
        print(f"  episode {i}: {a/fs:7.1f} s -> {b/fs:7.1f} s   "
              f"({dur:6.1f} s)")
    rec_dur = sig.size / fs
    print(f"\n  episodes            : {len(eps)}")
    print(f"  total episode time  : {total_dur:.1f} s "
          f"({100*total_dur/rec_dur:.1f}% of record)")

    # --- what lives inside the episodes ---
    ann_in = in_any(ann.sample, eps)
    inside = collections.Counter(np.array(ann.symbol)[ann_in])
    outside = collections.Counter(np.array(ann.symbol)[~ann_in])
    print(f"\n  annotation symbols INSIDE episodes  : {dict(inside)}")
    print(f"  annotation symbols OUTSIDE episodes : {dict(outside)}")

    n_bang_in = inside.get('!', 0)
    n_bang_tot = collections.Counter(ann.symbol).get('!', 0)
    if n_bang_tot:
        print(f"  '!' flutter waves: {n_bang_in}/{n_bang_tot} inside episodes "
              f"({100*n_bang_in/n_bang_tot:.1f}%)")

    # --- detect and match, exactly as 03 does ---
    # Deliberately reproduce convention (c) here, ignoring 03's exclusion
    # default, so the two conventions can be compared side by side.
    det, _ = qrs.detect_qrs(sig, fs)
    keep = np.array([s not in qrs.NON_BEAT for s in ann.symbol])
    ref = ann.sample[keep]
    tol = int(round(ev.TOLERANCE_S * fs))
    tp, fp, fn, err, claimed = ev.match_beats(ref, det, tol)

    print(f"\n=== Current scoring (convention (c), as in 03) ===")
    print(f"  reference beats : {ref.size}")
    print(f"  detected        : {det.size}")
    print(f"  TP / FP / FN    : {tp} / {fp} / {fn}")
    print(f"  sensitivity     : {100*tp/(tp+fn):.2f}%")
    print(f"  PPV             : {100*tp/(tp+fp):.2f}%")

    # --- THE TEST ---
    fps = det[~claimed]
    fp_in = in_any(fps, eps)
    print(f"\n=== THE TEST: where do the false positives fall? ===")
    print(f"  false positives total   : {fps.size}")
    print(f"  FP inside episodes      : {int(fp_in.sum())} "
          f"({100*fp_in.mean():.1f}%)" if fps.size else "  (no FPs)")
    print(f"  FP outside episodes     : {int((~fp_in).sum())}")

    if len(eps):
        print(f"\n  per-episode FP counts:")
        for i, (a, b) in enumerate(eps, 1):
            n = int(((fps >= a) & (fps <= b)).sum())
            print(f"    episode {i} ({(b-a)/fs:6.1f} s): {n:4d} FP")

    # --- what convention (a) would give ---
    ref_out = ref[~in_any(ref, eps)]
    det_out = det[~in_any(det, eps)]
    tp2, fp2, fn2, _, _ = ev.match_beats(ref_out, det_out, tol)
    se2 = tp2 / (tp2 + fn2) if (tp2 + fn2) else 0.0
    ppv2 = tp2 / (tp2 + fp2) if (tp2 + fp2) else 0.0

    print(f"\n=== Convention (a): exclude episode regions entirely ===")
    print(f"  reference beats : {ref_out.size}  (was {ref.size})")
    print(f"  detected        : {det_out.size}  (was {det.size})")
    print(f"  TP / FP / FN    : {tp2} / {fp2} / {fn2}")
    print(f"  sensitivity     : {100*se2:.2f}%  (was {100*tp/(tp+fn):.2f}%)")
    print(f"  PPV             : {100*ppv2:.2f}%  (was {100*tp/(tp+fp):.2f}%)")

    # --- knock-on effect for the 46-record headline ---
    csv = os.path.join(RESULTS_DIR, f'{DATASET}_e1_eval.csv')
    if os.path.exists(csv):
        df = pd.read_csv(csv, dtype={'record': str})
        head = df[~df['record'].isin(ev.NO_MLII)].copy()
        g_tp, g_fp, g_fn = head['TP'].sum(), head['FP'].sum(), head['FN'].sum()
        old_se = g_tp / (g_tp + g_fn)
        old_ppv = g_tp / (g_tp + g_fp)

        m = head['record'] == RECORD
        head.loc[m, ['TP', 'FP', 'FN']] = [tp2, fp2, fn2]
        n_tp, n_fp, n_fn = head['TP'].sum(), head['FP'].sum(), head['FN'].sum()
        new_se = n_tp / (n_tp + n_fn)
        new_ppv = n_tp / (n_tp + n_fp)

        print(f"\n=== Effect on the 46-record headline ===")
        print(f"  gross sensitivity : {100*old_se:.2f}%  ->  {100*new_se:.2f}%")
        print(f"  gross PPV         : {100*old_ppv:.2f}%  ->  {100*new_ppv:.2f}%")
        print(f"  total FP          : {int(g_fp)}  ->  {int(n_fp)}")
    else:
        print(f"\n  ({csv} not found: skipping headline recomputation)")

    print(f"\n=== Verdict ===")
    if fps.size:
        frac = fp_in.mean()
        if frac >= 0.5:
            print(f"  {100*frac:.0f}% of {RECORD}'s false positives fall inside")
            print("  annotated flutter episodes. The hypothesis HOLDS: these are")
            print("  largely an artifact of scoring convention (c), not detector")
            print("  failure. Convention (a) is justified and principled.")
        else:
            print(f"  Only {100*frac:.0f}% of false positives fall inside episodes.")
            print("  The hypothesis FAILS. The FPs are mostly elsewhere, so this")
            print("  is a real detector limitation, not a scoring artifact.")
            print("  Keep convention (c), report as-is, move on.")
    print("\n  Whichever way it lands, state the convention explicitly in the")
    print("  writeup. An unstated scoring convention is the thing that makes a")
    print("  number unreproducible.")


if __name__ == '__main__':
    main()
