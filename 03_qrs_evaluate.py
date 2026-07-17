"""
03_qrs_evaluate.py — E1: score the Pan-Tompkins detector against MIT-BIH truth.

Counts matching is not evidence. A detector can produce the right NUMBER of
beats in all the wrong PLACES. This script scores POSITIONS.

Metrics (ANSI/AAMI EC57):
    tolerance   : a detection counts as correct if within +/-150 ms of a
                  reference beat. Not exact equality: MIT-BIH annotations were
                  hand-corrected by cardiologists and mark the R-peak by human
                  convention, so a few ms of disagreement is expected.
    sensitivity : TP / (TP + FN)  -> of the REAL beats, how many did we find?
    PPV         : TP / (TP + FP)  -> of OUR detections, how many were real?
    DER         : (FP + FN) / reference beats  -> total error rate

Both sensitivity and PPV are required because either alone is gameable: detect
everything and sensitivity hits 100% while PPV collapses; detect only the one
obvious beat and PPV is perfect while sensitivity is nil.

Scoring convention for ventricular flutter (ANSI/AAMI EC57):
    Detections and reference beats falling inside an annotated flutter episode
    ('[' opens, ']' closes) are EXCLUDED from scoring. EC57 excludes
    flutter/fibrillation episodes from beat-by-beat comparison because "should a
    detector fire during VF?" has no agreed answer.

    Justified empirically by 04_flutter_audit.py, not assumed:
      - only record 207 carries brackets (6 episodes, 143 s, 7.9% of record)
      - all 472 '!' flutter waves fall inside them (100%)
      - 200 of 207's 202 false positives fall inside them (99%)
      - CRUCIALLY, the episodes contain NO beat labels at all: only '[', '+',
        '!', ']'. The MIT-BIH annotators declined to annotate beats during
        flutter. Scoring there means scoring against a deliberate blank and
        counting every detection as an error by default.
    So this exclusion removes no reference beats. It only stops penalising the
    detector for firing in a region the annotators refused to label.

    VERIFY the EC57 clause against the actual standard before citing it.

Record scope (decided BEFORE seeing any scores, so it cannot be motivated by
results):
    headline      : the 46 records that carry an MLII lead.
    supplementary : 102 and 104, which have no MLII (both are ['V5','V2']).
                    Excluded from the headline on LEAD-MISMATCH grounds: the
                    target hardware (AD8232) is a single-lead module in roughly
                    a Lead I/II placement, so MLII is the lead the deployed
                    system will actually see. Scoring a V5 chest lead measures
                    something the prototype will never produce. They are still
                    reported, separately, so nothing is hidden.
    NOTE: this is a different question from E2's paced-record exclusion
    (102, 104, 107, 217). Two decisions, two justifications. Keep them apart
    in the writeup.

The detector is IMPORTED from 02, never copied. If it were duplicated, tuning a
constant in 02 would silently leave 03 scoring a different algorithm than the
one shipped.
"""
import os
import importlib.util

import numpy as np
import pandas as pd

import wfdb

DETECTOR_FILE = '02_qrs_detect.py'
DATASET = 'mitdb'
DATA_DIR = os.path.join('data', DATASET)
RESULTS_DIR = 'results'

TOLERANCE_S = 0.150          # ANSI/AAMI EC57 matching window
NO_MLII = ['102', '104']     # confirmed by audit_leads.py

EXCLUDE_EPISODE_REGIONS = True   # convention (a); see module docstring
EPISODE_OPEN = '['
EPISODE_CLOSE = ']'


def find_episodes(ann, sig_len):
    """Pair '[' with the next ']' by walking annotations in order.

    Sequential walk rather than collecting starts and ends separately: that
    would mis-pair if an episode were unclosed. An unclosed '[' runs to the end
    of the record.

    Lives here, in the scorer, because it defines a SCORING CONVENTION. Analysis
    scripts import it from here rather than the reverse.
    """
    eps = []
    open_at = None
    for s, sym in zip(ann.sample, ann.symbol):
        if sym == EPISODE_OPEN and open_at is None:
            open_at = int(s)
        elif sym == EPISODE_CLOSE and open_at is not None:
            eps.append((open_at, int(s)))
            open_at = None
    if open_at is not None:
        eps.append((open_at, int(sig_len)))
    return eps


def in_any(positions, eps):
    """Boolean mask: which positions fall inside any episode."""
    positions = np.asarray(positions, dtype=np.int64)
    mask = np.zeros(positions.size, dtype=bool)
    for a, b in eps:
        mask |= (positions >= a) & (positions <= b)
    return mask


def load_detector():
    """Import 02_qrs_detect.py by path.

    A module name cannot start with a digit, so `import 02_qrs_detect` is a
    syntax error. Loading by file path keeps ONE definition of the detector.
    02 guards its main() behind __name__, so importing runs nothing.
    """
    spec = importlib.util.spec_from_file_location('qrs_detect', DETECTOR_FILE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def list_records(data_dir):
    return sorted(f[:-4] for f in os.listdir(data_dir) if f.endswith('.hea'))


def match_beats(ref, det, tol):
    """Greedy one-to-one matching of detections to reference beats.

    For each reference beat, claim the NEAREST detection within tolerance that
    no other reference beat has already claimed.

    The one-to-one constraint is load-bearing. Without it, a detector spraying
    five peaks around a single beat would score five true positives from one
    real beat, and PPV would flatter a detector that is actually misbehaving.

    Returns (tp, fp, fn, signed timing errors in samples, claimed mask).

    `claimed` is a boolean array over `det`: True where a reference beat matched
    that detection. So det[~claimed] IS the list of false positives. Exposed so
    downstream analysis can locate FPs without reimplementing this matching. Two
    copies of the matching logic would eventually diverge, and the reported
    metrics would stop describing the shipped detector.
    """
    ref = np.asarray(ref, dtype=np.int64)
    det = np.asarray(det, dtype=np.int64)
    if det.size == 0:
        return 0, 0, int(ref.size), np.array([]), np.zeros(0, dtype=bool)

    claimed = np.zeros(det.size, dtype=bool)
    tp = 0
    errors = []

    for r in ref:
        lo = np.searchsorted(det, r - tol, 'left')
        hi = np.searchsorted(det, r + tol, 'right')
        best, best_d = -1, tol + 1
        for j in range(lo, hi):
            if claimed[j]:
                continue
            d = abs(int(det[j]) - int(r))
            if d < best_d:
                best, best_d = j, d
        if best >= 0:
            claimed[best] = True
            tp += 1
            errors.append(int(det[best]) - int(r))

    fn = int(ref.size) - tp
    fp = int(det.size) - tp       # every unclaimed detection is a false positive
    return tp, fp, fn, np.array(errors), claimed


def evaluate_record(qrs, name, exclude_episodes=None):
    if exclude_episodes is None:
        exclude_episodes = EXCLUDE_EPISODE_REGIONS

    path = os.path.join(DATA_DIR, name)
    rec = wfdb.rdrecord(path)
    fs = rec.fs
    ch = qrs.find_channel(rec, 'MLII', 'II')
    lead = rec.sig_name[ch].strip()
    sig = rec.p_signal[:, ch]

    det, _ = qrs.detect_qrs(sig, fs)

    ann = wfdb.rdann(path, 'atr')
    keep = np.array([s not in qrs.NON_BEAT for s in ann.symbol])
    ref = ann.sample[keep]

    # Flutter-region exclusion. Drops detections AND references inside
    # annotated episodes. In practice only 207 is affected, and only its
    # detections are dropped, because the episodes carry no beat labels.
    excluded_s = 0.0
    if exclude_episodes:
        eps = find_episodes(ann, sig.size)
        if eps:
            excluded_s = sum(b - a for a, b in eps) / fs
            ref = ref[~in_any(ref, eps)]
            det = det[~in_any(det, eps)]

    tol = int(round(TOLERANCE_S * fs))
    tp, fp, fn, err, _ = match_beats(ref, det, tol)

    se = tp / (tp + fn) if (tp + fn) else 0.0
    ppv = tp / (tp + fp) if (tp + fp) else 0.0
    der = (fp + fn) / ref.size if ref.size else 0.0
    mae_ms = float(np.mean(np.abs(err)) / fs * 1000) if err.size else np.nan

    return {
        'record': name,
        'lead': lead,
        'ref_beats': int(ref.size),
        'detected': int(det.size),
        'TP': tp, 'FP': fp, 'FN': fn,
        'sensitivity': se,
        'PPV': ppv,
        'DER': der,
        'timing_mae_ms': mae_ms,
        'excluded_episode_s': excluded_s,
    }


def summarise(df, label):
    """Report gross AND average, because they answer different questions.

    gross   : pool every beat across records, then compute the ratio. This is
              what the literature quotes. Long/busy records pull it harder.
    average : mean of the per-record scores. Every record counts equally, so a
              single disastrous record is visible instead of being drowned out.
    A large gap between the two is itself a finding: it means the failures are
    concentrated in a few records rather than spread thin.
    """
    tp, fp, fn = df['TP'].sum(), df['FP'].sum(), df['FN'].sum()
    g_se = tp / (tp + fn) if (tp + fn) else 0.0
    g_ppv = tp / (tp + fp) if (tp + fp) else 0.0

    print(f"\n=== {label} ({len(df)} records) ===")
    print(f"  reference beats   : {int(df['ref_beats'].sum()):,}")
    print(f"  TP / FP / FN      : {int(tp):,} / {int(fp):,} / {int(fn):,}")
    print(f"  gross sensitivity : {100*g_se:.2f}%")
    print(f"  gross PPV         : {100*g_ppv:.2f}%")
    print(f"  gross DER         : {100*(fp+fn)/df['ref_beats'].sum():.2f}%")
    print(f"  avg  sensitivity  : {100*df['sensitivity'].mean():.2f}%")
    print(f"  avg  PPV          : {100*df['PPV'].mean():.2f}%")
    print(f"  median timing MAE : {df['timing_mae_ms'].median():.1f} ms")


def main():
    qrs = load_detector()
    all_recs = list_records(DATA_DIR)
    mlii_recs = [r for r in all_recs if r not in NO_MLII]

    print(f"Records total       : {len(all_recs)}")
    print(f"Headline (MLII)     : {len(mlii_recs)}")
    print(f"Supplementary (V5)  : {NO_MLII}")
    print(f"Tolerance           : +/-{TOLERANCE_S*1000:.0f} ms")
    print(f"Flutter regions     : "
          f"{'EXCLUDED (convention a)' if EXCLUDE_EPISODE_REGIONS else 'scored (convention c)'}\n")

    rows = []
    for i, name in enumerate(all_recs, 1):
        r = evaluate_record(qrs, name)
        rows.append(r)
        tag = '' if name not in NO_MLII else '  [supplementary]'
        print(f"  [{i:2d}/{len(all_recs)}] {name}  {r['lead']:>4s}  "
              f"Se {100*r['sensitivity']:6.2f}%  PPV {100*r['PPV']:6.2f}%  "
              f"FN {r['FN']:4d}  FP {r['FP']:4d}{tag}")

    df = pd.DataFrame(rows)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out = os.path.join(RESULTS_DIR, f'{DATASET}_e1_eval.csv')
    df.to_csv(out, index=False)

    head = df[~df['record'].isin(NO_MLII)]
    supp = df[df['record'].isin(NO_MLII)]

    summarise(head, 'E1 HEADLINE — MLII records')
    if len(supp):
        summarise(supp, 'E1 SUPPLEMENTARY — no MLII, V5 fallback')

    print("\n=== Worst 5 headline records by DER ===")
    worst = head.nlargest(5, 'DER')[
        ['record', 'ref_beats', 'sensitivity', 'PPV', 'FN', 'FP', 'DER']]
    for _, r in worst.iterrows():
        print(f"  {r['record']}: Se {100*r['sensitivity']:6.2f}%  "
              f"PPV {100*r['PPV']:6.2f}%  FN {int(r['FN']):4d}  "
              f"FP {int(r['FP']):4d}  DER {100*r['DER']:5.2f}%")

    print(f"\n  results -> {out}")
    print("\n  Expect a few records to score badly. That is the correct outcome,")
    print("  not a bug: 108 has severe baseline wander with tall P waves that")
    print("  mimic QRS, and 207 contains ventricular flutter where discrete")
    print("  beats barely exist. Published detectors struggle on these too.")
    print("  A CHARACTERISED detector beats a suspiciously perfect average.")


if __name__ == '__main__':
    main()
