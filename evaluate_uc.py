"""
evaluate_uc.py — score the frozen UC contraction detector against the
Romagnoli expert annotations.

DEV / TEST SPLIT
----------------
MIN_PEAK_SEP_S and MIN_DUR_S were set from percentiles of the annotated
duration and interval distributions, which were computed over the whole corpus.
That is mild leakage. It is handled the same way as the DS1/DS2 split in the
beat-classification work: records are split into disjoint DEV and TEST sets,
constants stay frozen, and the TEST number is the one that gets reported.

The split is deterministic and structural, not random: sorted record IDs,
alternating. Alternating rather than first-half/second-half because record IDs
in CTU-CHB may correlate with acquisition period, and a contiguous split would
confound the comparison with whatever changed over time in the delivery ward.

MATCHING
--------
Detections and references are time SPANS, not points, so "did it match" needs a
rule. Three are reported, because a result that only holds under one rule is a
result about the rule:

  any        any overlap at all                        (most lenient)
  half       overlap >= 50% of the SHORTER span        (PRIMARY)
  iou        intersection over union >= 0.5            (strictest)

`half` is primary because it asks the question that matters clinically — did
the detector find this contraction — without punishing boundary disagreement,
which is reported separately as a bias instead. Matching is greedy 1:1 by
descending overlap, so one detection can never claim two references. Same
discipline as match_beats in the QRS work.

WHAT IS NOT EXCLUDED
--------------------
Reference contractions that fall inside dropout regions are NOT dropped from
the denominator. The detector cannot find them, so they count as misses and
depress sensitivity. Removing them would flatter the result. They are counted
and reported separately so the size of that effect is visible.

FALSE POSITIVE DIAGNOSIS
------------------------
Record 1006 has obvious contractions nobody annotated, and record 1029 has
dozens of unmarked wobbles. Some false positives will therefore be real events
the annotator declined to mark, not detector errors. The script compares the
peak amplitude of matched against unmatched detections. If unmatched ones are
systematically smaller, the gap is an annotation-criterion mismatch and should
be reported as a finding, not tuned away. Structural diagnosis before any
attempt at a fix.
"""
import os
import csv
from collections import defaultdict

import numpy as np
import wfdb

from uc_detector import (
    detect_contractions, confidence, find_channel,
    DATA_DIR, RESULTS_DIR, FS,
)

PRIMARY = 'half'
CRITERIA = ('any', 'half', 'iou')


def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def accepts(a, b, criterion):
    ov = overlap(a, b)
    if ov <= 0:
        return False
    if criterion == 'any':
        return True
    if criterion == 'half':
        return ov >= 0.5 * min(a[1] - a[0], b[1] - b[0])
    union = (a[1] - a[0]) + (b[1] - b[0]) - ov
    return union > 0 and ov / union >= 0.5


def match_spans(ref, det, criterion):
    """Greedy 1:1 assignment by descending overlap.

    Returns (pairs, unmatched_ref_idx, unmatched_det_idx).
    """
    cands = []
    for i, r in enumerate(ref):
        for j, d in enumerate(det):
            if accepts(r, d, criterion):
                cands.append((overlap(r, d), i, j))
    cands.sort(key=lambda x: -x[0])

    used_r, used_d, pairs = set(), set(), []
    for _, i, j in cands:
        if i in used_r or j in used_d:
            continue
        used_r.add(i)
        used_d.add(j)
        pairs.append((i, j))
    return (pairs,
            [i for i in range(len(ref)) if i not in used_r],
            [j for j in range(len(det)) if j not in used_d])


def prf(tp, fp, fn):
    se = tp / (tp + fn) if tp + fn else float('nan')
    ppv = tp / (tp + fp) if tp + fp else float('nan')
    f1 = 2 * se * ppv / (se + ppv) if se + ppv else float('nan')
    return se, ppv, f1


def load_reference():
    path = os.path.join(RESULTS_DIR, 'ctu_contractions.csv')
    out = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f):
            out[r['record']].append((int(r['start_sample']), int(r['end_sample'])))
    return out


def load_excluded():
    path = os.path.join(RESULTS_DIR, 'ctu_ann_manifest.csv')
    with open(path) as f:
        return {r['record'] for r in csv.DictReader(f)
                if r['excluded'].lower() in ('true', '1')}


def main():
    ref_all = load_reference()
    excluded = load_excluded()
    available = sorted(f[:-4] for f in os.listdir(DATA_DIR) if f.endswith('.hea'))
    scorable = [r for r in available if r not in excluded and ref_all.get(r)]

    print(f"records on disk        : {len(available)}")
    print(f"excluded (no annotation): {len(available) - len(scorable)}")
    print(f"scorable                : {len(scorable)}")
    if len(scorable) < 100:
        print("\n  NOTE: only part of the corpus is downloaded. 473 records are\n"
              "  scorable in total. Pull the rest before quoting a headline number.")

    split = {r: ('DEV' if i % 2 == 0 else 'TEST')
             for i, r in enumerate(scorable)}

    per_record, detections = [], []
    boundary = {'DEV': [], 'TEST': []}
    peak_tp = {'DEV': [], 'TEST': []}
    peak_fp = {'DEV': [], 'TEST': []}
    fn_in_dropout = {'DEV': 0, 'TEST': 0}
    counts = {c: {'DEV': [0, 0, 0], 'TEST': [0, 0, 0]} for c in CRITERIA}

    for name in scorable:
        grp = split[name]
        rec = wfdb.rdrecord(os.path.join(DATA_DIR, name))
        uc = rec.p_signal[:, find_channel(rec, 'UC')].astype(float)
        det, d = detect_contractions(uc)
        conf, frac, rate = confidence(det, d['valid'])
        ref = ref_all[name]

        for crit in CRITERIA:
            pairs, ur, ud = match_spans(ref, det, crit)
            c = counts[crit][grp]
            c[0] += len(pairs)
            c[1] += len(ud)
            c[2] += len(ur)

        pairs, ur, ud = match_spans(ref, det, PRIMARY)
        tp, fp, fn = len(pairs), len(ud), len(ur)
        se, ppv, f1 = prf(tp, fp, fn)

        for i, j in pairs:
            boundary[grp].append((
                (det[j][0] - ref[i][0]) / FS,
                (det[j][1] - ref[i][1]) / FS,
                ((det[j][1] - det[j][0]) - (ref[i][1] - ref[i][0])) / FS,
            ))

        # amplitude of matched vs unmatched detections, relative to the record
        amp = d.get('amp', np.nan)
        detr = d.get('detr')
        if detr is not None and np.isfinite(amp) and amp > 0:
            matched_j = {j for _, j in pairs}
            for j, (s, e) in enumerate(det):
                rel = float(detr[s:e].max() / amp)
                (peak_tp if j in matched_j else peak_fp)[grp].append(rel)
                detections.append(dict(
                    record=name, split=grp, start_sample=s, end_sample=e,
                    duration_s=(e - s) / FS, peak_rel=rel,
                    matched=int(j in matched_j)))

        dead = ~d['valid']
        for i in ur:
            s, e = ref[i]
            if dead[s:e].mean() > 0.5:
                fn_in_dropout[grp] += 1

        per_record.append(dict(
            record=name, split=grp, confidence=conf,
            analysable_frac=frac, rate_per_10min=rate,
            n_ref=len(ref), n_det=len(det), tp=tp, fp=fp, fn=fn,
            sensitivity=se, ppv=ppv, f1=f1))

    with open(os.path.join(RESULTS_DIR, 'uc_evaluation.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(per_record[0].keys()))
        w.writeheader()
        w.writerows(per_record)
    if detections:
        with open(os.path.join(RESULTS_DIR, 'uc_detections.csv'), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(detections[0].keys()))
            w.writeheader()
            w.writerows(detections)

    # ---------------------------------------------------------------- report
    bar = "=" * 68
    print(f"\n{bar}\nUC CONTRACTION DETECTION\n{bar}")
    print(f"  split: alternating sorted record IDs")
    for g in ('DEV', 'TEST'):
        n = sum(1 for r in per_record if r['split'] == g)
        print(f"    {g:<5}: {n} records")

    print(f"\n{bar}\nPOOLED, criterion '{PRIMARY}' (overlap >= 50% of shorter span)\n{bar}")
    print(f"  {'':<8}{'ref':>7}{'TP':>7}{'FP':>7}{'FN':>7}{'Se':>9}{'PPV':>9}{'F1':>9}")
    for g in ('DEV', 'TEST'):
        tp, fp, fn = counts[PRIMARY][g]
        se, ppv, f1 = prf(tp, fp, fn)
        print(f"  {g:<8}{tp+fn:>7}{tp:>7}{fp:>7}{fn:>7}"
              f"{100*se:>8.2f}%{100*ppv:>8.2f}%{100*f1:>8.2f}%")

    print(f"\n{bar}\nSENSITIVITY TO THE MATCHING RULE (TEST)\n{bar}")
    print(f"  {'criterion':<12}{'TP':>7}{'FP':>7}{'FN':>7}{'Se':>9}{'PPV':>9}")
    for crit in CRITERIA:
        tp, fp, fn = counts[crit]['TEST']
        se, ppv, _ = prf(tp, fp, fn)
        print(f"  {crit:<12}{tp:>7}{fp:>7}{fn:>7}{100*se:>8.2f}%{100*ppv:>8.2f}%")

    print(f"\n{bar}\nBY CONFIDENCE TIER (TEST, criterion '{PRIMARY}')\n{bar}")
    for tier in ('HIGH', 'LOW'):
        rs = [r for r in per_record if r['split'] == 'TEST' and r['confidence'] == tier]
        if not rs:
            continue
        tp = sum(r['tp'] for r in rs); fp = sum(r['fp'] for r in rs)
        fn = sum(r['fn'] for r in rs)
        se, ppv, f1 = prf(tp, fp, fn)
        print(f"  {tier:<6} {len(rs):>4} records  ref {tp+fn:>5}  "
              f"Se {100*se:>6.2f}%  PPV {100*ppv:>6.2f}%  F1 {100*f1:>6.2f}%")

    b = np.array(boundary['TEST'])
    if b.size:
        print(f"\n{bar}\nBOUNDARY AGREEMENT (TEST, matched pairs)\n{bar}")
        for k, label in enumerate(('onset error', 'offset error', 'duration error')):
            v = b[:, k]
            print(f"  {label:<16} median {np.median(v):>7.1f} s   "
                  f"IQR {np.percentile(v,25):>6.1f} to {np.percentile(v,75):>6.1f} s")
        print("  Positive onset error = detector starts late.")

    print(f"\n{bar}\nFALSE POSITIVE DIAGNOSIS (TEST)\n{bar}")
    tpv, fpv = np.array(peak_tp['TEST']), np.array(peak_fp['TEST'])
    if tpv.size and fpv.size:
        print(f"  peak amplitude relative to record scale:")
        print(f"    matched   (n={tpv.size:>5}) : median {np.median(tpv):.2f}  "
              f"IQR {np.percentile(tpv,25):.2f} to {np.percentile(tpv,75):.2f}")
        print(f"    unmatched (n={fpv.size:>5}) : median {np.median(fpv):.2f}  "
              f"IQR {np.percentile(fpv,25):.2f} to {np.percentile(fpv,75):.2f}")
        if np.median(fpv) < np.median(tpv):
            print("  Unmatched detections are systematically smaller. Consistent")
            print("  with the annotator declining to mark low-amplitude events,")
            print("  which caps achievable PPV independently of the detector.")
    print(f"\n  reference contractions missed inside dropout regions: "
          f"{fn_in_dropout['TEST']} (TEST), counted as FN, not excluded")

    print(f"\n  wrote {os.path.join(RESULTS_DIR, 'uc_evaluation.csv')}")
    print(f"  wrote {os.path.join(RESULTS_DIR, 'uc_detections.csv')}")


if __name__ == '__main__':
    main()
