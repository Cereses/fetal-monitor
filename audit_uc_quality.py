"""
audit_uc_quality.py — diagnose why explore_ctu.py's flat-fraction metric
condemns most of the UC channel, and replace it with one that survives contact
with the expert annotations.

THE PROBLEM
-----------
explore_ctu.py reported median flat fraction 60.9% and concluded 18/50 records
usable. Cross-referencing against the Romagnoli annotations shows 42/50 records
carry expert-marked contractions, including record 1017 with 25 of them. The
metric is wrong, not the data.

WHY IT IS WRONG
---------------
The UC channel is coarsely quantised: integer-valued over roughly 0-127. A
contraction ramps over ~60 s. Sampled at 4 Hz, a rise of 30 units across a
contraction is one integer increment every ~8 samples, so ~7 of every 8
consecutive pairs are IDENTICAL by construction. A per-sample "did it change"
test therefore measures the quantisation step, not signal health.

The original comment in explore_ctu.py had this backwards: it assumed a slow
ramp guarantees local variation. Slow ramp plus coarse quantisation guarantees
the opposite.

THE FIX
-------
Genuine sensor dropout does not look like quantisation. Quantisation produces
SHORT flat runs (a few samples between increments). Dropout produces LONG flat
runs (tens of seconds of a frozen or zeroed value, belt loosened or off).

So measure RUN LENGTH, not per-sample equality:

  quantisation  : flat runs of a handful of samples, spread evenly
  dropout       : flat runs of >= DROPOUT_MIN_S seconds

Reported per record:
  q_step        median non-zero |diff|, i.e. the actual quantisation step
  flat_frac     the OLD metric, kept so the two can be compared directly
  dropout_frac  fraction of samples inside flat runs >= DROPOUT_MIN_S
  zero_frac     fraction of samples exactly zero
  zerorun_frac  fraction inside ZERO runs >= DROPOUT_MIN_S (dropout sentinel)
  max_run_s     longest single flat run, in seconds

VALIDATION, NOT INVENTION
-------------------------
A new threshold picked by eye is no better than the old one. This script scores
every candidate rule against the annotation manifest and reports the confusion
directly: how many annotated records a rule would discard (the costly error)
and how many unannotated ones it would keep. Pick the rule from that table,
then FREEZE it before the detector runs.

Run parse_ctu_annotations.py first: this reads results/ctu_ann_manifest.csv.
"""
import os
import csv

import numpy as np
import wfdb

DATASET = 'ctu-uhb-ctgdb'
DATA_DIR = os.path.join('data', DATASET)
RESULTS_DIR = 'results'
MANIFEST = os.path.join(RESULTS_DIR, 'ctu_ann_manifest.csv')

FS = 4.0
DROPOUT_MIN_S = 10.0                       # a flat stretch this long is not physiology
DROPOUT_MIN_N = int(DROPOUT_MIN_S * FS)    # 40 samples


def find_channel(rec, *targets):
    cleaned = [n.strip() for n in rec.sig_name]
    for t in targets:
        if t in cleaned:
            return cleaned.index(t)
    return 0


def list_records(data_dir):
    return sorted(f[:-4] for f in os.listdir(data_dir) if f.endswith('.hea'))


def flat_runs(x):
    """(start, length) of every maximal run of identical consecutive values.

    Working in run space is the whole point: it is what separates a quantised
    signal (many short runs) from a dead sensor (one enormous run).
    """
    if x.size == 0:
        return np.empty(0, int), np.empty(0, int)
    change = np.flatnonzero(np.diff(x) != 0) + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [x.size]))
    return starts, ends - starts


def main():
    records = list_records(DATA_DIR)
    print(f"Auditing UC quality on {len(records)} records\n")

    ann_n = {}
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as f:
            for row in csv.DictReader(f):
                ann_n[row['record']] = int(row['n_contractions'])
    else:
        print(f"  NOTE: {MANIFEST} not found. Run parse_ctu_annotations.py to\n"
              f"  enable annotation-based validation. Continuing without it.\n")

    rows = []
    for name in records:
        rec = wfdb.rdrecord(os.path.join(DATA_DIR, name))
        uc = rec.p_signal[:, find_channel(rec, 'UC')].astype(float)
        finite = np.isfinite(uc)
        nan_frac = float((~finite).mean())
        uc = np.where(finite, uc, np.nan)
        v = uc[finite]

        d = np.abs(np.diff(v))
        nz = d[d > 0]
        q_step = float(np.median(nz)) if nz.size else 0.0

        starts, lens = flat_runs(v)
        long_mask = lens >= DROPOUT_MIN_N
        dropout_frac = float(lens[long_mask].sum() / v.size) if v.size else 1.0
        max_run_s = float(lens.max() / FS) if lens.size else 0.0

        zero_frac = float((v == 0).mean()) if v.size else 1.0
        zero_run = sum(
            L for s, L in zip(starts, lens) if L >= DROPOUT_MIN_N and v[s] == 0
        )
        zerorun_frac = float(zero_run / v.size) if v.size else 1.0

        rows.append(dict(
            record=name,
            q_step=q_step,
            flat_frac=float((d == 0).mean()) if d.size else 1.0,   # OLD metric
            dropout_frac=dropout_frac,
            zero_frac=zero_frac,
            zerorun_frac=zerorun_frac,
            max_run_s=max_run_s,
            nan_frac=nan_frac,
            n_contractions=ann_n.get(name, -1),
        ))

    out = os.path.join(RESULTS_DIR, 'ctu_uc_quality.csv')
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    q = np.array([r['q_step'] for r in rows])
    old = np.array([r['flat_frac'] for r in rows])
    new = np.array([r['dropout_frac'] for r in rows])
    zr = np.array([r['zerorun_frac'] for r in rows])
    mx = np.array([r['max_run_s'] for r in rows])

    print("=" * 66)
    print("QUANTISATION")
    print("=" * 66)
    print(f"  median quantisation step : {np.median(q):.3f}")
    print(f"  range                    : {q.min():.3f} - {q.max():.3f}")
    print("  If this is ~1.0 the channel is integer-valued, which is the whole")
    print("  explanation for the old flat-fraction reading.\n")

    print("=" * 66)
    print("OLD METRIC vs NEW METRIC")
    print("=" * 66)
    print(f"  flat_frac    (per-sample) : median {np.median(old)*100:5.1f}%  max {old.max()*100:5.1f}%")
    print(f"  dropout_frac (runs >={DROPOUT_MIN_S:.0f}s) : median {np.median(new)*100:5.1f}%  max {new.max()*100:5.1f}%")
    print(f"  zerorun_frac              : median {np.median(zr)*100:5.1f}%  max {zr.max()*100:5.1f}%")
    print(f"  longest flat run          : median {np.median(mx):.0f}s  max {mx.max():.0f}s\n")

    annotated = [r for r in rows if r['n_contractions'] > 0]
    if not annotated:
        print("  No annotation data loaded, stopping before rule validation.")
        print(f"  wrote {out}")
        return

    unann = [r for r in rows if r['n_contractions'] == 0]
    print("=" * 66)
    print(f"RULE VALIDATION  ({len(annotated)} annotated, {len(unann)} unannotated)")
    print("=" * 66)
    print("  A rule DISCARDS a record when the value exceeds the threshold.")
    print("  'lost' = annotated records discarded. That is the expensive error:")
    print("  it throws away ground truth you already paid for.\n")
    print(f"  {'rule':<34}{'lost':<14}{'unann. caught'}")

    def report(label, key, thr):
        lost = sum(1 for r in annotated if r[key] > thr)
        caught = sum(1 for r in unann if r[key] > thr)
        print(f"  {label:<34}{lost:>3}/{len(annotated):<10}{caught:>3}/{len(unann)}")

    report("OLD flat_frac > 0.50", 'flat_frac', 0.50)
    for t in (0.20, 0.30, 0.50):
        report(f"dropout_frac > {t:.2f}", 'dropout_frac', t)
    for t in (0.20, 0.30, 0.50):
        report(f"zerorun_frac > {t:.2f}", 'zerorun_frac', t)

    print("\n  Note what the annotations already do for free: the no_annotation")
    print("  exclusion removes every record with no ground truth. A quality rule")
    print("  only earns its place if it catches something that rule misses")
    print("  WITHOUT discarding annotated records.")
    print(f"\n  wrote {out}")


if __name__ == '__main__':
    main()
