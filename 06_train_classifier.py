"""
06_train_classifier.py — E2: Random Forest beat classification, inter-patient.

THE SPLIT IS THE WHOLE METHODOLOGY.
DS1 trains, DS2 tests, no patient in both. Record lists are de Chazal et al.
(2004), the standard inter-patient partition, verified against five independent
published sources. Using the published split rather than inventing one means
these numbers are directly comparable to a large literature, and it means the
split was not chosen by us after seeing results.

Why not a shuffled split: heartbeats from one person are near-identical to each
other. Shuffle them and the model learns "I have seen this exact patient's
ventricular beats 300 times, here is the 301st". It memorises the patient, and
the test set rewards it for that, typically ~98% accuracy that means nothing.
The deployed question is "a new person puts on the device, it has never seen
her heart, can it classify her beats?". That is inter-patient by construction.

DS2 IS OPENED ONCE.
Any tuning happens inside DS1 by cross-validation across its patients. The
moment you look at DS2 performance and adjust anything, DS2 stops being a held
out test set in spirit even though the records never mixed.

CLASS Q IS DROPPED FROM TRAINING.
Only 15 Q beats survive the paced-record exclusion (8 in DS1, 7 in DS2). With
class_weight='balanced' each would carry a weight near 1250, so the model would
contort around statistical noise. The published comparison tables report
N, S, V, F for the same reason. The 15 beats are reported as excluded, not
quietly dropped.

TWO IMBALANCE STRATEGIES
  headline : class_weight='balanced'. Reweights the split criterion. Invents no
             data, so nothing has to be defended beyond the weighting itself.
  ablation : SMOTE oversampling of the training set only. Note the specific
             criticism in an inter-patient setting: SMOTE interpolates between
             a minority beat and its nearest neighbours, which often belong to
             DIFFERENT PATIENTS, producing synthetic beats that correspond to
             no real heart. Since patient-to-patient variation is exactly what
             makes this task hard, this may inject noise rather than signal.
             If SMOTE scores worse, that is a finding, not a failure.
"""
import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (classification_report, confusion_matrix,
                             precision_recall_fscore_support)

DATASET = 'mitdb'
RESULTS_DIR = 'results'
FEATURES = os.path.join(RESULTS_DIR, f'{DATASET}_e2_features.npz')

# de Chazal et al. 2004 inter-patient split. 22 records each, 44 total.
# The four paced records (102, 104, 107, 217) are already absent from both.
DS1 = ['101', '106', '108', '109', '112', '114', '115', '116', '118', '119',
       '122', '124', '201', '203', '205', '207', '208', '209', '215', '220',
       '223', '230']
DS2 = ['100', '103', '105', '111', '113', '117', '121', '123', '200', '202',
       '210', '212', '213', '214', '219', '221', '222', '228', '231', '232',
       '233', '234']

CLASSES = ['N', 'S', 'V', 'F']      # Q excluded, see module docstring
DROP_CLASSES = ['Q']

N_TREES = 100
RANDOM_STATE = 42                   # fixed so the run is reproducible


def load_features():
    d = np.load(FEATURES, allow_pickle=True)
    X = np.hstack([d['X_morph'], d['rr']])
    y = d['y'].astype(str)
    rec = d['record'].astype(str)

    # Names come FROM the feature file. They used to be hardcoded here, which
    # meant adding a feature in 05 would silently misalign every name in the
    # importance report without raising anything.
    if 'feat_names' in d:
        names = [str(n) for n in d['feat_names']]
    else:
        names = ([f'morph_{i:02d}' for i in range(d['X_morph'].shape[1])] +
                 [f'rr_{i:02d}' for i in range(d['rr'].shape[1])])
    assert len(names) == X.shape[1], (
        f"{len(names)} names for {X.shape[1]} feature columns: re-run 05")
    return X, y, rec, names


def check_split(rec):
    """Fail loudly if the split is malformed. A silent split bug would
    invalidate every number downstream, so this is worth the ten lines."""
    present = set(rec)
    overlap = set(DS1) & set(DS2)
    missing = present - set(DS1) - set(DS2)
    unused = (set(DS1) | set(DS2)) - present

    assert not overlap, f"records in BOTH DS1 and DS2: {sorted(overlap)}"
    assert len(DS1) == 22 and len(DS2) == 22, "DS1/DS2 should be 22 records each"
    if missing:
        print(f"  WARNING: records in features but in neither set: {sorted(missing)}")
    if unused:
        print(f"  WARNING: records in the split but not in features: {sorted(unused)}")
    return True


def report(y_true, y_pred, label):
    """Per-class precision/recall/F1 plus confusion matrix.

    No headline accuracy. N is ~90% of beats, so a model that always guesses N
    scores ~90% accuracy while catching zero arrhythmias. Same reason E1 uses
    sensitivity and PPV rather than accuracy: both deliberately ignore the
    abundant easy case and measure performance on the rare hard one.
    """
    print(f"\n{'='*62}")
    print(f"  {label}")
    print(f"{'='*62}")

    p, r, f, sup = precision_recall_fscore_support(
        y_true, y_pred, labels=CLASSES, zero_division=0)

    print(f"\n  {'class':>6}  {'precision':>10}  {'recall':>8}  "
          f"{'F1':>7}  {'support':>8}")
    for i, c in enumerate(CLASSES):
        print(f"  {c:>6}  {100*p[i]:9.2f}%  {100*r[i]:7.2f}%  "
              f"{100*f[i]:6.2f}%  {sup[i]:8,}")

    macro_f1 = float(np.mean(f))
    print(f"\n  macro F1 (unweighted mean across classes) : {100*macro_f1:.2f}%")
    print("  Macro, not weighted: weighting by support would let N drown out")
    print("  the classes that actually matter clinically.")

    cm = confusion_matrix(y_true, y_pred, labels=CLASSES)
    print(f"\n  Confusion matrix (rows = true, cols = predicted):")
    print(f"      {'':>6}" + ''.join(f"{c:>9}" for c in CLASSES))
    for i, c in enumerate(CLASSES):
        row = ''.join(f"{cm[i, j]:>9,}" for j in range(len(CLASSES)))
        print(f"      {c:>6}{row}")

    print("\n  Read the S row: where do the missed S beats GO? If they land in")
    print("  the N column, the model is calling early beats normal, which is")
    print("  the expected failure. S beats look near-normal and differ mainly")
    print("  by TIMING, which varies from person to person.")
    return {'precision': p, 'recall': r, 'f1': f, 'support': sup,
            'macro_f1': macro_f1, 'cm': cm}


def train_eval(Xtr, ytr, Xte, yte, label, use_smote=False):
    if use_smote:
        try:
            from imblearn.over_sampling import SMOTE
        except ImportError:
            print(f"\n  SKIPPING {label}: imbalanced-learn not installed.")
            print("  Install with:  python -m pip install imbalanced-learn")
            return None
        # SMOTE on the TRAINING SET ONLY. Resampling before the split, or
        # touching the test set at all, would leak synthetic neighbours of test
        # beats into training and inflate every number.
        n_min = min(np.bincount(pd.factorize(ytr)[0]))
        k = max(1, min(5, n_min - 1))
        Xtr, ytr = SMOTE(random_state=RANDOM_STATE, k_neighbors=k).fit_resample(Xtr, ytr)
        print(f"\n  SMOTE resampled training set to {len(ytr):,} beats "
              f"(k_neighbors={k})")
        for c in CLASSES:
            print(f"    {c}: {int((ytr == c).sum()):,}")

    clf = RandomForestClassifier(
        n_estimators=N_TREES,
        class_weight=None if use_smote else 'balanced',
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )
    # No feature scaling: Random Forest splits on thresholds per feature, so it
    # is scale invariant. Morphology is already z-normalised per record and RR
    # features are in seconds. Mixing units is fine here.
    clf.fit(Xtr, ytr)
    pred = clf.predict(Xte)
    res = report(yte, pred, label)
    res['clf'] = clf
    return res


def main():
    if not os.path.exists(FEATURES):
        print(f"Missing {FEATURES}. Run 05_extract_features.py first.")
        return

    X, y, rec, names = load_features()
    print(f"Loaded {X.shape[0]:,} beats x {X.shape[1]} features")

    check_split(rec)

    n_dropped = int(np.isin(y, DROP_CLASSES).sum())
    mask = ~np.isin(y, DROP_CLASSES)
    X, y, rec = X[mask], y[mask], rec[mask]
    print(f"Dropped class {DROP_CLASSES}: {n_dropped} beats "
          f"(too few to learn or evaluate; see docstring)")

    tr = np.isin(rec, DS1)
    te = np.isin(rec, DS2)
    Xtr, ytr = X[tr], y[tr]
    Xte, yte = X[te], y[te]

    print(f"\n=== Inter-patient split (de Chazal et al., 2004) ===")
    print(f"  DS1 train : {len(DS1)} records, {len(ytr):,} beats")
    print(f"  DS2 test  : {len(DS2)} records, {len(yte):,} beats")
    print(f"  patient overlap: none by construction\n")

    print(f"  {'class':>6}  {'DS1 train':>11}  {'DS2 test':>10}")
    for c in CLASSES:
        ntr, nte = int((ytr == c).sum()), int((yte == c).sum())
        print(f"  {c:>6}  {ntr:>11,}  {nte:>10,}")

    n_s_tr, n_s_te = int((ytr == 'S').sum()), int((yte == 'S').sum())
    if n_s_te > n_s_tr:
        print(f"\n  Note: there are MORE S beats in test ({n_s_te:,}) than in")
        print(f"  train ({n_s_tr:,}). This is a known property of the published")
        print("  split, not an error. Record 232 alone holds about half of all")
        print("  S beats in the database and sits in DS2. It is also the")
        print("  realistic case: you rarely have many examples of a rare")
        print("  arrhythmia at training time, but a deployed device meets them.")

    head = train_eval(Xtr, ytr, Xte, yte,
                      "HEADLINE — RF + class_weight='balanced'")

    if head:
        imp = pd.Series(head['clf'].feature_importances_, index=names)
        print(f"\n  Top 8 features by importance:")
        for n, v in imp.nlargest(8).items():
            print(f"    {n:<16} {v:.4f}")
        rr_share = imp[[n for n in names if n.startswith('rr_')]].sum()
        n_rr = len([n for n in names if n.startswith('rr_')])
        print(f"\n  RR features account for {100*rr_share:.1f}% of total importance")
        print(f"  ({n_rr} of {len(names)} features). Published work finds normalised RR")
        print("  intervals among the most discriminative features, so this is a")
        print("  useful sanity check on the pipeline rather than a new result.")

    abl = train_eval(Xtr, ytr, Xte, yte,
                     "ABLATION — RF + SMOTE oversampling", use_smote=True)

    if head and abl:
        print(f"\n{'='*62}")
        print("  COMPARISON")
        print(f"{'='*62}")
        print(f"\n  {'class':>6}  {'balanced F1':>12}  {'SMOTE F1':>10}  {'delta':>8}")
        for i, c in enumerate(CLASSES):
            d = 100 * (abl['f1'][i] - head['f1'][i])
            print(f"  {c:>6}  {100*head['f1'][i]:11.2f}%  "
                  f"{100*abl['f1'][i]:9.2f}%  {d:+7.2f}")
        d = 100 * (abl['macro_f1'] - head['macro_f1'])
        print(f"  {'macro':>6}  {100*head['macro_f1']:11.2f}%  "
              f"{100*abl['macro_f1']:9.2f}%  {d:+7.2f}")

        if abl['macro_f1'] < head['macro_f1']:
            print("\n  SMOTE performed WORSE. Report this, do not bury it. The")
            print("  likely reason is specific and worth stating: SMOTE builds")
            print("  synthetic beats by interpolating between a minority beat")
            print("  and its nearest neighbours, which here often belong to")
            print("  DIFFERENT PATIENTS. The result is a beat matching no real")
            print("  heart, in a task whose entire difficulty is patient to")
            print("  patient variation.")
        else:
            print("\n  SMOTE improved macro F1. Still report the caveat above:")
            print("  synthetic between-patient beats are a real methodological")
            print("  concern even when the number goes up.")

    print(f"\n{'='*62}")
    print("  This is ONE split, so there are no error bars on any of these")
    print("  numbers. The published split is fixed by design, which is what")
    print("  makes results comparable, but it means a single unusual test")
    print("  patient moves the metrics. Say so rather than implying precision")
    print("  the evaluation does not have.")


if __name__ == '__main__':
    main()
