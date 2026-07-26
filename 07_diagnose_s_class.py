"""
07_diagnose_s_class.py — why is S recall 7.41%?

===========================================================================
READ THIS FIRST — PROVENANCE (added after the feature fix)
This script was written to diagnose the ORIGINAL 5-RR-feature model, whose
headline S recall was 7.41%. Its H2 section proposed adding record-level RR
ratios and a baseline-free ratio to 05. That fix was IMPLEMENTED: 05 now emits
9 RR features and 06 was re-run against them.

Consequence for anyone reading the output:
  - H1 (per-record S recall breakdown) is STILL CURRENT and is the finding to
    cite. It shows S recall is dominated by record 232 (~75% of DS2 S beats,
    ~0% recall) while ordinary records score 30-100%. That is the real result.
  - H2 (local-window vs record-level AUC) is NOW HISTORICAL. It analyses whether
    the fix WAS worth doing, and its verdict text still says "add record-level
    RR ratios to 05 and re-run" even though that has already happened. Do NOT
    quote H2's verdict as a live recommendation. The features it recommends are
    already in the model being evaluated.
  - The headline numbers in this docstring (7.41%, 11.67%) are the PRE-FIX
    values. Post-fix the headline S recall is 7.19% / F1 12.14%: the fix
    improved separability and precision but could not move 232's recall, which
    is structurally capped. See 05's docstring for why.

Kept as-is rather than rewritten: it is the record of the diagnostic reasoning,
and re-running it reproduces the H1 table that belongs in the writeup.
===========================================================================

E2's headline (PRE-FIX): N F1 95.73%, V F1 72.37%, S F1 11.67%. Published RF
work on this same inter-patient split reports S F1 around 73%. So S is not
merely "hard here", it is underperforming what the split allows. Before changing
any features, find out where the failure actually lives.

WHAT THE SMOTE ABLATION ALREADY TOLD US
  SMOTE quadrupled the training data for S and V equally.
    V: F1 72.37% -> 81.31%  (+8.94)
    S: F1 11.67% -> 12.01%  (+0.34)
  If S were separable but merely outnumbered, rebalancing would have moved it
  the way it moved V. It did not. So the problem is not imbalance. S beats are
  not DISTINGUISHABLE from N with the current features.

TWO HYPOTHESES
  H1  Record 232 dominates. It holds 1,381 of DS2's 1,836 S beats (75%), so
      overall S recall is close to 232's S recall. If 232 is pathological and
      other records are respectable, the headline number is one patient.

  H2  rr_local is self-defeating during sustained arrhythmia. The current
      feature is rr_pre / mean(last 10 RR). That isolates an early beat only
      when the surrounding rhythm is REGULAR. If a patient is in near
      continuous atrial arrhythmia, every RR is irregular, so rr_local tracks
      the irregularity itself and the ratio drifts back toward 1. The feature
      erases the thing it was built to detect.

  Supporting hint: rr_pre_s (absolute, 0.0954 importance) ranks nearly as high
  as rr_pre_ratio (0.1013). The model is leaning on raw interval length, which
  is patient specific and transfers poorly, which is what you would expect if
  the ratio had gone flat.

THE TEST FOR H2
  Compare how well two normalisations separate S from N WITHIN each record:
    local  : rr_pre / mean(last 10 RR)         <- current
    record : rr_pre / median(all RR in record) <- proposed
  Separation is measured as AUC, which needs no threshold: 0.5 means the
  feature carries no information, 1.0 means perfect separation. If the record
  level baseline gives a materially higher AUC, the fix is justified. If both
  are near 0.5, the timing information is simply absent and the answer lies in
  morphology instead, so we would think again rather than add features on faith.

  Note this validates the fix BEFORE implementing it, using data already in the
  .npz. No re-extraction needed to find out whether it is worth doing.
"""
import os
import importlib.util

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

TRAIN_FILE = '06_train_classifier.py'
DATASET = 'mitdb'
RESULTS_DIR = 'results'
FEATURES = os.path.join(RESULTS_DIR, f'{DATASET}_e2_features.npz')

MIN_S_FOR_AUC = 15      # below this, an AUC is noise rather than evidence


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    if not os.path.exists(FEATURES):
        print(f"Missing {FEATURES}. Run 05_extract_features.py first.")
        return

    # Split, classes and model config are imported from 06 rather than restated,
    # so this diagnostic cannot silently drift from the thing it is diagnosing.
    tr_mod = load_module(TRAIN_FILE, 'train_classifier')
    DS1, DS2 = tr_mod.DS1, tr_mod.DS2
    CLASSES = tr_mod.CLASSES

    d = np.load(FEATURES, allow_pickle=True)
    X = np.hstack([d['X_morph'], d['rr']])
    y = d['y'].astype(str)
    rec = d['record'].astype(str)

    # rr columns, in the order 05 wrote them
    rr_pre = d['rr'][:, 0]
    rr_local = d['rr'][:, 2]
    rr_pre_ratio = d['rr'][:, 3]

    keep = np.isin(y, CLASSES)
    X, y, rec = X[keep], y[keep], rec[keep]
    rr_pre, rr_local, rr_pre_ratio = rr_pre[keep], rr_local[keep], rr_pre_ratio[keep]

    tr, te = np.isin(rec, DS1), np.isin(rec, DS2)

    # --- reproduce the headline model exactly ---
    clf = RandomForestClassifier(
        n_estimators=tr_mod.N_TREES, class_weight='balanced',
        random_state=tr_mod.RANDOM_STATE, n_jobs=-1)
    clf.fit(X[tr], y[tr])
    pred = clf.predict(X[te])

    yte, rte = y[te], rec[te]

    # ---------------------------------------------------------------- H1
    print("=" * 68)
    print("  H1: is record 232 the whole story?")
    print("=" * 68)
    print(f"\n  {'record':>7}  {'S beats':>8}  {'S recall':>9}  "
          f"{'-> N':>7}  {'-> V':>6}  {'share of DS2 S':>15}")

    total_s = int((yte == 'S').sum())
    rows = []
    for r in DS2:
        m = (rte == r)
        s_mask = m & (yte == 'S')
        n_s = int(s_mask.sum())
        if n_s == 0:
            continue
        p = pred[s_mask]
        rc = float((p == 'S').mean())
        to_n = float((p == 'N').mean())
        to_v = float((p == 'V').mean())
        rows.append({'record': r, 'n_s': n_s, 'recall': rc,
                     'to_n': to_n, 'to_v': to_v})
        print(f"  {r:>7}  {n_s:>8,}  {100*rc:8.2f}%  {100*to_n:6.1f}%  "
              f"{100*to_v:5.1f}%  {100*n_s/total_s:14.1f}%")

    df = pd.DataFrame(rows)
    big = df[df['n_s'] >= 50]
    if len(big):
        print(f"\n  Records with >=50 S beats: {len(big)}")
        print(f"  Their mean S recall      : {100*big['recall'].mean():.2f}%")
        print(f"  Beat-weighted S recall   : "
              f"{100*(big['recall']*big['n_s']).sum()/big['n_s'].sum():.2f}%")
        worst = big.nsmallest(1, 'recall').iloc[0]
        best = big.nlargest(1, 'recall').iloc[0]
        print(f"  Worst: {worst['record']} at {100*worst['recall']:.1f}% "
              f"({int(worst['n_s'])} S beats)")
        print(f"  Best : {best['record']} at {100*best['recall']:.1f}% "
              f"({int(best['n_s'])} S beats)")
        spread = big['recall'].max() - big['recall'].min()
        if spread > 0.3:
            print(f"\n  Spread across records is {100*spread:.0f} points. S recall")
            print("  is dominated by WHICH patients are in the test set, not by")
            print("  a uniform model weakness. H1 has support.")
        else:
            print(f"\n  Spread is only {100*spread:.0f} points: the failure is")
            print("  fairly uniform across patients, so this is not one record.")

    # ---------------------------------------------------------------- H2
    print("\n" + "=" * 68)
    print("  H2: does the local RR window erase the earliness signal?")
    print("=" * 68)
    print("\n  AUC of each normalisation as an S-vs-N discriminator, per record.")
    print("  0.50 = no information. Higher = better separation.")
    print("  Both are computed on the SAME beats, so they are directly")
    print("  comparable; only the denominator differs.\n")
    print(f"  {'record':>7}  {'S':>6}  {'N':>7}  {'local AUC':>10}  "
          f"{'record AUC':>11}  {'delta':>7}")

    aucs = []
    for r in sorted(set(rec)):
        m = (rec == r)
        s_mask = m & (y == 'S')
        n_mask = m & (y == 'N')
        n_s, n_n = int(s_mask.sum()), int(n_mask.sum())
        if n_s < MIN_S_FOR_AUC or n_n < MIN_S_FOR_AUC:
            continue

        sub = s_mask | n_mask
        target = (y[sub] == 'S').astype(int)

        # Record-level baseline: median RR over the whole record. Median rather
        # than mean because ectopic beats and pauses skew a mean.
        rec_rr = float(np.median(rr_pre[m]))
        rec_ratio = rr_pre[sub] / rec_rr if rec_rr > 0 else rr_pre[sub]

        # Lower ratio means earlier means more likely S, so negate for AUC.
        try:
            a_local = roc_auc_score(target, -rr_pre_ratio[sub])
            a_rec = roc_auc_score(target, -rec_ratio)
        except ValueError:
            continue

        split = 'DS1' if r in DS1 else 'DS2'
        aucs.append({'record': r, 'split': split, 'n_s': n_s,
                     'local': a_local, 'rec': a_rec})
        print(f"  {r:>7}  {n_s:>6,}  {n_n:>7,}  {a_local:10.3f}  "
              f"{a_rec:11.3f}  {a_rec-a_local:+7.3f}")

    ad = pd.DataFrame(aucs)
    if len(ad):
        print(f"\n  Records analysed        : {len(ad)}")
        print(f"  Mean local-window AUC   : {ad['local'].mean():.3f}")
        print(f"  Mean record-level AUC   : {ad['rec'].mean():.3f}")
        print(f"  Mean improvement        : {ad['rec'].mean()-ad['local'].mean():+.3f}")
        better = int((ad['rec'] > ad['local']).sum())
        print(f"  Record-level wins in    : {better}/{len(ad)} records")

        # Weight by S beats: improving a record with 1,381 S beats matters far
        # more to the headline metric than improving one with 20.
        w = ad['n_s']
        wl = float((ad['local'] * w).sum() / w.sum())
        wr = float((ad['rec'] * w).sum() / w.sum())
        print(f"\n  S-beat-weighted local AUC  : {wl:.3f}")
        print(f"  S-beat-weighted record AUC : {wr:.3f}")
        print(f"  weighted improvement       : {wr-wl:+.3f}")

    # ---------------------------------------------------------------- verdict
    print("\n" + "=" * 68)
    print("  VERDICT")
    print("=" * 68)
    if len(ad):
        gain = wr - wl
        if gain >= 0.03:
            print(f"\n  Record-level normalisation separates S from N better")
            print(f"  ({wl:.3f} -> {wr:.3f} weighted AUC). H2 had support and the")
            print("  fix was ALREADY APPLIED: 05 now emits these record-level RR")
            print("  ratios and 06 was re-run against them. This section is kept as")
            print("  the historical record of that decision (see the banner at the")
            print("  top of this file). It improved separability but could not lift")
            print("  record 232's recall, which is structurally capped.")
        elif gain <= -0.01:
            print(f"\n  Record-level normalisation is WORSE ({wl:.3f} -> {wr:.3f}).")
            print("  H2 is wrong. Do not add the feature. The local window is")
            print("  doing its job and the problem lies elsewhere.")
        else:
            print(f"\n  Little difference ({wl:.3f} -> {wr:.3f}). H2 is not")
            print("  supported. If both AUCs sit near 0.5, timing information is")
            print("  largely ABSENT for these S beats, and no RR feature will")
            print("  rescue them. That would point at morphology, or at an")
            print("  honest report that S is not separable with this feature set.")

        near_chance = ad[ad['local'] < 0.6]
        if len(near_chance):
            print(f"\n  {len(near_chance)} record(s) have local AUC below 0.60,")
            print("  meaning timing barely distinguishes their S beats at all:")
            print(f"    {', '.join(near_chance['record'].tolist())}")

    print("\n  Whatever this shows, it is now a measured finding rather than a")
    print("  guess, and it belongs in the writeup either way. A diagnosed")
    print("  limitation reads as understanding; an undiagnosed one reads as a gap.")


if __name__ == '__main__':
    main()
