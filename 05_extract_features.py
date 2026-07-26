"""
05_extract_features.py — E2: turn annotated beats into a feature table.

BEAT POSITIONS COME FROM THE .atr ANNOTATIONS, NOT FROM OUR DETECTOR.
E1 asks "where are the beats?". E2 asks "given a beat, what type is it?".
Feeding E2 our own detections would fold E1's errors into E2's score and make
the two sub-problems impossible to tell apart. In deployment you chain E1 -> E2;
for evaluation you isolate them. This is standard practice and worth stating.

FILTERING IS DIFFERENT FROM E1, DELIBERATELY.
E1's 5-15 Hz bandpass is a DETECTION filter: it destroys P and T waves on
purpose, because they are obstacles to finding the QRS. E2 needs those waves,
because their shape carries the diagnosis. So morphology uses a gentle
0.5-40 Hz bandpass: baseline wander out, shape preserved.

RR FEATURES, AND WHY THERE ARE NOW NINE OF THEM
The first version carried five RR features, all normalised against a rolling
10-beat local mean. 07_diagnose_s_class.py showed why that was not enough.

  Record 232 holds 1,381 of the test set's 1,836 S beats and scored 0.07%
  recall. Yet the local RR ratio separated S from N in that record with an AUC
  of 0.992, near perfect. AUC is RANK based and threshold free; a Random Forest
  learns ABSOLUTE thresholds. So the feature ordered 232's beats correctly
  while sitting in completely the wrong place on the number line.

  The cause: record 232 contains 397 N beats and 1,381 S beats. The ectopic
  class is the MAJORITY for that patient. Any baseline computed from the
  patient's own rhythm is therefore contaminated by the arrhythmia itself.
    DS1 records (normal balance): baseline ~ the normal N-N interval, so an S
      beat lands near rr_pre_ratio 0.6 and the model learns "well below 1 = S".
    Record 232: the baseline is dominated by the premature beats, so its S
      beats land near 1.0, exactly where DS1's NORMAL beats sit.
  The model does what it was taught and calls them normal. 74.3% of 232's S
  beats were classified N, which is precisely this.

  So the additions target baseline contamination, not imbalance:
    record-level ratios : median RR over the whole record is a steadier
                          baseline than a 10-beat window. Justified
                          empirically: improves S-vs-N separability in 8 of 15
                          records, notably training records 207 (+0.431) and
                          209 (+0.197). Better separability in DS1 means a
                          better learned boundary.
    rr_post_pre_ratio   : THE IMPORTANT ONE. Uses no baseline at all, so
                          contamination cannot reach it. An S beat is
                          short-then-long (premature, then compensatory pause)
                          giving a ratio well above 1; a normal beat gives
                          about 1. That signature holds whether or not the
                          patient's dominant rhythm is arrhythmic.
    p80 baseline        : the 80th percentile of RR is a robust proxy for the
                          patient's non-premature rhythm, and unlike the median
                          it survives ectopics being the majority.

  HONEST CAVEAT: record-level statistics use the whole recording, including
  beats after the one being classified. That is fine for offline analysis and
  is what the published feature sets do, but a real-time device would need a
  trailing estimate instead. Say so rather than implying the pipeline is
  already streaming-ready.

  HONEST EXPECTATION: this is unlikely to rescue record 232, because the median
  and the local window are contaminated by the same mechanism. If 232 stays
  near zero, maximum achievable S recall is 455/1836 = 24.8%. The gain to look
  for is on the OTHER 455 S beats, which are currently also being missed.

RECORD ID IS RETAINED FOR EVERY BEAT.
Without it an inter-patient split is impossible: you cannot hold out a patient
you cannot identify.
"""
import os
import collections

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

import wfdb

DATASET = 'mitdb'
DATA_DIR = os.path.join('data', DATASET)
RESULTS_DIR = 'results'

# --- AAMI EC57 five-class mapping -------------------------------------------
# The 23 MIT-BIH beat symbols collapse to 5 clinically meaningful classes.
# Using the standard mapping (rather than inventing one) is what makes these
# results comparable to the published literature.
AAMI = {
    # N: normal + bundle branch block + escape. NOTE these are NOT all "healthy"
    # beats: L (left BBB) is a conduction defect that looks nothing like a
    # normal beat, yet carries the same class label. This is a large part of why
    # inter-patient generalisation is hard.
    'N': 'N', 'L': 'N', 'R': 'N', 'e': 'N', 'j': 'N',
    # S: supraventricular ectopic. The hard class: near-normal shape, early timing.
    'A': 'S', 'a': 'S', 'J': 'S', 'S': 'S',
    # V: ventricular ectopic. Wide, bizarre, usually easy.
    'V': 'V', 'E': 'V',
    # F: fusion of ventricular and normal. Genuinely ambiguous by nature.
    'F': 'F',
    # Q: paced / unclassifiable.
    '/': 'Q', 'f': 'Q', 'Q': 'Q',
}
CLASSES = ['N', 'S', 'V', 'F', 'Q']

# Paced records, excluded per AAMI convention: pacing is not arrhythmia, and a
# paced beat's shape is set by the device rather than the heart.
# NOTE this is a DIFFERENT decision from E1's lead-mismatch exclusion of
# 102/104. It happens that {102,104} is a subset of the paced set, so E2's 44
# records all carry MLII anyway. Two justifications, one convenient overlap.
PACED = ['102', '104', '107', '217']

# --- signal / feature constants ---
BP_LOW_HZ = 0.5              # kill baseline wander
BP_HIGH_HZ = 40.0            # kill muscle/mains noise, keep P/QRS/T shape
BP_ORDER = 3
WIN_PRE_S = 0.25             # before R: captures the P wave
WIN_POST_S = 0.45            # after R: captures the T wave
N_MORPH = 32                 # resampled morphology points
RR_LOCAL_N = 10              # beats in the local RR average
RR_BASELINE_PCTL = 80        # robust proxy for the non-premature rhythm

# Column order is FIXED. 07 indexes rr[:, 0], [:, 2], [:, 3] positionally, so
# new features are APPENDED, never inserted.
RR_NAMES = [
    'rr_pre_s', 'rr_post_s', 'rr_local_s',
    'rr_pre_ratio', 'rr_post_ratio',
    'rr_pre_rec_ratio', 'rr_post_rec_ratio',
    'rr_post_pre_ratio', 'rr_pre_p80_ratio',
]


def clean_nonfinite(x):
    """Interpolate over NaN/inf before filtering.

    One non-finite sample makes scipy's filtfilt return all-NaN and silently
    poison everything downstream.
    """
    x = np.asarray(x, dtype=float)
    bad = ~np.isfinite(x)
    if bad.any():
        idx = np.arange(x.size)
        x = x.copy()
        x[bad] = np.interp(idx[bad], idx[~bad], x[~bad])
    return x


def find_channel(rec, *targets):
    """Index of first matching lead by name; falls back to 0. See 02."""
    cleaned = [n.strip() for n in rec.sig_name]
    for t in targets:
        if t in cleaned:
            return cleaned.index(t)
    return 0


def preprocess(sig, fs):
    """0.5-40 Hz bandpass, then per-record z-normalisation.

    Per-RECORD normalisation, not per-beat and not global. Absolute mV amplitude
    depends on electrode placement and body habitus, so it differs between
    patients for reasons that have nothing to do with arrhythmia. Left raw, the
    model could learn "record 205 is loud" instead of "this beat is wide" and
    that lesson would not transfer to a new patient. Normalising per record
    removes the patient-specific scale while preserving relative shape and the
    amplitude differences BETWEEN beats within a patient, which are real signal.
    """
    x = clean_nonfinite(sig)
    nyq = 0.5 * fs
    high = min(BP_HIGH_HZ, 0.95 * nyq)
    sos = butter(BP_ORDER, [BP_LOW_HZ / nyq, high / nyq], btype='band', output='sos')
    x = sosfiltfilt(sos, x)
    sd = np.std(x)
    return (x - np.mean(x)) / sd if sd > 0 else x - np.mean(x)


def morphology(sig, r, fs):
    """Fixed-length beat window resampled to N_MORPH points.

    Resampled rather than raw-sliced so the feature vector has the same length
    regardless of sampling rate. Keeps the door open for the ESP32 data later,
    which will not be 360 Hz.
    """
    a = r - int(round(WIN_PRE_S * fs))
    b = r + int(round(WIN_POST_S * fs))
    if a < 0 or b >= sig.size:
        return None
    seg = sig[a:b]
    src = np.linspace(0, 1, seg.size)
    dst = np.linspace(0, 1, N_MORPH)
    return np.interp(dst, src, seg)


def extract_record(name):
    path = os.path.join(DATA_DIR, name)
    rec = wfdb.rdrecord(path)
    fs = rec.fs
    ch = find_channel(rec, 'MLII', 'II')
    sig = preprocess(rec.p_signal[:, ch], fs)

    ann = wfdb.rdann(path, 'atr')
    keep = np.array([s in AAMI for s in ann.symbol])
    samples = ann.sample[keep]
    labels = [AAMI[s] for s in np.array(ann.symbol)[keep]]

    # --- record-level RR baselines, computed once ---
    # Median: steadier than a rolling window, but note it IS contaminated when
    # ectopic beats are the majority (record 232). p80: more robust in exactly
    # that case, because premature beats sit in the lower tail.
    all_rr = np.diff(samples) / fs if len(samples) > 1 else np.array([1.0])
    rr_rec = float(np.median(all_rr))
    rr_p80 = float(np.percentile(all_rr, RR_BASELINE_PCTL))
    if rr_rec <= 0:
        rr_rec = 1.0
    if rr_p80 <= 0:
        rr_p80 = rr_rec

    rows, X = [], []
    for i in range(1, len(samples) - 1):        # need both neighbours for RR
        r = int(samples[i])
        m = morphology(sig, r, fs)
        if m is None:                            # window would overrun the record
            continue

        rr_pre = (samples[i] - samples[i - 1]) / fs
        rr_post = (samples[i + 1] - samples[i]) / fs
        lo = max(0, i - RR_LOCAL_N)
        local = np.diff(samples[lo:i + 1]) / fs
        rr_local = float(np.mean(local)) if local.size else rr_pre

        rows.append({
            'record': name,
            'sample': r,
            'label': labels[i],
            # original five
            'rr_pre_s': rr_pre,
            'rr_post_s': rr_post,
            'rr_local_s': rr_local,
            'rr_pre_ratio': rr_pre / rr_local if rr_local > 0 else 1.0,
            'rr_post_ratio': rr_post / rr_local if rr_local > 0 else 1.0,
            # record-level baseline: steadier than the rolling window
            'rr_pre_rec_ratio': rr_pre / rr_rec,
            'rr_post_rec_ratio': rr_post / rr_rec,
            # BASELINE-FREE. Immune to the contamination that breaks the others.
            # S beat: short pre, long compensatory post -> well above 1.
            # N beat: roughly equal -> about 1.
            'rr_post_pre_ratio': rr_post / rr_pre if rr_pre > 0 else 1.0,
            # robust baseline, survives ectopics being the majority rhythm
            'rr_pre_p80_ratio': rr_pre / rr_p80,
        })
        X.append(m)

    return rows, np.array(X, dtype=np.float32)


def main():
    all_recs = sorted(f[:-4] for f in os.listdir(DATA_DIR) if f.endswith('.hea'))
    recs = [r for r in all_recs if r not in PACED]

    print(f"Records total   : {len(all_recs)}")
    print(f"Paced (excluded): {PACED}")
    print(f"E2 records      : {len(recs)}\n")

    all_rows, all_X = [], []
    for i, name in enumerate(recs, 1):
        rows, X = extract_record(name)
        all_rows.extend(rows)
        all_X.append(X)
        c = collections.Counter(r['label'] for r in rows)
        print(f"  [{i:2d}/{len(recs)}] {name}: {len(rows):5d} beats  "
              f"{ {k: c.get(k, 0) for k in CLASSES} }")

    meta = pd.DataFrame(all_rows)
    X = np.vstack(all_X)
    assert len(meta) == X.shape[0], "metadata and feature rows out of sync"

    morph_names = [f'morph_{i:02d}' for i in range(X.shape[1])]
    feat_names = morph_names + RR_NAMES

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out = os.path.join(RESULTS_DIR, f'{DATASET}_e2_features.npz')
    np.savez_compressed(
        out,
        X_morph=X,
        rr=meta[RR_NAMES].to_numpy(np.float32),
        y=meta['label'].to_numpy(),
        record=meta['record'].to_numpy(),
        sample=meta['sample'].to_numpy(),
        # Names travel WITH the data. Previously 06 hardcoded five RR names, so
        # adding a feature here would have silently desynced them.
        feat_names=np.array(feat_names),
    )

    counts = collections.Counter(meta['label'])
    total = len(meta)
    print(f"\n=== Class distribution ({len(recs)} records, {total:,} beats) ===")
    for c in CLASSES:
        n = counts.get(c, 0)
        print(f"  {c} : {n:>7,}  ({100*n/total:5.2f}%)")

    n_maj = counts.most_common(1)[0][1]
    print(f"\n  Majority class is {100*n_maj/total:.1f}% of the data. A model that")
    print("  always guesses it scores that as ACCURACY while catching zero")
    print("  arrhythmias. Report per-class precision/recall and a confusion")
    print("  matrix instead. Same reason E1 excludes true negatives.")

    # Sanity check on the baseline-free feature, using record 232 as the case
    # the whole change was designed around.
    print(f"\n=== rr_post_pre_ratio: does the baseline-free feature separate? ===")
    print("  (S beats should sit clearly above N beats. This feature uses no")
    print("   baseline, so arrhythmia-dominated records cannot contaminate it.)")
    print(f"\n  {'record':>7}  {'N median':>9}  {'S median':>9}  {'separation':>11}")
    for r in ['100', '209', '222', '232']:
        sub = meta[meta['record'] == r]
        n_med = sub.loc[sub['label'] == 'N', 'rr_post_pre_ratio'].median()
        s_med = sub.loc[sub['label'] == 'S', 'rr_post_pre_ratio'].median()
        if pd.notna(n_med) and pd.notna(s_med):
            print(f"  {r:>7}  {n_med:9.3f}  {s_med:9.3f}  {s_med-n_med:+11.3f}")

    print(f"\n  features -> {out}")
    print(f"  shape: {X.shape[0]:,} beats x {len(feat_names)} features "
          f"({X.shape[1]} morphology + {len(RR_NAMES)} RR)")
    print("\n  Add results/*.npz to .gitignore: this file is large and")
    print("  regenerable from the raw data, like the PhysioNet downloads.")


if __name__ == '__main__':
    main()
