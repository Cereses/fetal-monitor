# Phase 1 MHR — Summary

End-of-phase notes for the maternal heart rate (MHR) detection pipeline.
Companion to `phase1_fhr_summary.md`.

## Framing: why this is easier than FHR

FHR is a separation problem (weak fetal heart *sounds* pulled out of a noisy
abdominal recording). MHR is a direct measurement: the MAX30102 sits on the
mother and her pulse is the dominant, clean signal at the sensor. There is
nothing to separate — the fetal signal does not meaningfully appear in a
maternal finger/wrist PPG. So MHR is standard adult PPG heart-rate estimation,
and the pipeline is correspondingly simpler (peak-finding, not autocorrelation).

Pregnancy does not change the signal-processing problem; it only shifts the
*normal range* used for flagging (maternal resting HR rises over gestation).
Those thresholds come from clinical literature, not from the validation data —
which is why a non-pregnant validation dataset is acceptable here.

## Pipeline

Raw PPG → bandpass 0.5-5 Hz (SOS form) → auto-orient peaks up →
per-window: z-normalize → find_peaks → de-double check → inter-beat
intervals → BPM + confidence.

Key design choices and why:
- **SOS filter (sosfiltfilt), not b/a + filtfilt.** A 0.5 Hz cutoff at 125 Hz
  is a very low normalized frequency (~0.008) where the transfer-function form
  is numerically fragile. Second-order-sections form is stable there.
- **Per-window z-normalization before peak-finding.** Makes the prominence
  threshold amplitude-agnostic, so one detector serves BIDMC (PLETH ~0.4) and
  aromring (raw IR ~130,000) without change.
- **Skew-based auto-orientation.** PPG polarity isn't guaranteed (raw MAX30102
  IR can be inverted). A signal with sharp systolic peaks is skewed toward
  them; if skew is negative, flip. Manual override available.
- **Confidence = inter-beat-interval regularity.** Coefficient of variation of
  the IBIs maps to a 0-1 confidence. Two-tier thresholds match FHR exactly
  (per-window >= 0.45; record reliable if >= 50% of windows clear it), so the
  evaluation scaffolding carries over.
- **Everything spatial expressed via fs** (peak distance, window length), so
  the pipeline works at 125 Hz and 25 Hz unchanged.

## Harmonic-doubling correction

Real PPG with a strong dicrotic (diastolic) hump can cause find_peaks to
detect two peaks per cardiac cycle, halving the interval and DOUBLING the
reported rate. Because the false peaks are regularly spaced, the CV stays low
and confidence stays high — so the confidence metric alone does NOT catch it
(it measures regularity, not correctness).

The correction (`_dedouble_peaks`) fires only when BOTH guardrails hold:
1. **Rate gate:** apparent rate > 150 BPM (implausible for resting adult).
2. **Prominence check:** after re-detecting peaks with spacing widened to
   ~1.5x the current median interval, the *dropped* peaks must be markedly
   smaller than the *kept* ones (dropped/kept prominence <= 0.60).

The conjunction is what makes it safe:
- Dicrotic doubling (e.g. bidmc47): high rate + small dropped humps -> corrected.
- Genuine tachycardia: high rate but uniform peaks -> left alone.
- Bigeminy: alternating peaks but normal rate -> left alone (rate gate fails).
- Clean/normal records: below the rate gate -> never even examined.

Implementation note: an earlier even/odd parity-subset version was fragile
(broke when dicrotic humps were missed or varied in height). The re-detection
approach is robust to both.

## Datasets

| Dataset  | Role        | Records | fs (Hz) | Duration | Ground truth          |
|----------|-------------|---------|---------|----------|-----------------------|
| BIDMC    | validation  | 53      | 125     | ~8 min   | HR (ECG) + PULSE (PPG)|
| aromring | unit test   | 1       | 25      | 4 s      | none (sanity only)    |

BIDMC = "BIDMC PPG and Respiration", PhysioNet slug `bidmc`. ICU patients;
noisy/pathological signals make it a hard, conservative test. The PPG channel
is `PLETH`; references are in separate 1 Hz numerics records (`bidmcXXn`).
Channel names carry trailing commas in the headers (handled by `find_channel`).

aromring = the single `ExpectedGoodQualitySignals.csv` from the
MAX30102_by_RF repo. One 4-second capture from the *actual deploy sensor*.

## Results

### BIDMC (validation)
- 46/53 reliable, 7 flagged LOW
- Reliable records, MAE vs HR  : mean 1.13 BPM, median 0.75
- Reliable records, MAE vs PULSE: mean 1.43 BPM, median 0.98
- **44/46 reliable records meet the MAE < 3 BPM target vs HR**
- Scored over high-confidence windows only (a low-confidence window is one the
  device would not display, so it shouldn't count toward displayed accuracy)
- Dual reference: PULSE (same modality, isolates algorithm error) and HR
  (independent ECG, the harder cross-modality test). MAE vs PULSE runs slightly
  lower than vs HR, as expected.

### aromring (unit test)
- 68.2 BPM, 0.92 confidence, single 4 s window, 4 peaks detected
- Physiological-range and expected-value assertions both PASS
- Validates the 25 Hz path, the single-window path, and the auto-orientation
  on real MAX30102 sensor data — all with no detector changes

## Claims defensible from these results

1. **The MHR pipeline meets its accuracy target.** 44/46 reliable BIDMC
   records achieve MAE < 3 BPM against independent ECG-derived ground truth,
   median 0.75 BPM.
2. **The detector generalizes across sample rates and sensors unchanged.**
   Developed on 125 Hz BIDMC; runs on 25 Hz MAX30102 data with no modification.
3. **The two-tier confidence framework transfers from FHR.** Same thresholds,
   applied to an easy signal, yield high reliability (most records 90%+) where
   FHR's hard signal yielded conservative flagging (~58%). The framework adapts.
4. **A known PPG failure mode (harmonic doubling) is detected and corrected**
   with guardrails that leave genuine fast rhythms and clean records untouched.

## Things to mention if asked

- Why ICU data for a maternal monitor: validation measures *algorithm accuracy
  against known truth*, which is independent of patient health and of
  pregnancy. Noisy ICU data is a conservative stress test — good performance
  there implies at-least-as-good on a calm subject. The flagging *thresholds*
  (normal maternal range) come from literature, not from BIDMC.
- bidmc45 (MAE 7.59, flagged OK) and bidmc24 (MAE 5.20, OK at 53%) are the two
  reliable records that miss the target. bidmc45 is genuine physiological
  irregularity (HR and PULSE legitimately diverge), not an algorithm error —
  arguably the *reference* is the unreliable party there. A regularity-based
  flag can pass an irregular-but-real rhythm while MAE-vs-smoothed-reference
  stays elevated. Honest limitation, not a bug.
- The mean-vs-median gap (mean 1.13 vs median 0.75 vs HR) reflects a few
  high-error records pulling the mean; the median is the "typical record".
- BIDMC channel headers contain trailing commas; PPG units are `NU`
  (normalized) and references step in 1 BPM integer increments.
- aromring's panels 3-4 show a single point because the 4 s clip is one
  window; that is correct, not a rendering fault.

## What is NOT claimed

- Performance on pregnant subjects specifically (no such open dataset with
  ground truth exists; the argument is that the signal-processing problem is
  pregnancy-independent).
- Performance during arrhythmia (HR/PULSE decoupling makes ground truth itself
  unreliable; scope is non-arrhythmic HR estimation).
- Per-beat HRV accuracy. This is a heart-rate estimator, not an HRV tool.
- SpO2. The MAX30102 supports it and BIDMC includes a SpO2 reference, but it
  is out of scope for the MHR rate pipeline.

## Next step

Maternal ECG pipeline (AD8232 + Pan-Tompkins + Random Forest arrhythmia
detection on MIT-BIH) — the primary individual component. The MHR PPG work
here is the PPG-based half of maternal monitoring; ECG is the other half and
provides the arrhythmia-detection capability PPG alone cannot.
