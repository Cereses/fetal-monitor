# Phase 1 FHR — Summary

End-of-phase notes for the fetal heart rate (FHR) detection pipeline,
covering both PhysioNet datasets used in Phase 1.

## Pipeline (unchanged across both datasets)

Raw single-channel PCG → Butterworth bandpass (25 Hz to 0.95×Nyquist,
order 4) → Shannon-energy envelope → sliding-window autocorrelation
peak detection → per-window BPM + confidence.

Two-tier confidence:
- **Per-window**: confidence ≥ 0.45 → window is "high confidence"
- **Per-record**: ≥ 50% of windows must be high-confidence → record
  flagged OK; otherwise LOW

The bandpass high-cutoff is clamped to 0.95×Nyquist so the same code
works at any sampling rate without modification (333 Hz on fpcgdb,
1000 Hz on simfpcgdb, and arbitrary future rates).

## Datasets

| Dataset    | Source     | Records | fs (Hz) | Duration | Annotations | Type          |
|------------|------------|---------|---------|----------|-------------|---------------|
| simfpcgdb  | PhysioNet  | 37      | 1000    | ~30 s    | None        | Synthetic     |
| fpcgdb     | PhysioNet  | 26      | 333     | ~20 min  | None        | Real clinical |

simfpcgdb is synthetic PCG (Cesarelli et al. simulator) with controlled
additive noise at known SNR levels from -4.4 dB to -26.7 dB.
fpcgdb is real transabdominal microphone recordings from 26 women in
the last trimester of healthy singleton pregnancies (gestational weeks
31–40), Fetaphon device, 2010.

## Results

### simfpcgdb
- 22/37 reliable (59%)
- All OK records produce estimates of 140.2–140.5 BPM, matching the
  Cesarelli simulator's documented ~140 BPM baseline
- Mean confidence drops monotonically with SNR: 0.81 at -4.4 dB to
  0.24 at -26.7 dB
- Reliability cutoff lands cleanly between SNR levels of approximately
  -22 dB and -22.6 dB
- One LOW-flagged outlier estimates 153.3 BPM (only 2% reliable) —
  correctly identified as untrustworthy

### fpcgdb
- 15/26 reliable (58%)
- All 15 OK records produce estimates in 124.1–159.8 BPM (most in
  130–145), within the physiologically normal fetal range
- Best record: fetal_PCG_p21_GW_39, 0.76 confidence, 97% reliable
- Worst record: fetal_PCG_p02_GW_31, 0.15 confidence, 9% reliable,
  estimates 180.0 BPM (physiologically implausible, correctly flagged
  LOW)
- No correlation between gestational week and reliability — reliability
  is driven by recording quality
- Real recording quality is approximately equivalent to simulated data
  at -23 to -24 dB SNR

## Four claims defensible from these results

1. **The detector generalizes from simulated to real data.**
   58% reliable rate on fpcgdb vs 59% on simfpcgdb, using identical
   thresholds with no re-tuning between datasets.

2. **All reliably-estimated records fall within the physiologically
   normal fetal range.** On both datasets, every OK-flagged record
   produces an estimate within 110–160 BPM.

3. **The two-tier confidence system catches its own failures.**
   The records that produce implausible estimates (153.3 BPM on
   simfpcgdb, 180.0 BPM on fpcgdb) are automatically flagged
   unreliable rather than reported as results.

4. **On real-world data, low confidence is the dominant failure mode
   rather than wrong-number failure.** Even most LOW-flagged fpcgdb
   records produce physiologically plausible estimates — the algorithm
   is conservative rather than fragile.

## Report figures

- `plots/fpcgdb/02_fhr_fetal_PCG_p21_GW_39.png` — best case,
  illustrates the algorithm working on clean real data
- `plots/fpcgdb/02_fhr_fetal_PCG_p02_GW_31.png` — worst case,
  illustrates the confidence system correctly refusing to certify a
  noisy recording
- `results/simfpcgdb_eval.csv` and `results/fpcgdb_eval.csv` — frozen
  full-precision evaluation tables

## Things to mention if asked

- The natural BPM wandering visible in p21's panel 4 is real fetal
  heart rate variability, not algorithm noise. Continuous baseline
  variability is a healthy sign; a flat trace would be clinically
  abnormal.
- The horizontal stripes at ~110 and ~180 BPM in p02's panel 4 are
  autocorrelation lock-points at the search-window edges (min_bpm and
  max_bpm). They appear when the signal has no periodic structure for
  the algorithm to find — expected behaviour for an autocorrelation
  detector on near-noise.
- A handful of blue (high-confidence) dots on p02 cluster at 180 BPM —
  these are windows where artifacts happened to repeat at the
  autocorrelation lock-point. The per-record reliability flag (50%
  threshold) is what catches this: the record still gets flagged LOW
  because the confidently-wrong windows can't push it over the bar.
  This is why the system has two tiers, not one.
- fpcgdb has no annotation files. Evaluation against per-beat ground
  truth was not possible; physiological plausibility, within-record
  self-consistency, and confidence-flag rates were used instead.
- The reported amplitude units in fpcgdb headers are "mV" but the
  values are clearly raw ADC counts. This does not affect processing
  (the pipeline uses relative amplitudes throughout) but is worth
  flagging.

## What is *not* claimed

- Per-beat accuracy. Without annotations on either dataset, no claim
  about beat-detection precision is made — only baseline FHR estimation.
- Performance on pathological recordings. Both datasets contain only
  healthy pregnancies.
- Detection of fetal distress or specific arrhythmias. This is a
  baseline FHR estimator, not a clinical decision tool.

## Next step

Phase 1 MHR (maternal heart rate from PPG/ECG). The script scaffolding
— DATASET constant, plots/<dataset>/, results/<dataset>_eval.csv —
carries over directly. Natural starting dataset: aromring MAX30102 CSV.
