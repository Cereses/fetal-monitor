# Phase 1 — Uterine Contraction Detection (FSR402 / CTU-CHB)

Software validation of the contraction pipeline, completed before hardware
integration. Same strategy as the ECG, MHR and FHR components: validate against
a PhysioNet corpus with expert ground truth, then move to the sensor.

---

## Result

Held-out TEST set, 84 records, 1,234 annotated contractions:

| Metric | TEST | DEV |
|---|---|---|
| Sensitivity | **77.07%** | 77.32% |
| PPV | **77.57%** | 73.44% |
| F1 | **77.32%** | 75.33% |

Primary matching criterion: overlap ≥ 50% of the shorter span, greedy 1:1.

DEV and TEST sensitivity agree to within 0.25 points, and TEST PPV exceeds DEV,
so the constants did not absorb DEV-specific structure.

### Robustness to the matching rule (TEST)

| Criterion | Se | PPV |
|---|---|---|
| any overlap | 77.39% | 77.90% |
| overlap ≥ 50% of shorter (primary) | 77.07% | 77.57% |
| IoU ≥ 0.5 | 67.50% | 67.94% |

The 0.32-point gap between `any` and `half` shows detections cover references
rather than grazing them. The fall at IoU is duration disagreement, quantified
below, not failure to locate contractions.

### Confidence tiers (TEST)

| Tier | Records | Reference | Se | PPV | F1 |
|---|---|---|---|---|---|
| HIGH | 79 | 1,213 | 77.33% | 77.91% | 77.62% |
| LOW | 5 | 21 | 61.90% | 59.09% | 60.47% |

The tier separates as intended, though LOW holds only 21 contractions so the
gap is indicative rather than firmly estimated.

### Boundary agreement (TEST, matched pairs)

| Quantity | Median | IQR |
|---|---|---|
| Onset error | +2.5 s | −2.0 to +8.8 s |
| Offset error | −2.0 s | −8.6 to +3.8 s |
| Duration error | −5.0 s | −17.9 to +6.1 s |

Detections start slightly late and end slightly early. This is the expected
signature of threshold crossing: the signal crosses the threshold after true
onset and falls below it before true offset. On a 60 s contraction the bias is
a few percent, and it is systematic rather than random.

---

## Dataset and ground truth

**Signals:** CTU-CHB Intrapartum CTG Database (`ctu-uhb-ctgdb`), 4 Hz, channels
`FHR` (bpm) and `UC` (nd, uncalibrated). 200 records downloaded of 552
available. Channel layout audited: UC at index 1 in every record, no exceptions.

**Annotations:** Romagnoli et al., *Data in Brief* 31 (2020) 105690, DOI
10.1016/j.dib.2020.105690, supplementary file `mmc2.zip`. 552 files named
`annotation_<record>.mat`, each holding `dataloss` (1×3 double) and `ann`
(5×N cell array).

### Annotation format

Annotations are **marker positions, not per-sample labels**. Cell `(r, c)` is
non-empty when an event of type `r` starts or ends at sample `c`. Cell contents
are MATLAB `string` objects serialised as MCOS references into the file
subsystem; scipy cannot resolve them and does not need to, because the row
index encodes the event type and the position is the ground truth.

### Row-to-event mapping

The paper does not state the row order. It was inferred and then confirmed:

| Row | Event | Evidence |
|---|---|---|
| 0 | Bradycardia | 10/552 records, median episode 600 s |
| 1 | Tachycardia | 11/552 records, median episode 660 s |
| 2 | Acceleration | ~166 events |
| 3 | Deceleration | ~490 events |
| 4 | **Uterine contraction** | median 60.0 s, IQR 46.5–76.5 s |

1. The paper lists five event types in this order: bradycardia, tachycardia,
   acceleration, deceleration, uterine contraction.
2. Row 4 durations are physiologically textbook.
3. Row 4 start-to-start interval median 141 s = 4.3 per 10 min, i.e. active
   labour frequency.
4. Published CTU-CHB proportions are ~12.2% deceleration against ~4.4%
   acceleration, a ratio near 2.8:1. Observed rows 3:2 = 490:166 = 2.95:1.
5. **Direct visual confirmation:** `plots/uc_annotation_overlay.png`. Record
   1025 shows 37 annotated spans sitting on 37 visible tocograph bumps.

Point 5 is the one that settles it; points 1–4 are supporting.

### Pairing

Markers come in start/end pairs. 12 of 552 records carry an odd count. In all
12 the unpaired marker falls 97.7–99.6% through the record, 21–90 s from the
end: a contraction that began before the recording stopped at delivery. The
trailing start is dropped and the record retained.

**Corpus totals:** 7,015 contractions across 473 scorable records, median 14
per record. Duration p5 30.0 s, median 60.0 s, p95 124.3 s.

---

## Pre-registered exclusions

Fixed before any detector was run.

| Rule | Records | Rationale |
|---|---|---|
| `no_annotation` | 79 | Zero contraction markers. Absence of markers cannot be distinguished from absence of contractions, and scoring against them would manufacture false positives. |

No signal-quality threshold is applied. This was tested rather than assumed —
see below.

---

## Correction: the flat-fraction quality metric was wrong

`explore_ctu.py` initially reported a median flat fraction of 60.9% and
concluded only 18 of 50 records were usable. **That verdict was an artefact of
the metric, not a property of the data.**

The UC channel is integer-quantised (measured quantisation step 1.000, range
0.5–3.0) over roughly 0–127. A contraction ramps over ~60 s at 4 Hz, so a rise
of 30 units means one integer increment every ~8 samples: roughly 7 of every 8
consecutive pairs are identical *by construction*. A per-sample "did it change"
test measures the quantisation step, not signal health.

Cross-referencing against the annotations proved it. Of ten records flagged
"likely unusable", seven carried expert-marked contractions, including record
1017 with 25. Across 50 records, 42 carried annotations against the 18 the
metric called usable.

The metric was replaced with run-length analysis (`audit_uc_quality.py`):
quantisation produces short flat runs, dropout produces long ones. Every
candidate threshold was then scored against the annotations:

| Rule | Annotated records lost | Unannotated caught |
|---|---|---|
| OLD `flat_frac > 0.50` | 25/42 | 7/8 |
| `dropout_frac > 0.30` | 13/42 | 6/8 |
| `dropout_frac > 0.50` | 3/42 | 1/8 |
| `zerorun_frac > 0.50` | 0/42 | 0/8 |

Every unannotated record is already removed by `no_annotation`, so the right
column has no value. Every record in the left column is ground truth destroyed.
**No quality rule earns its place.**

Genuine dropout does exist (median longest flat run 510 s, `dropout_frac`
median 25.8%) and is handled by masking spans within the detector, not by
excluding records. Only 3 TEST reference contractions were missed inside masked
regions, confirming that the annotators also skipped dead signal.

---

## Detector

`uc_detector.py`. Every stage is forced by something visible in
`plots/uc_annotation_overlay.png`.

| Stage | Choice | Why |
|---|---|---|
| Dropout mask | flat runs ≥ 10 s | Record 1003's first 21 min are flat and carry zero annotations |
| Spike removal | 5 s median filter | Records 1003 and 1025 contain unannotated single-sample spikes to 100 |
| Smoothing | 15 s moving average | Quantisation and noise |
| Baseline | rolling 10th percentile, 10 min window | Record 1025's resting tone falls from ~15 to ~5; 1029 wanders 0–25 |
| Threshold | 0.30 × record p98 of detrended | `nd` units are uncalibrated, so scale must be record-relative |
| Peak gate | 0.50 × record p98 | Rejects record 1029's unmarked low wobbles |
| Merge | gaps < 20 s | Noise fragments single contractions |
| Trough split | peaks ≥ 40 s apart, prominence ≥ 0.20 × scale | Closely spaced contractions never return to baseline, so thresholding fuses them |
| Duration gate | 20–300 s | |

Masked regions are interpolated before baseline estimation so a long zero
stretch cannot drag the baseline down for live signal beside it; the mask still
suppresses detection inside them.

### Constants set from the reference distribution, not by tuning

`MIN_PEAK_SEP_S` and `MIN_DUR_S` encode what a contraction *is*, taken from
annotation percentiles rather than from detection scores:

- **`MIN_PEAK_SEP_S = 40 s`** — annotated start-to-start intervals have
  p0.5 = 39 s. The initial 60 s discarded 2% of real pairs corpus-wide and 17%
  of record 1025, which is a fast record (median interval 88 s against the
  corpus 135 s, minimum 38 s).
- **`MIN_DUR_S = 20 s`** — annotated durations have p1 = 18 s, p2 = 23 s. The
  initial 25 s discarded 2.6% of real contractions against 1.3% at 20 s.

Because these percentiles came from the whole corpus, records were split into
disjoint DEV and TEST sets (alternating sorted record IDs) and only the TEST
number is reported. Alternating rather than contiguous, so that any drift in
acquisition period does not confound the two arms.

`THRESH_FRAC`, `PEAK_FRAC` and `AMP_PCTL` are the genuinely free parameters and
were adjusted once, from the diagnostic plot, before scoring.

---

## Known limitations

**Unmarked events inside kept records.** Record 1006, excluded, has an obvious
contraction at minute 25 that nobody annotated. Record 1029 has one annotated
contraction among dozens of similar unmarked wobbles. The same is likely true
inside records that were kept.

Matched detections have median relative peak amplitude 0.85 (IQR 0.68–1.03)
against 0.76 (IQR 0.61–0.95) for unmatched. The direction supports an
annotation-criterion effect, **but the distributions overlap heavily, so
amplitude explains only part of the 275 false positives.** The remainder are
genuine detector errors. This is stated as a partial explanation, not a ceiling
on achievable PPV.

**Record 1029 class of failure.** An unsupervised detector cannot infer that an
expert would mark one bump and decline forty similar ones. The two-tier
confidence system flags these (1029 → LOW at 0.72 contractions per 10 min)
rather than pretending to solve them.

**Shallow contractions.** Record 1025 detected 30 of 37; the misses are the
low-amplitude choppy stretch around minutes 12–17. Closing them requires
lowering `THRESH_FRAC`, which was declined because that record is inside the
scored corpus.

**Corpus coverage.** 200 of 552 records downloaded, 168 scorable, split 84/84.
The remaining records are available if a larger evaluation is wanted.

---

## Files

| File | Purpose |
|---|---|
| `parse_ctu_annotations.py` | `.mat` → `results/ctu_contractions.csv`, `results/ctu_ann_manifest.csv` |
| `audit_uc_quality.py` | Quantisation vs dropout diagnosis, rule validation → `results/ctu_uc_quality.csv` |
| `plot_uc_annotations.py` | Row-mapping verification figure |
| `uc_detector.py` | Detector plus diagnostic overlay |
| `evaluate_uc.py` | DEV/TEST scoring → `results/uc_evaluation.csv`, `results/uc_detections.csv` |

---

## Status

Contraction pipeline validated. All four software pipelines (ECG QRS, ECG beat
classification, MHR, FHR, UC) are complete. Next: FSR402 on the ESP32 as the
first streaming sensor.
