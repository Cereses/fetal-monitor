# PRE-REGISTRATION — NIFECGDB thoracic vs abdominal maternal QRS detection

**Written:** 14 August 2026
**Status:** to be committed BEFORE `04_nifecgdb_evaluate.py` is written or run.
**Phase:** ECG bring-up, offline placement study

---

## 0. The question

The AD8232 smoke test and Phase A used **chest/torso** electrodes on a
non-pregnant male. The device is intended for a pregnant woman, and in a belt
configuration the electrodes would sit on the **abdomen**.

Nothing in the project so far establishes what that change costs. This study
answers exactly one question:

> **What does moving maternal ECG electrodes from the thorax to the abdomen,
> in a pregnant woman, cost the frozen Pan-Tompkins QRS detector?**

NIFECGDB records both placements **simultaneously in every record**, so subject,
session, gestational age, equipment and amplifier are all held constant. Only
placement varies. That is the experiment.

It is a *placement* study, not a validation. What it cannot support is in §7.

---

## 1. Data, verified before this document was written

`data/nifecgdb/`, 55 records, one subject, gestational weeks 22+1 to 40+2.

| property | value |
|---|---|
| records | 55, each with `.edf` and `.edf.qrs` |
| sampling rate | **1000 Hz, all records, all channels** |
| annotation scale | 1.0 (annotations index the signal directly) |
| alignment | all 55 pass (`check_alignment`) |
| gestation parsed | 55/55, from the EDF `patientname` field |
| reference beats | **47,980** |
| total duration | 549.7 min (9.16 h) |
| record length | 114 s to 2780 s |
| channel sets | 2 — `Thorax_1..2 + Abdomen_1..3` (17 records), `+ Abdomen_4` (38 records) |
| units | thoracic **mV**, abdominal **µV** |
| prefiltering | `HP:0.01Hz LP:100Hz NF:50Hz` (applied at acquisition) |

**Channel-records to score: 313** — 110 thoracic, 203 abdominal.

### Loading is not with wfdb, and that is load-bearing

`wfdb.io.convert.read_edf` computes `fs` as the **GCD across all channels
including the `EDF Annotations` text channel**. Here the ECG channels carry
6000 samples per 6 s block (1000 Hz) and the TAL channel carries 600 (100 Hz),
so `gcd(1000, 100) = 100`. It then **averages every 10 consecutive samples**
into one output sample.

That is not a resolution loss. Prefiltering passes content to 100 Hz, so
decimating to 100 Hz puts Nyquist at 50 Hz and the 10-point boxcar (first null
at 100 Hz) attenuates but does not remove the 50–100 Hz band, which **folds
back**. The same code path also mis-scales physical units
(`mean(raw − baseline/gain)` instead of `(mean(raw) − baseline)/gain`).

Had this gone unnoticed, annotation indices would have run ~10× past the end of
the signal, producing near-zero sensitivity and a hunt for a detector bug that
does not exist.

All loading therefore goes through **`nifecgdb_io.py`** (pyedflib), which
returns native per-channel rates and never exposes the TAL channel as a signal.

---

## 2. Reference quality — measured, and it constrains what can be reported

The `.qrs` annotations are **maternal R-peaks** (median 92.7 BPM across 55
records; fetal would be ~140). This resolves a contradiction in the published
literature, where one paper describes them as maternal, another as fetal, and a
third as absent. Settled empirically, not by citation.

**They are NOT exhaustive.** Two independent estimates agree:

| method | estimate |
|---|---|
| RR intervals above 1.75× record median, clustering at ~2× | median **7.33%**, worst 9.18% |
| beats ÷ duration vs median RR rate, pooled over all records | **5.84%** |

53/55 records exceed 5% by the first method. The `~2x` column tracks the flagged
count almost exactly (e.g. `ecgca711`: 110 flagged, 107 between 1.8× and 2.2×
median), which is the unambiguous signature of a **skipped beat** rather than a
physiological pause. A uniform ~6–7% across 55 recordings is systematic
incompleteness, not physiology.

This is expected: NIFECGDB exists to support **fetal** ECG extraction research.
Maternal annotations are supplied so the maternal component can be located and
cancelled — they are a means, not the product. MIT-BIH's are exhaustive and
cardiologist-corrected; these are not, and must not be treated as equivalent.

15/55 records also contain intervals below the 200 ms refractory floor, which
are physiologically impossible and therefore spurious annotations. `ecgca998` is
the outlier at **10.12%**; the rest sit between 0.07% and 3.42%.

---

## 3. Reporting convention — fixed here, before any number is seen

**Sensitivity is the headline.** It asks: did the detector find the beats that
are annotated? Reference incompleteness does not affect it.

**PPV is reported as a LOWER BOUND.** Every real beat the reference omits is
counted as a false positive through no fault of the detector.

This convention is chosen now precisely so it cannot be chosen later in response
to a disappointing number.

### Predicted PPV ceiling

If the detector is perfect and the reference is ~5.84% incomplete, it finds
~50,957 beats of which 47,980 are annotated:

> **PPV ceiling ≈ 94.2%** (93.2% using the per-record RR-gap estimate instead).
> **Expect PPV in the range 92–94% even from a flawless detector.**

### The test that turns this from an excuse into a finding

For every unmatched detection, check whether it falls **inside a double-length
reference gap**. If the large majority do, the detector is demonstrably finding
beats the reference omitted, and that is a positive result rather than a caveat.
If they do not, the false positives are real and the excuse does not apply.

**This test is mandatory. The PPV lower-bound framing is only defensible if it
passes.**

---

## 4. Scoring

Imported from `03_qrs_evaluate.py`, never reimplemented:

- tolerance `TOLERANCE_S = 0.150` → **150 samples at 1000 Hz**
- `match_beats`: greedy, nearest-first, **one-to-one**
- Se = TP/(TP+FN), PPV = TP/(TP+FP), DER = (FP+FN)/reference

Detector imported from `02_qrs_detect.py`, unmodified. SHA-256 of both recorded
in the results file, as in Phase A.

**No detector parameter may be changed in response to these results.** Those two
files produced the MIT-BIH headline (Se 99.52%, PPV 99.35%) and the Phase A
device result (31/31). Editing them retroactively invalidates both.

---

## 5. Design

Every channel is scored **individually** against the same maternal reference for
its record. No "best channel" selection — that would be choosing after seeing
results.

- **Headline:** pooled Se and PPV per group (thoracic, abdominal) across all
  channel-records.
- **Paired:** per record, thoracic median Se vs abdominal median Se. Paired by
  record, so subject/session/gestation cancel.
- **Reported:** the *gap*, not either number alone. The gap is the cost of
  placement, and it is the deliverable.

### Exclusions

**None on loading grounds** — all 55 records pass every check.

`ecgca998` (10.12% sub-refractory annotations) is **retained in the headline**
and reported **separately** as a supplementary line, following the same
headline/supplementary pattern already used for mitdb records 102/104. Dropping
a single record on a threshold that excludes exactly that record would look like
tuning; reporting it twice hides nothing.

---

## 6. Predictions

Falsifiable, recorded before the detector runs.

**P1 — Thoracic sensitivity ≥ 98% pooled.** Thoracic maternal ECG is a clean,
dominant signal; this is close to the MIT-BIH condition.

**P2 — Abdominal sensitivity is lower than thoracic in ≥ 75% of records.**
Abdominal maternal QRS is smaller (µV vs mV ranges) with fetal ECG and uterine
EMG in band. Magnitude is deliberately not predicted — that gap is the thing
being measured.

**P3 — PPV lands in 92–94% for BOTH groups**, dominated by reference
incompleteness rather than by placement. If PPV differs markedly between groups,
the incompleteness explanation is insufficient and something else is happening.

**P4 — Unmatched detections fall predominantly inside double-length reference
gaps** (see §3). Threshold for "predominantly": ≥ 70%.

**Exploratory, with NO directional prediction: performance versus gestational
age.** Ages span 22.1–40.3 weeks. A directional prediction is deliberately
withheld because the mechanism is unclear (the growing uterus displaces the
heart, but whether that increases or decreases abdominal QRS amplitude is not
something this study can establish) and because **electrode positions were not
fixed across recordings and were sometimes moved to improve SNR**. Any trend
across gestation is confounded with placement changes and must not be attributed
to gestation alone. Several records share a gestational age (three at 22+1),
which permits partial separation of session from gestation.

---

## 7. Scope — what this can and cannot establish

### Can establish

The cost, in Se and PPV, of moving maternal ECG electrodes from thorax to
abdomen in a pregnant woman, with everything else held constant, for the frozen
Pan-Tompkins detector.

### Cannot establish

- **Anything about population variation.** ONE subject. 55 recordings of one
  woman are not 55 women. No claim about generalisation across women is
  supported by this data at any sample size.
- **Anything about fetal heart rate.** Annotations are maternal. FHR in this
  project comes from the MAX4466 (acoustic), not from ECG.
- **Anything about this project's hardware.** Different amplifier, different
  electrodes, 1000 Hz, and prefiltering (`HP:0.01 LP:100 NF:50`) applied at
  acquisition that the AD8232 capture chain does not apply.
- **A benchmark against MIT-BIH or against Phase A.** Different databases,
  subjects, rates and reference quality. Context only.

### Statistical note

47,980 reference beats. A perfect score would carry a 95% lower bound of
99.99%, against 88.78% for Phase A's 31 beats. **The confidence interval is no
longer the limiting factor — reference incompleteness is.** Report the
incompleteness, not the interval, as the dominant uncertainty.

---

## 8. Integrity

`04_nifecgdb_evaluate.py` will record, per run:

- SHA-256 of `02_qrs_detect.py` and `03_qrs_evaluate.py`
- SHA-256 of this file
- pyedflib and wfdb versions
- per-channel-record Se, PPV, DER, TP/FP/FN, timing error
- gestational age and channel group per record

**This document is committed before that script is written.** The commit
timestamp is what makes the ordering checkable by someone who was not here.
