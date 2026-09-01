# ECG bring-up — NIFECGDB: placement study, BPM accuracy, and a diagnosed failure mode

**Status:** complete
**Date:** 15 August 2026
**Dataset:** NIFECGDB (PhysioNet), 55 records, one pregnant subject, 22+1 to 40+2 weeks
**Pre-registration:** `notes/PREREGISTRATION_nifecgdb.md`, committed before `04` was written

---

## 0. Scope statement — read before the numbers

This phase answered one question the project could not otherwise answer:

> **What does moving maternal ECG electrodes from the chest to the abdomen, on
> a pregnant woman, cost the frozen Pan-Tompkins detector?**

Phase A validated the AD8232 chain on **chest** electrodes on a **non-pregnant
male**. The device is for a pregnant woman, and a belt would put electrodes on
the abdomen. Nothing in the project established what that change costs.

NIFECGDB records both placements **simultaneously in every recording**, so
subject, session, gestational age, amplifier and electrodes are held constant
and only placement varies.

Four limits define what the result means.

**One subject.** Fifty-five recordings of one woman are not fifty-five women.
Nothing here supports any claim about variation across women, at any sample
size. The large beat counts buy precision, not generality.

**Different hardware.** A research analog amplifier at 1000 Hz, 16-bit, with a
**50 Hz notch applied at acquisition** and a 0.01–100 Hz passband. This
project's chain is an AD8232 at 250 Hz, 12-bit, with clipping R peaks and **no
mains notch** — hum was fought down to 3.9% by electrode contact alone. The
electrodes are the same class (Ag/AgCl); the amplifier and digitiser are not.

**This measures the ALGORITHM, not the device.** Rate accuracy on pregnant
chest ECG. Not the accuracy of this project's hardware on a pregnant woman,
which remains untested and requires ethics approval.

**The reference is not exhaustive.** ~6–7% of beats are missing from the
annotations (§4). PPV is reported as a lower bound, justified by measurement
(§5.3), not asserted.

---

## 1. What this phase did

| script | role |
|---|---|
| `nifecgdb_io.py` | native-rate loader (pyedflib); wfdb cannot read this database correctly |
| `03_explore_nifecgdb.py` | established what the annotations ARE before anything was built |
| `04_nifecgdb_evaluate.py` | thoracic vs abdominal detection, resolves the pre-registered predictions |
| `05_mhr_from_ecg.py` | windowed BPM accuracy — the number the device would display |
| `06_diagnose_window_bpm.py` | diagnosed the one channel that failed |

The detector (`02_qrs_detect.py`) and scorer (`03_qrs_evaluate.py`) are
**imported by path, never copied**, and their SHA-256 is recorded per run.
Neither was modified at any point in this phase.

---

## 2. The loading problem — wfdb silently destroys this database

**This is the single most important technical finding of the phase**, because
every subsequent number depended on catching it.

`wfdb.io.convert.read_edf` computes the sampling rate as the **GCD across all
channels, including the `EDF Annotations` text channel**:

```python
sample_rate    = [int(i / block_duration) for i in samps_per_block]
fs             = functools.reduce(math.gcd, sample_rate)
samps_per_frame= [int(s / min(samps_per_block)) for s in samps_per_block]
```

From the header of `ecgca102.edf`:

```
Number of data records       : 45
Duration of each data record : 6.0 s
Number of Samples per Record : [6000, 6000, 6000, 6000, 6000, 600]
```

6000 / 6.0 = **1000 Hz** for the five ECG channels; 600 / 6.0 = **100 Hz** for
the TAL text channel. `gcd(1000, 100) = 100`. The reader then **averages every
10 consecutive samples into one**, returning 27,000 samples where the file
holds 270,000.

Three consequences:

**It aliases, it does not merely smooth.** Prefiltering passes content to
100 Hz. Decimating to 100 Hz puts Nyquist at 50 Hz, and a 10-point boxcar
(first null at 100 Hz) attenuates but does not remove the 50–100 Hz band, which
folds back.

**Physical units are wrong in that code path.** The mixed-rate branch computes
`mean(raw − baseline/gain)` where the equal-rate branch correctly computes
`(mean(raw) − baseline)/gain`. Irrelevant to Pan-Tompkins (scale-invariant via
squaring plus the adaptive threshold) but fatal to any amplitude claim — and
this database mixes units across the very groups being compared (**mV
thoracic, µV abdominal**).

**It would have produced a silent, plausible, wrong result.** Annotations are
indexed at 1000 Hz, so they ran ~10× past the end of the decimated signal
(269,292 against 27,000, ratio 9.97). Sensitivity would have come out near
zero and the obvious next move would have been to hunt for a detector bug that
does not exist.

`read_edf`'s own docstring documents the behaviour and points at
`smooth_frames=False` — **a parameter `read_edf` does not accept.**

**Resolution:** all loading goes through `nifecgdb_io.py` using `pyedflib`,
which returns native per-channel rates and does not expose the TAL channel as
a signal at all (`signals_in_file` is 5, where wfdb reported 6).
`check_alignment()` asserts that annotation indices land inside the signal, so
this class of failure becomes a named error rather than a bad number.

**Gestational age** is parsed from the EDF patient field. The EDF+ local
patient identification packs `code sex birthdate name`, and this database
stores `X F X Gestation_22+1` — so the age is in the **name** slot. The loader
scans every header string rather than naming a field.

---

## 3. The dataset, verified before scoring

| property | value |
|---|---|
| records | 55, all with `.edf` and `.edf.qrs` |
| sampling rate | 1000 Hz, all records, all channels |
| annotation scale | 1.0 (annotations index the signal directly) |
| alignment | all 55 pass |
| gestation | 22.1–40.3 weeks, 55/55 parsed |
| reference beats | **47,980** |
| total duration | 549.7 min (9.16 h) |
| record length | 114 s to 2780 s |
| channel sets | `(2 thoracic, 3 abdominal)` × 17, `(2, 4)` × 38 |
| units | thoracic **mV**, abdominal **µV** |
| prefiltering | `HP:0.01Hz LP:100Hz NF:50Hz` at acquisition |
| electrodes | Ag/AgCl |

**Channel-records scored: 313** — 110 thoracic, 203 abdominal.

---

## 4. Reference quality — measured, and it constrains reporting

**The `.qrs` annotations are MATERNAL R-peaks**, median 92.7 BPM across 55
records (fetal would be ~140). This settles a contradiction in the published
literature, where one source describes them as maternal, another as fetal, and
a third as absent. Settled empirically, not by citation.

**They are not exhaustive.** Two independent estimates agree:

| method | estimate |
|---|---|
| RR intervals above 1.75× record median, clustering at ~2× | median **7.33%** |
| beats ÷ duration versus median RR rate, pooled | **5.84%** |

53/55 records exceed 5% by the first method, and the count of intervals between
1.8× and 2.2× median tracks the flagged count almost exactly — the signature of
a **skipped beat** rather than a physiological pause.

This is expected. NIFECGDB exists to support **fetal** ECG research; maternal
annotations are supplied so the maternal component can be located and
cancelled. They are a means, not the product, and are not the exhaustive
cardiologist-corrected reference MIT-BIH provides.

15/55 records also contain intervals below the 200 ms refractory floor, which
are physiologically impossible. `ecgca998` is the outlier at **10.12%**; the
rest sit between 0.07% and 3.42%.

---

## 5. Placement study (`04`)

### 5.1 Headline

| group | channels | ref beats | TP | FP | FN | Se | PPV |
|---|---|---|---|---|---|---|---|
| thoracic | 110 | 95,960 | 95,661 | 7,962 | 299 | **99.69%** | 92.32% |
| abdominal | 203 | 183,365 | 168,185 | 17,541 | 15,180 | **91.72%** | 90.56% |

**Pooled placement cost: 7.97 percentage points of sensitivity.**

`ecgca998` is **retained in the headline** and reported separately, per the
pre-registration (§8.1 records a deviation that was found and corrected).
Supplementary: Se 57.47%, PPV 89.57% over 6 channels.

### 5.2 The pooled abdominal number is misleading — the distribution is bimodal

| abdominal channel Se | channels | share |
|---|---|---|
| ≥99% | 144 | **70.9%** |
| 95–99% | 21 | 10.3% |
| 90–95% | 2 | 1.0% |
| 50–90% | 11 | 5.4% |
| **<50%** | **25** | **12.3%** |

**Median 100.00%. Mean 90.11%.** Thoracic median is also 100%, with only 2 of
110 channels below 95%.

So "abdominal costs 7.97 points" is wrong as a description. **Seven in ten
abdominal channels are indistinguishable from thoracic; one in eight collapses
below 50%**, some to under 2%.

Per record:

- **38 records:** no abdominal channel below 90%
- **13 records:** some channels fine, others fail — the outcome depends on
  *which electrode you pick*
- **4 records:** every abdominal channel fails

`ecgca848` has `Abdomen_4` at 1.92% and `Abdomen_3` at 4.41% while other
channels in the same recording are fine. Same woman, same instant, same
amplifier — **electrode position alone.**

The per-record paired gap has median **+0.00 pp** with a range to +71.84. Most
recordings show no difference at all; a minority fail completely.

### 5.3 False positives, and the test that makes the PPV framing honest

Every unmatched detection was classified as falling inside a double-length
reference gap, outside the annotated span, or genuine.

| group | in gap | genuine | outside | % in gap |
|---|---|---|---|---|
| thoracic | 7,247 | 648 | 67 | 91.8% |
| abdominal | 13,405 | 4,022 | 114 | 76.9% |
| **overall** | | | | **81.6%** |

**P4 passes at a 70% threshold.** The detector is demonstrably finding beats
the reference omitted, so reporting PPV as a lower bound is justified *by
measurement* rather than asserted.

Genuine false positives per 1000 reference beats: **thoracic 6.75, abdominal
21.93** — 3.2× worse. Gap-corrected PPV: **thoracic 99.26%, abdominal 97.60%**.

### 5.4 Pre-registered predictions

| | prediction | result | |
|---|---|---|---|
| P1 | thoracic Se ≥ 98% | 99.69% | **HELD** |
| P2 | abdominal worse in ≥75% of records | 42% | **DID NOT HOLD** |
| P3 | PPV in 92–94% for both groups | thor 92.32%, abd 90.56% | **DID NOT HOLD** |
| P4 | ≥70% of FPs in reference gaps | 81.6% | **HELD** |

**P2 failed for an informative reason and that is the finding.** It predicted
uniform degradation; the truth is that most records show no difference and a
minority fail catastrophically. P2 failing is more useful than P2 holding would
have been, and it is exactly what pre-registration is for.

**P3 half-held, and the half that held matters.** Thoracic PPV landed inside
the predicted 92–94% band, independently confirming the reference-incompleteness
ceiling arithmetic. Abdominal fell below it because of the excess genuine false
positives in §5.3.

### 5.5 No incidental fetal detection

Abdominal channels carry ~118,700 channel-seconds, so a fetus at 140 BPM
contributes ~276,900 complexes. The detector produced **4,022** genuine false
positives — **1.45%**.

**The maternal QRS detector does not incidentally detect fetal beats.** The
fetal complex is 10–50 µV against maternal 100–1000 µV at the abdomen; the
adaptive threshold learns SPKI from maternal beats and pushes the decision
level far above anything fetal, while every rejected fetal peak *raises* NPKI.
Refractory and T-wave windows discard whatever survives.

This is the project's own evidence that **FHR cannot come from a single
abdominal ECG channel**, and that the MAX4466 acoustic path is the correct
architecture rather than a workaround.

### 5.6 Gestational age

No trend. Failures scatter across 23.9, 29.3, 30.7, 31.6, 32.4, 33.6 and 37.1
weeks. Consistent with electrode repositioning driving the variance —
**PhysioNet documents that electrode positions were varied to improve SNR** —
which the pre-registration named as an unavoidable confound. Reported
qualitatively; no statistic computed, because none could be interpreted.

---

## 6. Maternal BPM accuracy (`05`)

Everything above scores **beat positions**. This scores **the number the device
would display**, which is derived from RR intervals and is a different quantity.

**Method:** 8-second non-overlapping windows, **median** RR within window, ≥4
intervals required. Median rather than mean because ~7% of reference beats are
missing, doubling roughly one interval per window — a median is immune, a mean
is not. **The reference incompleteness that caps PPV at ~94% therefore does not
contaminate the rate.** 8 seconds matches the frozen PPG pipeline so the two
maternal heart rate paths are directly comparable.

Thoracic channels only. Abdominal excluded by design following §5.

### 6.1 Pooled result — 110 channels, 8,190 windows

| statistic | value |
|---|---|
| reference rate | 93.5 BPM (75.0–213.9) |
| **median \|error\|** | **0.146 BPM** |
| MAE | 0.939 BPM |
| bias | +0.689 BPM |
| SD of error | 8.877 BPM |
| 95% limits of agreement | [−16.71, +18.09] BPM |
| max \|error\| | 163.51 BPM |
| **within 2 BPM** | **98.47%** of windows |
| within 5 BPM | 99.38% of windows |
| window coverage | 99.7% |

Per channel: **median MAE 0.251 BPM, mean 1.451 BPM, max 119.39.**

| per-channel MAE | channels |
|---|---|
| <0.5 BPM | **100** |
| 0.5–1 | 3 |
| 1–2 | 3 |
| 2–5 | 3 |
| >5 | **1** |

**Excluding the single failing channel: 0.343 BPM over 8,149 windows.** That
channel is not excluded from the headline; the figure is given to show what a
working channel does.

**Device capture:** 3 windows from the 30 s AD8232 record — MAE 0.13 BPM, max
0.26 BPM, reference rate 62.2 BPM. Consistent, but three windows carry no
weight.

**Comparability:** the frozen PPG pipeline reports 1.13 BPM MAE on BIDMC at the
same window length with the same estimator. Different databases and subjects,
so this compares *paths*, not a controlled head-to-head.

### 6.2 How to report this

**Do not quote the limits of agreement alone.** An SD 60× the median error is
not a spread — it is a tight core plus a small number of catastrophic outliers,
and ±17 BPM misdescribes a system that is within 2 BPM 98.47% of the time.

Report the **median error and the percentage-within-bounds** as the headline,
with MAE, SD and LoA given in full alongside the failure mode in §7.

---

## 7. The failure mode: T-wave over-detection (`06`)

One channel, `ecgca699 Thorax_2`, produced 119.39 BPM MAE — 0.4% of the windows
carrying most of the pooled error.

### 7.1 The reference is sound

| | |
|---|---|
| median rate | 108.7 BPM |
| **maximum rate** | **119.0 BPM** |
| lag-1 RR correlation | −0.014 (no alternation) |
| half-median cluster | 0.00% |

A median of 108.7 BPM at 32w3d is **plausible physiology** — sinus tachycardia
begins above 100 and pregnancy adds 10–20. No stretch of this record is
doubled.

### 7.2 The detector is doubling — measured, not inferred

Same record, same instant, two chest channels:

| channel | MAE | mean ref/det ratio |
|---|---|---|
| `Thorax_1` | **0.22 BPM** | 1.00 |
| `Thorax_2` | **119.39 BPM** | 0.48 |

Worst window: 13 reference beats, **28 detections**, reference 104.3 BPM against
detector 267.9 BPM.

**The decisive measurement.** For every detection with no reference beat within
±150 ms, the offset from the preceding R peak was measured:

> **638 unmatched detections. Median offset 211 ms. IQR 6 ms. 38% of the RR
> interval.**

**Six milliseconds of spread across 638 events** is not scatter; it is a fixed
landmark in the cardiac cycle, and 38% of RR is where the T wave sits. The
detector is counting T waves as beats.

**Why the reported rate is 268 and not 209.** Intervals alternate ~211 ms
(R→T) and ~341 ms (T→R). With one more short than long, the *median* lands on
the short group: 60/0.224 = 268 BPM. Counting doubles; a median-RR rate
over-reports by more.

**Why this channel and not the other.** `Thorax_1` has R ≈ 3.4 mV against
T ≈ 1.3 mV (ratio 0.38); `Thorax_2` has R ≈ 2.5 mV against T ≈ 1.7 mV (**0.68**).
Stage 5 rejects a peak inside the 360 ms window only when its slope is under
half the previous beat's. On this channel that test loses.

### 7.3 Why this matters to the device

**The AD8232 is single-channel.** `Thorax_1` and `Thorax_2` are the same heart,
same instant, same amplifier — one gives 0.22 BPM MAE, the other reports double
the rate. The only difference is electrode position.

If an AD8232 lands in a `Thorax_2`-like position, **the device would report
~210 BPM to a woman whose heart rate is 105.** That is a false alarm capable of
sending someone to hospital, not an academic error rate.

### 7.4 Mitigation — and one approach that does NOT work

**Do not modify `02_qrs_detect.py`.** It produced the MIT-BIH headline and the
Phase A device result. The mitigation belongs on the **output**, not in the
detector.

**A lag-1 RR correlation test does not work at window scale.** In the five worst
windows the values are −0.46, −0.09, −0.11, −0.39, −0.06 — only two clear a
−0.30 threshold despite all five being near-complete doubling. On 8-second
windows the statistic is too noisy. This was proposed and then measured and
rejected.

**What is defensible: a plausibility ceiling.** A resting pregnant woman at home
sits at 60–110 BPM. Above ~140 BPM is either a genuine emergency or a detection
failure, and in both cases the correct output is *"reading unreliable —
reposition the electrodes; seek advice if it persists"*, never a number.

This is the **second independent argument** for a confidence check on reported
output, alongside §5.2. The MHR and FHR pipelines already carry two-tier
confidence systems; this is the same pattern with two measurements behind it.

---

## 8. Errors found in this phase's own tooling

Recorded because each produced a wrong number that looked reasonable.

### 8.1 The headline excluded a record the pre-registration retained

`04` implemented `head = [r for r in rows if not r['supplementary']]`, dropping
`ecgca998` from the headline where the pre-registration said to retain it. It
flattered the result — thoracic Se 99.77% instead of 99.69%, abdominal 92.18%
instead of 91.72%. Found by reconciling channel counts (108/199 against the
expected 110/203) and corrected. **The numbers in §5 are the corrected ones.**

### 8.2 A physiological filter hid the worst failures

`05` applied a 30–220 BPM bound to **both** the reference and the detector. The
failing channel reported ~268 BPM, above the ceiling, so **24 of its 41 windows
were silently discarded** — the filter removed precisely the windows where the
detector failed worst.

Corrected so the bound applies to the reference only: implausible *annotations*
cannot be scored against, but an implausible *detector output* is a wrong
answer and must be scored as one. Pooled MAE moved 0.547 → 0.939 BPM, and the
entire change (+3,224 BPM of summed error) reconciles to that one channel
(+3,226). **Median error moved 0.145 → 0.146** — the core never changed; only
the hidden tail.

### 8.3 The diagnostic inspected the healthy channel

`06` originally took `thor[0]` unconditionally, which on `ecgca699` is
`Thorax_1` at 0.22 BPM MAE — while the failure was on `Thorax_2`. It returned a
confident wrong answer. Corrected to evaluate every thoracic channel and
diagnose the worst.

### 8.4 A mechanism was claimed from a plot before it was measured

T-wave over-detection was asserted from marker positions in a rendered figure.
Those positions are misleading: `refine_to_r_peak` returns the index of maximum
|bandpassed| signal and the marker is drawn at the **raw** amplitude there, so
a detection on a genuine QRS can render well below the visible R peak. The
claim was challenged, withdrawn, and replaced with the phase measurement in
§7.2. The conclusion survived; the reasoning that produced it did not.

---

## 9. What this establishes, and what it does not

### Established

- Moving maternal ECG electrodes from chest to abdomen costs **7.97 points of
  pooled sensitivity**, but the effect is **bimodal**: 71% of abdominal
  channels are indistinguishable from thoracic and 12% collapse below 50%.
- On **chest** electrodes on a pregnant woman, the frozen detector achieves
  **99.69% sensitivity over 95,960 beats** and reports heart rate with a
  **median error of 0.146 BPM, within 2 BPM in 98.47% of 8-second windows**.
- A maternal QRS detector does **not** incidentally detect fetal beats (1.45%),
  confirming acoustic PCG as the correct FHR architecture.
- The detector has a real, characterised failure mode — **T-wave
  over-detection**, mechanism measured to 211 ms ± 3 ms — which is
  electrode-position dependent and therefore a live risk for a single-channel
  device.

### Not established

- **Anything about variation across women.** One subject.
- **Anything about fetal heart rate.** The annotations are maternal.
- **Anything about this project's hardware on a pregnant woman.** Different
  amplifier, different resolution, and a 50 Hz notch this project does not have.
- **A benchmark against MIT-BIH or Phase A.** Different databases, subjects,
  rates and reference quality. Context only.

### Statement for the writeup

> The frozen QRS detector was run unmodified against 55 multichannel recordings
> from a pregnant subject spanning 22 to 40 weeks. On chest electrodes it
> detected 99.69% of 95,960 annotated beats and reported heart rate with a
> median error of 0.146 BPM, within 2 BPM in 98.47% of 8-second windows. Moving
> to abdominal electrodes cost 7.97 points of pooled sensitivity, but the effect
> is bimodal: 71% of abdominal channels performed indistinguishably from
> thoracic while 12% failed below 50%, depending on electrode position. On this
> evidence the design adopted chest placement. One thoracic channel exhibited
> T-wave over-detection, reporting approximately double the true rate; the
> mechanism was measured and a plausibility check on reported output is
> recommended rather than any change to the frozen detector.

---

## 10. Design decisions taken on this evidence

**Chest placement adopted.** Not because abdominal placement offers nothing,
but because it carries a measured 12% catastrophic failure rate that chest
placement does not. Since FHR comes from the MAX4466, nothing requires
abdominal electrodes.

**MAX30102 dropped.** MHR is covered by the AD8232; the MAX30102 was only ever
a redundant MHR path or an SpO2 channel that cannot be calibrated without
induced-hypoxia reference data. Sensor set is now MAX4466 → FHR, AD8232 → MHR,
FSR402 → contraction timing.

**Plausibility check on reported rate — recommended, not yet built.** Two
independent measurements now support it (§5.2, §7.3).

---

## 11. Requirements

### Closed

- [x] Maternal QRS detection validated on **pregnant chest ECG**
- [x] Cost of abdominal placement measured, with the bimodal structure characterised
- [x] Maternal **BPM accuracy** quantified — the number the device displays
- [x] Fetal-ECG-via-AD8232 ruled out by measurement
- [x] The one catastrophic failure diagnosed to mechanism
- [x] Pre-registration committed before the evaluation script existed; all four predictions resolved

### Open

- [ ] **Plausibility ceiling on reported rate.** Highest-value software item.
      Reject and warn above ~140 BPM rather than displaying a number.
- [ ] **Digital 50 Hz notch** in the capture script, validated offline against
      the five existing captures. Zero cost, directly de-risks any live session.
- [ ] **Live-participant protocol** documented as designed-but-not-executed,
      with the ethics requirement named (§0 of the live protocol).
- [ ] **Non-pregnant second subject** on chest placement — broadens beyond one
      male subject without pregnancy-specific approval. Check with supervisor
      whether KNUST requires approval for this.
- [ ] `tolerance_samples` recorded in the device eval JSON (carried over from Phase A).

---

## 12. Files

```
hardware/ecg_bringup/
├── nifecgdb_io.py                     native-rate loader (pyedflib)
├── 03_explore_nifecgdb.py
├── 04_nifecgdb_evaluate.py
├── 05_mhr_from_ecg.py
├── 06_diagnose_window_bpm.py
├── quick_check.py
├── results/
│   ├── nifecgdb_channel_results.csv
│   ├── nifecgdb_summary.json
│   ├── mhr_ecg_channel_results.csv
│   ├── mhr_ecg_summary.json
│   └── 06_diagnose_ecgca699.json
├── plots/nifecgdb/
│   ├── 03_explore_ecgca102.png
│   ├── nifecgdb_placement_comparison.png
│   ├── mhr_ecg_agreement.png
│   └── 06_diagnose_ecgca699.png
└── notes/
    ├── notes.md                       Phase A
    ├── PREREGISTRATION_nifecgdb.md
    └── notes_nifecgdb.md              this file
```

`download_nifecgdb.py` sits at the repo root beside `download_data.py`, which
was left unmodified — it is correct for every WFDB database in the project and
cannot handle this EDF one.

Imported, never copied: `02_qrs_detect.py`, `03_qrs_evaluate.py`.

---

## 13. Next

**MAX4466 / fetal PCG bring-up.** Under the project's own framing — give a
pregnant woman a cheap way to check her own and the baby's heart rate — the
fetal channel is the entire differentiator, and it has not been started. It is
also the only channel that can fail at *acquisition*, where no pipeline work
rescues it. This is the critical path.

Then the multi-rate ISR (PCG ~1 kHz, ECG 250 Hz, UC ~4 Hz on a shared clock —
required because FHR and contraction timing must share a time base for
deceleration phase to mean anything), and the causal-filter decision, which
should be resolved on paper before it consumes a week.
