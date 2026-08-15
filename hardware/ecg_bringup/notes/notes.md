# ECG bring-up — Phase A: frozen Pan-Tompkins on AD8232 device data

**Status:** complete
**Date:** 13 August 2026
**Capture used:** `raw_ecg_capture_20260729_180936.csv` (recorded 29 July 2026)
**Freeze commit:** `4bdd372`

---

## 0. Scope statement — read this before the numbers

This phase ran the **already-frozen** Pan-Tompkins QRS detector against a
**pre-existing** recording from this project's own AD8232 hardware. It closes
the gap recorded since the smoke test: *"AD8232 not yet run against the frozen
Pan-Tompkins pipeline."*

Three limits define what the result means, and none of them are incidental.

**It reuses smoke-test data, not a purpose-built capture.** The recording is
from 29 July; only the reference and the detector run are new. Nothing was
re-wired or re-recorded. This validates the **pipeline against device data**.
It does not re-validate the circuit, the firmware or the capture script — those
rest on the smoke test, separately documented.

**The capture was selected after inspection.** Five captures exist. `180936`
was chosen because it was the cleanest. The other four were recorded while
electrode contact was still being improved (mains hum fell 90.7% -> 79.9% ->
3.9% across the session). So this result characterises the detector under
**good electrode contact**, not under typical or worst-case contact. Calling it
"the AD8232 result" without that qualifier overstates it. See §9 for the
supplementary runs that would remove the qualifier.

**The reference is not an independent measurement.** The plateau proposal
(§4) is an independent *algorithm*, but it reads the same ADC samples the
detector reads. If the AD8232 had dropped a beat through a momentary contact
glitch, the reference and the detector would miss it identically and the score
would still be 100%. Closing that requires a second instrument — pulse
oximeter, manual pulse count, or a second channel — and none was available.

What this phase **does** establish is stated in §8, and the boundary is not
blurred anywhere in this document.

---

## 1. What was already in place

From the AD8232 smoke test (29 July):

- 3-lead AD8232 breakout: output on **G34**, LO- on **G33**, LO+ on **G32**
  (all ADC1 — ADC2 is unusable with WiFi active)
- Clean P-QRS-T morphology confirmed on a 30-beat aligned overlay
- Mains hum reduced from 90.7% to **3.9%** of signal power purely through
  electrode contact improvement
- ~1.1% of samples clipping at the ADC ceiling
- ESP32 battery-powered whenever electrodes contact a person (the AD8232 has
  no patient isolation; a mains-connected USB host is not acceptable)

From the offline pipeline (E1, MIT-BIH):

- Pan-Tompkins implemented from scratch, five stages, zero-phase throughout
- Headline over 46 MLII records: **Se 99.52%, PPV 99.35%**
- Flutter-episode exclusion (record 207) justified empirically by
  `04_flutter_audit.py`, not assumed
- Records 102/104 excluded from the headline on lead-mismatch grounds,
  reported separately

---

## 2. What this phase did

1. Audited the frozen detector for hardcoded sample counts (§3)
2. Audited the capture chain and found three defects (§3)
3. Built a hand-reviewed R-peak reference from device data and froze it (§4, §5)
4. Ran the frozen detector against that reference and scored it (§6, §7)

Scripts, in `hardware/ecg_bringup/`:

| script | role |
|---|---|
| `01_annotate.py` | validation gate, candidate proposal, review plots, freeze + hash |
| `02_detect_device.py` | integrity gate, import frozen code, detect, score, report |

Neither script contains a copy of the detector or the scorer. Both are imported
by path from the pipeline root, resolved by walking up the directory tree.

---

## 3. Defects found before any scoring

### 3.1 Detector audit — PASSED

Every timing constant in `02_qrs_detect.py` is expressed in **seconds** and
converted using the record's own `fs`:

| constant | value | conversion |
|---|---|---|
| `MWI_WIDTH_S` | 0.150 s | `int(round(MWI_WIDTH_S * fs))` |
| `REFRACTORY_S` | 0.200 s | `int(round(REFRACTORY_S * fs))` |
| `REFINE_WIN_S` | 0.075 s | `int(round(REFINE_WIN_S * fs))` |
| `TWAVE_WIN_S` | 0.360 s | `int(round(TWAVE_WIN_S * fs))` |
| derivative kernel | 5-point | scaled by `fs / 8.0` |
| learning phase | 2 s | `int(2 * fs)` |

No raw sample counts anywhere, and the Nyquist clamp is present in
`stage1_bandpass`. This was checked **before** the first run, because a
hardcoded 360 Hz sample count would not crash at 250 Hz — it would degrade
silently, which is the worse failure.

### 3.2 Dropped serial byte at 230400 baud — capture `180817`

Row 6874 of `180817` reads `7510,1776`, sitting between `27506` and `27514`.
The ADC value is correct and fits the local baseline; the timestamp lost its
leading character. **`27510` -> `7510`, a single dropped byte.**

The smoke-test firmware opens at `Serial.begin(230400)`. The project's
established constraint is 115200, with 230400 known to garble on this setup.
This capture predates that discovery and is plausibly the evidence for it.

At 250 Hz and ~12 bytes per line the stream needs about 3 kB/s — roughly a
quarter of 115200's capacity. **230400 buys nothing and should be dropped.**

`01_annotate.py` now refuses any capture whose timestamps are not strictly
increasing, and exits non-zero. Verified against `180817`.

### 3.3 The `-1` lead-off sentinel — a live landmine

```c
int val = lo ? -1 : analogRead(ECG_PIN);
```

Never fired in any of the five captures (checked: zero occurrences, minima
987-1264). But if it fires it destroys a run **silently**:

- `clean_nonfinite()` interpolates NaN and inf only. **`-1` is finite** and
  passes through untouched.
- Baseline sits near 1800, so `-1` is a ~1800-count negative excursion —
  comparable in magnitude to an R wave, opposite in sign.
- The step down and back produces two large slope events, and Pan-Tompkins
  detects on slope. Guaranteed false positives, no warning anywhere.

**Fixes:** in firmware, move lead-off to a **separate third column** and keep
the real ADC reading — no data lost, diagnostic channel gained. In software,
`01_annotate.py` asserts no `-1` rows and refuses to proceed.

### 3.4 Timestamp resolution — a claim that had to be withdrawn

The sketch schedules with `micros()` (absolute, non-drifting — good) but
timestamps with `millis()` (1 ms resolution). Sub-millisecond jitter is
therefore **invisible by construction**.

An earlier note claimed "zero jitter, better than the toco firmware's 9.4 us."
That is not supportable and has been withdrawn. The toco firmware logged
`micros()` and could measure 9.4 us; this one cannot. The defensible claim is
**"no sample interval deviates by 1 ms or more."** Future ECG firmware should
log `micros()` if jitter is to be characterised.

---

## 4. The reference: method, and why it is not circular

Device data has no annotations. Without a reference, running the detector
yields "31 beats, median 61 BPM" — a plausibility check that cannot separate a
detector which found every beat from one that missed three and invented three.

### Proposal method

Candidates are the **centres of contiguous runs of `ADC == 4095`** — the
saturation plateaus. Human review then accepts, rejects or adds.

This is how MIT-BIH itself was built: a crude detector proposed beats,
cardiologists corrected them. The human remains the arbiter; the machine saves
the clicking. (Interactive click-to-mark was not available — matplotlib is
pinned to Agg by the Python 3.14 windowing bug on this setup.)

### Independence from the detector

| | proposal | detector |
|---|---|---|
| method | contiguous `ADC == 4095` runs | 5-15 Hz bandpass -> derivative -> square -> 150 ms integration -> adaptive SPKI/NPKI threshold |
| shared parameters | **none** | |

The proposal exploits ADC saturation, a property of *this* hardware that
MIT-BIH does not have. If Pan-Tompkins were broken, the proposal would be
unaffected. That is the requirement for an independent reference — subject to
the shared-signal caveat in §0.

### The blind spot, and how it was closed

The method cannot see a QRS that fails to saturate. Tested directly: masking
±200 ms around every marked beat, **the tallest surviving sample in the entire
record is 1990**, against a ceiling of 4095. Every unmasked sample sits at
least 51% below the ceiling. No unmarked QRS-sized feature exists anywhere.

This also settles the record edges numerically rather than by eye:

- lead-in 248 ms vs median RR 987 ms — nothing cut off at the start
- tail 808 ms vs shortest observed RR 748 ms — **a beat would have fitted on
  timing**, and is ruled out only because the tail lies inside that same
  1990-maximum region

### Amendment made before the detector ran

The first version used `int(round(mean))` for the plateau centre. For an
even-width plateau the mean lands on `x.5`, and Python's `round()` is
**banker's rounding** — it goes to the nearest *even* integer:

```
round(3000.5) = 3000     round(3001.5) = 3002
```

So the reference position depended on the **parity of the sample index**, an
accident of when the capture started. This affected **10 of 31 plateaus** (all
width 2); five rounded up, five down, net 0.000 ms — by coincidence, not
design.

Irrelevant to Se/PPV (tolerance is 37.5 samples). But the pre-registered
prediction concerned the signed timing error at about ±1 sample, and
parity-dependent noise of ±0.5 sample is the same order as the effect being
measured.

**Resolution.** The frozen file now carries two columns:

- `sample` = `floor(centre)` — deterministic, stated, fed to `match_beats()`
- `sample_exact` = the true float centre — used for timing error only

`03_qrs_evaluate.py` was **not modified**. Mean offset: 0.16129 samples
(0.645 ms), recorded in the manifest so the correction travels with the data.

---

## 5. Pre-registration

Fixed and committed at `4bdd372` **before** `02_detect_device.py` existed.

| item | value |
|---|---|
| capture | `raw_ecg_capture_20260729_180936.csv`, full 30 s, no exclusions |
| reference | plateau centres, human-reviewed, frozen |
| scoring | ±150 ms tolerance, greedy one-to-one matching, **imported** from `03_qrs_evaluate.py` |
| prediction | signed timing error small, within about ±1 sample (4 ms) |

**The prediction was weakened before the run, and that matters.** An earlier
version claimed a systematic *early* bias, reasoning that `np.argmax` returns
the first index of a flat plateau. That mechanism is weak: `refine_to_r_peak()`
operates on the **bandpassed** signal, and a 3-sample plateau (12 ms) is far
shorter than the 5-15 Hz passband period (67-200 ms), so the flat top is
smoothed into a single rounded peak. The weaker claim is what went on record.

The result (§6) shows the original mechanism was not merely weak but **wrong in
sign** — every error is zero or positive. Had the stronger claim stayed on the
record it would have been falsified. This is the pre-registration doing its job
and should be reported as such.

### Integrity chain

`02_detect_device.py` recomputes SHA-256 of the capture and the annotations and
refuses to score on any mismatch, so "annotated blind, then scored" is
verifiable rather than asserted.

| artifact | sha256 (first 16) |
|---|---|
| capture | `1caf85070d638014` |
| annotations | `db909bfafb92bbd4` |
| `02_qrs_detect.py` | `71eb60e950f64edd` |
| `03_qrs_evaluate.py` | `a768b352f89b7efd` |

The detector and scorer hashes were reproduced on an independent machine,
confirming the imported code is byte-identical to the file that produced the
MIT-BIH headline.

---

## 6. Results

**Capture:** 7500 samples, 250.0 Hz, 30.0 s, 31 reference beats (0 manual).
**Applied tolerance:** 38 samples = 152 ms (`int(round(0.150 * 250))` = 38;
banker's rounding again, this time upward — the JSON records the nominal 150 ms
and should also record the applied 38 samples).

| metric | value |
|---|---|
| TP / FP / FN | **31 / 0 / 0** |
| sensitivity | 100.00% — 95% CI **[88.78%, 100%]** |
| PPV | 100.00% — 95% CI **[88.78%, 100%]** |
| DER | 0.00% |

**Quote the interval, never the point estimate.** A perfect score on 31 beats
is compatible with a true sensitivity near 89%. For scale, a perfect score
gives a lower bound of 96.4% at n=100 and 99.3% at n=500. One beat here is
worth 3.23 percentage points.

### Timing

| statistic | samples | ms |
|---|---|---|
| signed mean (corrected) | +0.290 | **+1.16** |
| signed median | 0.000 | 0.00 |
| MAE | 0.290 | 1.16 |
| worst single beat | +1.0 | +4.0 |
| uncorrected signed mean | +0.452 | +1.81 |

Error distribution across all 31 beats:

| offset | beats |
|---|---|
| 0.0 samples | 17 |
| +0.5 samples | 10 |
| +1.0 samples | 4 |

The ten beats at exactly +0.5 are the even-width plateaus: the detector lands
on a real sample while the true centre sits between two. Correction removes
exactly 0.645 ms, matching the manifest.

**MAE equals the signed mean exactly (0.290 = 0.290).** That is only possible
if no error is negative, and it independently confirms the direction: every
detection lands at or after its reference, never before.

**Prediction: HELD.** |+0.290| ≤ 1 sample.

### Heart rate

Reference median 60.8 BPM (range 55.4-80.2); detected 60.9 BPM (identical
range). The spread is respiratory sinus arrhythmia — rate rises on inhalation,
falls on exhalation — not instability. The 0.1 BPM difference is the 0-1 sample
offset propagating into RR intervals.

---

## 7. Detection margin — why the perfect score is not luck

A count of 31/31 says the detector succeeded. It does not say by how much,
which is what indicates whether it survives a worse recording. From the
integration stage:

| quantity | value |
|---|---|
| weakest QRS blob | 7.22e8 |
| strongest non-beat feature | 5.52e7 |
| **separation ratio** | **13.1x** |
| weakest blob vs adaptive threshold | 3.22x |
| strongest noise vs adaptive threshold | 0.246x |

Better than an order of magnitude of clear air between signal and noise. The
decision was never close.

The stages plot confirms the mains-hum reasoning empirically: the 50 Hz ripple
plainly visible in the raw trace is **completely absent** after the 5-15 Hz
bandpass, as predicted — 50 Hz sits deep in the stopband, and at 250 Hz
(Nyquist 125 Hz) there is no aliasing.

---

## 8. What this establishes, and what it does not

### Established

The frozen Pan-Tompkins detector, **unmodified**, correctly locates R peaks in
a real AD8232/ESP32 recording at 250 Hz with saturating R waves, under good
electrode contact, on one subject at rest, with a 13.1x detection margin and
sub-sample timing agreement.

The sampling-rate transfer (360 Hz -> 250 Hz) works because every constant is
expressed as a duration. R-peak *timing* is what transfers; this is the same
principle that let the MHR detector move from 125 Hz BIDMC to 25 Hz MAX30102
data unchanged, and the same principle that means **SpO2 calibration will not
transfer**, being a photometric ratio rather than a timing measurement.

### Not established

- **Performance under poor contact, motion or lower SNR.** The best of five
  captures was selected after inspection (§0).
- **Agreement with an independent instrument.** The reference shares the
  signal with the detector (§0).
- **Anything about today's hardware setup.** No re-wiring or re-capture; the
  hardware claim rests on the smoke test.
- **Anything about arrhythmia.** Ectopic beats, flutter and baseline wander are
  characterised by the MIT-BIH result, not by 30 seconds of sinus rhythm.
- **Beat classification (E2).** Out of scope — no beat-type labels exist for
  device data, and R-wave clipping degrades the morphology features E2 depends
  on. Detection is slope-based and unaffected; classification is not.

### Honest statement for the writeup

> The frozen QRS detector was run without modification against a 30-second
> AD8232 recording at 250 Hz, scored against 31 hand-reviewed R peaks frozen
> and committed before the detector was run. It detected every beat with no
> false positives (95% CI [88.8%, 100%]) and a mean timing error of +1.16 ms.
> This demonstrates that the pipeline transfers to the prototype hardware. It
> is not a validation of the detector, which rests on the MIT-BIH evaluation,
> and it characterises performance under good electrode contact only.

---

## 9. Requirements

### Closed

- [x] AD8232 device data run against the frozen Pan-Tompkins pipeline
- [x] Reference created and frozen before detection, with a verifiable hash chain
- [x] Detector confirmed sampling-rate agnostic (audit + 250 Hz run)
- [x] Timing agreement quantified and the pre-registered prediction resolved
- [x] Capture-chain defects identified and gated against

### Open

- [ ] **Supplementary captures.** All five captures saturate (0.99-1.36%
      clipping), so the plateau proposal works on every one and annotation is
      cheap. Running `180512` (90.7% mains hum) alongside `180936` (3.9%) would
      convert a single best-case point into a **characterised range across the
      contact quality actually observed**. `180817` is excluded on integrity
      grounds (§3.2) — pre-register that exclusion before running the others.
      **Zero electrode cost. Highest value-per-hour item currently open.**
- [ ] **Independent reference instrument.** A fingertip pulse oximeter gives a
      simultaneous HR reference and would remove the shared-signal caveat. Also
      the highest-value purchase for the MAX30102 phase.
- [ ] **Firmware revision:** `Serial.begin(115200)`; lead-off as a third column
      rather than a `-1` sentinel; `micros()` timestamps if jitter matters.
- [ ] **Fresh purpose-built capture (Phase B).** Now optional rather than
      load-bearing. Six gel pads available — enough for two 3-lead captures.
      **Pads are single-use**: reuse raises contact impedance, which is the
      exact mechanism behind the 90.7% -> 3.9% hum reduction.
- [ ] `tolerance_samples: 38` to be added to the device eval JSON.

---

## 10. Files

```
hardware/ecg_bringup/
├── 01_annotate.py
├── 02_detect_device.py
├── captures/
│   └── raw_ecg_capture_20260729_180936.csv
├── plots/
│   ├── *_annot_overview.png
│   ├── *_annot_panels.png
│   ├── *_device_stages.png
│   └── *_device_scoring.png
├── results/
│   ├── *_annotations_draft.csv
│   ├── *_annotations_frozen.csv        [committed 4bdd372]
│   ├── *_annotations_manifest.json     [committed 4bdd372]
│   └── *_device_eval.json
└── notes/
    └── notes.md
```

Imported, never copied: `02_qrs_detect.py`, `03_qrs_evaluate.py` (pipeline root).

---

## 11. Next

1. **MAX4466 / fetal PCG bring-up** — the schedule risk. Unlike the other three
   channels it can fail at *acquisition*, and no amount of pipeline work
   rescues a signal that was never captured. Start before the multi-rate ISR.
2. **MAX30102** — deferred pending board identification (the 1.8V pull-up trap
   on generic GY-MAX30102 breakouts can prevent I2C detection entirely) and the
   reference oximeter. Both are lead-time items; order now, work in parallel.
3. **Causal filter conversion — resolve on paper before spending a week.**
   `filtfilt` is acausal across a record, but applied to a *completed* window it
   is block processing with one window of latency. If the device reports a
   windowed estimate rather than responding within a beat, the existing frozen
   filters run unchanged and the revalidation disappears. Nothing in a
   reassurance device needs beat-latency response. Decide deliberately; a
   documented design decision is defensible in a viva, an undocumented shortcut
   is not.
4. **Multi-rate ISR** — single high-rate timer with per-channel decimation
   (PCG ~1 kHz, ECG 250 Hz, UC ~4 Hz). Note G36/G39 are the SENSOR_VP/VN pins
   with a documented ADC glitch erratum; LO+/LO- currently occupy ADC1 pins
   (G32/G33) but are read digitally and could be relocated to free them.
