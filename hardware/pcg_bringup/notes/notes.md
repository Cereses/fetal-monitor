# MAX4466 PCG hardware bring-up

Phase 4 of the multi-sensor home fetal monitor. Fetal phonocardiography channel:
wiring, firmware, capture, and a frozen-pipeline run against external ground
truth.

Completed 2026-08-31 / 2026-09-01. Commit `98efe36`.

---

## 0. Scope, and what this phase does not establish

The headline result is that the **acquisition chain works end to end**: capsule
→ MAX4466 → ESP32 ADC → serial → capture script → frozen pipeline recovers a
known fetal heart rate from a played recording, matching the file-only
reference to the detector's own resolution.

Four limits, stated first because they bound every number below.

 **Airborne playback is the easy case, but surface conduction was subsequently tested.** 
 Run C plays a recording at a bare electret capsule, which an electret is designed to transduce. 
 It says nothing about fetal heart sounds arriving as surface vibration through tissue. 
 §14 addresses that separately and establishes that surface-conducted heart sound does reach the
 capsule — contact-dependent, with S1/S2 doublet structure and a
 physiologically correct systolic fraction — on an **adult**, with the
 absolute rate unresolved (§14.5) and the chamber-versus-bare comparison
 invalidated by a hardware fault (§14.6).

**The stimulus is the best case.** `fetal_PCG_p21_GW_39` is the highest-scoring
record in fpcgdb (0.7631 confidence, 97.41% of windows reliable over the full
20 minutes). A recorded signal also carries none of the maternal heart sounds,
bowel sounds, movement artefact or coupling loss a live acquisition faces.

**No on-body validation, and none planned.** No pregnant participant, and there
will not be one without institutional ethics approval, informed consent and
clinical supervision. That protocol is documented as designed-but-not-executed.
The 3D-printed belt enclosure was descoped on cost.

**The frozen FHR pipeline has no per-beat ground truth.** Neither simfpcgdb nor
fpcgdb carries beat annotations, so the pipeline has no Se, no PPV, no MAE — 
unlike the ECG channel (99.69% Se over 95,960 annotated beats) or the UC channel
(F1 77.32% against Romagnoli annotations). Its reliable-record rate on real data
is 58%. Every comparison in this phase is against the pipeline's own output on
file data, not against truth.

---

## 1. Hardware configuration

| Item | Value |
|---|---|
| Board | GY-MAX4466 clone, silkscreen revision **V440**, purple PCB |
| Supply | 3V3 rail (**not** VIN/5V — see below) |
| Signal pin | **GPIO39**, silkscreened **`SN`** |
| ADC | 12-bit, 11 dB attenuation |
| Sample rate | 500 Hz |
| Baud | 115200 |
| Gain trimpot | factory position, never moved |
| Speaker volume | 40% |
| Mic position | a few cm from the subwoofer port, not resting on the cabinet |

**Finding: GPIO36 and GPIO39 are silkscreened `SP` and `SN` on this devkit.**
They are the ESP32's dedicated sensor/pre-amp inputs (SENSOR_VP / SENSOR_VN) and
are labelled by function rather than GPIO number. Working from a pin table alone
leads to the conclusion that the board has no G36 — which is what happened here.
`SP`/`SN` sit between `G34` and `EN`/`3V3`; a miscount toward `EN` holds the chip
in reset, and toward `3V3` connects the op-amp output to a stiff rail.

**Supply constraint.** The MAX4466 accepts 2.4–5.5 V and biases its output at
VCC/2. On 5 V that bias is 2.5 V with swing toward 5 V, exceeding the ESP32
ADC's absolute maximum. The board would work and quietly damage the pin.

**Pin allocation after this phase:**

| Pin | Silkscreen | Channel |
|---|---|---|
| GPIO34 | `G34` | AD8232 output |
| GPIO33 | `G33` | AD8232 LO− (digital) |
| GPIO32 | `G32` | AD8232 LO+ (digital) |
| GPIO35 | `G35` | FSR402 |
| GPIO39 | `SN` | **MAX4466** |
| GPIO36 | `SP` | free |

No passive components. Unlike the FSR402, which needed a fixed 10 kΩ to form a
divider, the MAX4466 presents a buffered op-amp output already biased at VCC/2
and drives the ADC pin directly.

**Assembly note.** The header ships unsoldered. The pin/PCB junction is a loose
friction fit until soldered, and an intermittent contact there produces baseline
*steps* that pass every validation gate — finite values, monotonic timestamps,
no sentinel — while looking like signal. It also has the same signature as the
SENSOR_VP/VN erratum, which would have made §6 unattributable. Header soldered
before any characterisation, pins protruding from the trimpot face so the
capsule face stays flat.

---

## 2. Firmware

`pcg_capture/pcg_capture.ino`, derived from the measured-good
`toco_capture.ino` with **three changes only**: pin 35 → 39, header string, and
settling delay 200 ms → 2000 ms. The timing loop, the signed-cast wraparound
guard, the ADC configuration and the line format are untouched. Rewriting would
have discarded the toco phase's measured 500.000 Hz / 9.4 µs evidence.

`pcg_smoke/pcg_smoke.ino` is a first-power diagnostic that aggregates per-second
statistics on-device. **No claim in this document derives from it.**

**Why 500 Hz.** Not a compromise. The frozen pipeline's bandpass evaluates its
upper corner as `min(200, 0.95 × Nyquist)`; at 500 Hz that is 25–200 Hz — the
*identical* passband used on 1000 Hz simfpcgdb data, and wider than the
25–158.2 Hz used on 333 Hz fpcgdb. 1 kHz would also not fit the serial link:
15 bytes/line × 1000 = 15,000 B/s against 11,520 B/s available at 115200 8N1.
At 500 Hz it is 7,500 B/s, ~65% utilisation.

**Lessons carried from previous phases.** `micros()` not `millis()`, because the
AD8232 firmware timestamped at 1 ms resolution while scheduling on microseconds,
making sub-millisecond jitter invisible by construction and forcing a zero-jitter
claim to be withdrawn. No in-band sentinel, because the AD8232's `−1` lead-off
marker is finite and passed straight through `clean_nonfinite()` as a huge
negative spike; a microphone has no lead-off condition, so that mode cannot
recur here.

---

## 3. Timing and integrity

Four 90 s characterisation captures plus the 90 s run C, ~225,000 samples total.

| | base 1 | base 2 | base 3 | run C |
|---|---|---|---|---|
| samples | 45,005 | 45,005 | 45,004 | 45,001 |
| dt mean (µs) | 2000.000 | 2000.000 | 2000.000 | 2000.000 |
| dt std (µs) | 0.61 | 0.62 | **0.14** | — |
| dt range (µs) | 1998–2001 | 1998–2001 | 1999–2001 | — |
| gaps >2× | 0 | 0 | 0 | 0 |
| malformed lines | 0 | 0 | 0 | 0 |

Jitter of 0.6 µs is 0.03% of the 2000 µs period. Run 3's last 80 s recorded
40,003 consecutive intervals of **exactly** 2000 µs (std 0.0) — the busy-wait
phase-locking to the microsecond counter. Not something to rely on, but it
confirms the ±1 µs elsewhere is quantisation rather than scheduling slip.

Zero gaps and zero malformed lines across ~225k samples settles the serial
budget: 65% utilisation at 115200 is comfortable, and the 230400-baud dropped-byte
failure from the ECG phase does not recur.

**Finding: the achieved-rate measurement is self-referential.** `fs_achieved`
came out at 500.0000111 Hz, a rate error of 2×10⁻⁶ %. That number is too good
because the ESP32 timestamps with `micros()` and schedules from `micros()` — the
chip timing its own scheduler. It establishes that **no scheduled slot was
missed or slipped**, which is the thing that would have broken. It does **not**
verify crystal accuracy, which needs an external reference. At a typical ±50 ppm
that is ±0.025 Hz, or ±0.007 BPM at 140 — negligible, so no external reference
is needed. The notes must say "no scheduler slip", not "sample rate verified".

**Settling: measured, not assumed.** Rather than discarding a lead-in on faith,
full-run and last-80-second statistics were compared. They agree to under 1% in
every run. No settling exclusion is needed, and that is now evidence.

---

## 4. Noise floor

Three quiet-room baselines, factory trimpot:

| | base 1 | base 2 | base 3 |
|---|---|---|---|
| DC bias (counts) | 1911.05 | 1910.16 | 1911.21 |
| std (counts) | 3.82 | 4.11 | 4.00 |
| p2p (counts) | 120 | 131 | 122 |

Bias of ~1910 at 11 dB attenuation is where a VCC/2 bias of 1.65 V should land,
and it reproduced across power cycles, re-uploads, pin changes and physical
handling.

**The channel is ADC-noise-limited, not amplifier-noise-limited.** The toco
phase measured the bare ADC floor at 3.92 counts std at the same rate and
attenuation. With the MAX4466 connected and running the floor is 3.82–4.11. The
amplifier contributes essentially nothing above the converter, which is the
measured justification for raising gain: extra gain lifts signal without lifting
total noise proportionally, so SNR improves nearly one-for-one until amplifier
noise finally appears.

---

## 5. The 47.8 Hz line

Every baseline capture shows exactly one narrow line clearing 6 dB over local
background in the 25–200 Hz passband:

| | base 1 | base 2 | base 3 |
|---|---|---|---|
| frequency (Hz) | 47.80 | 47.80 | 47.80 |
| strength | +12.1 dB | +8.6 dB | +12.4 dB |
| share of in-band power | 4.07% | 2.58% | 4.43% |

The strength row is `quick_check_lines.py` stdout for the 47.80 Hz line and is
not persisted. The `db_over_bg` values in `results/*.json` for these same three
captures are much lower — `spectrum_full.mains[].db_over_bg` at 50/100/150 Hz
reads +1.55/+0.62/+1.76 (base 1, `pcg_char_baseline_20260831_183354.json`),
+0.11/+1.03/+0.74 (base 2, `pcg_char_baseline_20260831_183912.json`) and
+1.07/+1.60/+1.92 (base 3, `pcg_char_baseline_20260831_185145.json`) — because
those are measured at 50/100/150 Hz and 47.80 Hz falls outside every mains
window (§11.1). Different quantities, not a discrepancy.

Base 1/2/3 are mapped to these files by matching `signal_full.mean` against
§4's DC-bias row (exact to 2 dp for all three, and consistent with capture
timestamp order). No results file records the base labels, so the mapping is
inferred, not recorded. The nine dB values are read directly from the files and
hold regardless of the mapping.

**47.80 ± 0.00 Hz** across three runs spanning eighteen minutes — a persistent
source, not a transient.

**It is not mains.** Mains was later injected directly and unambiguously (§7)
and measures **50.40 Hz** with a harmonic ladder. 47.80 and 50.40 are 2.6 Hz
apart, thirteen resolution bins at 0.2 Hz. Different sources.

The source was not identified. Candidates include a motor or fan near 2870 rpm,
a mechanical resonance, or an alias from above 250 Hz (there is no analog
anti-aliasing filter). Not pursued further because at 12 dB it does not fool the
detector (§7) and is not on the critical path.

**No anti-aliasing filter.** The MAX4466 passes audio to ~20 kHz and everything
above 250 Hz folds into the sampled band. An RC low-pass was considered and
deliberately deferred so the baseline could be measured without it; it was never
added. This is an open limitation, not a solved problem.

---

## 6. SENSOR_VN erratum test

GPIO36/39 carry a documented ADC glitch erratum; GPIO34 does not. The mic was
brought up on G34 first so that anything anomalous would be attributable to the
sensor, the gain or the wiring rather than to a suspect converter, then moved.

**Metric and prediction registered before the measurement.** A glitch injects
isolated outliers, which inflates peak-to-peak far more than standard deviation.
The discriminator is therefore the p2p/std *ratio*, not either number alone.

| | G34 (control) | **G39 (SN)** |
|---|---|---|
| mean | 1907.3 | 1907.2 |
| mean std | 3.92 | 3.43 |
| mean p2p | 57.0 | 51.0 |
| **mean p2p/std** | **14.5** | **14.9** |
| quietest second, std | 1.70 | 1.24 |
| quietest second, p2p | 20 | 12 |

**The erratum did not manifest.** Both p2p and the ratio are *lower* on G39, and
the ratio moved by under 3% — inside the spread between two G34 runs. The
quietest-second figures, least contaminated by ambient sound and therefore
closest to pure converter behaviour, also favour G39.

Supporting structure: elevated-p2p seconds arrive in runs of two or three
*consecutive* seconds on both pins. Sample-level glitches would scatter randomly;
second-scale clusters are room events.

**Limitation, recorded before the result was seen.** This used single-channel
continuous sampling. The erratum is most often reported when the ADC alternates
between channels, which nothing here exercised. **This does not clear GPIO39 for
the multi-rate ISR.** It establishes a single-channel baseline against which the
multi-channel case must be tested separately.

---

## 7. Detector behaviour on non-cardiac input

Two findings, both accidental, both more important than the 47.8 Hz line.

### 7.1 On pure noise: correctly LOW, but the BPM is always plausible

The frozen pipeline was run on all three quiet-room baselines. Prediction
registered first: 90 s of quiet room contains no heart sound, so confidence
should sit far below the 0.45 threshold and the verdict should be LOW.

| | base 1 | base 2 | base 3 |
|---|---|---|---|
| windows | 87 | 87 | 87 |
| mean confidence | 0.064 | 0.043 | 0.066 |
| windows ≥ 0.45 | 0.0% | 0.0% | 0.0% |
| verdict | LOW | LOW | LOW |
| BPM mean / std | 138.2 / 24.0 | 140.2 / 23.7 | 150.5 / 22.4 |

Held in all three. Mean confidence around 0.05 against a 0.45 threshold — an
order of magnitude below. **An OK verdict later therefore carries meaning.**

**But the BPM values are the finding.** With no periodicity, `argmax` selects
essentially arbitrarily within the lag search range (166–271 samples at 500 Hz),
so BPM = 30000/lag over 110.7–180.7. For uniformly selected lags that
distribution has mean **140.1 BPM** and std **20.0**. Observed means were 138.2,
140.2 and 150.5, with stds 22.4–24.0.

**Pure noise produces a BPM clustered around 140 — dead centre of the normal
fetal range.** A number from this device will look physiologically plausible
whether or not there is a fetus. This is direct measured evidence that a BPM must
never be displayed without its confidence flag, and it feeds the plausibility
ceiling still outstanding.

Compounding it: `03_evaluate.py` falls back to the median of *all* windows when
no window clears threshold. It still returns a number, flagged LOW.

### 7.2 On a strong pure tone: confident, reliable, and wrong

During an occlusion test (§11.3) a fingertip on the capsule injected mains hum
capacitively. The frozen detector on that capture:

```
mean confidence  : 0.809  (max 0.863)
windows >= 0.45  : 97.7%
BPM mean/std     : 175.7 / 4.1
VERDICT          : OK
```

**On a pure 50 Hz sinusoid.** Not a marginal pass. Higher confidence than
`fetal_PCG_p21_GW_39` scores over its full record (0.7631), with a BPM stable to
±4.1 across 87 windows.

**Hypothesised mechanism, not verified.** `shannon_energy_envelope` squares the
signal and applies a 50 ms moving average (25 samples at 500 Hz). Squaring moves
a 50 Hz carrier to 100 Hz, where the smoother's response is roughly −27 dB, so a
steady tone is heavily suppressed — but *relatively*, not absolutely. The tone
here sat ~47 dB above the noise floor, so the residual still dominated the
envelope. That would explain why the 47.8 Hz line at 12 dB is harmless while
50.4 Hz at 57 dB is catastrophic: **the failure depends on interference
amplitude, not merely its presence.** This should be verified before being
stated as mechanism.

**Consequence.** Combined with §7.1, the device can produce a physiologically
plausible number at high confidence from a signal containing no heart.
**Confidence alone is not a sufficient gate.** The measured separation is stark
— dominant-line share of in-band power was 97.91% for interference against
2.58–4.43% for the quiet floor — which suggests a signal-quality gate ahead of
the detector. That would be new code, never a modification of frozen code.

(2.58–4.43% is `quick_check_lines.py` stdout from the §5 baselines and appears
in no results file — the same provenance gap §11.5 describes.)

The 97.91% is `results/pcg_char_occluded_20260831_191155.json`,
`spectrum_full`, 50 Hz `pct_of_inband`. The strongest in-band line in that
capture was at 50.40 Hz, inside the ±1 Hz window around 50, so the
fixed-frequency and strongest-line figures coincide there.

**Both ends of this separation are non-cardiac.** The share was measured on an
occluded capture and on the quiet floor, and never on a capture containing
fetal heart sounds. An earlier version of this paragraph cited 4.07% as "real
fetal PCG": 4.07% is baseline 1 of the quiet-floor row in §5, a 47.80 Hz line
on a capture with no acoustic stimulus. See §11.5.

### 7.3 Body contact injects mains at 235× the noise floor

The occlusion capture measured std **939.01** counts against a 3.82–4.11 baseline,
with 50.40 Hz at +57.1 dB carrying **97.91%** of in-band power
(`results/pcg_char_occluded_20260831_191155.json`, `spectrum_full`, 50 Hz
`pct_of_inband`) and harmonics at 100.60, 151.00 and 198.20 Hz. The strongest
in-band line in that capture was at 50.40 Hz, inside the ±1 Hz window around
50, so the fixed-frequency and strongest-line figures coincide there.

The body is capacitively coupled to mains wiring and acts as an antenna; a
fingertip on a high-impedance electret input couples that in directly. Recorded
as a **characterised bench interference mode**. It is not a deployment
requirement — there is no belt and no on-body validation in this project — but it
corroborates the existing battery-power protocol documented for the on-body
procedure.

---

## 8. Playback chain characterisation

Tone bursts at 25/50/75/100/150/200 Hz, 3 s each with 1 s gaps, played through a
2.1 speaker system with the capsule at the subwoofer. Six bursts detected at
3.98 ± 0.04 s spacing against a manifest 4.00 s; alignment independently
confirmed by the 10 s settling split, which removed the 50 Hz burst's energy
entirely (3.88% → 0.0011%) while leaving 100 and 150 Hz intact.

| Hz | rel dB | harmonic % | ramp dB |
|---|---|---|---|
| 25 | −25.7 | **39.3** | +0.6 |
| 50 | −11.3 | 1.3 | −0.1 |
| 75 | −2.0 | 0.2 | **+2.0** |
| 100 | 0.0 | 0.2 | +0.1 |
| 150 | −8.9 | 0.0 | +0.0 |
| 200 | −18.5 | 0.0 | −0.3 |

**The 25 Hz burst is 39.3% harmonic content**, with 75 Hz only 6.6 dB below the
fundamental. The driver is not reproducing 25 Hz, it is distorting. Its true
response there is worse than −25.7 dB because part of what the capsule receives
is manufactured harmonics. Measuring narrowband at each fundamental and reporting
harmonics separately is what caught this; a broadband measurement would have
scored 25 Hz as present.

The 75 Hz ramp is real at +2.0 dB, an order of magnitude above every other burst.
Cause not established — room mode, cabinet resonance, or slow amplifier dynamics.

**Weighted transmission of the reference segment: 7.3% (−11.4 dB).** Band shares
invert:

| Band | Share of original | Share of what survives |
|---|---|---|
| 25–50 Hz | 80.7% | 24.3% |
| 50–100 Hz | 18.8% | 73.5% |
| 100–158 Hz | 0.2% | 1.9% |

What reaches the capsule is not p21 — it is p21 heard through a bandpass centred
near 100 Hz.

**This was explicitly not treated as a gate on run C**, and the reasoning was
recorded before the measurement: `estimate_fhr_autocorr` runs on the Shannon
energy envelope and detects the *timing* of S1/S2 bursts, not their spectral
content. A transient bandpassed to 50–100 Hz is still a transient at the same
instant.

---

## 9. The reference segment

`fetal_PCG_p21_GW_39`, seconds 60–120, 19,980 samples at 333 Hz. Zero non-finite
samples, so `clean_nonfinite()` never fires and run A carries no interpolated
samples that run C lacks.

**Quantisation.** `wfdb`'s `p_signal` is scaled by the header gain: values are
quantised in steps of **200**, so the true raw ADC span is **−107 to +97, 204
counts, about 7.7 bits**. The reference data on which the entire loopback rests
is itself coarsely quantised, which means the ESP32's 12-bit ADC is comfortably
not the bottleneck. (A step of 100 was predicted from two round endpoints; the
direction was right and the magnitude was a guess. The measured value is what
stands.)

**Band energy**, against the pipeline's own 25–158.2 Hz passband at 333 Hz:

| Band | Share |
|---|---|
| 0–25 Hz | 20.7% (discarded by the highpass) |
| **25–50 Hz** | **64.1%** |
| 50–100 Hz | 14.7% |
| 100–158.2 Hz | 0.2% |

In-band centroid 42.0 Hz, median 38.5 Hz, and **99.8% of in-band energy below
100 Hz**. Crest factor 6.27 (15.9 dB). This is the physical basis for the
acoustic chamber argument: fetal PCG lives at the bottom of the audio band, where
airborne coupling from a vibrating surface is worst.

---

## 10. Runs A, B, C — acquisition validation

### 10.1 Scoring rule validated before use

`03_evaluate.py`'s record-level rule (`estimated = median(bpms[high_mask])`;
`reliable = high_frac >= 0.5`) could not be called directly because
`evaluate_record()` reads through `wfdb.rdrecord`. It was reimplemented, and the
reimplementation was checked against the row `03_evaluate.py` had already written
into `results/fpcgdb_eval.csv`:

| field | reimplemented | stored | match |
|---|---|---|---|
| estimated_bpm | 136.8493 | 136.8493 | yes |
| mean_confidence | 0.7631 | 0.7631 | yes |
| high_conf_fraction | 0.9741 | 0.9741 | yes |
| reliable | True | True | yes |

Exact to four decimals. `04_loopback.py` then imports that same validated
function rather than holding a second copy.

### 10.2 Results

Predictions registered before each run.

| | **A** (file, 333 Hz) | **B** (resampled 500 Hz) | **C** (device) | C tail |
|---|---|---|---|---|
| windows | 56 | 56 | 56 | 56 |
| reliable | 100% | 100% | **100%** | 100% |
| BPM | 134.09 | 133.93 | **133.93** | 134.23 |
| confidence | 0.839 | 0.840 | **0.808** | 0.820 |
| BPM std | 2.42 | 2.46 | 2.18 | 2.20 |
| verdict | OK | OK | **OK** | OK |

**A → B: resampling costs essentially nothing.** 0.17 BPM and 0.002 confidence.
Predicted ≤1.0 BPM on the grounds that lag quantisation *improves* with rate
(0.98 BPM/step at 333 Hz, 0.65 at 500 Hz). Run B's passband is 25–200 Hz against
run A's 25–158.2 Hz, and the answers are identical to 0.17 BPM — nothing the
pipeline needs lives above 158 Hz, consistent with §9.

**B → C: all three predictions held.** Verdict OK, BPM difference 0.00 against a
3.0 tolerance, confidence lower as expected (0.808 vs 0.840).

**On the 0.00 BPM difference.** `peak_lag` is an integer, so BPM is quantised.
30000/224 = 133.92857…, which both runs report: **they landed on the same median
lag of 224 samples**, adjacent lags being 134.53 and 133.33. The tail's
134.2288597 is not 30000/n for any integer — with 56 windows the median averages
two values, and (30000/223 + 30000/224)/2 gives exactly that. The honest
statement is **"identical median lag"**, not agreement to fourteen decimals.

**Level.** Volume 40%, trimpot factory, first attempt GOOD: peak excursion 575
counts (target 400–900), RMS 66.4, SNR 24.4 dB, no flat-topping. Run C measured
26.8 dB — the 2.4 dB rise is content, not drift: the capture covers a louder
passage of the record, 7 samples sat within 90% of peak, and peak 642 is nowhere
near either rail. Crest factor fell 8.66 → 7.36 over the same interval, which in
isolation is the compression signature and is ruled out by those two checks.

Crest factor *above* the source value of 6.27 is itself consistent with §8: the
playback chain acts as a high-pass, removing steady low-frequency energy while
leaving transients intact.

**The seam did not matter.** Primary and tail windows agree to 0.30 BPM, both OK.

**The result in one line:** the acquisition chain recovers the reference heart
rate to the detector's own resolution, at 96% of the reference confidence,
**despite 92.7% of in-band energy being lost in playback**. That is direct
evidence for the claim registered beforehand: the pipeline needs burst timing,
which survives filtering, not spectral fidelity, which does not.

---

## 11. Errors made in this phase

All five are analysis or process errors, recorded rather than quietly fixed.
Three share one shape.

**11.1 Fixed-frequency mains window.** The mains test integrated 50/100/150 Hz
± 1 Hz and returned "not detected" (< 2.1 dB) on all nine checks, while a
+12 dB line sat at 47.8 Hz — eleven bins outside the window. The test was not
wrong about 50 Hz; it was **blind by construction**, because it pre-registered
*where* to look instead of *what* to look for. Fixed with a general narrow-line
detector, not a widened mains window; widening would have been tuning the test
until it found a thing already known to be there.

**11.2 Coefficient of variation without a reference.** The line's Hilbert-envelope
CV was framed as "low = steady tone, high = modulated". The envelope of
narrowband Gaussian noise is Rayleigh-distributed with a **fixed CV of
√(4/π − 1) = 0.523** regardless of amplitude. Observed values were 0.572, 0.802,
0.548 — two sitting almost exactly on the noise value. The metric cannot separate
a modulated tone from narrowband noise without that reference, and was presented
as though it could.

**11.3 One-directional invalidity condition.** The occlusion test pre-registered
"if the broadband floor *collapses*, the comparison is invalid". The floor
**exploded** instead, 235×. The condition should have read "changes materially in
either direction". The test was invalid, but produced §7.2 and §7.3, which are
worth more than the test would have been.

**Pattern.** All three specified one direction of failure and were blind to the
other. Recorded as a pattern in test design, not three unrelated incidents.

**11.4 Scope drift.** The §7.3 result was initially written up as a *deployment
requirement for belt operation*, despite the belt having been descoped in this
phase's own opening brief. The measurement stands; the framing was wrong and was
corrected to a characterised bench interference mode. Worse than a retrieval
failure, since the constraint was in context throughout.

**11.5 A figure that outlived its provenance.** The dominant-line share quoted
in §7.2 and §13 was produced by `quick_check_lines.py`, which prints
`pct_inband` to stdout and **persists nothing** — its only file output is a
plot. With no results file to check against, the number survived in the notes
past the run that made it and was later re-attributed: 4.07%, baseline 1 of the
§5 quiet-floor row, was written up as the share for "real fetal PCG", a capture
class never measured for this quantity at all. The 97.94% quoted alongside it
was not a `spectrum_full` value either; the persisted figure is 97.91%
(`results/pcg_char_occluded_20260831_191155.json`, 50 Hz `pct_of_inband`).

**This is not the 11.1–11.3 pattern.** Those are metrics that specified
something *adjacent* to what mattered — a fixed window beside the real line, a
CV without its noise reference, a collapse condition beside an explosion. The
metric was pointed slightly wrong. Here the metric is **correctly designed and
measures exactly the right thing**; what failed is that its output was never
written to disk, so the number lost its link to the capture and the run that
produced it and drifted onto a different claim. Distinct failure, distinct fix:
persist what gets cited. Found by provenance audit against `results/` — which
only works for quantities that are in `results/`.

---

## 12. Established / not established

**Established by measurement:**

- The acquisition chain reproduces a known FHR end to end (§10).
- Firmware holds 500 Hz with 0.6 µs jitter, zero gaps, zero malformed lines over
  ~225k samples (§3).
- Noise floor 3.82–4.11 counts, at or below the bare-ADC floor; the channel is
  converter-limited, not amplifier-limited (§4).
- DC bias 1907–1911 counts, stable across power cycles and handling (§4).
- No settling exclusion is needed (§3).
- The SENSOR_VN erratum does not manifest under single-channel continuous
  sampling (§6).
- The frozen detector returns LOW at ~0.05 confidence on quiet-room noise (§7.1).
- The frozen detector returns OK at 0.809 confidence on a pure 50 Hz tone (§7.2).
- On noise, reported BPM clusters near 140 — always physiologically plausible
  (§7.1).
- Body contact near the capsule injects mains at 235× the noise floor (§7.3).
- The 47.8 Hz line is persistent, is not mains, and does not fool the detector
  (§5).
- Resampling 333 → 500 Hz costs 0.17 BPM and 0.002 confidence (§10.2).
- The reference segment is quantised to ~7.7 usable bits; the ESP32 ADC is not
  the bottleneck (§9).
- Surface-conducted heart sound reaches the capsule: contact-dependent, air
  control at exactly the noise floor, alternating S1/S2 doublet passing all four
  pre-registered criteria, systolic fraction 0.365 (§14.3, §14.4).
- Body contact plus USB power drives the channel to 95–97% mains; unplugging the
  charger collapses it by ~40 dB (§14.2).

**Not established:**

- **How much of the ~30 dB tissue/air impedance loss the chamber recovers.**
  Surface conduction is established (§14.3, §14.4), but no quantitative
  transmission figure was measured — and the comparison that would have
  produced one was invalidated (§14.6).
- Absolute sample-rate accuracy. Only "no scheduler slip" (§3).
- The mechanism behind the pure-tone false positive — hypothesised, not verified
  (§7.2).
- Erratum behaviour under multi-channel alternating sampling (§6).
- The identity of the 47.8 Hz source (§5).
- Aliasing behaviour: there is no analog anti-alias filter and none was added
  (§5).
- Performance on anything but the best record in fpcgdb, played back cleanly
  (§0).
- The **absolute rate** of the contact signal: 87.0 BPM by three converging
  routes against two careful 60-second counts of 71 and 70 (§14.5).
- **Chamber versus bare capsule.** Attempted, invalidated by intermittent
  connection reaching both ADC rails (§14.6).

---

## 13. Outstanding

- **Signal-quality gate** ahead of the detector, motivated by §7.2. New code.
  **The measured separation cannot supply the threshold.** Both ends of it are
  non-cardiac — an occluded capture and the quiet floor — and no capture
  containing fetal heart sounds has ever been measured for this quantity
  (§11.5). A gate needs a pass-side measurement first, and the quantity has to
  be persisted before it can be one.
- **Plausibility ceiling** on reported heart rate — §7.1 is measured evidence for
  it.
- ~~**Multi-rate ISR.**~~ **DONE** (2026-09-06), two-channel: PCG 500 Hz and
  UC 4 Hz from one timer, so FHR and contraction timing share a time base. The
  erratum was re-tested under channel alternation and **cleared** by a paired
  same-session control, closing §6's deferred item. See
  `hardware/mrisr/notes/notes.md`. The ECG channel was excluded on a measured
  serial-budget constraint, not an omission.
- **Causal filter conversion.** All frozen pipelines use `filtfilt`, which is
  acausal and incompatible with live streaming.
- **RC anti-alias filter**, deferred and never added.
- **Resolve the §14.5 rate discrepancy**, ideally by simultaneous AD8232
  capture using the validated QRS detector as the measuring instrument.
- **Re-run the chamber comparison** if time allows, after adding mechanical
  strain relief to the capsule leads.
