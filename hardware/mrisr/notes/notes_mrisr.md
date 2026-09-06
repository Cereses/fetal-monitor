# Multi-rate ISR — two-channel (PCG + UC)

Phase 5. MAX4466 and FSR402 on a single 500 Hz time base, with per-channel
decimation. Closes the SENSOR_VN erratum question that the PCG phase explicitly
deferred to this point.

Completed 2026-09-06.

---

## 0. Why this phase exists, and what it does not cover

A fetal heart rate deceleration is interpreted by its **phase relative to the
contraction**. Early decelerations — nadir coincident with the contraction peak
— are benign head compression. Late decelerations — nadir *after* the peak —
indicate uteroplacental insufficiency. Same deceleration shape, opposite
clinical meaning, and the only thing separating them is timing.

Measure the two channels on independent clocks and that phase relationship is
unrecoverable. **The shared time base is not an optimisation; it is what makes
the FHR/UC pairing — cardiotocography — mean anything.** Without it the FSR402
channel has no purpose in this system.

**Two channels, not three.** The AD8232 is not in the concurrent demonstration.
That was a scoping decision under time pressure, and the serial arithmetic
below records why it is also a measured constraint rather than an omission. The
ECG channel remains individually validated: 31/31 on device data at Se 100%,
and 99.69% Se over 95,960 annotated beats on NIFECGDB.

**Not covered here:** no scripted contraction protocol was run, so the frozen
`uc_detector` has not been exercised against concurrent device data. The
captures are constant-preload characterisation runs — there are no contractions
in them to find, and running the detector would be meaningless.

---

## 1. Design

One hardware-paced 500 Hz schedule. Every tick reads **both** channels; PCG
emits immediately, UC accumulates and emits every 125th tick.

| Channel | Pin | Divisor | Rate |
|---|---|---|---|
| MAX4466 (PCG) | GPIO39 (`SN`) | 1 | 500 Hz |
| FSR402 (UC) | GPIO35 | 125 | 4 Hz |

`mrisr_capture.ino` derives its schedule exactly as `pcg_capture.ino` did — the
same `next_us += PERIOD_US` accumulation and the same signed-cast wraparound
guard, both previously measured at 500.0000 Hz with 0.6 µs jitter over 225k
samples.

### 1.1 Averaging, not selection — a compatibility requirement

The UC channel emits the **sum of 125 raw samples**; the host divides by 125.0.

This is not a preference. `04_detect.py` records that device UC data is
"float-valued (means of 125 raw samples), so exact-equality flat runs
essentially never occur" — which makes `uc_detector`'s `flat_run_mask`
**inert on device data**, a property the toco phase documented honestly rather
than hiding.

**Taking every 125th sample instead would restore integer values and silently
reactivate that mask.** Same frozen detector, different behaviour, no error
message. Averaging reproduces exactly what the frozen detector was validated
against.

Two incidental benefits. A 125-tap boxcar is a real (if weak) anti-alias filter
with nulls at multiples of 4 Hz. And sending the integer sum rather than a
formatted float avoids on-device rounding, reproducing the toco means
bit-for-bit.

### 1.2 UC timestamp is the window centre

`(t_first + t_last) / 2`. A mean represents the centre of its window; tagging
it with the emitting tick would introduce a systematic **125 ms lag** between
channels. Cross-channel phase is the entire reason this firmware exists, so
that lag would defeat the purpose.

Confirmed in the raw stream: a `U` line at t = 273,126,000 µs against a
neighbouring PCG tick at 273,250,000 — a 124 ms offset, which is the half
window.

### 1.3 Serial budget, and why the third channel was excluded

```
PCG   500 lines/s x 15 bytes  =  7,500 B/s
UC      4 lines/s x 19 bytes  =     76 B/s
                                 ---------
                                  7,576 B/s  =  65.8% of 11,520
```

0.1 points above `pcg_capture.ino`, which ran clean across 225k samples.

Adding ECG at 250 Hz with per-channel tagged lines gives **11,818 B/s, 102.6%
— it does not fit.** A shared-timestamp variable-field format would come to
roughly 8,770 B/s (76.1%), which fits but exceeds the 74% the toco phase
flagged as uncomfortable, and requires the parser to handle four line shapes.

**The two-channel scope is therefore a measured constraint, documented, not a
gap left unexplained.**

Line formats: PCG `t_us,adc` (starts with a digit), UC `U,t_us_centre,sum`.
**No sentinel values anywhere** — the AD8232's `−1` lead-off marker was finite
and passed straight through `clean_nonfinite()` as a huge negative spike.

---

## 2. Wiring

FSR402 returned to the board alongside the MAX4466. The ESP32's 3V3 and GND
feed the breadboard rails; both sensors tap those rails, so no ESP32 pin
needed reassigning.

```
3V3 ──── FSR402 ──┬── GPIO35
                  │
                1 kΩ
                  │
                 GND
```

Fixed resistor on the **ground** side so force up gives reading up, matching
the CTU-CHB UC channel sign convention with no software inversion. Rows 15/16
on the a–e side; row 16 is the divider node. MAX4466 unchanged on GPIO39 with
its OUT at row 5.

Resting UC read **~1510 counts** against the toco phase's ~1085. Same 1 kΩ,
same spring clip, different preload tension after months apart. Well clear of
both the ADC dead zone and saturation.

**Wiring was verified by measurement, not by photograph.** A low-angle photo
could not resolve which holes were occupied, and reading it as confirmation
would have been the eyeballing error this project keeps catching. Instead: the
circuit has no low-impedance failure path (every branch is limited by 1–3 kΩ),
so powering on and reading both channels *is* the verification. PCG at 1909 and
UC at 1510 simultaneously confirm supply, both rails, both pins, divider
orientation, clip preload, and the absence of a row-5 collision.

---

## 3. Integrity and timing

Four 90 s captures, ~180,000 PCG samples.

| | fan on | fan off | paired multi | paired single |
|---|---|---|---|---|
| PCG samples | 45,006 | 45,003 | 45,006 | 45,001 |
| achieved fs | 499.9889 | 499.9778 | 499.9889 | 499.9778 |
| UC fs | 4.0000 | 4.0000 | 4.0000 | — |
| PCG:UC ratio | 125.02 | 125.01 | 125.02 | — |
| dt std (µs) | 9.45 | 18.86 | 9.45 | 13.33 |
| malformed lines | 0 | 1 | 0 | 2 |

Decimation is exact and UC lands on 4.0000 Hz. Largest UC gap 0.250 s, so
`04_detect.py`'s `find_trim_point()` — which trims everything before any gap
exceeding `MIN_GAP_S = 5.0` — trims nothing, as intended.

### 3.1 The dt std figures are single events, and the cause is not the ISR

Each elevated `dt std` is **one interval**, and the arithmetic proves it: one
6000 µs interval among 45,002 gives `4000/√45002 = 18.86` (reported 18.857);
one 4000 µs interval among 45,005 gives `2000/√45005 = 9.43` (reported 9.4466).

I hypothesised that the UC burst every 125th tick — 34 bytes into a 2 ms window
— was blocking `Serial.print` and slipping the tick.

**The paired test refuted it.** The single-channel firmware, which emits no UC
line at all, slipped **twice** against the multi-channel run's once, and had two
malformed lines against zero.

The real cause is visible in what `02_capture.py` printed:

```
'♦♦♦♦♦...♦♦000,1911'
'45958000♦♦♦...♦♦1'
```

Long runs of replacement characters — **serial corruption**, the same family as
the dropped byte that forced 230400 → 115200 in the ECG phase, reappearing at
115200 at ~66% utilisation on this USB/CP2102 combination. The 4000 µs "slip"
is the *consequence*: a corrupted line is unparseable, the parser drops it, and
the surrounding timestamps show a 2 ms hole. **The ESP32 never missed a tick.**

Rate ~1 in 20,000 lines. Host-side, independent of channel count, and
timestamps survive it — so cross-channel alignment is unaffected. Documented,
not fixed.

---

## 4. SENSOR_VN erratum under channel alternation — CLEARED

PCG notes §6 cleared GPIO39 using **single-channel continuous** sampling and
stated explicitly that this did **not** clear it for the multi-rate ISR,
because the erratum is most reported when the ADC alternates channels. This
firmware alternates on every tick.

### 4.1 The first attempt failed — my metric, not the ADC

The pre-registered discriminator was the **p2p/std ratio**, on the reasoning
that a glitch injects isolated outliers and so inflates p2p far more than std.
It returned ERRATUM ACTIVE on every multi-channel capture:

| | measured | §6 baseline |
|---|---|---|
| PCG std | 3.23–3.28 | 3.82–4.11 |
| PCG p2p | 88–93 | 120–131 |
| ratio | 26.9–28.7 | 14.5–14.9 |

**Both absolute inputs went DOWN.** Std below the quietest single-channel run,
p2p below the smallest. The ratio doubled because std fell proportionally
further than p2p — arithmetic, not glitching. Std is set by the bulk of the
distribution; p2p by two samples. A quieter capture shrinks the bulk faster than
the extremes, so the ratio rises.

**The metric has a floor-level dependence I did not account for**, which makes
it uninterpretable whenever the noise floor itself moves.

Definitive proof came from the paired single-channel control: same firmware,
same pin, **no alternation**, and the ratio measured **26.26** against the
14.5–14.9 baseline recorded six days earlier. The ratio tracks the noise floor,
not glitches.

**Fifth instance of the §11 pattern** (see §6 below).

### 4.2 Paired same-session test — the one that decides

The August baselines were taken in an unknown room state; the ISR captures were
taken in a room deliberately quietened. Two sessions, not a comparison — the
same confound that invalidated the chamber comparison in PCG §14.6.

**Pre-registered:** erratum active if multi-channel p2p exceeds
**same-session** single-channel p2p by more than 30%. Absolute p2p, no ratio,
no stale baseline.

**Sequencing chosen to handicap the expected result.** Multi ran first, at
00:58, in the noisier part of the night; single ran at 01:02. Since the
expectation was that the erratum is inactive, giving the multi-channel run the
noisier window makes any "multi is not worse" finding conservative.

| | multi (00:58) | single (01:02) | change |
|---|---|---|---|
| PCG std | 3.277 | 3.466 | −5.5% |
| **PCG p2p** | **88** | **91** | **−3.3%** |
| PCG bias | 1909.1 | 1910.8 | −1.7 |

**Multi is 3.3% LOWER, against a +30% threshold, while carrying the handicap.**

**VERDICT: the SENSOR_VN erratum does not manifest under channel alternation at
500 Hz with two ADC1 channels.** GPIO39 is cleared for the multi-rate ISR. No
move to GPIO36, no settling delay, no restructured read order.

**Limitation:** two channels at 500 Hz, ADC1 only, 11 dB attenuation. Three
channels or a higher rate would need re-testing.

---

## 5. Cross-channel settling — no contamination

The ESP32 has one SAR converter behind a multiplexer. A read too soon after a
channel switch can return a value contaminated by the previous channel. PCG
sits near 1909 and UC near 1510, so contamination would pull each toward the
other and raise both floors.

**Pre-registered:** contaminated if either channel's std exceeds 1.5× its
single-channel baseline.

| | measured | baseline | limit | |
|---|---|---|---|---|
| PCG std | 3.23–3.28 | 3.82–4.11 | 6.17 | ok |
| PCG bias | 1909.1–1910.6 | 1907–1911 | — | ok |
| UC std | 0.68–1.40 | 3.92 | 5.88 | ok |

Every measurement is **below** its single-channel baseline. The bias moved 1.7
counts between paired captures, where contamination from a 1510-count channel
would pull it down by far more.

**`analogRead()` settles adequately between channels at this rate. No
inter-read delay is required, and that is now measured rather than assumed.**

**The UC figure also confirms the averaging.** Averaging N independent samples
reduces noise by √N: 3.92/√125 = 0.35 expected. Measured 0.68–1.40, above that
because consecutive samples are not fully independent, but the direction and
magnitude confirm the boxcar is working. UC drift over 90 s: −0.7 to +1.7
counts on a ~1510 baseline.

---

## 6. Errors in this phase

**(h) The ratio metric had an unaccounted floor dependence.** Registered
p2p/std as the erratum discriminator; it returned a false ERRATUM ACTIVE on
every multi-channel capture while both of its absolute inputs were *below*
baseline. Disproved by the paired single-channel control, which produced the
same elevated ratio with no alternation at all.

**Fifth instance of the pattern in PCG notes §11**: the fixed-frequency mains
window, the CV without a Rayleigh reference, the one-directional occlusion
invalidity condition, the unimplemented half of a pre-registered criterion, and
now this. Each specified something adjacent to the quantity that mattered. The
common failure is **not validating the metric against a control that should
produce a null result** before trusting it on the condition of interest.

**(i) The UC-burst timing hypothesis was wrong.** Proposed that the 34-byte
burst every 125th tick was blocking `Serial.print` and slipping the tick. The
paired single-channel run — no UC line at all — slipped *more*. Withdrawn, and
the real cause identified as serial corruption.

---

## 7. Established / not established

**Established by measurement:**

- Two channels share one 500 Hz time base; decimation exact at 125.02, UC at
  4.0000 Hz (§3).
- SENSOR_VN erratum does not manifest under two-channel alternation at 500 Hz,
  by same-session paired control with adverse sequencing (§4.2).
- No cross-channel settling contamination; both floors below their
  single-channel values (§5).
- 125-sample averaging reduces UC noise 3.92 → 0.68–1.40 counts (§5).
- UC output is byte-compatible with `04_detect.py`; `flat_run_mask` stays inert
  as documented; nothing is trimmed (§1.1, §3).
- Serial corruption ~1 in 20,000 lines, host-side, independent of channel
  count, timestamps preserved (§3.1).
- Three ADC channels do not fit the serial budget in a per-channel line format
  (§1.3).

**Not established:**

- Frozen `uc_detector` against concurrent device data — no contraction protocol
  was run (§0).
- Erratum behaviour with three channels or at rates above 500 Hz (§4.2).
- Any FHR/UC phase measurement. The time base is in place; nothing has yet been
  measured *through* it.

---

## 8. Outstanding

- **Causal filter conversion.** All frozen pipelines use `filtfilt`, which is
  acausal and incompatible with live streaming.
- **Plausibility ceiling** on reported heart rate — PCG §7.1 is the measured
  evidence for it.
- **Signal-quality gate** ahead of the detector — PCG §7.2 is the measured
  evidence for it.
- **Scripted contraction protocol** through the ISR, if time allows, to
  exercise `uc_detector` on concurrent data.
- Resolve the PCG §14.5 rate discrepancy (87.0 BPM vs 70–71 counted).
- Final writeup and demonstration.
