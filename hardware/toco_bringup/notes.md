# FSR402 Tocodynamometer: Hardware Bring-Up Notes

Phase: first sensor of the hardware integration stage, following completion of
all five offline software validation pipelines.
Dates: 2026-08-07 to 2026-08-08.

---

## 0. Scope statement

**The belt enclosure is out of scope for this project.** Fabrication cost placed
a 3D-printed wearable housing beyond the project's means. The system was
designed throughout with belt mounting as the intended deployment form, and the
requirements that form imposes are recorded in section 10 for any future work.

This is a functional descope, not only a packaging one, and the distinction is
recorded deliberately rather than glossed. The FSR402 requires a defined static
preload to operate at all: unloaded, it reads exactly zero and sits inside the
ESP32 ADC's dead zone (section 4). Belt tension was the intended source of that
preload. **In its absence, preload was supplied by a spring-clip bench fixture**,
and every measurement in these notes was taken under that fixture. Where a
figure depends on the preload source, this is stated.

The device was therefore demonstrated, not deployed. No measurement here was
taken on a human abdomen, and none is claimed.

---

## 1. Circuit

FSR402 in a voltage divider against a fixed pull-down resistor, read by the
ESP32 ADC.

```
3V3 ──── FSR402 ──┬── GPIO35 (ADC1_CH7)
                  │
              R_fixed
                  │
                 GND
```

**Pin choice.** GPIO35 is on ADC1. ADC2 is unusable whenever WiFi is active, so
all analog sensors must sit on ADC1 (GPIO 32-39). GPIO34/33/32 were already
allocated to the AD8232 during its smoke test, making GPIO35 the next free
channel. GPIO34-39 are input-only with no internal pull-up or pull-down, which
is desirable here: the external resistor is the only thing setting the node's
resting voltage, with no internal pull-up in parallel to shift readings.

**Orientation.** Fixed resistor on the ground side, so increasing force lowers
FSR resistance and raises the measured voltage. Pressure up means reading up,
matching the sign convention of the CTU-CHB UC channel. No inversion needed in
software.

**Mounting.** The sensor is on flying leads (female-to-male jumpers onto the
tail tabs) rather than seated directly in the breadboard. Adopted to avoid
bending the flat tabs; it also matches the configuration a belt-mounted sensor
would need. Flex tail taped down for strain relief so tugs land on the tape
rather than the tab junction.

**Breadboard layout.** Whole circuit confined to rows 1-15 of one terminal
block, with 3V3 and GND on the adjacent rail strips. Rationale: split rails, if
present, break at the board midpoint around row 30. Keeping every connection
within ten rows means a split rail cannot affect the circuit whether or not one
exists. No multimeter was available to test for it directly.

---

## 2. Firmware

`toco_capture/toco_capture.ino`, 500 Hz, CSV over serial at 115200 baud,
columns `t_us,adc`.

**Why 500 Hz for characterisation.** Chosen to resolve 50 Hz mains hum and its
harmonics clearly, on the expectation that hum would need to be measured before
choosing a decimation ratio. See section 5: this expectation proved wrong.

**Timestamp is relative to start, not raw `micros()`.** Eight digits rather than
ten. At 500 Hz and 115200 baud the link runs near 65% utilisation; ten-digit
timestamps push it past 74%. When the TX buffer fills, `Serial.print` blocks and
corrupts the sample timing being measured. 115200 was retained rather than
raised because 230400 garbled intermittently during the AD8232 smoke test.

**Wraparound handling.** The pacing loop casts the `micros()` difference to
`int32_t` before comparing. `micros()` wraps every ~71 minutes; comparing raw
unsigned values would stall the loop for most of an hour at that point.

**Measured timing.** 500.000 Hz over 90 s, interval std 9.4 us, one late slot in
45,000 samples.

---

## 3. Resistor selection

| R_fixed | condition | resting counts | implied R_FSR | verdict |
|---|---|---|---|---|
| 10 k | unloaded | 0 | > 1 M | ADC dead zone |
| 10 k | spring clip | ~3,170 | ~2,920 ohm | too little headroom |
| 1 k | spring clip | ~1,085 | ~2,780 ohm | **selected** |

**10 k was correct for finger-press testing** (a light press landed mid-scale at
~1,500 counts, implying ~17 kohm) but wrong once static preload was applied. With
a clip holding the sensor at ~2,900 ohm, 10 k puts the operating point in the
compressed upper region of the divider curve, leaving under 1,000 counts of
headroom and reduced sensitivity.

**Swapped to 1 k on 2026-08-08.** Predicted resting level 1,045 counts,
measured 1,085. Counterintuitively this is *more* sensitive at this operating
point, not less: for a squeeze taking R_FSR from 3,650 to 2,000 ohm, the
excursion is 485 counts at 1 k versus 412 counts at 10 k.

**1 k is the final value for this project.** It was initially recorded as a
bench expedient pending belt-derived preload measurements. With the enclosure
descoped (section 0), there is no belt operating point to select against, so the
value chosen for the bench fixture stands. A future belt implementation must
re-derive it: the correct fixed resistor is a function of the preload the
housing applies, and the method for choosing it is `02_preload.py` plus the
target window in section 4.

---

## 4. ADC dead zone, and why preload is a functional requirement

An unloaded FSR402 reads **exactly 0**, with zero variance across 9,750 samples.

The divider puts the node at roughly 33 mV unloaded (R_FSR > 1 Mohm against
10 k), and the ESP32 ADC reports 0 for anything below approximately 150 mV.
Everything under that is clamped flat.

**Consequences:**

- The noise floor is unmeasurable without preload. The first capture returned a
  baseline standard deviation of 0.000, which is not a low-noise result, it is
  a clipped one.
- A device resting at hard zero would be pathological for the detector. The
  contraction detector uses a rolling baseline with a record-relative threshold
  and needs a resting level that exists and varies slightly. Additionally, the
  pipeline's dropout mask treats long flat runs as dead sensor, so most of such
  a record would be masked out as invalid.
- Therefore **static preload is a functional requirement of the toco channel,
  not a convenience.** This is the finding that makes the enclosure a functional
  dependency rather than packaging, and it is why section 0 names the substitute
  rather than treating the belt as cosmetic.

**Target resting window: 600 to 2,000 counts.** Below 600 the dead zone distorts
readings; above 2,000 there is insufficient headroom for a contraction to rise
into. Any future preload mechanism should be validated against this window.

**Substitute used.** A spring clip (clothes peg) clamped across the sensing pad,
supplying constant opposed force. This is mechanically the correct analogue of a
belt, which compresses the sensor between strap and abdomen. It is *not* the
correct analogue of a weight resting on the sensor: see section 10.

---

## 5. Mains hum: prediction wrong, empirically

The 500 Hz sample rate was chosen on the prediction that 50 Hz hum would alias
into the contraction band if sampled directly at 4 Hz. Measured share of AC
power at 50 Hz, across three independent segments:

| segment | 50 Hz share |
|---|---|
| finger press 2 | 0.03 % |
| finger press 3 | 0.01 % |
| clip static load | 0.74 % |

Hum is absent, not merely small.

**Cause.** Hum couples capacitively, so susceptibility scales with source
impedance. The FSR divider node sits between roughly 800 ohm and 17 kohm, whereas
an AD8232 electrode input looks into megohms. The ECG experience did not
transfer; impedance should have been reasoned about before predicting.

**Where the AC power actually sits** (clip static load, detrended):

| band | share |
|---|---|
| below 0.5 Hz | 4.07 % |
| 0.5-2 Hz | 0.86 % |
| 2-10 Hz | 3.36 % |
| 10-45 Hz | 12.36 % |
| 51-250 Hz | 76.53 % |

The high-frequency bulk is SAR converter noise and averages away. The low-
frequency residual does not.

**Consequence for the multi-rate ISR.** The toco channel does not need 500 Hz.
Roughly 6-11 % of power sits in 1-10 Hz and would fold into band if sampled
directly at 4 Hz, so some oversampling is still warranted, but approximately
50 Hz with a 12-sample box-car is sufficient. This reduces the toco's share of
the timer budget by about 90 %, which matters when the PCG channel needs 1 kHz
on the same ISR.

---

## 6. Noise floor and decimation gain

Measured on the settled portion of the clip static-load capture (t > 20 s,
R_fixed = 10 k):

| quantity | value |
|---|---|
| detrended std at 500 Hz | 9.70 counts |
| detrended std at 4 Hz (125-sample box-car) | 3.64 counts |
| improvement factor | 2.67x |
| theoretical white-noise factor | 11.18x (sqrt 125) |

The shortfall is expected, not a fault: averaging removes the high-frequency
converter noise but the 4 % of power below 0.5 Hz passes through unchanged, and
that residual sets the 4 Hz floor.

With R_fixed = 1 k the baseline std fell to **3.92 counts at 500 Hz**, since
lower source impedance picks up less. Contraction amplitudes of 1,000-1,570
counts give a signal-to-noise ratio around 275:1. **Detection is not
noise-limited; it is drift-limited.**

---

## 7. Creep

Static load, spring clip, 83 s of settled data:

- **+0.615 counts/s**, rising (resistance falls under sustained load)
- **+51 counts, +1.50 %** across the segment
- Approximately 37 counts over a 60 s window, roughly ten times the 4 Hz noise
  floor

Creep is monotonic and slow, whereas a contraction rises and falls. The rolling
baseline is the right mechanism for it. The constraint this imposes is that the
baseline window must be long compared to a contraction and short compared to the
creep timescale.

**Note on scope.** A spring clip applies constant force. A belt's tension would
vary as the abdomen moves. This measurement isolates the sensor's own polymer
relaxation and must not be written up as representative of belt behaviour.

---

## 8. Baseline drift and hysteresis

### 8.1 Initial finding (5-minute capture) and its correction

The 5-minute capture (3 cycles, 40 s rests) showed baseline climbing +57, +97,
then **+410** counts across successive rest phases, and this was initially
recorded as an accumulating hysteresis ratchet: the most serious threat to the
detector on long recordings.

**The 20-minute capture disproved that interpretation.** With 7 cycles and 100 s
rests, drift does not accumulate:

| rest phase | baseline | drift from first |
|---|---|---|
| REST | 1,531.8 | — |
| REST 1 | 1,565.4 | +33.7 |
| REST 2 | 1,558.0 | +26.2 |
| REST 3 | 1,567.1 | +35.4 |
| REST 4 | 1,554.7 | +22.9 |
| REST 5 | 1,565.5 | +33.7 |
| REST 6 | 1,548.6 | +16.9 |
| REST 7 | 1,553.3 | +21.6 |

Maximum drift **+35.4 counts over seven cycles**, versus +410 over three. The
baseline steps once after the first contraction, then oscillates about that
level without growing.

**What was actually happening.** The 5-minute run's 40 s rests were too short for
the polymer to finish relaxing, so each baseline measurement caught the sensor
mid-decay and the apparent drift compounded. The +410 figure was **incomplete
recovery, not permanent hysteresis.** Measured residual recovery slope inside
each 100 s rest is +0.3 to +0.7 counts/s, still decaying but nearly settled by
the end of the window.

The single ~+34 step after the first cycle is real and is the standard FSR
conditioning effect: the polymer settles into a repeatable hysteresis loop after
its first load cycle and then stays there.

### 8.2 Consequences

- **Drift does not accumulate over a long recording.** The threat model is much
  smaller than initially recorded. The detector does not need special defences
  against a growing baseline.
- **Preconditioning is a real requirement.** Any future belt implementation
  should load the sensor once before recording begins, so the first genuine
  contraction is not the one that absorbs the conditioning step.
- **Protocol design lesson.** Rest phases must be long enough for recovery, or
  the recovery itself is misread as drift. 40 s was insufficient, 100 s is
  adequate. This applies to any future FSR characterisation.

### 8.3 Trough margins

Detrended inter-contraction troughs against the detection threshold:

| capture | trough peak | threshold | margin |
|---|---|---|---|
| 5 min, 40 s rests | 423.9 | 476.5 | 11 % |
| 20 min, 100 s rests | 40.9-83.4 | 380.9 | 78-89 % |

The near-miss false positive recorded in the 5-minute run was an artifact of the
short rest phases, not an inherent property of the sensor.

---

## 9. Detector run on device data

`04_detect.py` imports `uc_detector.py` unmodified from the repo root and prints
the resolved path, so the run is traceable to a specific copy of the algorithm.
No constants redefined.

**Result (5-minute capture):** 3 detected, 3 matched, 0 false positives. IoU
0.72 / 0.82 / 0.89. Amplitude reference 1,588.2 counts, threshold 476.5.
Confidence HIGH, rate 6.03 per 10 minutes, analysable 100 %.

**Result (20-minute capture):** 7 detected, 7 matched, 0 false positives. IoU
0.70 / 0.78 / 0.79 / 0.80 / 0.78 / 0.80 / 0.81. Amplitude reference 1,269.6
counts, threshold 380.9. Confidence HIGH, rate 3.46 per 10 minutes, analysable
100 %. IoU tightens to a consistent 0.78-0.81 after the first cycle, which
matches the conditioning step described in section 8.

**This is not a performance figure.** Three events on a scripted hand stimulus
scored against a protocol file, not expert annotations. It establishes that the
detector recognises device-produced morphology. Nothing more. The 77.32 % F1
belongs to CTU-CHB clinical data and cannot be claimed for the device: no
annotated device-specific ground truth exists, and with the enclosure descoped,
none will be produced within this project.

**Two behavioural differences from CTU-CHB, both recorded rather than fixed:**

1. **Dropout mask is inert.** `flat_run_mask` uses exact equality, which fires
   on CTU-CHB because that channel is integer-quantised. Device 4 Hz samples are
   means of 125 raw values and are never bitwise identical, so the mask reports
   0.00 %. Correct here (no dropout occurred) but the mechanism is silently
   absent. A real sensor disconnection would not be caught. Needs replacing with
   a variance-floor test.

2. **Detected spans run short.** Coverage of 72 % and 82 % on cycles 1 and 2,
   because a threshold crossing clips the low tails of the ramp and fall. This
   is definitional rather than error, but device-derived durations will read
   short if duration is ever reported.

---

## 10. Scope boundaries and requirements passed forward

These are **not pending work.** They are closed as out of scope, with the
requirement recorded so that future work does not have to rediscover it.

### Closed: belt enclosure

Descoped on cost (section 0). Requirements a future housing must satisfy:

- **Force concentrator, mandatory.** The FSR402's sensing area is only 12.7 mm
  across. A soft strap laid flat spreads force across the whole footprint rather
  than concentrating it on the sensing circle, which destroys sensitivity. The
  housing needs a small rigid puck or dome bearing directly on that circle.
- **Opposed force, not resting weight.** Confirmed empirically: seven books
  resting on the bare sensor produced a reading identical to no books, because
  the load bridged over the 0.4 mm sensor onto the desk surface. A belt
  compresses the sensor between strap and abdomen, which is mechanically a
  pinch. Any bench substitute must reproduce that geometry, which is why a
  spring clip was used rather than a weight.
- **Preload target 600-2,000 counts** at the chosen fixed resistor (section 4),
  re-derived per housing using `02_preload.py`.
- **No load on the flex tail.** The conductive traces run along it; sustained
  load or a sharp crease cracks them.

### Closed: on-body validation

Never in scope. No measurement was taken on a human abdomen. All contraction
morphology was hand-produced against a scripted protocol.

### Closed: fixed resistor selection

1 k, final for this build (section 3).

### Open: dropout detection on device

The exact-equality mask does not transfer to float device data (section 9). A
variance-floor test is the likely replacement. Unblocked and small, but not yet
designed.

### Open: causal filter conversion

Blocked on the multi-rate ISR. The toco channel tolerates latency well:
contractions last 30-60 s and arrive minutes apart, so reporting one 30 s late
is clinically meaningless for a reassurance device. A buffered sliding window
with deliberate output lag is therefore preferable to a strict causal rewrite,
and preserves the validated behaviour. This is a genuine asymmetry with the ECG
channel, which needs low latency for a live heart rate and does require causal
conversion.

---

## 11. Method notes worth keeping

- **Pinch, not weight.** Load-transfer geometry must match the intended device.
  The book test failed because resting an object on a 0.4 mm sensor is not the
  same mechanism as a belt compressing it.
- **Serial buffer must be flushed before timed captures.** The board streams
  continuously while Python blocks on `input()`. The OS buffer fills and
  overflows, so a capture begins with stale backlog and then loses samples.
  Observed as a 55.47 s gap and an 86.35 s offset between protocol cues and
  device timestamps. Fixed by `ser.reset_input_buffer()` at the top of
  `read_for` in both `02_preload.py` and `03_simulated_contractions.py`.
- **Detect event onset from data, not from cue timing.** Host-clock cues lag
  buffered serial data. `02_preload.py` originally assumed load onset at t=0 and
  consequently reported creep 14x too high and noise 60x too high, because the
  onset step sat inside the analysis window. Cues are navigation aids; the
  signal is the arbiter. Same principle as the Romagnoli annotations.
- **Check the operating point before running a long protocol.** Two seconds in
  Serial Monitor tells you what 95 seconds of capture would.
- **Record length must exceed the analysis window.** `BASELINE_WIN_S` is 600 s;
  a 300 s record makes `percentile_filter` degenerate to a global percentile and
  the rolling baseline cannot roll. Check this before interpreting any baseline
  result.

---

## 12. 20-minute capture: result

Protocol: 60 s baseline, then 7 cycles of 20 s ramp, 25 s hold, 20 s fall, 100 s
rest. Total 1,215 s (20.25 min), rate 3.46 contractions per 10 min.

**Purpose.** The 5-minute record was shorter than `BASELINE_WIN_S` (600 s), so
`percentile_filter` degenerated to a global percentile and the rolling baseline
could not roll: it tracked only 14 % of the apparent drift. This record is 2.02x
the window, making it the first honest test of the baseline mechanism.

### Capture quality

| quantity | value |
|---|---|
| samples | 607,508 |
| duration | 1,215.0 s |
| effective rate | 500.001 Hz |
| interval std | 0.7 us |
| gaps over 0.5 s | 0 |
| saturated / at zero | 0.00 % / 0.00 % |

The serial-flush fix held over twenty minutes of sustained streaming at ~65 %
link utilisation. No stale prefix, no overflow, no cue offset.

### Findings

1. **Baseline drift does not accumulate.** Maximum +35.4 counts across seven
   cycles. The 5-minute run's +410 was incomplete recovery caused by 40 s rest
   phases. See section 8, which has been corrected accordingly.

2. **The rolling baseline worked and could be observed working.** The estimator
   tracked 57 % of the residual drift, against 14 % on the short record. Drift
   is small enough that this matters little in practice, but the mechanism is
   now demonstrated rather than assumed.

3. **Trough margins are comfortable.** 78-89 % below threshold, against 11 % on
   the short record.

4. **Detection is consistent across cycles.** 7/7 matched, 0 false positives,
   IoU converging to 0.78-0.81 after the conditioning cycle. Hand-produced ramps
   varied noticeably in amplitude (HOLD means 2,613 to 2,819) without affecting
   detection, which is the expected behaviour of a record-relative threshold.

5. **Baseline noise rises within each rest phase.** Detrended std goes from 4.02
   counts in the opening rest to 8-11 counts in later ones, with a residual
   slope of +0.3 to +0.7 counts/s. This is ongoing recovery inside the rest
   window, not added noise, and it is the mechanism that produced the false
   ratchet reading on the short record.

### Status

This closes the last substantive open question on the toco channel. The
remaining open items (section 10) are the dropout mask and causal conversion,
neither of which is blocked on further capture work.
