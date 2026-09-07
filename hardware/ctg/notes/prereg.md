# Pre-registration — concurrent CTG capture

Phase 6. `hardware/ctg/`.

**This document is written and committed before any capture in this phase
exists.** Everything below is fixed at the commit that introduces this file.
Nothing in it may be revised in response to a capture, a detector output or a
figure. If something here turns out to be wrong, it gets a correction section
at the end of the phase notes, not an edit here.

---

## 0. What this phase is for

Every channel in this system is individually validated. The shared time base is
validated. **The system has never produced a cardiotocograph:** fetal heart rate
and contraction timing on one time axis from one concurrent capture. That figure
is what the project exists to produce and it does not exist.

This phase produces it once, from a single capture, using the frozen detectors
unmodified.

### 0.1 What the resulting figure will and will not show

**Will show.** Two channels acquired concurrently on one hardware-paced clock,
each scored by the detector already validated for it, plotted against a shared
time axis.

**Will not show.** Any physiological relationship between the two traces. The
FHR channel is a recorded fetal PCG played through a speaker at a microphone.
The UC channel is a hand-pressed force sensor. **There is no causal link between
them and none may be inferred.** No deceleration will appear, and no early or
late classification can be drawn from this figure.

This distinction goes in the figure caption itself, not only in surrounding
text. The shared time base is the *prerequisite* for deceleration analysis. It
is not deceleration analysis.

---

## 1. Frozen inputs

Imported by path with `importlib`, never copied. The resolved path of each is
printed at runtime and recorded in the results JSON.

| Component | File | Role |
|---|---|---|
| FHR detector | `02_fhr_detector.py` | `detect_fhr(raw, fs)` |
| FHR record rule | `03_evaluate.py` | record-level verdict |
| UC detector | `uc_detector.py` | `detect_contractions`, `confidence` |
| UC device runner | `hardware/toco_bringup/04_detect.py` | device-format UC scoring |
| Firmware | `hardware/mrisr/mrisr_capture/mrisr_capture.ino` | **unmodified** |

The firmware is not edited for this phase. The build that passed the paired
SENSOR_VN erratum control is the build that produces the headline figure, with
no intervening change to account for.

### 1.1 Detector parameters, recorded because they constrain the result

`02_fhr_detector.py`: bandpass 25–200 Hz; `win_sec=4.0`, `hop_sec=1.0`;
`min_bpm=110`, `max_bpm=180`; `CONFIDENCE_THRESHOLD = 0.45`.

`uc_detector.py`: `FS=4.0`; `BASELINE_WIN_S=600.0` (**module-level constant, not
an argument**); `BASELINE_PCTL=10`; `AMP_PCTL=98`; `THRESH_FRAC=0.30`;
`PEAK_FRAC=0.50`; `MERGE_GAP_S=20.0`; `MIN_DUR_S=20.0`; `MAX_DUR_S=300.0`;
`MIN_ANALYSABLE_FRAC=0.40`; `PLAUSIBLE_RATE=(1.0, 8.0)` per 10 analysable
minutes.

None of these is changed. If a result is disappointing, it is reported.

---

## 2. Hardware and physical setup, fixed before the capture

- ESP32 devkit, 115200 baud. PCG on GPIO39 at 500 Hz, UC on GPIO35 at 4 Hz,
  both from the single 500 Hz timer schedule, decimation 125.
- **Acoustic geometry identical to run C**: bare capsule, no chamber, same
  speaker-to-capsule distance, speaker volume 40%, MAX4466 trimpot untouched at
  factory position. Run C's 133.93 BPM / 0.808 is the reference for C2 and
  ceases to be a valid reference if the geometry changes.
- WAV looping, playback started **before** the capture begins, so every second
  of the record contains real audio. Same rule as run C.
- FSR402 through the 1 kΩ divider, fixed resistor on the ground side, spring
  clip providing static preload. Unloaded the sensor reads exactly 0.
- **Laptop unplugged from mains for the entire capture.** USB power with the
  charger connected was measured to raise mains content on the MAX4466 by
  roughly 40 dB; unplugging is the validated mitigation and it applies for the
  whole run, not only during electrode contact. Audio plays from the laptop's
  own output so nothing else near the bench is on mains.
- **The FSR sits on a physically separate surface from the microphone.** Not the
  same bench. See C4.

---

## 3. Protocol

**The 20-minute toco protocol, timings unchanged.**

| Phase | Duration | Repeats |
|---|---|---|
| Baseline, no load beyond preload | 60 s | 1 |
| Ramp | 20 s | 7 |
| Hold | 25 s | 7 |
| Fall | 20 s | 7 |
| Rest | 100 s | 7 |

Total 1,215 s. Expected rate 3.46 contractions per 10 minutes, inside the
detector's 1.0–8.0 plausibility gate. Record length is 2.02× `BASELINE_WIN_S`,
which is the condition under which the rolling baseline can actually roll.

**Reusing the exact protocol is deliberate.** The single-channel run of this
protocol returned 7/7 with 0 false positives. Repeating it with the PCG channel
live makes this a replication with a known prior answer, so any deviation is
attributable to the presence of the second channel rather than to a protocol
invented for this phase.

Cues are printed by the capture script but are **navigation aids only**. Onset
is determined from the signal, never from cue timing. Host-clock cues lag
buffered serial data; this was measured at 86.35 s in the toco phase.

---

## 4. Analysis rules, fixed before the data exists

**4.1 The FHR detector is called once on the whole record.** Not segmented.

`shannon_energy_envelope` normalises by the maximum of the whole array it is
given, and a 1,215 s capture of a looped 60 s WAV contains around 20 loop-seam
transients. I predicted this would force segmentation into 60 s chunks. **That
prediction was tested and it was wrong.** On a synthetic five-loop record with a
single transient inserted mid-record, whole-record and 60 s-segmented calls were
compared window by window on matched centre times:

| transient, × signal peak | median \|ΔBPM\| | max \|ΔBPM\| | max \|Δconf\| | windows crossing 0.45 |
|---|---|---|---|---|
| 2× | 0.0000 | 0.0000 | 0.010 | 0.0% |
| 20× | 0.0000 | 0.0000 | 0.016 | 0.0% |
| 100× | 0.0000 | 0.0000 | 0.020 | 0.0% |
| 1000× | 0.0000 | 0.0000 | 0.026 | 0.0% |

No window's BPM changed at any amplitude. Rescaling the normaliser by K turns
`−u·ln u` into `(1/K²)[original + 2·lnK·u]`; the autocorrelation removes the mean
and divides by `ac[0]`, so the `1/K²` cancels and only a small change in the
mix between envelope and squared signal survives. Synthetic, so it establishes a
property of the code rather than of any capture, but it holds across three
orders of magnitude.

Whole-record it is. Expected output on a 1,215 s record: 1,211 windows, first
centre at 2.0 s, last at 1,212.0 s.

**4.2 No trim.** The whole record is analysed. Run C needed a trim rule because
a 90 s capture of a 60 s WAV had to be reduced to a comparable 60 s; that does
not apply here.

**4.3 UC scoring is `04_detect.py` unmodified**, pointed at the
`contractions_<stamp>_uc4hz.csv` the capture script writes. That file is already
in `04_detect.py`'s format, float-valued, so `flat_run_mask` stays inert as
recorded in the toco phase. No new UC scoring code is written.

**4.4 One capture, one analysis.** The gold capture is scored once. If a script
bug is found after scoring, the fix is recorded, the rescore is recorded, and
**both results are reported**.

---

## 5. The null, measured before this document was written

The frozen FHR detector was run on three already-committed `mrisr` captures.
These are 90 s each of this exact hardware with no acoustic stimulus.

| capture | windows | median BPM | BPM IQR | median conf | max conf | % ≥ 0.45 |
|---|---|---|---|---|---|---|
| `mrisr_baseline_20260906_003840` | 87 | 146.34 | 41.93 | 0.058 | 0.298 | 0.0% |
| `mrisr_baseline_nofan_20260906_005014` | 87 | 129.30 | 46.37 | 0.036 | 0.244 | 0.0% |
| `mrisr_paired_multi_20260906_005816` | 87 | 150.00 | 53.36 | 0.056 | 0.240 | 0.0% |

**0 of 261 windows clear 0.45. The highest single window anywhere is 0.298.**

This sets no threshold: 0.45 was already frozen. It establishes that C1 below is
a criterion that discriminates on this hardware, rather than a formality. Null
0%, run C 100%.

It also bounds an existing open finding. The recorded finding is that the
detector returns OK at 0.809 on a pure 50 Hz tone and clusters BPM near 140 on
noise. Both stand. What is now measured is that **this device's own broadband
noise floor is rejected outright**, with wide BPM scatter (IQR 41.93–53.36)
rather than a tight false rate. The failure mode is therefore specific to
periodicity, which is what an autocorrelation detector is built to reward. That
makes the pure-tone result more pointed, not less.

---

## 6. Pre-registered criteria

### C1 — the FHR channel survives concurrency

Verdict OK: at least 50% of windows clear confidence 0.45. This is
`03_evaluate.py`'s record-level rule, unchanged.

Null 0% of 261 windows. Run C 100% of windows.

### C2 — the rate is preserved

Record median BPM within **3.0** of 133.93.

Same tolerance and same reasoning as run C: the capture is not phase-aligned to
the WAV, spans many loops, and run B's own window-to-window std was 2.46.

### C3 — UC replicates the single-channel result

**7 contractions matched, 0 false positives**, and `confidence()` returns HIGH.

Descriptive only. There is no annotated ground truth for a hand-pressed bench
capture, exactly as recorded in the toco phase. These are not performance
figures and must not be reported as sensitivity or PPV.

### C4 — mechanical cross-talk between channels

Seven 65 s presses will be performed roughly a metre from an electret capsule
whose passband is 25–200 Hz. Hand and bench transients live in that band. If
they couple in, they will corrupt FHR **precisely during the contractions**,
which is the only part of a cardiotocograph anyone reads.

**Measurement:** per-window FHR confidence during HOLD phases versus during REST
phases, reported as measured, with medians and distributions.

**Prediction: no systematic difference.**

**No threshold, and no pass/fail.** Inventing one now would be inventing a
threshold. The useful quantity is the measured difference in either direction.
If confidence drops during holds, that is a finding and it goes in the report.

### C5 — capture integrity

All must hold for the capture to be accepted. These are objective and
signal-independent.

- All PCG and UC timestamps strictly increasing.
- Achieved PCG fs within 0.1% of 500 Hz.
- Achieved UC fs within 0.1% of 4 Hz.
- Largest UC gap under 5.0 s, so `04_detect.py` trims nothing.
- The full 1,215 s protocol was completed.

---

## 7. Numerical predictions

Recorded so they can be wrong.

| Quantity | Prediction | Basis |
|---|---|---|
| PCG samples | ~607,500 | 1,215 s × 500 Hz |
| UC samples | ~4,860 | 1,215 s × 4 Hz |
| FHR windows | 1,211 | measured on an array of that length |
| Malformed lines | ~30 | ~1 in 20,000, measured across earlier runs |
| Raw PCG file size | ~10.3 MB | 14.88 B/sample measured, +2 B for 10-digit `t_us` |
| Decimation ratio | ~125.0 | 125.02 measured previously |
| Serial utilisation | 65.8% | measured, unchanged by this phase |
| Strongest in-band line | 47.80 Hz | present in all three `mrisr` captures |
| Dominant-line share | **below 8.21%** | share is relative to in-band power, and playback adds broadband in-band energy |

The last row is the one most likely to be wrong and is the most informative if
it is.

---

## 8. Capture acceptance and the re-run rule

**This is the rule that keeps the phase honest, so it is stated in full.**

The gold capture is **the first capture that satisfies C5**. Not the best one.

A capture may be **voided and repeated only for an operational reason**,
declared before any detector is run on it:

- The capture script crashed or the USB link dropped.
- The protocol was not completed, or a press cue was missed or badly mistimed.
- A C5 integrity gate failed.
- The speaker stopped, the clip slipped, or the charger was found plugged in.

**A capture may never be voided because of a detector result.** C1, C2, C3 and
C4 outputs are not grounds for a re-run. Every voided capture is recorded in the
notes with its reason and its files are kept.

If two captures are accepted under C5 through operator error, **the earlier one
is the gold capture.**

---

## 9. Reported, not deciding

Computed and recorded in the results JSON. None of these gates anything.

**9.1 Dominant in-band line and its share.** The strongest line between 25 and
200 Hz is **located**, not assumed to sit at 50 Hz, and its share of in-band
power is reported alongside its frequency. On this hardware the assumption is
known to be wrong: the strongest line is at 47.80 Hz, while 50 Hz measures only
+1.6 to +4.1 dB over background, below the +6 dB "detected" threshold.

Reference values, all with their conditions attached:

| Condition | Share | Source |
|---|---|---|
| Occluded capture, line at 50.40 Hz | 97.91% | `pcg_char_occluded_20260831_191155.json`, `spectrum_full` |
| Quiet floor, `pcg_bringup` baselines 1–3 | 2.58 / 4.07 / 4.43% | `notes.md:196` |
| Quiet floor, three `mrisr` captures at 47.80 Hz | 8.21 / 8.38 / 8.59% | measured for this pre-registration |

**This will be the first persisted measurement of this quantity on a capture
containing fetal heart sounds.** Precisely stated, that is speaker playback into
a bare capsule, not tissue-conducted fetal sound, so it does not close the gap
identified in §11. It converts "never measured on signal" into "measured once,
descriptively, on the loopback condition."

**9.2 Malformed line count and rate**, against the ~1 in 20,000 prediction.

**9.3 Cue-to-detection offset** for each of the seven contractions, so the lag
between host-clock cues and device timestamps is visible rather than asserted.

---

## 10. Scope boundaries — not established by this phase

**10.1 The FHR detector cannot represent a deceleration.**
`estimate_fhr_autocorr` has `min_bpm=110`. The lag search cannot return a value
below 110 BPM. Clinically significant decelerations go below that. So the frozen
detector, by construction, cannot represent the event that the FHR/UC pairing
exists to interpret. Widening the range would require full revalidation and is
**not** done here.

**10.2 Acausal filtering.** `bandpass_filter` uses `filtfilt`, which is acausal
and incompatible with live streaming. This prototype does not stream and has
never claimed to. Documented, not built. The toco channel additionally tolerates
latency well, since contractions last 30–60 s and arrive minutes apart, so a
buffered sliding window with deliberate output lag would preserve validated
behaviour better than a strict causal rewrite. The ECG channel does not share
that tolerance.

**10.3 The ECG channel is not in the concurrent capture.** Excluded on a
measured serial-budget constraint of 102.6% at 115200 baud, not as an omission.
It remains individually validated at 31/31 on device data and 99.69% Se over
95,960 annotated NIFECGDB beats.

**10.4 Signal-quality gate.** Not built. See §11.

**10.5 No physiological signal.** Restated because it matters most: no pregnant
participant is involved, the fetal PCG is a played-back recording, and the
contractions are a hand press.

---

## 11. Corrections carried in from the provenance audit

A repository audit was run before this document was written, tracing four
figures cited in `notes.md` back to the results files that produced them. Three
findings, all of which change what this phase does.

**11.1 A misattribution in `notes.md`.** The figure 4.07%, cited at `notes.md:319`
and `notes.md:578` as the dominant-line share for **real fetal PCG**, is
baseline 1 of the quiet-floor row at `notes.md:196`. It is a 47.80 Hz line on a
capture with no stimulus. The range 2.58–4.43% quoted as the noise floor is
baselines 2 and 3 of that same row. **4.07% sits inside 2.58–4.43% because it is
one of the three values defining it.** No dominant-line share for fetal PCG
exists anywhere in the repository; the loopback and reference outputs contain no
spectral share of any kind.

What the evidence supports is a two-point comparison between interference and
silence, which is real and large. It is not a comparison between interference
and signal. `notes.md:319` and `:578` are corrected accordingly, and the finding
is kept in reduced form rather than deleted.

**11.2 The signal-quality gate is closed as documented-only.** Three reasons,
all evidential rather than schedule-driven:

1. No fetal-PCG measurement of this metric exists, so a gate threshold has
   nothing to be calibrated against on the pass side.
2. The quiet floor of the same metric is 2.58–4.43% in `pcg_bringup` and
   8.21–8.59% in `mrisr`. Same hardware, no stimulus, a factor of two apart.
   A threshold fitted to one phase would misclassify the other.
3. The metric is not persisted, so it is not auditable.

**11.3 An unpersisted diagnostic detached from its condition.**
`quick_check_lines.py` prints its `pct_inband` and writes nothing. The figure
survived in the notes while its condition did not, and was later re-attributed
to a different condition entirely. This is **not** an instance of the §11
metric-design pattern in the PCG notes: that pattern is a metric specifying
something adjacent to the quantity that mattered, whereas this is a correctly
designed metric whose output lost its provenance. Recorded as a distinct
finding so the existing pattern stays legible.

The figure previously cited as 97.94% is replaced by **97.91%** from
`pcg_char_occluded_20260831_191155.json`, `spectrum_full`, 50 Hz
`pct_of_inband`, with a note that the strongest in-band line in that capture was
at 50.40 Hz and therefore falls inside the ±1 Hz window around 50, so the
fixed-frequency and strongest-line figures coincide there. Auditable from a
committed file, and the 0.03 difference does no work in the argument.

---

## 12. Demo configuration, declared here so it is not mistaken for a result

The live demonstration is a **90 s** capture, same rig, same firmware, same
scripts.

- The FHR detector runs and is fully valid at that length.
- **The UC detector does not run.** The UC channel is displayed as the
  conditioned trace only.

The reason, measured: `BASELINE_WIN_S` is 600 s, so `percentile_filter` receives
a size of 2,400 samples against a 90 s record of 360. The rolling baseline
collapses to a **single constant value** across the whole record. On a synthetic
90 s trace the baseline had exactly 1 distinct value, and `confidence()` still
returned **HIGH** at 6.67 contractions per 10 min, because it gates on
analysable fraction and rate and knows nothing about whether the record is long
enough for the mechanism those depend on to function.

The 90 s output is not wrong so much as produced by a configuration nothing
validated, while being labelled HIGH. `BASELINE_WIN_S` is **not** reduced for
the demo. The gold figure is open as the fallback.

---

## 13. Deliverables

```
hardware/ctg/
  01_capture.py                 1,215 s capture with cue schedule
  02_ctg.py                     both frozen detectors, CTG figure
  captures/
    mrisr_ctg_gold_<stamp>_pcg500.csv     ~10.3 MB, committed by exception
    contractions_<stamp>_uc4hz.csv
  plots/
    ctg_gold_<stamp>.png                  the headline figure
  results/
    ctg_gold_<stamp>.json
  notes/
    prereg.md                             this file, committed first
    notes.md                              written at the end of the phase
```

`.gitignore` gains a targeted exception for the gold PCG capture only, following
the precedent set for the run C CSV and the two chest captures. Without it the
headline figure cannot be regenerated from the repository.
