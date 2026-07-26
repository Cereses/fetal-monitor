# Home Fetal Monitoring System

A research/educational prototype for non-invasive home monitoring of three
physiological signals during pregnancy — **Fetal Heart Rate (FHR)**,
**Maternal Heart Rate (MHR)**, and **uterine contractions** — with IoT
integration on an ESP32.

> ⚠️ **Not a clinical device.** This is a university capstone prototype for
> research and educational purposes only. It is not validated, certified, or
> intended for medical diagnosis, monitoring, or patient care(at least not yet).

## Signal → Sensor mapping

| Signal | Sensor | Method |
|---|---|---|
| Fetal Heart Rate (FHR) | MAX4466 electret mic | Fetal phonocardiography (PCG) |
| Maternal Heart Rate (MHR) | AD8232 ECG module | Maternal ECG (+ artifact cancellation) |
| Maternal Heart Rate (MHR) | MAX30102 | PPG + SpO₂ |
| Uterine contractions | FSR402 force sensor | Tocodynamometry |

## Pipelines

**FHR (phonocardiography)** — Butterworth bandpass (25–200 Hz) →
Shannon-energy envelope → sliding-window autocorrelation → BPM estimation →
two-tier confidence (per-window threshold 0.45; per-record reliability flag
requires ≥50% of windows to clear the threshold).

**MHR — PPG (MAX30102)** — Bandpass 0.5–5 Hz → scipy peak-finding.
Target MAE < 3 BPM.

**MHR — ECG (AD8232)** — Pan–Tompkins (bandpass 5–15 Hz → differentiate →
square → moving-window integrate → adaptive threshold) → Random Forest
arrhythmia detection trained on MIT-BIH.

**Contractions (FSR402)** — Low-pass 1 Hz → rolling-median baseline
subtraction → duration-based event detection.

## Datasets (downloaded separately — not committed to git)

Raw waveform data is **not** stored in this repo. Download from PhysioNet:

| Dataset | PhysioNet slug | Use |
|---|---|---|
| Simulated fetal PCG | `simfpcgdb` | FHR validation |
| Fetal PCG | `fpcgdb` | FHR validation (real transabdominal) |
| BIDMC PPG & Respiration | `bidmc` | MHR (PPG) validation |
| MIT-BIH Arrhythmia | `mitdb` | ECG arrhythmia training |
| CTU-UHB Intrapartum CTG | `ctu-uhb-ctgdb` | Contraction validation |

Each script's `DATASET` constant parameterizes all paths — change one line to
switch datasets.

## Hardware constraints

- **AD8232 has no isolation.** The ESP32 must be **battery-powered** (never
  USB-to-mains) whenever electrodes are on a person.
- **ADC2 pins are unusable when WiFi is active.** All analog sensors must use
  ADC1 pins (GPIO 32–39).

## Setup

```bash
pip install -r requirements.txt
```

## Status

PC-side algorithm prototyping and dataset validation precede any hardware
integration. See `notes/` for phase summaries and `results/` for frozen
evaluation outputs.

## Team
The JAD Trio


## License


