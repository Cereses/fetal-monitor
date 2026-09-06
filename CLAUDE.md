# CLAUDE.md

Rules for working in this repository. Read before doing anything.

---

## What this project is

A multi-sensor non-invasive fetal monitoring prototype. Final-year Computer
Science capstone, KNUST. **Research/educational prototype, not a clinical
device.** No pregnant participant has been or will be involved without
institutional ethics approval, informed consent and clinical supervision; that
protocol is documented as designed-but-not-executed.

Three channels: fetal heart rate by acoustic phonocardiography (MAX4466),
maternal ECG (AD8232), uterine contraction timing (FSR402). The FHR/UC pairing
is the point — a deceleration is interpreted by its phase relative to the
contraction, so the two must share a time base or the output is meaningless.

---

## Your remit in this repository

**You are here for the writeup and for consistency auditing.** Nothing else was
delegated.

**In scope**

- Drafting and editing the final report from the notes, results and plots.
- Verifying that every number cited in a `notes.md` appears in the
  `results/*.json` or CSV that produced it.
- Cross-phase consistency: do the "established" lists match the measurements?
  Are the outstanding lists current? Do cross-references resolve?
- Checking that no script has copied a frozen file instead of importing it.
- Git hygiene: staging, ignore rules, commit messages.

**Out of scope unless explicitly asked in the session**

- Hardware bring-up. That is done step-by-step in chat with confirmation
  between stages, and that rhythm has caught real faults.
- Running detectors or generating new results.
- Design decisions on the causal filter conversion, the plausibility ceiling or
  the signal-quality gate. Implementation is fine once the approach is settled
  elsewhere.

If a task feels like it is drifting out of scope, say so rather than proceeding.

---

## Hard rules
### Git: do not commit, and never add attribution

**Do not run `git commit`, `git push`, `git add`, or any other command that
writes to the repository's history.** The author makes all commits personally
from his own account. If work is ready to commit, say so and stop.

If a commit is ever explicitly requested in-session, it must carry **no
attribution to Claude, Claude Code, or any AI tool** — no `Co-Authored-By`
trailer, no "Generated with" line, no emoji marker, nothing in the body. The
commit must be indistinguishable from one the author wrote by hand.

This is not a preference about credit. It is a final-year individually assessed
capstone, and the commit history is part of what is assessed.

### Frozen pipeline files are immutable

```
02_fhr_detector.py     FHR: bandpass -> Shannon envelope -> autocorrelation
03_evaluate.py         FHR record-level scoring rule
02_qrs_detect.py       Pan-Tompkins QRS
03_qrs_evaluate.py     QRS scoring
uc_detector.py         CTU-CHB contraction detector
```

**Never edit these. Never edit them in response to a disappointing number.**
They are marked read-only at the filesystem level; if you hit a permission
error on one of them, that is the guard working — stop and report it, do not
work around it.

**Import by path, never copy.** A copy silently forks the moment either file is
edited, which voids every comparison built on it. The established pattern is a
walk-up path resolver plus `importlib.util.spec_from_file_location` (the
filenames start with digits, so plain `import` will not work). See
`hardware/pcg_bringup/03_reference.py`.

If a frozen file genuinely needs to change, that is a new frozen version with a
new name and a re-run of everything downstream — not an edit.

### Pre-register before any detector run

Exclusions, thresholds and predictions are stated in writing **before** the run,
not after seeing the output. This is not ceremony; five metrics in this project
have been wrong, and pre-registration is what made each failure visible instead
of absorbing it into the conclusion.

### Verify by measurement, not by reading a plot

Several wrong conclusions here came from eyeballing a figure instead of
computing the number. If a claim is not backed by a value in a results file,
either compute it or do not make it.

### External ground truth as arbiter; never invent a threshold

Validate against established references: MIT-BIH, CTU-CHB, BIDMC, NIFECGDB,
fpcgdb, simfpcgdb. Where a threshold is needed, derive it from a measured
separation and say where the number came from.

### Structural diagnosis before parameter tuning

Understand the failure mechanism first. Do not tune away an honest failure.

### Report failure modes honestly

Named findings, explicit established/not-established splits, own errors
recorded. The notes contain five documented instances of the author's own
metrics being wrong. **Do not soften or remove these when drafting the report.**
They are among the strongest content in it.

---

## Repository layout

```
fetal-monitor/
├── 01_explore.py  02_fhr_detector.py  03_evaluate.py     FROZEN, offline FHR
├── 02_qrs_detect.py  03_qrs_evaluate.py                  FROZEN, offline ECG
├── uc_detector.py                                        FROZEN, offline UC
├── data/                     gitignored, PhysioNet
├── results/  plots/  notes/                              offline pipeline
└── hardware/
    ├── toco_bringup/         FSR402   — 7/7, 0 FP
    ├── ecg_bringup/          AD8232   — Se 100% device; 99.69% / 95,960 beats
    ├── pcg_bringup/          MAX4466  — run C 133.93 BPM, 100% reliable
    │                                    + notes_chamber_amendment.md (§14)
    └── mrisr/                two-channel ISR — erratum cleared
```

Each hardware phase: numbered scripts `01_`, `02_`, … plus `quick_check_*.py`
diagnostics, with `captures/ plots/ results/ notes/` beneath. Notes are written
at the **end** of a phase, not inline.

---

## Writing conventions for the report

- Scope limits and uncomfortable findings **first**, not buried.
- Measured numbers throughout, with their source file.
- Explicit established / not-established split.
- The author's own errors recorded as named findings.
- Claims stated at the strength the evidence supports. "Surface-conducted heart
  sound reaches the capsule" is supported; a rate measurement is not, because
  87.0 BPM by three routes against two counts of 71 and 70 is unresolved.

---

## Environment

Windows, PowerShell, `.venv`. `pip install` needs no special flags here.
`matplotlib.use('Agg')` must be set **before** importing pyplot or any frozen
module that imports it. Python 3.14: prefer `np.trapezoid` with a `np.trapz`
fallback.

Hardware, for reference only — do not attempt to drive it from here:
ESP32 devkit, 500 Hz, 115200 baud. GPIO39 (`SN`) MAX4466, GPIO35 FSR402,
GPIO34/33/32 AD8232. ADC1 only; ADC2 is unusable with WiFi active.
