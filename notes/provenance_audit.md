# Provenance audit — every cited figure against the committed results

**Date:** 2026-09-07
**Git HEAD:** `bf564f924b275c614850ee75dddcbb66788c18ff` — "Revised final note
claims for pcg and ctg" (working tree also carried one untracked file,
`hardware/ctg/01_capture.py`, outside any `results/` directory and outside
this audit's scope)
**Method, one line:** every notes file read in full and every numeric claim in
it checked by hand against the committed `results/` file and JSON/CSV key path
that should hold it, with aggregates recomputed from the raw rows rather than
trusted from prose.

**Scope:** all eleven notes documents under `hardware/*/notes/`,
`hardware/toco_bringup/` and `notes/`, audited against every committed file
under any `results/` directory (50 files, JSON and CSV).

Read-only pass. No notes file was edited to produce this, no pipeline script was
run, nothing was committed. Aggregates were recomputed by reading the committed
CSVs directly.

---

## 1. What was audited

Every figure presented as **a measurement, or as an aggregate over
measurements.**

Excluded as not figures in that sense: dates, commit hashes, GPIO numbers,
section cross-references, byte counts quoted inside design arithmetic, and
pre-registered thresholds that are stipulations rather than results. Where a
pre-registered constant was justified *by* a measurement — `MIN_DUR_S`,
`MIN_PEAK_SEP_S` — the justifying measurement is audited.

### Status meanings

| Status | Meaning |
|---|---|
| `EXACT` | The figure appears at full precision in a committed results file. |
| `ROUNDED` | It is the committed value, rounded or truncated as printed. |
| `DERIVED` | It is arithmetic over committed values. The arithmetic is stated. |
| `NOT FOUND` | It appears in no committed `results/` file of any extension. |

Figures whose mapping to a committed value is genuinely uncertain are marked
`[AMB]` and described rather than resolved.

### Totals

| Status | Count |
|---|---|
| EXACT | 59 |
| ROUNDED | 88 |
| DERIVED | 115 |
| NOT FOUND | 160 |
| **Total figures** | **422** |

36 of the 422 rows additionally carry `[AMB]`.

### File tags used in the ledger

| Tag | File |
|---|---|
| `fhr` | `notes/phase1_fhr_summary.md` |
| `mhr` | `notes/phase1_mhr_summary.md` |
| `uc` | `notes/phase1_uc_summary.md` |
| `toco` | `hardware/toco_bringup/notes.md` |
| `ecga` | `hardware/ecg_bringup/notes/notes.md` |
| `ecgn` | `hardware/ecg_bringup/notes/notes_nifecgdb.md` |
| `ecgp` | `hardware/ecg_bringup/notes/PREREGISTRATION_nifecgdb.md` |
| `pcg` | `hardware/pcg_bringup/notes/notes.md` |
| `pcgc` | `hardware/pcg_bringup/notes/notes_chamber_amendment.md` |
| `mri` | `hardware/mrisr/notes/notes_mrisr.md` |
| `ctg` | `hardware/ctg/notes/prereg.md` |

`pcg` line numbers reflect the file **after** the 2026-09-07 provenance edits
(§5 note, §7.2 correction, §7.3 citation, §11.5, §13 bullet).

---

## 2. The ledger

| Figure | Notes | Results source and key path | Status |
|---|---|---|---|
| 22/37 reliable (59%) | fhr:37 | `results/simfpcgdb_eval.csv` · `reliable` — file holds 36 records, 21 reliable = 58.3%; no `SNR-4_4` row exists | NOT FOUND `[AMB]` |
| 140.2–140.5 BPM, all OK records | fhr:38 | `results/simfpcgdb_eval.csv` · `estimated_bpm` — reliable rows take exactly 140.18691588785046 / 140.5152224824356 | EXACT |
| mean confidence 0.81 at −4.4 dB | fhr:40 | `results/simfpcgdb_eval.csv` · `mean_confidence` — no −4.4 dB record; highest committed is 0.8000111796 at SNR-6_6dB | NOT FOUND `[AMB]` |
| 0.24 at −26.7 dB | fhr:41 | `results/simfpcgdb_eval.csv` · `mean_confidence` [SNR-26_7dB] = 0.23675428021147069 | ROUNDED |
| cutoff between −22 dB and −22.6 dB | fhr:42 | `results/simfpcgdb_eval.csv` · `reliable` — SNR-22dB True, SNR-22_6dB False | EXACT |
| LOW outlier 153.3 BPM, 2% reliable | fhr:44 | `results/simfpcgdb_eval.csv` [SNR-26_2dB] — 153.2569549558954; `high_conf_fraction` 0.02109704641350211 | ROUNDED |
| 15/26 reliable (58%) | fhr:48 | `results/fpcgdb_eval.csv` · `reliable` — count True = 15; 15/26 = 57.7% | DERIVED |
| estimates in 124.1–159.8 BPM | fhr:49 | `results/fpcgdb_eval.csv` · `estimated_bpm` — min reliable 124.09937888198758 (p11), max 159.84 (p19) | ROUNDED |
| best p21: 0.76 conf, 97% reliable | fhr:51 | `results/fpcgdb_eval.csv` [fetal_PCG_p21_GW_39] — 0.7630745554720815; 0.9741235392320534 | ROUNDED |
| worst p02: 0.15 conf, 9% reliable, 180.0 BPM | fhr:52 | `results/fpcgdb_eval.csv` [fetal_PCG_p02_GW_31] — 0.15221715962960486; 0.09136212624584718; 180.0 | ROUNDED |
| 58% vs 59% across datasets | fhr:63 | both eval CSVs — fpcgdb 57.7%→58% holds; the 59% depends on the missing simfpcgdb record | NOT FOUND `[AMB]` |
| every OK record within 110–160 BPM | fhr:68 | both eval CSVs · `estimated_bpm` where reliable — range 124.10–159.84 and 140.19–140.52; holds | DERIVED |
| real quality ≈ −23 to −24 dB SNR | fhr:58 | — no committed file records this mapping | NOT FOUND `[AMB]` |
| 46/53 reliable, 7 LOW | mhr:85 | `results/bidmc_eval.csv` · `reliable` — 7 rows at 0 (bidmc17,26,28,33,41,49,53); 53−7 = 46 | DERIVED |
| MAE vs HR mean 1.13, median 0.75 | mhr:86 | `results/bidmc_eval.csv` · `mae_hr` where reliable==1 — Σ = 52.142 over 46 rows → 1.1335; sorted 23rd/24th = 0.747, 0.752 → 0.7495 | DERIVED |
| MAE vs PULSE mean 1.43, median 0.98 | mhr:87 | `results/bidmc_eval.csv` · `mae_pulse` where reliable==1 — Σ = 65.556 → 1.4251; 23rd/24th = 0.977, 0.990 → 0.9835 | DERIVED |
| 44/46 meet MAE < 3 BPM | mhr:88 | `results/bidmc_eval.csv` · `mae_hr` — only bidmc24 (5.195) and bidmc45 (7.588) exceed 3; 46−2 = 44 | DERIVED |
| aromring 68.2 BPM, 0.92 conf, 4 peaks | mhr:96 | — no aromring results file exists anywhere in `results/` | NOT FOUND |
| bidmc45 MAE 7.59; bidmc24 MAE 5.20 at 53% | mhr:121 | `results/bidmc_eval.csv` — 7.588; 5.195; coverage 0.529 | ROUNDED |
| median 0.75 vs mean 1.13 gap | mhr:127 | `results/bidmc_eval.csv` · `mae_hr` — same aggregation as above | DERIVED |
| most records 90%+ reliability | mhr:109 | `results/bidmc_eval.csv` · `high_conf_fraction` — 33 of 53 rows ≥ 0.90 | DERIVED |
| TEST: 84 records, 1,234 contractions | uc:11 | `results/uc_evaluation.csv` · split=TEST — 84 rows, Σ `n_ref` = 1234 | DERIVED |
| TEST Sensitivity 77.07% | uc:15 | `results/uc_evaluation.csv` · Σtp/(Σtp+Σfn) — 951/1234 = 0.770665 | DERIVED |
| TEST PPV 77.57% | uc:16 | `results/uc_evaluation.csv` · Σtp/(Σtp+Σfp) — 951/1226 = 0.775693 | DERIVED |
| TEST F1 77.32% | uc:17 | harmonic mean of the two above = 0.773197 | DERIVED |
| DEV Se 77.32 / PPV 73.44 / F1 75.33 | uc:15-17 | `results/uc_evaluation.csv` · split=DEV — TP 951, FP 344, FN 279 → 0.773171 / 0.734363 / 0.753265 | DERIVED |
| DEV and TEST Se agree within 0.25 points | uc:21 | `results/uc_evaluation.csv` — 77.3171 − 77.0665 = 0.2506 pp | DERIVED |
| any overlap: Se 77.39%, PPV 77.90% | uc:28 | — `uc_evaluation.csv` carries one matching criterion only | NOT FOUND |
| IoU ≥ 0.5: Se 67.50%, PPV 67.94% | uc:30 | — no IoU-criterion columns in any committed file | NOT FOUND |
| 0.32-point gap between any and half | uc:32 | — depends on the uncommitted any-overlap figures | NOT FOUND |
| HIGH: 79 records, 1,213 ref, 77.33 / 77.91 / 77.62 | uc:40 | `results/uc_evaluation.csv` · TEST, confidence=HIGH — TP 938, FP 266, FN 275 | DERIVED |
| LOW: 5 records, 21 ref, 61.90 / 59.09 / 60.47 | uc:41 | `results/uc_evaluation.csv` · TEST, confidence=LOW — TP 13, FP 9, FN 8 | DERIVED |
| Onset error median +2.5 s, IQR −2.0 to +8.8 | uc:50 | — `uc_detections.csv` holds no reference boundaries | NOT FOUND |
| Offset error median −2.0 s, IQR −8.6 to +3.8 | uc:51 | — same, no `ref_start`/`ref_end` columns | NOT FOUND |
| Duration error median −5.0 s, IQR −17.9 to +6.1 | uc:52 | — same | NOT FOUND |
| 200 records downloaded of 552 | uc:65 | `results/ctu_ann_manifest.csv` — 552 rows EXACT; the 200-download figure is not a column | NOT FOUND `[AMB]` |
| Bradycardia 10/552, median episode 600 s | uc:86 | — manifest has no per-event-type columns | NOT FOUND |
| Tachycardia 11/552, median 660 s | uc:87 | — same | NOT FOUND |
| Acceleration ~166 events | uc:88 | — same | NOT FOUND |
| Deceleration ~490 events | uc:89 | — same | NOT FOUND |
| rows 3:2 = 490:166 = 2.95:1 | uc:98 | — ratio checks (2.9518) but neither input is committed | NOT FOUND |
| Row 4 median 60.0 s, IQR 46.5–76.5 s | uc:90 | `results/ctu_contractions.csv` · `duration_s` — median 60.0 EXACT; quartiles compute to 46.0 and 76.0 | NOT FOUND `[AMB]` |
| start-to-start median 141 s = 4.3 per 10 min | uc:95 | `results/ctu_contractions.csv` · `start_s` — within-record intervals give 135.125 s (4.44 per 10 min) | NOT FOUND `[AMB]` |
| Record 1025: 37 annotated spans | uc:100 | `results/ctu_ann_manifest.csv` [1025] · `n_contractions` = 37 | EXACT |
| 12 of 552 records carry an odd count | uc:106 | `results/ctu_ann_manifest.csv` · `n_truncated` > 0 on 12 rows | EXACT |
| unpaired marker 97.7–99.6% through, 21–90 s from end | uc:107 | — manifest records the count, not the position | NOT FOUND |
| 7,015 contractions | uc:111 | `results/ctu_contractions.csv` — row count 7,015 | EXACT |
| 473 scorable records | uc:111 | `results/ctu_ann_manifest.csv` — 552 − 79 excluded = 473 | DERIVED |
| median 14 contractions per record | uc:111 | `results/ctu_ann_manifest.csv` · `n_contractions` over 473 records = 14 | DERIVED |
| Duration p5 30.0 / median 60.0 / p95 124.3 s | uc:112 | `results/ctu_contractions.csv` · `duration_s` — linear percentiles 30.0000, 60.0000, 124.3250 | DERIVED |
| no_annotation excludes 79 records | uc:122 | `results/ctu_ann_manifest.csv` · `exclusion_reason` — 79 rows | EXACT |
| median flat fraction 60.9% | uc:131 | `results/ctu_uc_quality.csv` · `flat_frac` — median 0.609471 | ROUNDED |
| only 18 of 50 records usable | uc:132 | `results/ctu_uc_quality.csv` — 32 rows exceed flat_frac 0.50; 50 − 32 = 18 | DERIVED |
| quantisation step 1.000, range 0.5–3.0 | uc:135 | `results/ctu_uc_quality.csv` · `q_step` — 0.5×11, 1.0×31, 1.5×2, 2.0×5, 3.0×1 | EXACT |
| record 1017 carries 25 contractions | uc:143 | `results/ctu_ann_manifest.csv` [1017] · `n_contractions` = 25 | EXACT |
| 42 of 50 carried annotations vs 18 called usable | uc:143 | `results/ctu_uc_quality.csv` · `n_contractions` > 0 on 42 rows | DERIVED |
| OLD flat_frac > 0.50: 25/42 lost, 7/8 caught | uc:152 | `results/ctu_uc_quality.csv` — cross-tab against `n_contractions` > 0 | DERIVED |
| dropout_frac > 0.30: 13/42, 6/8 | uc:153 | `results/ctu_uc_quality.csv` — same cross-tab | DERIVED |
| dropout_frac > 0.50: 3/42, 1/8 | uc:154 | `results/ctu_uc_quality.csv` — same | DERIVED |
| zerorun_frac > 0.50: 0/42, 0/8 | uc:155 | `results/ctu_uc_quality.csv` — no record exceeds the threshold | DERIVED |
| median longest flat run 510 s | uc:161 | `results/ctu_uc_quality.csv` · `max_run_s` — median 509.625 | ROUNDED |
| dropout_frac median 25.8% | uc:162 | `results/ctu_uc_quality.csv` · `dropout_frac` — median 0.257885 | ROUNDED |
| only 3 TEST references missed inside masked regions | uc:163 | — no committed file records mask-interior misses | NOT FOUND |
| MIN_PEAK_SEP_S: start-to-start p0.5 = 39 s | uc:194 | `results/ctu_contractions.csv` — linear p0.5 over 6,542 intervals = 39.029 s | DERIVED |
| 60 s discarded 2% of real pairs corpus-wide | uc:195 | `results/ctu_contractions.csv` — 131 of 6,542 intervals < 60 s = 2.00% | DERIVED |
| and 17% of record 1025 | uc:196 | `results/ctu_contractions.csv` [1025] — 36 intervals; not reproduced by a <60 s count | NOT FOUND `[AMB]` |
| record 1025 median interval 88 s, minimum 38 s | uc:197 | `results/ctu_contractions.csv` [1025] — median 88.125, min 38.00 | ROUNDED |
| corpus median interval 135 s | uc:197 | `results/ctu_contractions.csv` — median 135.125 s | ROUNDED |
| MIN_DUR_S: durations p1 = 18 s, p2 = 23 s | uc:198 | `results/ctu_contractions.csv` · `duration_s` — linear p1 18.035, p2 23.000 | DERIVED |
| 25 s discarded 2.6% against 1.3% at 20 s | uc:199 | `results/ctu_contractions.csv` — 180/7015 = 2.57%; 88/7015 = 1.25% | DERIVED |
| record 1006 unmarked contraction at minute 25 | uc:213 | `results/ctu_ann_manifest.csv` [1006] — `n_contractions` 0, excluded no_annotation; the event itself is not recorded | DERIVED |
| matched peak amplitude median 0.85, IQR 0.68–1.03 | uc:218 | `results/uc_detections.csv` · `peak_rel`, matched=1, TEST (n=951) — q25 0.6813, med 0.8526, q75 1.0281 | DERIVED |
| unmatched 0.76, IQR 0.61–0.95 | uc:219 | `results/uc_detections.csv` · matched=0, TEST (n=275) — q25 0.6093, med 0.7592, q75 0.9505 | DERIVED |
| 275 false positives | uc:221 | `results/uc_evaluation.csv` · Σfp TEST = 275; 275 unmatched rows in `uc_detections.csv` | EXACT |
| record 1029 → LOW at 0.72 per 10 min | uc:227 | `results/uc_evaluation.csv` [1029] · `rate_per_10min` 0.7172207695181172 | ROUNDED |
| record 1025 detected 30 of 37 | uc:230 | `results/uc_evaluation.csv` [1025] · `n_det` 30, `n_ref` 37 | EXACT |
| 168 scorable, split 84/84 | uc:235 | `results/uc_evaluation.csv` · `split` — 168 rows, 84 DEV / 84 TEST | EXACT |
| 500.000 Hz over 90 s, std 9.4 µs, one late slot in 45,000 | toco:87 | — no toco results file records capture timing | NOT FOUND |
| 10 k unloaded: 0 counts, R > 1 M | toco:96 | `results/preload_20260807_161207.csv` — `empty_mean`/`empty_std` 0.0; `r_ohms` inf | DERIVED |
| 10 k spring clip ~3,170 counts, ~2,920 Ω | toco:97 | — no committed preload CSV holds a 3,170-count row | NOT FOUND |
| 1 k spring clip ~1,085 counts, ~2,780 Ω | toco:98 | — not committed | NOT FOUND |
| finger press ~1,500 counts, ~17 kΩ | toco:100 | — not committed | NOT FOUND |
| predicted 1,045, measured 1,085 | toco:106 | — not committed | NOT FOUND |
| excursion 485 counts at 1 k vs 412 at 10 k | toco:108 | — not committed | NOT FOUND |
| unloaded FSR reads exactly 0, zero variance | toco:123 | `results/preload_20260807_161207.csv` — both rows read 0.0 / 0.0 | EXACT |
| across 9,750 samples | toco:123 | — sample count is not a column in the preload schema | NOT FOUND |
| node ~33 mV unloaded; ADC reports 0 below ~150 mV | toco:125 | — not committed | NOT FOUND |
| baseline std 0.000 on first capture | toco:132 | `results/preload_20260807_160807.csv` · `settled_std` 0.0, `empty_std` 0.0 | EXACT |
| target resting window 600–2,000 counts | toco:144 | — stipulated design window, not a measurement | NOT FOUND |
| 50 Hz share: finger press 2 = 0.03% | toco:162 | — no toco script computes a spectrum | NOT FOUND |
| finger press 3 = 0.01% | toco:163 | — same | NOT FOUND |
| clip static load = 0.74% | toco:165 | — same | NOT FOUND |
| divider node 800 Ω to 17 kΩ | toco:170 | `results/preload_20260807_170647.csv` · `r_ohms` = 1958.36 for the one real row | NOT FOUND `[AMB]` |
| below 0.5 Hz = 4.07% | toco:178 | — no spectral producer in the toco phase | NOT FOUND |
| 0.5–2 Hz = 0.86% | toco:179 | — same | NOT FOUND |
| 2–10 Hz = 3.36% | toco:180 | — same | NOT FOUND |
| 10–45 Hz = 12.36% | toco:181 | — same | NOT FOUND |
| 51–250 Hz = 76.53% | toco:182 | — same | NOT FOUND |
| 6–11% of power in 1–10 Hz | toco:188 | — same | NOT FOUND |
| detrended std 9.70 counts at 500 Hz | toco:203 | `results/preload_*.csv` · `noise_std` — the only committed non-zero row reads 614.678, the pre-fix run | NOT FOUND `[AMB]` |
| detrended std 3.64 counts at 4 Hz | toco:204 | — no script decimates and re-measures | NOT FOUND |
| improvement factor 2.67× | toco:205 | 9.70/3.64 = 2.665 — arithmetic holds, neither input committed | DERIVED |
| theoretical factor 11.18× (√125) | toco:206 | √125 = 11.1803 | DERIVED |
| 4% of power below 0.5 Hz | toco:209 | — restates the uncommitted band table | NOT FOUND |
| baseline std 3.92 counts at 500 Hz with 1 k | toco:212 | `hardware/mrisr/results/mrisr_*.json` · `baselines.uc_std` = 3.92 | EXACT |
| contraction amplitudes 1,000–1,570 counts | toco:213 | — not committed | NOT FOUND |
| SNR around 275:1 | toco:214 | 1570/3.92 = 400; 1000/3.92 = 255; 275 corresponds to ~1078 counts | NOT FOUND `[AMB]` |
| creep +0.615 counts/s | toco:222 | `results/preload_*.csv` · `creep_cps` — committed value is 8.7835, the pre-fix run | NOT FOUND `[AMB]` |
| creep +51 counts, +1.50% | toco:223 | `results/preload_*.csv` · `creep_total`, `creep_pct` — committed 862.948 and 25.200 | NOT FOUND `[AMB]` |
| ~37 counts over 60 s window | toco:225 | 0.615 × 60 = 36.9 — arithmetic on an uncommitted input | DERIVED |
| 83 s of settled data | toco:221 | — not a column in the preload schema | NOT FOUND |
| 5-min rests climbed +57, +97, +410 counts | toco:243 | — rest-phase baselines are printed, not written | NOT FOUND |
| 20-min rest table 1,531.8 → 1,553.3, drifts +16.9 to +35.4 | toco:254-260 | — no baseline column in any committed file | NOT FOUND |
| maximum drift +35.4 counts over seven cycles | toco:262 | — same | NOT FOUND |
| residual recovery slope +0.3 to +0.7 counts/s | toco:270 | — same | NOT FOUND |
| single ~+34 step after first cycle | toco:273 | — same | NOT FOUND |
| 5-min trough 423.9, threshold 476.5, margin 11% | toco:295 | — printed by `04_detect.py`, not in the detect CSV | NOT FOUND |
| 20-min troughs 40.9–83.4, threshold 380.9, margin 78–89% | toco:296 | — same | NOT FOUND |
| 5-min: 3 detected, 3 matched, 0 FP | toco:309 | `results/contractions_20260808_151139_detect.csv` — three rows, all `matched` 1 | EXACT |
| 5-min IoU 0.72 / 0.82 / 0.89 | toco:309 | `results/contractions_20260808_151139_detect.csv` · `iou` 0.7200, 0.8200, 0.8929 | ROUNDED |
| amplitude reference 1,588.2, threshold 476.5 | toco:310 | — printed only | NOT FOUND |
| rate 6.03 per 10 min, analysable 100% | toco:311 | — printed only | NOT FOUND |
| 20-min: 7 detected, 7 matched, 0 FP | toco:313 | `results/contractions_20260808_193404_detect.csv` — seven rows, all `matched` 1 | EXACT |
| 20-min IoU 0.70/0.78/0.79/0.80/0.78/0.80/0.81 | toco:313 | `results/contractions_20260808_193404_detect.csv` · `iou` 0.7038 … 0.8077 | ROUNDED |
| amplitude reference 1,269.6, threshold 380.9 | toco:314 | — printed only | NOT FOUND |
| rate 3.46 per 10 min | toco:315 | 7 contractions over 1,215 s → 3.4568; neither input in a results file | DERIVED |
| IoU tightens to 0.78–0.81 after first cycle | toco:316 | `results/contractions_20260808_193404_detect.csv` — cycles 2–7 span 0.7769–0.8077 | DERIVED |
| 77.32% F1 belongs to CTU-CHB | toco:321 | `results/uc_evaluation.csv` · TEST micro F1 — cross-phase citation, matches | DERIVED |
| dropout mask reports 0.00% | toco:331 | — printed only | NOT FOUND |
| coverage 72% and 82% on cycles 1 and 2 | toco:335 | `results/contractions_20260808_151139_detect.csv` · `cover` 0.7200, 0.8200 | EXACT |
| seven books read identical to no books | toco:356 | `results/preload_20260807_161207.csv` [books_x7] — all columns 0.0, identical to phone_x1 | EXACT |
| sensor 12.7 mm across, 0.4 mm thick | toco:351 | — datasheet values, not measurements from this project | NOT FOUND |
| 55.47 s gap and 86.35 s cue offset | toco:401 | — not committed | NOT FOUND |
| creep 14× and noise 60× too high | toco:406 | `results/preload_20260807_170647.csv` — 8.7835/0.615 = 14.3×; 614.678/9.70 = 63.4× | DERIVED |
| BASELINE_WIN_S 600 s; record 2.02× the window | toco:425 | 1,215/600 = 2.025; record duration not committed | DERIVED |
| 607,508 samples, 1,215.0 s, 500.001 Hz, 0.7 µs, 0 gaps | toco:432-436 | — no capture-quality columns in any toco results file | NOT FOUND |
| rolling baseline tracked 57% vs 14% of drift | toco:449 | — printed only | NOT FOUND |
| HOLD means 2,613 to 2,819 | toco:458 | — printed only | NOT FOUND |
| detrended std 4.02 → 8–11 counts across rests | toco:461 | — printed only | NOT FOUND |
| mains hum 90.7% → 79.9% → 3.9% | ecga:27 | — no committed file records hum share | NOT FOUND |
| ~1.1% of samples clipping | ecga:54 | — not committed | NOT FOUND |
| Se 99.52% over 46 MLII records | ecga:64 | `results/mitdb_e1_eval.csv` · Σtp/(Σtp+Σfn), lead=MLII — TP 104,570, FN 508 → 99.5165% | DERIVED |
| PPV 99.35% | ecga:64 | `results/mitdb_e1_eval.csv` — micro MLII 99.7215%, macro 99.6993%, all-48 99.7219%; no subset yields 99.35% | NOT FOUND `[AMB]` |
| records 102/104 excluded on lead grounds | ecga:66 | `results/mitdb_e1_eval.csv` · `lead` — the only two non-MLII rows, both V5 | EXACT |
| MWI 0.150 / refractory 0.200 / refine 0.075 / T-wave 0.360 s | ecga:97-100 | — source-code constants, not results | NOT FOUND |
| row 6874 of 180817 reads 7510,1776 | ecga:111 | — raw capture inspection; `captures/` is not `results/` | NOT FOUND |
| minima 987–1264 across five captures | ecga:132 | — not committed | NOT FOUND |
| tallest surviving sample 1990 against ceiling 4095 | ecga:190 | — the masking test is not persisted | NOT FOUND |
| every unmasked sample ≥51% below ceiling | ecga:192 | 1990/4095 = 48.6%, so 51.4% below; input uncommitted | DERIVED |
| lead-in 248 ms vs median RR 987 ms | ecga:196 | — not committed | NOT FOUND |
| tail 808 ms vs shortest RR 748 ms | ecga:197 | — not committed | NOT FOUND |
| 10 of 31 plateaus affected, net 0.000 ms | ecga:213 | `results/*_annotations_manifest.json` · `centre_convention.n_offset_beats` = 10 | EXACT |
| mean offset 0.16129 samples (0.645 ms) | ecga:226 | `results/*_annotations_manifest.json` · `mean_offset_samples` 0.16129032258064516, `mean_offset_ms` 0.6451612903225806 | EXACT |
| capture sha 1caf85070d638014 | ecga:262 | `results/*_device_eval.json` · `capture_sha256` | EXACT |
| annotations sha db909bfafb92bbd4 | ecga:263 | `results/*_device_eval.json` · `annotation_sha256` | EXACT |
| detector sha 71eb60e950f64edd | ecga:264 | `results/*_device_eval.json` · `detector_sha256`; also in both nifecgdb summaries | EXACT |
| evaluator sha a768b352f89b7efd | ecga:265 | `results/*_device_eval.json` · `evaluator_sha256` | EXACT |
| 7500 samples, 250.0 Hz, 30.0 s | ecga:275 | `results/*_annotations_manifest.json` · `fs_hz` 250.0; 30.0 × 250 = 7,500 | DERIVED |
| 31 reference beats, 0 manual | ecga:275 | `results/*_annotations_manifest.json` · `n_reference_beats` 31, `n_manual_entries` 0 | EXACT |
| applied tolerance 38 samples = 152 ms | ecga:276 | `results/*_device_eval.json` · `tolerance_ms` records the nominal 150.0 only; the notes flag this at line 422 | NOT FOUND |
| TP/FP/FN 31 / 0 / 0 | ecga:282 | `results/*_device_eval.json` · `TP`, `FP`, `FN` | EXACT |
| sensitivity 100.00%, PPV 100.00%, DER 0.00% | ecga:283-285 | `results/*_device_eval.json` · `sensitivity` 1.0, `PPV` 1.0, `DER` 0.0 | EXACT |
| 95% CI [88.78%, 100%] | ecga:283 | Clopper–Pearson exact lower bound for 31/31: 0.025^(1/31) = 0.88784 | DERIVED |
| lower bound 96.4% at n=100, 99.3% at n=500 | ecga:289 | 0.025^(1/100) = 0.96382; 0.025^(1/500) = 0.99265 | DERIVED |
| one beat worth 3.23 percentage points | ecga:290 | 1/31 = 3.226% | DERIVED |
| signed mean +0.290 samples / +1.16 ms | ecga:296 | `results/*_device_eval.json` · `timing_signed_mean_ms` 1.1612903225806452; ÷4 ms = 0.29032 | ROUNDED |
| MAE 0.290 samples / 1.16 ms | ecga:298 | `results/*_device_eval.json` · `timing_mae_ms` 1.1612903225806452 | ROUNDED |
| uncorrected signed mean +0.452 / +1.81 ms | ecga:300 | 0.29032 + 0.16129 = 0.45161 samples; × 4 = 1.8065 ms | DERIVED |
| worst single beat +1.0 sample / +4.0 ms | ecga:299 | — per-beat offsets are not written to any committed file | NOT FOUND |
| offset distribution 17 / 10 / 4 beats | ecga:305-307 | only the 10 even-width plateaus are recorded; 17·0 + 10·0.5 + 4·1 = 9.0 reconciles to the committed mean | NOT FOUND `[AMB]` |
| correction removes exactly 0.645 ms | ecga:312 | `results/*_annotations_manifest.json` · `mean_offset_ms` 0.6451612903225806 | ROUNDED |
| reference median 60.8 BPM (55.4–80.2); detected 60.9 | ecga:322 | — `device_eval.json` holds no rate fields; `mhr_ecg_summary` device_capture is a different quantity (62.2 BPM) | NOT FOUND `[AMB]` |
| weakest QRS blob 7.22e8 | ecga:337 | — integration-stage margins are printed | NOT FOUND |
| strongest non-beat feature 5.52e7 | ecga:338 | — same | NOT FOUND |
| separation ratio 13.1× | ecga:339 | 7.22e8/5.52e7 = 13.08; neither input committed | DERIVED |
| weakest blob 3.22× threshold; noise 0.246× | ecga:340-341 | — printed only | NOT FOUND |
| five captures at 0.99–1.36% clipping | ecga:406 | — not committed | NOT FOUND |
| 55 records, 22+1 to 40+2 weeks | ecgn:5 | `results/nifecgdb_summary.json` · `records_scored` 55 | EXACT |
| annotations maternal, median 92.7 BPM | ecgn:149 | — `03_explore_nifecgdb.py` prints only | NOT FOUND |
| RR method: median 7.33% missing | ecgn:158 | — same | NOT FOUND |
| beats ÷ duration method: 5.84% | ecgn:159 | — same | NOT FOUND |
| 53/55 records exceed 5% | ecgn:161 | — same | NOT FOUND |
| 15/55 sub-200 ms; ecgca998 10.12%; rest 0.07–3.42% | ecgn:170-172 | — same | NOT FOUND |
| reference beats 47,980 | ecgn:136 | `results/nifecgdb_summary.json` — (95,661 + 299) over 110 thoracic channels = 2 per record → 47,980 | DERIVED |
| total duration 549.7 min (9.16 h) | ecgn:136 | — no duration column in `nifecgdb_channel_results.csv` | NOT FOUND |
| record length 114 s to 2780 s | ecgn:137 | — same | NOT FOUND |
| channel sets (2,3)×17 and (2,4)×38 | ecgn:138 | `results/nifecgdb_channel_results.csv` — 17·3 + 38·4 = 203 abdominal | DERIVED |
| 313 channel-records — 110 thoracic, 203 abdominal | ecgn:143 | `results/nifecgdb_summary.json` · `channel_records`, `thoracic.n`, `abdominal.n` | EXACT |
| thoracic ref 95,960 / TP 95,661 / FP 7,962 / FN 299 | ecgn:182 | `results/nifecgdb_summary.json` · `thoracic` — TP/FP/FN exact; ref = TP + FN | EXACT |
| thoracic Se 99.69%, PPV 92.32% | ecgn:182 | `results/nifecgdb_summary.json` — 0.9968841183826594, 0.9231637763816913 | ROUNDED |
| abdominal ref 183,365 / TP 168,185 / FP 17,541 / FN 15,180 | ecgn:183 | `results/nifecgdb_summary.json` · `abdominal` | EXACT |
| abdominal Se 91.72%, PPV 90.56% | ecgn:183 | `results/nifecgdb_summary.json` — 0.9172142993482943, 0.9055544188751171 | ROUNDED |
| pooled placement cost 7.97 pp | ecgn:185 | `results/nifecgdb_summary.json` · `placement_cost_se_pp` 7.96698190343651 | ROUNDED |
| supplementary Se 57.47%, PPV 89.57% over 6 channels | ecgn:189 | `results/nifecgdb_channel_results.csv` · supplementary=True — 57.4713% / 89.5713% | DERIVED |
| abdominal Se bands 144 / 21 / 2 / 11 / 25 | ecgn:195-199 | `results/nifecgdb_channel_results.csv` · `sensitivity` at ≥99 / 95–99 / 90–95 / 50–90 / <50 | DERIVED |
| shares 70.9 / 10.3 / 1.0 / 5.4 / 12.3% | ecgn:195-199 | each band count ÷ 203 | DERIVED |
| abdominal median 100.00%, mean 90.11% | ecgn:201 | `results/nifecgdb_channel_results.csv` — median 1.000000, mean 0.901122 | DERIVED |
| thoracic median 100%, 2 of 110 below 95% | ecgn:201 | `results/nifecgdb_channel_results.csv` — median 1.000000, two rows below 0.95 | DERIVED |
| 38 records none below 90%; 13 mixed; 4 all fail | ecgn:210-213 | `results/nifecgdb_channel_results.csv` — per-record grouping; 38 + 13 + 4 = 55 | DERIVED |
| ecgca848 Abdomen_4 1.92%, Abdomen_3 4.41% | ecgn:215 | `results/nifecgdb_channel_results.csv` — 1.9157%, 4.4061% | ROUNDED |
| paired gap median +0.00 pp, range to +71.84 | ecgn:219 | — not a committed column; `abdominal_worse_records` 23 is | NOT FOUND `[AMB]` |
| thoracic FPs 7,247 / 648 / 67 | ecgn:228 | `results/nifecgdb_channel_results.csv` · `fp_in_gap`, `fp_genuine`, `fp_outside_span` — sums to 7,962 | DERIVED |
| abdominal FPs 13,405 / 4,022 / 114 | ecgn:229 | same columns — sums to 17,541 | DERIVED |
| % in gap 91.8 / 76.9 / overall 81.6 | ecgn:228-231 | in_gap/(in_gap+genuine): 7247/7895 = 91.79%; 13405/17427 = 76.92%; pooled 20652/25322 = 81.56% = committed `p4_frac_in_gap` | DERIVED |
| genuine FP per 1000: thoracic 6.75, abdominal 21.93 | ecgn:237 | 648/95,960 × 1000 = 6.753; 4,022/183,365 × 1000 = 21.934 | DERIVED |
| 3.2× worse | ecgn:238 | 21.934/6.753 = 3.248 | DERIVED |
| gap-corrected PPV 99.26% / 97.60% | ecgn:238 | TP/(TP+genuine+outside): 95661/96376 = 99.258%; 168185/172321 = 97.600% | DERIVED |
| P1 thoracic Se ≥98% → 99.69% HELD | ecgn:244 | `results/nifecgdb_summary.json` · `predictions.P1` = true | EXACT |
| P2 ≥75% of records → 42% DID NOT HOLD | ecgn:245 | `abdominal_worse_records` 23 / `paired_records` 55 = 41.8%; `predictions.P2` false | DERIVED |
| P3 PPV 92–94% → 92.32 / 90.56 DID NOT HOLD | ecgn:246 | `results/nifecgdb_summary.json` · `predictions.P3` = false | EXACT |
| P4 ≥70% FPs in gaps → 81.6% HELD | ecgn:247 | `predictions.P4` true; `p4_frac_in_gap` 0.8155753889898112 | EXACT |
| abdominal ~118,700 channel-seconds | ecgn:261 | — requires per-record durations, not committed | NOT FOUND |
| ~276,900 fetal complexes at 140 BPM | ecgn:261 | 118,700 × 140/60 = 276,967; input uncommitted | DERIVED |
| 4,022 genuine FPs = 1.45% | ecgn:263 | `fp_genuine` 4,022 committed; 4022/276,900 = 1.4525% | DERIVED |
| failures at 23.9 … 37.1 weeks | ecgn:277 | `results/nifecgdb_channel_results.csv` · `gestation_weeks`, joined to failing channels | DERIVED |
| pooled 110 channels, 8,190 windows | ecgn:299 | `results/mhr_ecg_summary.json` · `nifecgdb_thoracic.n` = 8190 | EXACT |
| reference rate 93.5 BPM (75.0–213.9) | ecgn:301 | `results/mhr_ecg_summary.json` · `ref_bpm_mean` 93.53486634190939, `ref_bpm_range` | ROUNDED |
| median \|error\| 0.146 BPM | ecgn:302 | `results/mhr_ecg_summary.json` · `median_abs_bpm` 0.1458002245323513 | ROUNDED |
| MAE 0.939 BPM | ecgn:303 | `results/mhr_ecg_summary.json` · `mae_bpm` 0.9393003382263891 | ROUNDED |
| bias +0.689 BPM | ecgn:304 | `results/mhr_ecg_summary.json` · `bias_bpm` 0.6887781336591334 | ROUNDED |
| SD of error 8.877 BPM | ecgn:305 | `results/mhr_ecg_summary.json` · `sd_bpm` 8.877130763298913 | ROUNDED |
| 95% LoA [−16.71, +18.09] | ecgn:306 | `results/mhr_ecg_summary.json` · `loa_lower`, `loa_upper` | ROUNDED |
| max \|error\| 163.51 BPM | ecgn:307 | `results/mhr_ecg_summary.json` · `max_abs_bpm` 163.5093167701863 | ROUNDED |
| within 2 BPM 98.47% | ecgn:308 | `results/mhr_ecg_summary.json` · `within_2_bpm` 0.9847374847374848 | ROUNDED |
| within 5 BPM 99.38% | ecgn:309 | `results/mhr_ecg_summary.json` · `within_5_bpm` 0.9937728937728938 | ROUNDED |
| window coverage 99.7% | ecgn:310 | `results/mhr_ecg_summary.json` · `mean_window_coverage` 0.9972592149062738 | ROUNDED |
| per-channel MAE median 0.251, mean 1.451, max 119.39 | ecgn:312 | `results/mhr_ecg_channel_results.csv` · `mae_bpm` over 110 rows — 0.2512, 1.4510, 119.3950 | DERIVED |
| per-channel bands 100 / 3 / 3 / 3 / 1 | ecgn:316-320 | `results/mhr_ecg_channel_results.csv` · `mae_bpm` — <0.5: 100, 0.5–1: 3, 1–2: 3, 2–5: 3, >5: 1 | DERIVED |
| excluding failing channel: 0.343 BPM over 8,149 windows | ecgn:322 | `results/mhr_ecg_channel_results.csv` — window-weighted mean excluding the 119.395 row; 8,190 − 41 | DERIVED |
| device capture 3 windows, MAE 0.13, max 0.26, ref 62.2 | ecgn:326 | `results/mhr_ecg_summary.json` · `device_capture` — 3; 0.13365580325229112; 0.25826112782634425; 62.175465492218414 | ROUNDED |
| PPG pipeline 1.13 BPM MAE on BIDMC | ecgn:330 | `results/bidmc_eval.csv` — cross-phase citation, matches the recomputed mean | DERIVED |
| SD 60× the median error | ecgn:336 | 8.877130763298913 / 0.1458002245323513 = 60.9 | DERIVED |
| ecgca699 Thorax_2 MAE 119.39 BPM | ecgn:347 | `results/06_diagnose_ecgca699.json` · `window_mae_bpm` 119.39495908562915 | ROUNDED |
| 0.4% of windows carrying most of the pooled error | ecgn:348 | 41 affected windows / 8,190 = 0.50%; the 0.4% denominator is not identified | NOT FOUND `[AMB]` |
| reference median rate 108.7 BPM | ecgn:354 | `results/06_diagnose_ecgca699.json` · `median_bpm` 108.69565217391303 | ROUNDED |
| maximum rate 119.0 BPM | ecgn:355 | `results/06_diagnose_ecgca699.json` · `bpm_max` 119.04761904761905 | ROUNDED |
| lag-1 RR correlation −0.014 | ecgn:356 | `results/06_diagnose_ecgca699.json` · `record_lag1_corr` −0.01442102981547793 | ROUNDED |
| half-median cluster 0.00% | ecgn:357 | `results/06_diagnose_ecgca699.json` · `record_half_cluster` 0.0 | EXACT |
| Thorax_1 MAE 0.22 BPM | ecgn:369 | `results/mhr_ecg_channel_results.csv` [ecgca699.edf, Thorax_1] = 0.2216 | ROUNDED |
| mean ref/det ratio 1.00 and 0.48 | ecgn:369-370 | per-window ratios exist for the five worst windows only (0.389–0.409); channel means not committed | NOT FOUND `[AMB]` |
| worst window 13 ref, 28 det, 104.3 vs 267.9 BPM | ecgn:372 | `results/06_diagnose_ecgca699.json` · `worst_windows[0]` | ROUNDED |
| 638 unmatched detections | ecgn:378 | `results/06_diagnose_ecgca699.json` · `unmatched_detections` 638 | EXACT |
| median offset 211 ms | ecgn:378 | `results/06_diagnose_ecgca699.json` · `unmatched_offset_median_s` 0.211 | EXACT |
| IQR 6 ms | ecgn:378 | `results/06_diagnose_ecgca699.json` · `unmatched_offset_iqr_s` 0.006000000000000005 | ROUNDED |
| 38% of the RR interval | ecgn:378 | `results/06_diagnose_ecgca699.json` · `unmatched_frac_of_rr_median` 0.3826247689463955 | ROUNDED |
| intervals alternate ~211 ms and ~341 ms | ecgn:385 | `median_rr_s` 0.552 − 0.211 = 0.341 s | DERIVED |
| 60/0.224 = 268 BPM | ecgn:387 | matches the committed `det_bpm` 267.85714285714283 | DERIVED |
| Thorax_1 R≈3.4/T≈1.3 (0.38); Thorax_2 R≈2.5/T≈1.7 (0.68) | ecgn:390-391 | — amplitude ratios are printed, not written to the diagnostic JSON | NOT FOUND |
| five worst lag-1 −0.46, −0.09, −0.11, −0.39, −0.06 | ecgn:412 | `results/06_diagnose_ecgca699.json` · `worst_windows[].lag1_corr` | ROUNDED |
| 8.1 bug: thoracic Se 99.77%, abdominal 92.18% | ecgn:437 | `results/nifecgdb_channel_results.csv` · supplementary=False — 108 → 99.7667%; 199 → 92.1751% | DERIVED |
| channel counts 108/199 against 110/203 | ecgn:438 | `results/nifecgdb_channel_results.csv` · `supplementary` — 307 False, 6 True | DERIVED |
| 24 of 41 windows silently discarded | ecgn:443 | `windows_paired` 41 is committed; the 24 discarded is not | NOT FOUND `[AMB]` |
| pooled MAE moved 0.547 → 0.939 | ecgn:449 | 0.939 committed; the pre-correction 0.547 appears in no file | NOT FOUND `[AMB]` |
| summed error +3,224 reconciles to +3,226 | ecgn:450 | — intermediate reconciliation, not committed | NOT FOUND |
| median error moved 0.145 → 0.146 | ecgn:451 | 0.146 committed; the pre-correction 0.145 is not | NOT FOUND `[AMB]` |
| mechanism measured to 211 ms ± 3 ms | ecgn:486 | IQR 0.006 s → ±3 ms about the 0.211 s median | DERIVED |
| 55 records, 47,980 beats, 313 channel-records | ecgp:35-47 | `results/nifecgdb_summary.json` — same values as the notes | DERIVED |
| total duration 549.7 min; length 114–2780 s | ecgp:41-42 | — no duration column committed | NOT FOUND |
| annotations maternal, median 92.7 BPM | ecgp:74 | — `03_explore_nifecgdb.py` prints only | NOT FOUND |
| median 7.33% missing, worst 9.18% | ecgp:83 | — same | NOT FOUND |
| 5.84% pooled | ecgp:84 | — same | NOT FOUND |
| ecgca711: 110 flagged, 107 between 1.8× and 2.2× | ecgp:87 | — same | NOT FOUND |
| 15/55 sub-200 ms; ecgca998 10.12%; rest 0.07–3.42% | ecgp:97-99 | — same | NOT FOUND |
| p21 best in fpcgdb: 0.7631 conf, 97.41% reliable | pcg:30 | `results/fpcgdb_eval.csv` — 0.7630745554720815; 0.9741235392320534 | ROUNDED |
| reliable-record rate 58% on real data | pcg:43 | `results/fpcgdb_eval.csv` · `reliable` — 15/26 = 57.7% | DERIVED |
| ECG 99.69% Se over 95,960 beats | pcg:41 | `results/nifecgdb_summary.json` · `thoracic` — cross-phase citation | ROUNDED |
| UC F1 77.32% | pcg:42 | `results/uc_evaluation.csv` · TEST — cross-phase citation | DERIVED |
| 1 kHz = 15,000 B/s against 11,520 available | pcg:113 | 15 × 1000 = 15,000; 115200/10 = 11,520 | DERIVED |
| 500 Hz = 7,500 B/s, ~65% utilisation | pcg:114 | 7,500/11,520 = 65.1% | DERIVED |
| samples 45,005 / 45,005 / 45,004 | pcg:132 | `results/pcg_char_baseline_*.json` · `timing_full.n` | EXACT |
| run C 45,001 samples | pcg:132 | `results/04_loopback_20260901_130246.json` records fs and malformed but no sample count | NOT FOUND `[AMB]` |
| dt mean 2000.000 µs | pcg:133 | `results/pcg_char_baseline_*.json` · `timing_full.dt_mean_us` — 1999.999955559506 ×2, 1999.999977779259 | ROUNDED |
| dt std 0.61 / 0.62 / 0.14 µs | pcg:134 | `timing_full.dt_std_us` — 0.6053930570409967 / 0.6220834520622113 / 0.13863960592992833 | ROUNDED |
| dt range 1998–2001 / 1998–2001 / 1999–2001 | pcg:135 | `timing_full.dt_min_us`, `dt_max_us` | EXACT |
| gaps >2× = 0; malformed = 0 | pcg:136-137 | `n_gaps_2x`, `malformed_lines` — all zero | EXACT |
| 40,003 intervals of exactly 2000 µs, std 0.0 | pcg:140 | `results/pcg_char_baseline_20260831_185145.json` · `timing_tail` — n 40003, dt_std 0.0 | EXACT |
| fs_achieved 500.0000111 Hz, error 2×10⁻⁶ % | pcg:149 | `results/pcg_char_baseline_20260831_183354.json` — 500.0000111101238; 2.2220247615223343e-06 | EXACT |
| ±50 ppm → ±0.025 Hz, ±0.007 BPM at 140 | pcg:153 | 500 × 50e-6 = 0.025; 140 × 50e-6 = 0.007 | DERIVED |
| full vs last-80 s agree to under 1% | pcg:158 | `signal_full` vs `signal_tail` — baseline 1 std 3.8185 vs 3.7837 = 0.91% | DERIVED |
| DC bias 1911.05 / 1910.16 / 1911.21 | pcg:169 | `signal_full.mean` — 1911.0510387734696 / 1910.164581713143 / 1911.2082037152254 | ROUNDED |
| std 3.82 / 4.11 / 4.00 counts | pcg:170 | `signal_full.std` — 3.8185071514431517 / 4.107222897533617 / 3.995106138754428 | ROUNDED |
| p2p 120 / 131 / 122 counts | pcg:171 | `signal_full.p2p` | EXACT |
| bare ADC floor 3.92 counts (toco) | pcg:178 | `hardware/mrisr/results/mrisr_*.json` · `baselines.uc_std` = 3.92 | EXACT |
| 47.80 Hz in all three baselines | pcg:194 | `spectrum_full.peak_inband_hz` — 47.80000106212784 ×2, 47.80000053107571 | ROUNDED |
| strength +12.1 / +8.6 / +12.4 dB | pcg:195 | `quick_check_lines.py` stdout; committed `db_over_bg` is +0.11 to +1.92 at 50/100/150 Hz — a different quantity | NOT FOUND |
| share of in-band power 4.07 / 2.58 / 4.43% | pcg:196 | `quick_check_lines.py` stdout; persists nothing (recorded as §11.5) | NOT FOUND |
| 50 Hz measures +1.6 to +4.1 dB | ctg:297 | `spectrum_full.mains[].db_over_bg` — committed 50 Hz values are +1.55 / +0.11 / +1.07; mrisr files carry no spectral block | NOT FOUND `[AMB]` |
| G34 vs G39 mean 1907.3 / 1907.2 | pcg:244 | `02_capture.py` persists `signal_full.mean`, but no G34 capture JSON was committed | NOT FOUND `[AMB]` |
| mean std 3.92 / 3.43 | pcg:245 | — same, no committed G34 capture | NOT FOUND `[AMB]` |
| mean p2p 57.0 / 51.0 | pcg:246 | — same | NOT FOUND `[AMB]` |
| mean p2p/std 14.5 / 14.9 | pcg:247 | `hardware/mrisr/results/mrisr_*.json` · `baselines.pcg_ratio` = [14.5, 14.9] | EXACT |
| quietest second std 1.70 / 1.24; p2p 20 / 12 | pcg:248-249 | — per-second statistics are not among the `02_capture.py` outputs | NOT FOUND |
| ratio moved by under 3% | pcg:251 | (14.9 − 14.5)/14.5 = 2.76% | DERIVED |
| noise runs 87 windows each | pcg:280 | — `quick_check_lines.py` runs the frozen detector and prints; persists nothing | NOT FOUND |
| mean confidence 0.064 / 0.043 / 0.066 | pcg:281 | — same | NOT FOUND |
| windows ≥ 0.45: 0.0% all three | pcg:282 | — same | NOT FOUND |
| BPM mean/std 138.2/24.0, 140.2/23.7, 150.5/22.4 | pcg:284 | — same | NOT FOUND |
| lag search 166–271 samples → 110.7–180.7 BPM | pcg:290-291 | 30000/271 = 110.70; 30000/166 = 180.72 | DERIVED |
| uniform-lag mean 140.1 BPM, std 20.0 | pcg:292 | — analytic result over the lag range, not a measurement | NOT FOUND |
| occlusion detector conf 0.809, 97.7%, BPM 175.7/4.1 | pcg:295-297 | — `quick_check_lines.py` stdout | NOT FOUND |
| dominant-line share 97.91% | pcg:333 | `results/pcg_char_occluded_20260831_191155.json` · `spectrum_full.mains[0].pct_of_inband` 97.91427611706176 | ROUNDED |
| quiet floor 2.58–4.43% | pcg:334 | — `quick_check_lines.py` stdout; the gap §11.5 describes | NOT FOUND |
| occlusion 50.40 Hz at +57.1 dB | pcg:354 | `pcg_char_occluded_*.json` · `peak_inband_hz` 50.40000055996267; `db_over_bg` 57.05785277778271 | ROUNDED |
| occlusion std 939.01 counts | pcg:353 | `pcg_char_occluded_*.json` — 936.8395155719737 is the tail; the cited value is the full-record std | ROUNDED `[AMB]` |
| harmonics at 100.60, 151.00, 198.20 Hz | pcg:356 | — `02_capture.py` records nominal 50/100/150 Hz centres, never measured peak positions | NOT FOUND |
| 235× the noise floor | pcg:352 | 939.01 / 3.99 (mean baseline std) = 235.3 | DERIVED |
| six bursts at 3.98 ± 0.04 s against manifest 4.00 s | pcg:372 | `results/quick_check_sweep.json` · `bursts[].t0` spacings 3.9, 4.0, 4.0, 4.0, 4.0; `stimulus_manifest.json` · `tones.layout` start_s 1,5,9,13,17,21 | DERIVED |
| 50 Hz burst 3.88% → 0.0011% on settling split | pcg:375 | `pcg_char_sweep_sub_*.json` · `spectrum_full`/`spectrum_tail` `mains[0].pct_of_inband` — 3.882229109342931, 0.0011155828814281815 | ROUNDED |
| 25 Hz: −25.7 dB, 39.3% harmonic, +0.6 ramp | pcg:379 | `results/quick_check_sweep.json` · `bursts[0]` | ROUNDED |
| 50 Hz: −11.3 dB, 1.3%, −0.1 | pcg:380 | `results/quick_check_sweep.json` · `bursts[1]` | ROUNDED |
| 75 Hz: −2.0 dB, 0.2%, +2.0 | pcg:381 | `results/quick_check_sweep.json` · `bursts[2]` | ROUNDED |
| 100 Hz: 0.0 dB, 0.2%, +0.1 | pcg:382 | `results/quick_check_sweep.json` · `bursts[3]` | ROUNDED |
| 150 Hz: −8.9 dB, 0.0%, +0.0 | pcg:383 | `results/quick_check_sweep.json` · `bursts[4]` | ROUNDED |
| 200 Hz: −18.5 dB, 0.0%, −0.3 | pcg:384 | `results/quick_check_sweep.json` · `bursts[5]` | ROUNDED |
| 75 Hz only 6.6 dB below the fundamental | pcg:387 | — harmonic-level breakdown is not a sweep JSON field | NOT FOUND |
| weighted transmission 7.3% (−11.4 dB) | pcg:396 | a literal in `04_loopback.py`'s docstring, not a results file; 10·log10(0.073) = −11.37 | NOT FOUND |
| band shares of original 80.7 / 18.8 / 0.2% | pcg:401-403 | — no committed file holds the segment's band decomposition | NOT FOUND |
| band shares surviving 24.3 / 73.5 / 1.9% | pcg:401-403 | — same | NOT FOUND |
| reference segment 19,980 samples, seconds 60–120 | pcg:418 | `stimulus/stimulus_manifest.json` · `fpcg.n_samples` 19980, `start_s` 60.0, `dur_s` 60.0 | EXACT |
| zero non-finite samples | pcg:418 | `stimulus/stimulus_manifest.json` · `fpcg.n_nonfinite` 0 | EXACT |
| quantised in steps of 200 | pcg:423 | — the manifest carries no quantisation fields | NOT FOUND |
| raw span −107 to +97, 204 counts, ~7.7 bits | pcg:424 | — same; log2(204) = 7.67 checks internally | NOT FOUND |
| band energy 20.7 / 64.1 / 14.7 / 0.2% | pcg:434-437 | — not in the manifest or any results file | NOT FOUND |
| centroid 42.0 Hz, median 38.5 Hz, 99.8% below 100 Hz | pcg:439 | — same | NOT FOUND |
| crest factor 6.27 (15.9 dB) | pcg:440 | — same; 20·log10(6.27) = 15.94 checks internally | NOT FOUND |
| reimplementation matches stored 136.8493 / 0.7631 / 0.9741 / True | pcg:457-460 | `results/fpcgdb_eval.csv` [fetal_PCG_p21_GW_39] | EXACT |
| run A: 56 windows, 100%, 134.09 BPM, 0.839, std 2.42 | pcg:472-477 | `results/03_reference_*.json` · `run_A` | ROUNDED |
| run B: 133.93 BPM, 0.840, std 2.46 | pcg:472-477 | `results/03_reference_*.json` · `run_B` | ROUNDED |
| run C: 133.93 BPM, 0.808, std 2.18 | pcg:472-477 | `results/04_loopback_*.json` · `run_C_primary` | ROUNDED |
| run C tail: 134.23 BPM, 0.820, std 2.20 | pcg:472-477 | `results/04_loopback_*.json` · `run_C_secondary` | ROUNDED |
| A→B costs 0.17 BPM and 0.002 confidence | pcg:479 | 134.09395973 − 133.92857143 = 0.16539; 0.84046103 − 0.83868957 = 0.00177 | DERIVED |
| 0.98 BPM/step at 333 Hz, 0.65 at 500 Hz | pcg:481 | lag-quantisation arithmetic from the detector's lag grid | DERIVED |
| B→C BPM difference 0.00 against 3.0 tolerance | pcg:485 | both report 133.92857142857142 exactly | DERIVED |
| 30000/224 = 133.92857; adjacent 134.53 and 133.33 | pcg:488-490 | 30000/223 = 134.529; 30000/225 = 133.333 | DERIVED |
| (30000/223 + 30000/224)/2 = 134.2288597 | pcg:491 | reproduces the committed `run_C_secondary.estimated_bpm` 134.2288597053171 | DERIVED |
| peak 575 counts (target 400–900), RMS 66.4, SNR 24.4 dB | pcg:495 | `results/level_finding_log.csv` · `peak_excursion`, `rms`, `snr_db`; verdict GOOD | EXACT |
| crest factor 8.66 → 7.36 | pcg:499 | `level_finding_log.csv` · `crest_factor` 8.66; `04_loopback_*.json` · `level.crest` 7.362697162928811 | ROUNDED |
| run C SNR 26.8 dB, a 2.4 dB rise | pcg:496 | `04_loopback.json` has no SNR field; 26.8 − 24.4 = 2.4 checks internally | NOT FOUND `[AMB]` |
| 7 samples within 90% of peak; peak 642 | pcg:497-498 | `results/04_loopback_*.json` · `level.n_within_90pct_peak` 7, `level.peak` 642.0562431945957 | EXACT |
| 96% of reference confidence | pcg:509 | 0.8076714135146965 / 0.8404610290857341 = 96.1% | DERIVED |
| 92.7% of in-band energy lost | pcg:510 | — a literal in `04_loopback.py`'s print, complement of the uncommitted 7.3% | NOT FOUND |
| seam: primary and tail agree to 0.30 BPM | pcg:506 | 134.2288597053171 − 133.92857142857142 = 0.30029 | DERIVED |
| CV values 0.572, 0.802, 0.548 against Rayleigh 0.523 | pcg:501 | — `quick_check_lines.py` stdout; √(4/π−1) = 0.5227 is analytic | NOT FOUND |
| bare on chest USB: +46.8 dB, 97.4%, std 210.4 | pcgc:41 | `results/pcg_char_chest_bare_p68_*.json` — 46.750248413951496; 97.40777421574914; 210.36782439674286 | ROUNDED |
| chamber on chest USB: +39.8 dB, 94.9%, std 260.9 | pcgc:42 | `results/pcg_char_chest_cap_p66_*.json` — 39.77924232719606; 94.91052512119644; 260.90766240629677 | ROUNDED |
| chamber battery: +5.0 dB, 1.86%, std 54.7 | pcgc:43 | `results/pcg_char_chest_cap_batt_p70_*.json` — 4.991607303949451; std 54.710042372898066. The 1.86% is `mains_total_pct_inband` (1.8551), not the 50 Hz `pct_of_inband` (1.5446) used in the two rows above | ROUNDED `[AMB]` |
| collapsed by roughly 40 dB | pcgc:48 | 46.75 − 4.99 = 41.8 dB (bare); 39.78 − 4.99 = 34.8 dB (chamber) | DERIVED |
| ~95–97% mains | pcgc:46 | `pct_of_inband` 94.91% and 97.41% | ROUNDED |
| quiet-room baselines 3.82 / 4.11 / 4.00 | pcgc:60 | `results/pcg_char_baseline_*.json` · `signal_full.std` | ROUNDED |
| near chest, no contact: 3.92 counts | pcgc:61 | `results/pcg_char_chest_air_batt_p74_*.json` · `signal_full.std` 3.918854349682361 | ROUNDED |
| bare capsule on chest: 50.1 counts | pcgc:62 | `results/pcg_char_chest_bare_batt_p72_*.json` · `signal_full.std` 50.13436438182181 | ROUNDED |
| chamber on chest: 54.7 counts | pcgc:63 | `results/pcg_char_chest_cap_batt_p70_*.json` · `signal_full.std` 54.710042372898066 | ROUNDED |
| envelope peaks ~0.68 s, CV 0.05, 100% reliable | pcgc:77 | — `quick_check_chest.py` prints only | NOT FOUND |
| peak finder distance 0.4 s; systole ~0.3 s | pcgc:79 | — script constant plus a physiological reference | NOT FOUND |
| D1 count ratio 2.96 / 1.85 / 2.74 | pcgc:94 | — `quick_check_doublet.py` prints only | NOT FOUND |
| D2 lag-1 autocorr −0.082 / −0.402 / −0.157 | pcgc:95 | — same | NOT FOUND |
| D3 short fraction 0.373 / 0.365 / 0.344 | pcgc:96 | — same | NOT FOUND |
| D4 cycle vs coarse Δ29.8 / Δ0.8 / Δ12.3 | pcgc:97 | — same | NOT FOUND |
| short 0.252 s, long 0.438 s, cycle 0.690 s | pcgc:104 | — same | NOT FOUND |
| 148 consecutive intervals | pcgc:105 | — same | NOT FOUND |
| systolic fraction 0.365 | pcgc:111 | — same. §12 calls this the strongest single number in the stage | NOT FOUND |
| coarse peak intervals 87.7 BPM | pcgc:124 | — same | NOT FOUND |
| envelope autocorrelation 88.2 BPM | pcgc:125 | — same | NOT FOUND |
| doublet cycle 87.0 BPM | pcgc:126 | — same | NOT FOUND |
| two radial counts of 71 and 70 BPM | pcgc:129 | — manual counts; no file records them | NOT FOUND |
| a 22% gap | pcgc:132 | (87.0 − 71)/71 = 22.5%; both inputs uncommitted | DERIVED |
| ctrl_cap_B both rails: min 0, max 4095, p2p 4095 | pcgc:159 | `results/pcg_char_ctrl_cap_B_p71_*.json` · `signal_full` | EXACT |
| collapses at ~0.4, 6, 12, 16–22, 28, 31 s | pcgc:161 | — time positions are not a `02_capture.py` output | NOT FOUND |
| ctrl_cap_A one 4095 spike | pcgc:163 | `results/pcg_char_ctrl_cap_A_p70_*.json` · `signal_full.max` 4095 with p2p 2908 | EXACT |
| bias 1857.7 → 1911.8 → 1987.0 → 1813.4 | pcgc:166 | `results/pcg_char_ctrl_*.json` · `signal_full.mean` | ROUNDED |
| against 1907–1911 stable | pcgc:167 | `hardware/mrisr/results/mrisr_*.json` · `baselines.pcg_bias` = [1907.0, 1911.0] | EXACT |
| verify_bench bias 1911.1, no rails, zero gaps/malformed | pcgc:170 | `results/pcg_char_verify_bench_*.json` — mean 1911.0993200906546; min 1712, max 2230; 0; 0 | ROUNDED |
| verify capture std 35.3 against a 4-count baseline | pcgc:176 | `results/pcg_char_verify_bench_*.json` · `signal_full.std` 35.29654215480487 | ROUNDED |
| energy piled below 30 Hz | pcgc:177 | `results/pcg_char_verify_bench_*.json` · `peak_inband_hz` 26.4; `inband_pct_of_total` 30.16% | DERIVED |
| chamber printed CONSISTENT on 21.2% of windows | pcgc:200 | — `quick_check_chest.py` prints only | NOT FOUND |
| ECG 31/31 at Se 100% | mri:28 | `hardware/ecg_bringup/results/*_device_eval.json` — 31, 31, 1.0 | EXACT |
| 99.69% Se over 95,960 beats | mri:28 | `hardware/ecg_bringup/results/nifecgdb_summary.json` · `thoracic` | ROUNDED |
| 500.0000 Hz, 0.6 µs jitter over 225k samples | mri:49 | `hardware/pcg_bringup/results/pcg_char_baseline_*.json` — fs 500.0000111; dt_std 0.605/0.622/0.139 µs; Σn ≈ 225k | DERIVED |
| UC timestamp offset 124 ms at t = 273,126,000 µs | mri:80 | — raw stream inspection; `captures/` is not `results/` | NOT FOUND |
| serial budget 7,500 + 76 = 7,576 B/s = 65.8% | mri:86-89 | 7,576/11,520 = 65.76% | DERIVED |
| adding ECG gives 11,818 B/s, 102.6% | mri:94 | 11,818/11,520 = 102.59% | DERIVED |
| shared-timestamp format ~8,770 B/s, 76.1% | mri:96 | 8,770/11,520 = 76.13% | DERIVED |
| resting UC ~1510 against toco's ~1085 | mri:127 | `mrisr_*.json` · `uc.mean` 1509.62 / 1507.89 / 1510.00; the ~1085 toco figure is not committed anywhere | ROUNDED `[AMB]` |
| PCG at 1909 and UC at 1510 simultaneously | mri:135 | `mrisr_paired_multi_*.json` · `pcg.mean` 1909.101097631427, `uc.mean` 1509.995577777778 | ROUNDED |
| four 90 s captures, ~180,000 PCG samples | mri:143 | 45,006 + 45,003 + 45,006 + 45,001 = 180,016 | DERIVED |
| PCG samples 45,006 / 45,003 / 45,006 / 45,001 | mri:147 | `mrisr_*.json` · `pcg.n`; `pcg_char_paired_single_*.json` · `timing_full.n`. The paired-single capture is stored under `hardware/pcg_bringup/results/` | EXACT |
| achieved fs 499.9889 / 499.9778 / 499.9889 / 499.9778 | mri:148 | `mrisr_*.json` · `fs_pcg`; paired_single · `fs_achieved` | ROUNDED |
| UC fs 4.0000 across three captures | mri:149 | `mrisr_*.json` · `fs_uc` — 3.9999999554317553 / 4.0 / 4.0 | ROUNDED |
| PCG:UC ratio 125.02 / 125.01 / 125.02 | mri:150 | `pcg.n` ÷ `uc.n` — 45,006/360 = 125.0167; 45,003/360 = 125.0083 | DERIVED |
| dt std 9.45 / 18.86 / 9.45 / 13.33 µs | mri:151 | `mrisr_*.json` · `dt_std_us`; paired_single · `timing_full.dt_std_us` | ROUNDED |
| malformed lines 0 / 1 / 0 / 2 | mri:152 | `mrisr_*.json` · `malformed`; paired_single · `malformed_lines` | EXACT |
| largest UC gap 0.250 s | mri:154 | `mrisr_*.json` · `largest_uc_gap_s` 0.2500010000000046 | ROUNDED |
| 4000/√45002 = 18.86 (reported 18.857) | mri:161 | √45,002 = 212.14; 4000/212.14 = 18.855 against the committed 18.8579 | DERIVED |
| 2000/√45005 = 9.43 (reported 9.4466) | mri:162 | 2000/212.14 = 9.428 against the committed 9.4466 | DERIVED |
| corruption rate ~1 in 20,000 lines | mri:184 | `mrisr_*.json` · `malformed` — 3 across ~180,000 lines ≈ 1 in 60,000; the notes' figure is a rounder upper estimate | DERIVED `[AMB]` |
| PCG std 3.23–3.28 measured | mri:205 | `mrisr_*.json` · `pcg.std` — 3.5733 (fan on), 3.2349, 3.2771. The range excludes the fan-on capture | NOT FOUND `[AMB]` |
| PCG p2p 88–93 | mri:206 | `mrisr_*.json` · `pcg.p2p` — 101, 93, 88. The range excludes the fan-on capture | NOT FOUND `[AMB]` |
| ratio 26.9–28.7 | mri:207 | `mrisr_*.json` · `pcg.ratio` — 26.853 / 28.265 / 28.749; this range does cover all three | ROUNDED |
| §6 baseline std 3.82–4.11, ratio 14.5–14.9 | mri:205-207 | `mrisr_*.json` · `baselines.pcg_std`, `.pcg_ratio` | EXACT |
| paired single ratio 26.26 | mri:219 | `pcg_char_paired_single_*.json` — p2p 91 ÷ std 3.4657526809555876 = 26.256 | DERIVED |
| multi std 3.277 vs single 3.466, −5.5% | mri:242 | 3.2770704079313018 and 3.4657526809555876; (3.277−3.466)/3.466 = −5.45% | DERIVED |
| multi p2p 88 vs single 91, −3.3% | mri:243 | 88 and 91; (88−91)/91 = −3.30% | DERIVED |
| bias 1909.1 vs 1910.8, −1.7 | mri:244 | 1909.101097631427 and 1910.7519832892604; difference −1.651 | ROUNDED |
| contamination limit 6.17 (PCG) | mri:269 | 1.5 × 4.11 = 6.165, from `baselines.pcg_std` and `preregistered.std_tolerance` | DERIVED |
| PCG bias range 1909.1–1910.6 | mri:270 | `mrisr_*.json` · `pcg.mean` — 1909.101 / 1910.294 / 1910.636; covers all three | ROUNDED |
| UC std 0.68–1.40 against baseline 3.92, limit 5.88 | mri:271 | `uc.std` 0.6787 / 0.8277 / 1.3951; 1.5 × 3.92 = 5.88 | ROUNDED |
| 3.92/√125 = 0.35 expected | mri:281 | 3.92/11.1803 = 0.3506 | DERIVED |
| UC drift over 90 s: −0.7 to +1.7 counts | mri:283 | the JSONs carry `uc.p2p` (4.97–6.83) but no drift measure | NOT FOUND `[AMB]` |
| run C reference 133.93 BPM / 0.808 | ctg:79 | `hardware/pcg_bringup/results/04_loopback_*.json` · `run_C_primary` | ROUNDED |
| cue lag measured at 86.35 s in the toco phase | ctg:119 | — not committed in any toco results file | NOT FOUND |
| amplitude-scaling table 2× / 20× / 100× / 1000× | ctg:136-139 | — computed for this pre-registration; no producing script writes it | NOT FOUND |
| null runs: 87 windows, BPM 146.34 / 129.30 / 150.00 | ctg:173-175 | — the mrisr result files hold no detector fields | NOT FOUND |
| std 41.93 / 46.37 / 53.36; conf 0.058 / 0.036 / 0.056 | ctg:173-175 | — same | NOT FOUND |
| max window confidence 0.298 / 0.244 / 0.240 | ctg:173-175 | — same | NOT FOUND |
| 0 of 261 windows clear 0.45 | ctg:177 | 261 = 3 × 87; neither input committed | DERIVED |
| serial utilisation 65.8% | ctg:258 | 7,576/11,520 = 65.76% — cross-phase | DERIVED |
| 47.80 Hz present in all three mrisr captures | ctg:259 | — mrisr result files contain no `spectrum_full`, no `peak_inband_hz`, no `db_over_bg` | NOT FOUND `[AMB]` |
| mrisr quiet floor 8.21 / 8.38 / 8.59% | ctg:306 | — marked in the prereg as "measured for this pre-registration"; no script persists it | NOT FOUND |
| occluded 97.91% | ctg:304 | `hardware/pcg_bringup/results/pcg_char_occluded_*.json` · `spectrum_full.mains[0].pct_of_inband` | ROUNDED |
| pcg_bringup quiet floor 2.58 / 4.07 / 4.43% | ctg:305 | — labelled in the prereg as `quick_check_lines.py` stdout, unpersisted | NOT FOUND |
| ECG serial constraint 102.6% | ctg:339 | 11,818/11,520 — cross-phase | DERIVED |
| 31/31 and 99.69% Se over 95,960 | ctg:340 | ECG bring-up results — cross-phase, both verified above | ROUNDED |

---

## 3. Not found, grouped by producing script

Every `NOT FOUND` figure above, with the script that computed it and whether
that script writes its output to disk. The gaps are not evenly distributed:
they cluster almost entirely in the diagnostic scripts that persist nothing.

### Prints only — no file output at all

| Figures | Cited at | Script |
|---|---|---|
| 47.80 Hz line strengths and in-band shares; §7.1 noise-run detector table; occlusion detector output; the Rayleigh CV comparison | pcg:195-196, 280-297, 501 | `quick_check_lines.py` — prints `pct_inband`, per-line dB, confidence and verdict. Only file output is a plot (line 302). |
| Doublet criteria D1–D4; interval clusters 0.252 / 0.438 / 0.690 s; 148 intervals; systolic fraction 0.365; the three rate routes 87.7 / 88.2 / 87.0 BPM | pcgc:94-126 | `quick_check_doublet.py` — prints only. |
| Coarse envelope peaks at ~0.68 s; interval CV 0.05; CONSISTENT at 21.2% of windows | pcgc:77, 200 | `quick_check_chest.py` — prints only. |
| Annotations are maternal at median 92.7 BPM; 7.33% and 5.84% incompleteness; 53/55; 15/55; ecgca998 at 10.12%; the 0.07–3.42% spread | ecgn:149-172, ecgp:74-99 | `03_explore_nifecgdb.py` — the only script in the ECG phase with no write call. |
| aromring unit test: 68.2 BPM, 0.92 confidence, 4 peaks | mhr:96 | `unit_test_aromring.py` — asserts and prints; no results file for this dataset. |

### Persists, but not these fields

| Figures | Cited at | Script |
|---|---|---|
| Toco creep +0.615 counts/s, +51 counts, +1.50%; noise floor 9.70 counts | toco:203-225 | `02_preload.py` — computes and writes `creep_cps`, `creep_total`, `creep_pct`, `noise_std` to `results/preload_<stamp>.csv`. **The committed rows are the pre-fix run**; see §4. |
| Detection thresholds, amplitude references, per-capture rates, dropout percentage | toco:310-331 | `04_detect.py` — writes only `cycle, ref_start_s, ref_end_s, det_start_s, det_end_s, cover, iou, matched`. |
| Capture timing, rest-phase baselines, HOLD means, trough margins | toco:87, 243-296, 432-465 | `01_capture.py` and `03_simulated_contractions.py` — write raw captures and cue files, no derived statistics. |
| Reference-segment quantisation, band energy, centroid, crest factor | pcg:423-440 | `01_stimulus.py` — persists `stimulus_manifest.json` (record identity, indices, resampling ratios), no spectral or quantisation fields. |
| Playback transmission 7.3% / 92.7% and the band-share inversion | pcg:396-403, 510 | `04_loopback.py` — the two percentages are hardcoded literals in its docstring and a print string. The script computes no spectral share. |
| G34-vs-G39 erratum table (mean, std, p2p, quietest second) | pcg:244-249 | `02_capture.py` — computes and persists exactly `signal_full.mean/std/p2p`, but **no G34 capture JSON was committed**. Only the ratio pair survives, carried into the mrisr baselines. |
| Occlusion harmonic frequencies 100.60 / 151.00 / 198.20 Hz | pcg:356 | `02_capture.py` — records `hz` as the nominal 50/100/150 Hz centre, never the measured peak position. |
| Smoke-test hum shares; clipping rates; detection margins 7.22e8 / 5.52e7; per-beat offsets; the 60.8 BPM rate summary; tolerance 38 samples | ecga:27-54, 190-197, 299-341, 406 | `02_detect_device.py` persists the scoring block; `01_annotate.py` persists the annotation set and manifest. The rest is printed. |
| NIFECGDB record durations and lengths | ecgn:136-137, ecgp:41-42 | `04_nifecgdb_evaluate.py` — writes 22 per-channel columns; duration is not among them. |
| Per-channel R/T amplitude ratios | ecgn:390-391 | `06_diagnose_window_bpm.py` — persists the phase measurement and worst-window array; ratios are printed. |
| Alternative matching criteria (any overlap, IoU ≥ 0.5); boundary-agreement medians | uc:28-52 | `evaluate_uc.py` — writes `uc_evaluation.csv` and `uc_detections.csv`; neither carries alternative-criterion scores or reference boundaries. |
| Per-event-type counts: bradycardia, tachycardia, acceleration, deceleration | uc:86-98 | `parse_ctu_annotations.py` — writes contraction rows only. The other four annotation rows are read and discarded. |

### No producer in the repository

| Figures | Cited at | Note |
|---|---|---|
| Every toco mains and band share: 0.03 / 0.01 / 0.74%; 4.07 / 0.86 / 3.36 / 12.36 / 76.53%; the 6–11% in 1–10 Hz; the 4 Hz decimated std 3.64 | toco:162-209 | No script in the toco phase computes a spectrum. `grep` for `welch`, `periodogram`, `fft` over `hardware/toco_bringup/*.py` returns nothing. |
| CTG null-run detector output (87 windows, BPM, std, confidence) and the mrisr line shares 8.21 / 8.38 / 8.59% | ctg:173-177, 306 | The prereg states these were "measured for this pre-registration". No committed script produces them, and `hardware/mrisr/results/*.json` contains no spectral or detector fields. |

---

## 4. Ambiguous, contradictory, or otherwise flagged

Described, not resolved.

**4.1 simfpcgdb is one record short of the notes.**
`results/simfpcgdb_eval.csv` holds 36 rows, 21 reliable. The FHR summary cites
37 records, 22 reliable (59%), and a mean confidence of 0.81 at −4.4 dB. No
`SNR-4_4dB` row exists in the file, and the highest committed confidence is
0.8000 at −6.6 dB. Three separate figures depend on the missing row. Whether
the record was dropped from the evaluation or never downloaded is not recorded
anywhere readable.

**4.2 MIT-BIH PPV 99.35% is not reproducible from the committed CSV.**
`hardware/ecg_bringup/notes/notes.md:64` cites Se 99.52%, PPV 99.35% over 46
MLII records. Sensitivity reproduces exactly (micro over MLII = 99.5165%). PPV
does not: micro 99.7215%, macro 99.6993%, all-48 micro 99.7219%, MLII minus
record 207 99.7185%. No subset tested lands on 99.35%. Not inferred.

**4.3 The committed toco preload CSV holds the numbers the notes retracted.**
`results/preload_20260807_170647.csv` records `creep_cps` 8.7835 and
`noise_std` 614.678. Toco §11 records that `02_preload.py` originally
"reported creep 14× too high and noise 60× too high" — and 8.7835/0.615 = 14.3,
614.678/9.70 = 63.4. So the one committed non-zero preload row **is** the
erroneous run, and the corrected §6/§7 figures were never written to disk. The
other two preload CSVs are all zeros.

**4.4 The mrisr §4.1 and §5 std/p2p ranges exclude one of their own captures.**
Both tables state PCG std 3.23–3.28 and p2p 88–93 as the multi-channel
measurements. The three committed mrisr JSONs read std 3.5733 / 3.2349 / 3.2771
and p2p 101 / 93 / 88 — the fan-on baseline falls outside both ranges. The ratio
range (26.9–28.7) and the bias range (1909.1–1910.6) in the same tables do cover
all three. Either the ranges describe a two-capture subset or the fan-on run was
meant to be excluded; the notes do not say which.

**4.5 One column in §14.2 mixes two different JSON keys.**
The "mains share of in-band" column reads 97.4% and 94.9% for the two USB rows —
those are `spectrum_full.mains[0].pct_of_inband`, the 50 Hz line alone. The
battery row reads 1.86%, which is `mains_total_pct_inband`, the 50+100+150 Hz
sum; its 50 Hz value is 1.54%. The three rows are not the same quantity. This is
the same class of defect the phase recorded as §11.5.

**4.6 The CTG pre-registration's +1.6 to +4.1 dB has no committed source.**
`hardware/ctg/notes/prereg.md:297-298` states that 50 Hz measures only +1.6 to
+4.1 dB over background on this hardware. The committed 50 Hz `db_over_bg` for
the three pcg baselines is +1.55 / +0.11 / +1.07 — range +0.11 to +1.55. If the
sentence refers to the mrisr captures instead, those result files carry no
spectral block at all: no `db_over_bg`, no `pct_of_inband`, no
`peak_inband_hz`. The neighbouring claim at line 259 that 47.80 Hz is "present
in all three mrisr captures" is unverifiable for the same reason.

**4.7 Two UC figures disagree with the committed contraction table.**
The row-mapping table cites a start-to-start median of 141 s (4.3 per 10 min);
`ctu_contractions.csv` gives 135.125 s (4.44 per 10 min), and the notes' own
line 197 cites 135 s for the same quantity. The same table cites an IQR of
46.5–76.5 s; linear-interpolated quartiles over all 7,015 durations give 46.0
and 76.0. Both may predate the final parse; nothing in the repository says.

**4.8 Superseded ECG figures are cited alongside current ones.**
§8.1 and §8.2 of the NIFECGDB notes deliberately record pre-correction values
(Se 99.77%, abdominal 92.18%, MAE 0.547, median 0.145). The first two are
recoverable from the committed CSV by reproducing the bug's filter; 0.547 and
0.145 are in no file. This is intended — the notes are recording their own
errors — but a reader cannot distinguish a recoverable superseded figure from an
unrecoverable one without doing this audit.

**4.9 Run C's sample count is the only §3 cell without a source.**
`hardware/pcg_bringup/notes/notes.md:132` gives run C 45,001 samples.
`04_loopback_20260901_130246.json` records `fs_achieved` 500.0 and `malformed` 0
but carries no sample count or timing block. The identical figure 45,001 appears
in `pcg_char_paired_single_20260906_010212.json`, a different capture five days
later. That is not treated as the source.

---

## 5. Method

- Notes read in full: eleven files, 3,783 lines.
- Results files examined: 50, across `results/`, `hardware/*/results/`.
- Aggregates (`bidmc_eval`, `mitdb_e1_eval`, `uc_evaluation`, `uc_detections`,
  `ctu_contractions`, `ctu_ann_manifest`, `ctu_uc_quality`,
  `nifecgdb_channel_results`, `mhr_ecg_channel_results`) recomputed by reading
  the committed CSVs directly. Percentiles use linear interpolation, matching
  numpy's default, except where noted.
- Script persistence classified by locating every `json.dump`, `to_csv`,
  `savetxt`, `csv.writer` and file-open-for-write call across all project
  scripts.
- No pipeline script was executed. No frozen file was read except to confirm
  constants. Nothing was modified.
