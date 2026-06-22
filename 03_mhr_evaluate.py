"""
03_mhr_evaluate.py — Evaluate the MHR detector against BIDMC ground truth.

Unlike the FHR evaluator (no ground truth -> physiological plausibility +
self-consistency), BIDMC gives us per-second reference heart rate, so here we
compute true error (MAE).

Two references are scored independently:
  - PULSE : pulse rate derived from the PPG (same modality -> isolates the
            algorithm's own error)
  - HR    : heart rate derived from the ECG (independent modality -> the
            harder cross-sensor test; can diverge during arrhythmia)

Alignment: detector estimates are at window centres (~every HOP_S seconds);
references are at 1 Hz. Each window centre is matched to the nearest 1 Hz
reference sample.

Scoring: MAE is computed over HIGH-CONFIDENCE windows only, where the chosen
reference is finite. Rationale: a low-confidence window is one the device
would NOT display, so it shouldn't count toward "when we show a number, how
accurate is it." PULSE has occasional NaN gaps, so its window count can be
slightly lower than HR's.
"""
import os
import csv
import importlib.util

import numpy as np
import wfdb


# ----------------------------------------------------------------------
# Import the detector (filename starts with a digit -> load via importlib)
# ----------------------------------------------------------------------
def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_HERE = os.path.dirname(os.path.abspath(__file__))
_mhr = _load_module(os.path.join(_HERE, '02_mhr_detector.py'), 'mhr_detector')
detect_mhr = _mhr.detect_mhr
find_channel = _mhr.find_channel
CONF_THRESH = _mhr.CONF_THRESH
PPG_CHANNEL = _mhr.PPG_CHANNEL


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
DATASET = 'bidmc'
DATA_DIR = os.path.join('data', DATASET)
RESULTS_DIR = 'results'
os.makedirs(RESULTS_DIR, exist_ok=True)
RESULTS_PATH = os.path.join(RESULTS_DIR, f'{DATASET}_eval.csv')

MAE_TARGET = 3.0             # project target: MAE < 3 BPM


# ----------------------------------------------------------------------
# Evaluation
# ----------------------------------------------------------------------
def ref_at(window_times, ref, ref_fs):
    """Sample the 1 Hz reference at each window-centre time (nearest sample)."""
    idx = np.clip(np.round(window_times * ref_fs).astype(int), 0, len(ref) - 1)
    return ref[idx]


def evaluate_record(record_name):
    """Run the detector on one record and score it against HR and PULSE."""
    rec = wfdb.rdrecord(os.path.join(DATA_DIR, record_name))
    pidx = find_channel(rec, PPG_CHANNEL)
    if pidx is None:
        return None

    raw = rec.p_signal[:, pidx]
    fs = rec.fs
    result = detect_mhr(raw, fs)

    times = result['window_times']
    bpms = result['bpms']
    high = result['confidence'] >= CONF_THRESH

    row = {
        'record': record_name,
        'fs': int(fs),
        'n_windows': int(len(bpms)),
        'est_median_bpm': float(np.nanmedian(bpms)) if np.isfinite(bpms).any() else np.nan,
        'high_conf_fraction': round(float(result['high_conf_fraction']), 3),
        'reliable': int(result['reliable']),
    }

    # Load the matching numerics record (XXn) for the references.
    num_path = os.path.join(DATA_DIR, record_name + 'n')
    if not os.path.exists(num_path + '.hea'):
        for key in ('mae_hr', 'mae_hr_n', 'mae_pulse', 'mae_pulse_n', 'coverage'):
            row[key] = np.nan
        return row
    num = wfdb.rdrecord(num_path)

    for ref_name, mae_key, n_key in (('HR', 'mae_hr', 'mae_hr_n'),
                                     ('PULSE', 'mae_pulse', 'mae_pulse_n')):
        ridx = find_channel(num, ref_name)
        if ridx is None:
            row[mae_key], row[n_key] = np.nan, 0
            continue
        ref_vals = ref_at(times, num.p_signal[:, ridx], num.fs)
        mask = high & np.isfinite(bpms) & np.isfinite(ref_vals)
        if mask.sum() == 0:
            row[mae_key], row[n_key] = np.nan, 0
        else:
            row[mae_key] = round(float(np.mean(np.abs(bpms[mask] - ref_vals[mask]))), 3)
            row[n_key] = int(mask.sum())

    # Coverage: fraction of all windows that were scorable against HR.
    row['coverage'] = round(row.get('mae_hr_n', 0) / len(bpms), 3) if len(bpms) else 0.0
    return row


def _fmt(x, width, prec=2):
    """Format a possibly-NaN float for the table."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return f"{'--':>{width}}"
    return f"{x:>{width}.{prec}f}"


def main():
    headers = sorted(f[:-4] for f in os.listdir(DATA_DIR) if f.endswith('.hea'))
    sig_recs = [h for h in headers if not h.endswith('n')]
    if not sig_recs:
        print(f"No signal records in {DATA_DIR}/. Run the BIDMC download first.")
        return

    rows = []
    hdr = (f"{'Record':<10}{'fs':>5}{'n_win':>7}{'est BPM':>9}"
           f"{'MAE_HR':>8}{'MAE_PUL':>9}{'cov':>6}{'flag':>6}")
    print(hdr)
    print('-' * len(hdr))

    for name in sig_recs:
        res = evaluate_record(name)
        if res is None:
            continue
        rows.append(res)
        print(f"{res['record']:<10}{res['fs']:>5}{res['n_windows']:>7}"
              f"{_fmt(res['est_median_bpm'], 9, 1)}"
              f"{_fmt(res['mae_hr'], 8)}{_fmt(res['mae_pulse'], 9)}"
              f"{100*res['coverage']:>5.0f}%"
              f"{'OK' if res['reliable'] else 'LOW':>6}")

    print('-' * len(hdr))

    # ---- Aggregates ----
    reliable = [r for r in rows if r['reliable']]
    mae_hr = [r['mae_hr'] for r in reliable if np.isfinite(r['mae_hr'])]
    mae_pul = [r['mae_pulse'] for r in reliable if np.isfinite(r['mae_pulse'])]

    print(f"Records processed : {len(rows)}  "
          f"({len(reliable)} reliable, {len(rows) - len(reliable)} flagged LOW)")
    if mae_pul:
        print(f"Reliable records, mean MAE vs PULSE : {np.mean(mae_pul):.2f} BPM "
              f"(median {np.median(mae_pul):.2f})")
    if mae_hr:
        print(f"Reliable records, mean MAE vs HR    : {np.mean(mae_hr):.2f} BPM "
              f"(median {np.median(mae_hr):.2f})")
        hit = sum(v < MAE_TARGET for v in mae_hr)
        print(f"Reliable records meeting MAE < {MAE_TARGET:.0f} vs HR : "
              f"{hit}/{len(mae_hr)}")

    # ---- Write CSV ----
    fieldnames = ['record', 'fs', 'n_windows', 'est_median_bpm',
                  'mae_hr', 'mae_hr_n', 'mae_pulse', 'mae_pulse_n',
                  'coverage', 'high_conf_fraction', 'reliable']
    with open(RESULTS_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, '') for k in fieldnames})
    print(f"\nResults written to {RESULTS_PATH}")


if __name__ == '__main__':
    main()
