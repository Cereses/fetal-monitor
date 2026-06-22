"""03_evaluate.py
Run the FHR detector on every record in the dataset and compute simple
quality metrics. Use this to tune filter cutoffs, window sizes, etc.

If the dataset has annotation files with beat positions, ground-truth BPM is
computed from those. Otherwise we fall back to the median estimate as a
self-consistency check.
"""
import os
import numpy as np
import wfdb
import csv

from importlib import import_module
fhr_mod = import_module('02_fhr_detector')

DATASET = 'fpcgdb'
DATA_DIR = os.path.join('data', DATASET)

RESULTS_DIR = 'results'
os.makedirs(RESULTS_DIR, exist_ok=True)
RESULTS_PATH = os.path.join(RESULTS_DIR, f'{DATASET}_eval.csv')

def ground_truth_bpm(record_path, fs):
    """Try common annotation extensions; return mean BPM or None."""
    for ext in ('qrs', 'fqrs', 'atr'):
        try:
            ann = wfdb.rdann(record_path, ext)
        except Exception:
            continue
        if len(ann.sample) < 2:
            continue
        rr_intervals_sec = np.diff(ann.sample) / fs
        bpm = 60.0 / rr_intervals_sec
        # Trim physiologically implausible RR intervals
        bpm = bpm[(bpm > 80) & (bpm < 220)]
        if len(bpm):
            return float(np.mean(bpm))
    return None


def evaluate_record(data_dir, record_name,
                    conf_threshold=fhr_mod.CONFIDENCE_THRESHOLD,
                    min_reliable_fraction=0.5):
    path = os.path.join(data_dir, record_name)
    record = wfdb.rdrecord(path)
    raw = record.p_signal[:, 0]
    fs = record.fs

    result = fhr_mod.detect_fhr(raw, fs)
    if not len(result['bpms']):
        return None

    bpms = np.asarray(result['bpms'])
    conf = np.asarray(result['confidence'])

    high_mask = conf >= conf_threshold
    n_total = len(conf)
    n_high = int(high_mask.sum())
    high_frac = n_high / n_total if n_total else 0.0

    if high_mask.any():
        estimated = float(np.median(bpms[high_mask]))
    else:
        estimated = float(np.median(bpms))

    truth = ground_truth_bpm(path, fs)
    return {
        'record': record_name,
        'fs': fs,
        'duration_s': len(raw) / fs,
        'estimated_bpm': estimated,
        'truth_bpm': truth,
        'mean_confidence': float(np.mean(conf)),
        'high_conf_fraction': high_frac,
        'reliable': high_frac >= min_reliable_fraction,  
    }


def main():
    headers = sorted(f for f in os.listdir(DATA_DIR) if f.endswith('.hea'))
    if not headers:
        print(f"No records in {DATA_DIR}/. Run download_data.py first.")
        return

    rows = []
    for h in headers:
        name = os.path.splitext(h)[0]
        try:
            r = evaluate_record(DATA_DIR, name)
            if r:
                rows.append(r)
        except Exception as e:
            print(f"  [{name}] error: {e}")

    # Header row
    print(f"\n{'Record':<35} {'fs':>5} {'dur(s)':>7} {'est BPM':>8} "
          f"{'true BPM':>10} {'conf':>5} {'reliable%':>10} {'flag':>5}")
    print('-' * 90)

    abs_errors = []
    abs_errors_reliable = []
    for r in rows:
        truth_str = f"{r['truth_bpm']:.1f}" if r['truth_bpm'] is not None else '--'
        flag = 'OK' if r['reliable'] else 'LOW'
        print(f"{r['record']:<35} {r['fs']:>5} {r['duration_s']:>7.1f} "
              f"{r['estimated_bpm']:>8.1f} {truth_str:>10} "
              f"{r['mean_confidence']:>5.2f} {r['high_conf_fraction']:>9.0%} "
              f"{flag:>5}")
        if r['truth_bpm'] is not None:
            err = abs(r['estimated_bpm'] - r['truth_bpm'])
            abs_errors.append(err)
            if r['reliable']:
                abs_errors_reliable.append(err)

    print('-' * 90)
    if abs_errors:
        print(f"Mean absolute error (all records):       "
              f"{np.mean(abs_errors):.2f} BPM (n = {len(abs_errors)})")
        print(f"Max absolute error  (all records):       "
              f"{np.max(abs_errors):.2f} BPM")
        if abs_errors_reliable:
            print(f"Mean absolute error (reliable only):     "
                  f"{np.mean(abs_errors_reliable):.2f} BPM "
                  f"(n = {len(abs_errors_reliable)})")
            print(f"Max absolute error  (reliable only):     "
                  f"{np.max(abs_errors_reliable):.2f} BPM")
    else:
        print("No annotation files found in dataset; only self-consistency reported.")

    n_reliable = sum(1 for r in rows if r['reliable'])
    print(f"Records processed: {len(rows)} / {len(headers)}  "
          f"({n_reliable} reliable, {len(rows) - n_reliable} flagged LOW)")
    
    # Write results to CSV
    if rows:
        fieldnames = ['record', 'fs', 'duration_s', 'estimated_bpm',
                      'truth_bpm', 'mean_confidence',
                      'high_conf_fraction', 'reliable']
        with open(RESULTS_PATH, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in rows:
                writer.writerow(r)
        print(f"\nResults written to {RESULTS_PATH}")

if __name__ == '__main__':
    main()
