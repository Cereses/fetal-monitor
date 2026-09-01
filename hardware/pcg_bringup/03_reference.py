"""
03_reference.py -- runs A and B: the numbers run C gets judged against.

    RUN A   frozen pipeline on the p21 segment at native 333 Hz
    RUN B   frozen pipeline on the SAME segment resampled to 500 Hz

Neither touches hardware. A -> B isolates what the rate change alone costs, so
that when run C (device capture of that segment played acoustically) is scored,
any further loss is attributable to the acquisition chain rather than to
resampling.

WHY RUN A IS NOT THE PUBLISHED p21 FIGURE
The published 0.76 confidence / 97% reliable is for the FULL ~20 minute record.
shannon_energy_envelope normalises by the GLOBAL maximum of whatever array it
is handed, so scoring a 60 s excerpt and scoring the whole record are different
operations on different data. Run A must be computed on the exact segment that
will be played. Comparing run C to a whole-record figure would compare two
different things.

WHY THE SCORING RULE IS VALIDATED, NOT TRUSTED
The record-level verdict lives in 03_evaluate.py, but its evaluate_record()
reads through wfdb.rdrecord and cannot be pointed at a NumPy array. The rule
therefore has to be reimplemented here, and a reimplementation that drifts
would corrupt every comparison built on it.

So before either run, this script runs the frozen detector on the FULL p21
record, scores it with the reimplementation, and checks the result against the
row 03_evaluate.py already wrote into results/fpcgdb_eval.csv. Mismatch is a
hard stop. That converts "I think I copied the rule correctly" into a check
against the frozen phase's own output.

The rule, quoted from 03_evaluate.py:
    high_mask  = conf >= fhr_mod.CONFIDENCE_THRESHOLD
    high_frac  = n_high / n_total
    estimated  = median(bpms[high_mask]) if any else median(bpms)
    reliable   = high_frac >= 0.5
Note the fallback: with NO window above threshold it still returns a number,
flagged LOW. Given that pure noise was measured to yield BPM clustered near
140 BPM, that fallback deserves its own line in the notes.

PRE-REGISTERED PREDICTIONS -- recorded before either run

  P1  Run B agrees with run A to within 1.0 BPM.
      peak_lag is an integer, so BPM is quantised. Near 140 BPM at 333 Hz the
      lag is 142.7 and adjacent integer lags give 140.7 / 139.7, a step of
      0.98 BPM. At 500 Hz the lag is 214.3, giving 140.19 / 139.53, a step of
      0.65 BPM. Resolution IMPROVES with rate, so a disagreement larger than
      ~1 BPM is not explained by quantisation and means resampling changed the
      signal.

  P2  Run B confidence within 0.10 of run A. Polyphase resampling preserves
      the passband; the periodicity structure should survive it.

  P3  Both runs produce 56 windows (60 s, 4 s window, 1 s hop, at either rate).

  P4  If run A comes back LOW, the segment is unusable as a reference and a
      different one must be chosen ON STATED CRITERIA BEFORE run C -- never
      after seeing run C's result.

The frozen detector is IMPORTED BY PATH, never copied.

USAGE
    python 03_reference.py
    python 03_reference.py --npy stimulus/<other>_native333.npy
"""
import os
import sys
import csv
import glob
import json
import argparse
import importlib.util

import numpy as np
from scipy.signal import resample_poly

import matplotlib
matplotlib.use('Agg')          # MUST precede the frozen module's pyplot import
import matplotlib.pyplot as plt

import wfdb

# ------------------------------------------------------------------ config
STIM_DIR = 'stimulus'
RES_DIR = 'results'
PLOT_DIR = 'plots'

FROZEN_NAME = '02_fhr_detector.py'
EVAL_CSV_REL = os.path.join('results', 'fpcgdb_eval.csv')
DATA_REL = os.path.join('data', 'fpcgdb')
SEARCH_DEPTH = 5

FS_NATIVE = 333.0
FS_DEVICE = 500.0              # the achieved device rate, measured at
                               # 500.0000111 Hz over 45,005 samples

RECORD_OK_FRACTION = 0.50      # from 03_evaluate.py min_reliable_fraction

PRED_BPM_TOL = 1.0             # P1
PRED_CONF_TOL = 0.10           # P2


# ------------------------------------------------------------- resolution
def resolve_upward(relpath, depth=SEARCH_DEPTH):
    here = os.path.dirname(os.path.abspath(__file__))
    tried = []
    for _ in range(depth + 1):
        cand = os.path.join(here, relpath)
        tried.append(cand)
        if os.path.exists(cand):
            return cand
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    return None, tried


def need(relpath, what):
    r = resolve_upward(relpath)
    if isinstance(r, tuple) or r is None:
        print(f"could not find {what} ({relpath}). Looked in:")
        for t in (r[1] if isinstance(r, tuple) else []):
            print("   ", t)
        sys.exit(1)
    return r


def load_frozen():
    path = need(FROZEN_NAME, 'the frozen detector')
    spec = importlib.util.spec_from_file_location('frozen_fhr', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)      # side effect: creates plots/fpcgdb/
    print(f"frozen detector : {path}")
    return mod


# ------------------------------------------------------------------ rule
def score(bpms, conf, thr):
    """Reimplementation of 03_evaluate.py's record-level rule. Validated below."""
    bpms = np.asarray(bpms)
    conf = np.asarray(conf)
    n_total = len(conf)
    if n_total == 0:
        return None
    high = conf >= thr
    n_high = int(high.sum())
    frac = n_high / n_total
    est = float(np.median(bpms[high])) if high.any() else float(np.median(bpms))
    return {
        'n_windows': n_total, 'n_high': n_high,
        'high_conf_fraction': frac,
        'estimated_bpm': est,
        'mean_confidence': float(np.mean(conf)),
        'max_confidence': float(np.max(conf)),
        'bpm_std_high': float(np.std(bpms[high])) if high.any() else None,
        'reliable': bool(frac >= RECORD_OK_FRACTION),
        'verdict': 'OK' if frac >= RECORD_OK_FRACTION else 'LOW',
        'used_fallback_median': not bool(high.any()),
    }


def validate_scoring(frozen, thr, record):
    """Reproduce 03_evaluate.py's own stored row for this record, or stop."""
    print("\n" + "=" * 68)
    print("SCORING VALIDATION -- reimplemented rule vs 03_evaluate.py's output")
    print("=" * 68)

    csv_path = need(EVAL_CSV_REL, "03_evaluate.py's results CSV")
    row = None
    with open(csv_path, newline='') as fh:
        for r in csv.DictReader(fh):
            if r['record'] == record:
                row = r
                break
    if row is None:
        print(f"  [FAIL] {record} not in {csv_path}")
        print("         Run 03_evaluate.py from the repo root first.")
        sys.exit(1)
    print(f"  stored row from : {csv_path}")

    data_dir = need(DATA_REL, 'the fpcgdb data directory')
    rec = wfdb.rdrecord(os.path.join(data_dir, record))
    raw = rec.p_signal[:, 0]
    print(f"  full record     : {len(raw)} samples at {rec.fs} Hz "
          f"({len(raw)/rec.fs:.1f} s)")

    res = frozen.detect_fhr(raw, rec.fs)
    mine = score(res['bpms'], res['confidence'], thr)

    checks = [
        ('estimated_bpm', mine['estimated_bpm'], float(row['estimated_bpm']), 0.05),
        ('mean_confidence', mine['mean_confidence'], float(row['mean_confidence']), 0.005),
        ('high_conf_fraction', mine['high_conf_fraction'],
         float(row['high_conf_fraction']), 0.005),
    ]
    ok = True
    print(f"\n  {'field':<22} {'reimplemented':>14} {'stored':>12}  match")
    for name, a, b, tol in checks:
        good = abs(a - b) <= tol
        ok &= good
        print(f"  {name:<22} {a:>14.4f} {b:>12.4f}  "
              f"{'yes' if good else 'NO'}")

    stored_rel = row['reliable'].strip().lower() in ('true', '1', 'yes')
    good = (mine['reliable'] == stored_rel)
    ok &= good
    print(f"  {'reliable':<22} {str(mine['reliable']):>14} "
          f"{str(stored_rel):>12}  {'yes' if good else 'NO'}")

    if not ok:
        print("\n  [FAIL] the reimplemented rule does not reproduce")
        print("         03_evaluate.py's own output. Every comparison built on")
        print("         it would be wrong. Stopping.")
        sys.exit(1)
    print("\n  [ok] rule validated against the frozen phase's stored result")
    return mine


# ------------------------------------------------------------------ runs
def run_reference(frozen, x, fs, label, stem, thr):
    print(f"\n--- {label}: fs={fs:g} Hz, {len(x)} samples, "
          f"{len(x)/fs:.3f} s ---")
    res = frozen.detect_fhr(x, fs)
    s = score(res['bpms'], res['confidence'], thr)
    if s is None:
        print("  no windows produced")
        return None

    print(f"  windows          : {s['n_windows']}")
    print(f"  mean confidence  : {s['mean_confidence']:.3f} "
          f"(max {s['max_confidence']:.3f})")
    print(f"  windows >= {thr}  : {s['high_conf_fraction']:.1%} "
          f"({s['n_high']}/{s['n_windows']})")
    print(f"  estimated BPM    : {s['estimated_bpm']:.2f}"
          + ("   [FALLBACK median-of-all: no window cleared threshold]"
             if s['used_fallback_median'] else "   (median of high-conf windows)"))
    if s['bpm_std_high'] is not None:
        print(f"  BPM std (high)   : {s['bpm_std_high']:.2f}")
    print(f"  VERDICT          : {s['verdict']}")

    out = os.path.join(PLOT_DIR, f'03_{stem}.png')
    frozen.plot_diagnostic(x, res, fs, label, out, conf_threshold=thr)
    print(f"  plot -> {out}")
    s['label'] = label
    s['fs'] = fs
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--npy', default=None)
    args = ap.parse_args()
    for d in (RES_DIR, PLOT_DIR):
        os.makedirs(d, exist_ok=True)

    path = args.npy
    if path is None:
        hits = sorted(glob.glob(os.path.join(STIM_DIR, '*_native333.npy')))
        if not hits:
            print(f"no *_native333.npy in {STIM_DIR}/ -- run 01_stimulus.py first")
            sys.exit(1)
        path = hits[0]

    base = os.path.basename(path).replace('_native333.npy', '')
    record = base.rsplit('_', 2)[0]      # strip the _60s_60s suffix
    seg = np.load(path)
    print(f"segment file    : {path}")
    print(f"source record   : {record}")

    frozen = load_frozen()
    thr = getattr(frozen, 'CONFIDENCE_THRESHOLD', 0.45)
    print(f"threshold       : {thr} (imported)")
    print(f"record OK if    : >= {RECORD_OK_FRACTION:.0%} of windows clear it")

    full = validate_scoring(frozen, thr, record)

    print("\n" + "=" * 68)
    print("PRE-REGISTERED: B within 1.0 BPM and 0.10 confidence of A; "
          "56 windows each")
    print("=" * 68)

    a = run_reference(frozen, seg, FS_NATIVE, 'RUN A (native 333 Hz)',
                      f'runA_{base}', thr)

    # Exact-ratio resample to the device rate.
    up, down = int(FS_DEVICE), int(FS_NATIVE)     # gcd(500,333)=1
    resamp = resample_poly(seg, up, down)
    d_in, d_out = len(seg) / FS_NATIVE, len(resamp) / FS_DEVICE
    print(f"\nresample 333 -> 500: up={up} down={down}, "
          f"{len(seg)} -> {len(resamp)} samples")
    print(f"  duration {d_in:.6f} s -> {d_out:.6f} s "
          f"(delta {abs(d_out-d_in)*1000:.4f} ms)")
    if abs(d_out - d_in) > 1e-3:
        print("  [FAIL] duration changed; playback speed and BPM would shift.")
        sys.exit(1)
    print("  [ok] duration preserved")

    b = run_reference(frozen, resamp, FS_DEVICE, 'RUN B (resampled 500 Hz)',
                      f'runB_{base}', thr)

    # ---- verdicts on the predictions ----
    print("\n" + "=" * 68)
    print("PREDICTIONS")
    print("=" * 68)
    if a and b:
        d_bpm = abs(b['estimated_bpm'] - a['estimated_bpm'])
        d_cnf = abs(b['mean_confidence'] - a['mean_confidence'])
        print(f"  P1 BPM  |{b['estimated_bpm']:.2f} - {a['estimated_bpm']:.2f}| "
              f"= {d_bpm:.2f}  (tol {PRED_BPM_TOL})   "
              f"{'HELD' if d_bpm <= PRED_BPM_TOL else 'FAILED'}")
        print(f"  P2 conf |{b['mean_confidence']:.3f} - "
              f"{a['mean_confidence']:.3f}| = {d_cnf:.3f}  "
              f"(tol {PRED_CONF_TOL})   "
              f"{'HELD' if d_cnf <= PRED_CONF_TOL else 'FAILED'}")
        p3 = (a['n_windows'] == 56 and b['n_windows'] == 56)
        print(f"  P3 windows {a['n_windows']} / {b['n_windows']}   "
              f"{'HELD' if p3 else 'FAILED'}")
        print(f"  P4 run A verdict {a['verdict']}   "
              + ("-> usable as the run-C reference"
                 if a['verdict'] == 'OK' else
                 "-> LOW. Pick a different segment on stated criteria NOW,"
                 " before run C is captured."))

        print(f"\n  run C target: BPM {b['estimated_bpm']:.2f}, "
              f"confidence {b['mean_confidence']:.3f}, "
              f"{b['high_conf_fraction']:.1%} reliable")

    out = {'segment_file': path, 'record': record,
           'full_record_validation': full, 'run_A': a, 'run_B': b,
           'threshold': thr, 'record_ok_fraction': RECORD_OK_FRACTION,
           'predictions': {'bpm_tol': PRED_BPM_TOL,
                           'conf_tol': PRED_CONF_TOL, 'windows': 56}}
    rp = os.path.join(RES_DIR, f'03_reference_{base}.json')
    with open(rp, 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f"\nresults -> {rp}")


if __name__ == '__main__':
    main()
