"""
04_loopback.py -- RUN C: the acoustic loopback.

Captures p21 played at the microphone, scores it with the frozen pipeline, and
compares against run B. This is the claim the phase exists to test.

WHAT A PASS WOULD ESTABLISH
That the acquisition chain -- capsule, MAX4466, ESP32 ADC, serial link,
capture script, frozen pipeline -- recovers a KNOWN fetal heart rate end to
end against external ground truth.

WHAT IT WOULD NOT ESTABLISH
Anything about surface-conducted sound reaching the capsule. Airborne playback
into a bare electret is what an electret is designed for. The acoustic chamber
addresses a different problem (~30 dB tissue-to-air impedance loss) and remains
untested.

PRE-REGISTERED BEFORE THE CAPTURE EXISTS

  Reference, from run B: BPM 133.93, confidence 0.840, 100% reliable.

  TRIM: the FIRST 60 s of the capture. Fixed now because the capture is not
  phase-aligned to the WAV and a 90 s capture spans more than one loop; picking
  the best-scoring window afterwards would be choosing the result.

  C1  Verdict OK (>= 50% of windows clear 0.45).
  C2  BPM within 3.0 of 133.93. Wider than the 1.0 used for A->B because the
      capture starts at an arbitrary point in the record and may straddle the
      loop seam, so it covers a different stretch. Run B's own window-to-window
      std was 2.46.
  C3  Confidence BELOW 0.840. Only 7.3% of p21's in-band energy survives the
      measured playback response (-11.4 dB), and SNR is 24.4 dB against a file
      with effectively infinite SNR. Colouration and noise should cost
      something; if confidence came back equal or higher, that would need
      explaining rather than celebrating.

  Interpretation, fixed in advance:
    C1 and C2 hold  -> the chain reproduces a known FHR. Headline result.
    C1 holds, C2 fails -> detection works, rate is wrong. Timing or rate fault;
                          the loop seam is the first suspect.
    C1 fails -> the chain does not carry it. The measured playback response
                attributes the loss instead of leaving speaker vs microphone
                open.

SECONDARY, REPORTED NOT DECIDING
The same analysis on the LAST 60 s of the capture. If the two windows agree,
the loop seam did not matter. The primary result remains the first 60 s
regardless of which scores better -- that is the whole point of fixing the
trim rule in advance.

THE SCORING RULE IS NOT REIMPLEMENTED HERE
03_reference.py's score() was validated to four decimals against the row
03_evaluate.py wrote into results/fpcgdb_eval.csv. It is imported by path and
called. Two copies of a rule is two chances to drift.

USAGE
    python 04_loopback.py
    python 04_loopback.py --seconds 90 --delay 5 --note "volume 40%, factory"
"""
import os
import re
import sys
import glob
import json
import time
import argparse
import importlib.util
from datetime import datetime

import numpy as np

import matplotlib
matplotlib.use('Agg')          # MUST precede the frozen module's pyplot import
import matplotlib.pyplot as plt

import serial

# ------------------------------------------------------------------ config
PORT = 'COM8'
BAUD = 115200
SECONDS = 90.0
DELAY_S = 5.0

EXPECT_PIN = 39
EXPECT_FS = 500
NOMINAL_FS = 500.0

TRIM_S = 60.0                  # pre-registered

CAP_DIR = 'captures'
RES_DIR = 'results'
PLOT_DIR = 'plots'

REF_GLOB = os.path.join(RES_DIR, '03_reference_*.json')
REFERENCE_MODULE = '03_reference.py'

PRED_BPM_TOL = 3.0             # C2
BASELINE_NOISE_STD = 4.0

HEADER_TIMEOUT_S = 25.0


# --------------------------------------------------------------- imports
def import_by_path(filename, modname):
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, filename)
    if not os.path.isfile(path):
        print(f"[FAIL] {filename} not found beside this script ({here})")
        sys.exit(1)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, path


# ---------------------------------------------------------------- serial
def open_and_sync(port, baud):
    print(f"opening {port}...")
    ser = serial.Serial(port, baud, timeout=1)
    time.sleep(0.2)
    ser.reset_input_buffer()
    t0, header, prompted = time.time(), None, False
    while time.time() - t0 < HEADER_TIMEOUT_S:
        raw = ser.readline()
        if not raw:
            if not prompted and time.time() - t0 > 6.0:
                print("  no header yet -- press EN/RST on the ESP32")
                prompted = True
            continue
        line = raw.decode('utf-8', errors='replace').strip()
        if line.startswith('# pcg_capture'):
            header = line
            break
    if header is None:
        ser.close()
        print("[FAIL] no header. Is the Arduino Serial Monitor holding the port?")
        sys.exit(1)

    m_pin = re.search(r'pin=(\d+)', header)
    m_fs = re.search(r'fs=(\d+)', header)
    if not m_pin or not m_fs or int(m_pin.group(1)) != EXPECT_PIN \
            or int(m_fs.group(1)) != EXPECT_FS:
        ser.close()
        print(f"[FAIL] header mismatch: {header}")
        sys.exit(1)
    print(f"  {header}")
    ser.readline()                     # consume 't_us,adc'
    return ser, header


def capture(ser, seconds):
    t_us, adc, bad = [], [], 0
    t0, nxt = time.time(), 15.0
    while time.time() - t0 < seconds:
        raw = ser.readline()
        if not raw:
            continue
        line = raw.decode('utf-8', errors='replace').strip()
        if not line or line.startswith('#') or line.startswith('t_us'):
            continue
        try:
            a, b = line.split(',')
            t_us.append(int(a))
            adc.append(int(b))
        except ValueError:
            bad += 1
            continue
        el = time.time() - t0
        if el >= nxt:
            print(f"    {el:5.0f} s  {len(t_us):7d} samples")
            nxt += 15.0
    return (np.array(t_us, dtype=np.int64),
            np.array(adc, dtype=np.int64), bad)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', default=PORT)
    ap.add_argument('--seconds', type=float, default=SECONDS)
    ap.add_argument('--delay', type=float, default=DELAY_S)
    ap.add_argument('--note', default='')
    args = ap.parse_args()
    for d in (CAP_DIR, RES_DIR, PLOT_DIR):
        os.makedirs(d, exist_ok=True)

    # ---- reference ----
    hits = sorted(glob.glob(REF_GLOB))
    if not hits:
        print(f"[FAIL] no reference at {REF_GLOB}. Run 03_reference.py first.")
        sys.exit(1)
    with open(hits[-1]) as fh:
        ref = json.load(fh)
    B = ref['run_B']
    print(f"reference : {hits[-1]}")
    print(f"  run B   : BPM {B['estimated_bpm']:.2f}, "
          f"conf {B['mean_confidence']:.3f}, "
          f"{B['high_conf_fraction']:.1%} reliable")

    refmod, refpath = import_by_path(REFERENCE_MODULE, 'ref_mod')
    print(f"scoring   : {refpath} (validated rule, imported not copied)")
    frozen = refmod.load_frozen()
    thr = getattr(frozen, 'CONFIDENCE_THRESHOLD', 0.45)

    print("\n" + "=" * 68)
    print(f"PRE-REGISTERED  trim: first {TRIM_S:.0f} s")
    print(f"  C1 verdict OK | C2 BPM within {PRED_BPM_TOL:g} of "
          f"{B['estimated_bpm']:.2f} | C3 confidence below "
          f"{B['mean_confidence']:.3f}")
    print("=" * 68)

    ser, header = open_and_sync(args.port, BAUD)
    if args.delay > 0:
        print(f"\n  >>> START THE WAV NOW -- LOOPING ON <<<")
        for i in range(int(args.delay), 0, -1):
            print(f"      capture starts in {i}...")
            time.sleep(1.0)
        ser.reset_input_buffer()       # discard the countdown's backlog
        ser.readline()                 # drop the partial line at the seam
    print(f"\n  capturing {args.seconds:.0f} s")
    try:
        t_us, adc, bad = capture(ser, args.seconds)
    finally:
        ser.close()

    if len(adc) < 1000:
        print(f"[FAIL] only {len(adc)} samples")
        sys.exit(1)

    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    stem = f'pcg_runC_{stamp}'
    cap = os.path.join(CAP_DIR, stem + '.csv')
    with open(cap, 'w') as fh:
        fh.write(header + '\n')
        fh.write(f'# runC {stamp} note={args.note}\n')
        fh.write('t_us,adc\n')
        for a, b in zip(t_us, adc):
            fh.write(f'{a},{b}\n')
    print(f"\nraw -> {cap}  ({len(adc)} samples)")

    # ---- gates ----
    d = np.diff(t_us)
    print("\n=== integrity ===")
    print(f"  malformed lines : {bad}")
    nb = int(np.sum(d <= 0))
    if nb:
        print(f"  [FAIL] {nb} non-increasing timestamps. Raw CSV kept; "
              f"analysis refused.")
        sys.exit(1)
    fs = (len(t_us) - 1) / ((t_us[-1] - t_us[0]) / 1e6)
    print(f"  timestamps      : strictly increasing")
    print(f"  achieved fs     : {fs:.4f} Hz")
    print(f"  gaps >2x        : {int(np.sum(d > 2 * 1e6 / NOMINAL_FS))}")

    # ---- level, confirming nothing drifted since level-finding ----
    x = adc.astype(np.float64)
    xc = x - np.mean(x)
    rms, peak = float(np.sqrt(np.mean(xc**2))), float(np.max(np.abs(xc)))
    print("\n=== level ===")
    print(f"  bias   : {np.mean(x):.1f}   rms : {rms:.1f}")
    print(f"  peak   : {peak:.0f}         crest : {peak/(rms+1e-12):.2f}")
    print(f"  SNR    : {20*np.log10(rms/BASELINE_NOISE_STD):.1f} dB "
          f"(level-finding measured 24.4 dB)")

    # ---- seam / isolated-transient check ----
    near = int(np.sum(np.abs(xc) > 0.9 * peak))
    print(f"\n=== envelope-normalisation risk ===")
    print(f"  samples within 90% of peak : {near}")
    if near <= 5:
        print("  [!] the peak is ISOLATED. shannon_energy_envelope normalises by")
        print("      the global max, and -u*ln(u) turns over at 0.607, so one")
        print("      stray transient (the loop seam?) both scales everything")
        print("      else down AND contributes nothing itself.")
    else:
        print("  -> peak is one of many similar transients; normalisation safe")

    # ---- score ----
    n_trim = int(TRIM_S * NOMINAL_FS)

    def evaluate(sig, label):
        print(f"\n--- {label}: {len(sig)} samples, {len(sig)/fs:.2f} s ---")
        res = frozen.detect_fhr(sig, fs)
        s = refmod.score(res['bpms'], res['confidence'], thr)
        if s is None:
            print("  no windows")
            return None, None
        print(f"  windows          : {s['n_windows']}")
        print(f"  mean confidence  : {s['mean_confidence']:.3f} "
              f"(max {s['max_confidence']:.3f})")
        print(f"  windows >= {thr}  : {s['high_conf_fraction']:.1%} "
              f"({s['n_high']}/{s['n_windows']})")
        print(f"  estimated BPM    : {s['estimated_bpm']:.2f}"
              + ("   [FALLBACK: no window cleared threshold]"
                 if s['used_fallback_median'] else ""))
        if s['bpm_std_high'] is not None:
            print(f"  BPM std (high)   : {s['bpm_std_high']:.2f}")
        print(f"  VERDICT          : {s['verdict']}")
        return s, res

    if len(x) < n_trim:
        print(f"\n[FAIL] capture shorter than the {TRIM_S:.0f} s trim window")
        sys.exit(1)

    print("\n" + "=" * 68)
    print("PRIMARY -- first 60 s (pre-registered)")
    print("=" * 68)
    primary, res_p = evaluate(x[:n_trim], 'RUN C')

    print("\n" + "=" * 68)
    print("SECONDARY -- last 60 s (reported, does not decide)")
    print("=" * 68)
    secondary, _ = evaluate(x[-n_trim:], 'RUN C tail')

    # ---- predictions ----
    print("\n" + "=" * 68)
    print("PREDICTIONS")
    print("=" * 68)
    if primary:
        c1 = primary['verdict'] == 'OK'
        dbpm = abs(primary['estimated_bpm'] - B['estimated_bpm'])
        c2 = dbpm <= PRED_BPM_TOL
        c3 = primary['mean_confidence'] < B['mean_confidence']
        print(f"  C1 verdict {primary['verdict']:<4}   "
              f"{'HELD' if c1 else 'FAILED'}")
        print(f"  C2 BPM |{primary['estimated_bpm']:.2f} - "
              f"{B['estimated_bpm']:.2f}| = {dbpm:.2f} "
              f"(tol {PRED_BPM_TOL:g})   {'HELD' if c2 else 'FAILED'}")
        print(f"  C3 conf {primary['mean_confidence']:.3f} vs "
              f"{B['mean_confidence']:.3f}   {'HELD' if c3 else 'FAILED'}")

        print()
        if c1 and c2:
            print("  ACQUISITION CHAIN VALIDATED. The device recovers a known")
            print("  fetal heart rate end to end against external ground truth,")
            print("  despite 92.7% of in-band energy lost in playback.")
            print("  This says NOTHING about surface-conducted sound; the")
            print("  acoustic chamber claim remains untested.")
        elif c1:
            print("  Detection works, rate is wrong. Timing or rate fault --")
            print("  check the loop seam and the achieved fs first.")
        else:
            print("  The chain does not carry it. The measured playback response")
            print("  (-11.4 dB weighted, 39.3% distortion at 25 Hz) is available")
            print("  to attribute the loss.")

        if secondary:
            db = abs(secondary['estimated_bpm'] - primary['estimated_bpm'])
            print(f"\n  seam check: tail BPM {secondary['estimated_bpm']:.2f} "
                  f"vs primary {primary['estimated_bpm']:.2f} (delta {db:.2f}),"
                  f" tail verdict {secondary['verdict']}")

    if res_p is not None:
        out = os.path.join(PLOT_DIR, f'04_{stem}.png')
        frozen.plot_diagnostic(x[:n_trim], res_p, fs, 'RUN C (device)',
                               out, conf_threshold=thr)
        print(f"\nplot -> {out}")

    rp = os.path.join(RES_DIR, f'04_loopback_{stamp}.json')
    with open(rp, 'w') as fh:
        json.dump({'capture': cap, 'note': args.note, 'header': header,
                   'fs_achieved': fs, 'malformed': bad,
                   'trim_s': TRIM_S, 'reference_run_B': B,
                   'run_C_primary': primary, 'run_C_secondary': secondary,
                   'level': {'rms': rms, 'peak': peak,
                             'crest': peak / (rms + 1e-12),
                             'n_within_90pct_peak': near}},
                  fh, indent=2)
    print(f"results -> {rp}")


if __name__ == '__main__':
    main()
