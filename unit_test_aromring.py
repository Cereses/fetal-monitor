"""
unit_test_aromring.py — Sanity check the MHR detector on real MAX30102 data.

This is a UNIT TEST, not a validation dataset. The aromring repo ships one
4-second "expected good quality" capture (100 samples @ 25 Hz, RED + IR
columns) used to verify a PPG algorithm offline before flashing to hardware.

Purpose here: confirm the detector produces a physiologically sensible BPM on
data from the *actual sensor we will deploy* (MAX30102), and that it works at
25 Hz without modification. With only ~4-5 beats, the estimate is necessarily
coarse — this checks "does it produce a sane number", not accuracy.

CSV: download from
  https://github.com/aromring/MAX30102_by_RF/blob/master/ExpectedGoodQualitySignals.csv
and place at  data/aromring/ExpectedGoodQualitySignals.csv
"""
import os
import importlib.util

import numpy as np

# Load the detector (filename starts with a digit -> import via importlib).
_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    'mhr_detector', os.path.join(_HERE, '02_mhr_detector.py'))
_mhr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mhr)
detect_mhr = _mhr.detect_mhr
plot_diagnostic = _mhr.plot_diagnostic

CSV_PATH = os.path.join('data', 'aromring', 'ExpectedGoodQualitySignals.csv')
FS = 25                      # MAX30102 sample rate for this capture (#define FS 25)
PLOT_OUT = os.path.join('plots', 'aromring', 'unit_test_aromring.png')

# Expected result (computed offline against this exact file).
EXPECTED_BPM = 68.2
BPM_TOLERANCE = 3.0          # allow small drift from library/version differences


def main():
    if not os.path.exists(CSV_PATH):
        print(f"CSV not found at {CSV_PATH}.")
        print("Download ExpectedGoodQualitySignals.csv from the aromring repo "
              "and place it there.")
        return

    data = np.genfromtxt(CSV_PATH, delimiter=',', names=True)
    ir = data['IR'].astype(float)    # IR channel is the standard one for HR
    print(f"Loaded {len(ir)} samples ({len(ir)/FS:.1f} s @ {FS} Hz)")
    print(f"IR raw range: {ir.min():.0f}-{ir.max():.0f} (ADC counts)")

    result = detect_mhr(ir, FS)
    bpms = result['bpms']
    valid = bpms[np.isfinite(bpms)]

    print(f"\nWindows analyzed : {len(bpms)}")
    print(f"Peaks found      : {len(result['peaks'])} at {result['peaks'].tolist()}")
    print(f"Per-window BPM   : {np.round(bpms, 1).tolist()}")
    print(f"Per-window conf  : {np.round(result['confidence'], 2).tolist()}")

    if valid.size == 0:
        print("\nFAIL: no valid BPM estimate produced.")
        return

    est = float(np.median(valid))
    print(f"\nEstimated MHR    : {est:.1f} BPM")

    # ---- assertions ----
    ok_range = 40 <= est <= 200
    ok_match = abs(est - EXPECTED_BPM) <= BPM_TOLERANCE
    print(f"Physiological (40-200 BPM) : {'PASS' if ok_range else 'FAIL'}")
    print(f"Matches expected {EXPECTED_BPM} BPM : "
          f"{'PASS' if ok_match else 'FAIL'} (tol +/-{BPM_TOLERANCE})")

    os.makedirs(os.path.dirname(PLOT_OUT), exist_ok=True)
    plot_diagnostic(ir, result, FS, 'aromring_MAX30102', PLOT_OUT, zoom_s=4)
    print(f"\nDiagnostic plot saved to: {PLOT_OUT}")

    print("\nUNIT TEST PASSED" if (ok_range and ok_match) else "\nUNIT TEST FAILED")


if __name__ == '__main__':
    main()
