"""
quick_check_lines.py -- what is in the noise floor, and does the detector bite on it?

WHY THIS EXISTS
02_capture.py pre-registered a mains test at 50/100/150 Hz +/- 1 Hz. It returned
"not detected" on all three harmonics in all three runs, while every run
independently reported peak_inband_hz = 47.8 Hz -- a sharp line eleven
resolution bins outside the window it was looking through.

The test was not wrong about 50 Hz. It was BLIND BY CONSTRUCTION, because it
pre-registered WHERE to look instead of WHAT to look for. Same failure shape as
the flat-fraction metric that confused quantisation with dropout.

The fix is a general line detector, NOT a widened mains window. Widening the
window after seeing the peak would be tuning the test until it found the thing
already known to be there.

TWO TESTS, both on captures already recorded. No new hardware runs.

  TEST 1 -- LINE DETECTION
    Every narrow peak in 25-200 Hz that stands above its own local background,
    wherever it sits. Run across all captures so cross-run consistency
    distinguishes a stable environmental source from a transient one.

  TEST 2 -- FROZEN DETECTOR ON PURE NOISE   [PRE-REGISTERED]
    90 s of quiet room contains no heart sound. The frozen pipeline should
    therefore report mean confidence well below CONFIDENCE_THRESHOLD and a
    record-level LOW verdict (<50% of windows reliable).

    If it reports OK, the detector is finding structure in noise. That is a
    more serious finding than the 47.8 Hz line itself, and it would apply to
    every FHR capture this project ever takes.

    Prediction recorded BEFORE the run. Whatever comes back is the result.

WHY A STEADY TONE MAY NOT MATTER, AND WHEN IT WOULD
shannon_energy_envelope squares the signal, then applies a 50 ms moving average
(25 samples at 500 Hz). Squaring moves a 47.8 Hz carrier to 95.6 Hz, where the
smoother's response is about -27 dB. A STEADY tone therefore yields a nearly
constant envelope and contributes little to the autocorrelation. An
AMPLITUDE-MODULATED tone does not get suppressed, because detecting amplitude
modulation is exactly what the envelope stage is for. So the modulation depth
of the strongest line is measured, not assumed.

The frozen detector is IMPORTED BY PATH, never copied. Its filename starts with
a digit, so importlib is required.

USAGE
    python quick_check_lines.py
    python quick_check_lines.py --glob "captures/pcg_char_baseline_*.csv"
"""
import os
import sys
import glob as globmod
import argparse
import importlib.util

import numpy as np
from scipy.signal import welch, butter, filtfilt, hilbert

import matplotlib
matplotlib.use('Agg')          # MUST precede the frozen module's pyplot import
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ config
CAP_GLOB = os.path.join('captures', 'pcg_char_*.csv')
PLOT_DIR = 'plots'

BAND_LOW, BAND_HIGH = 25.0, 200.0        # pipeline passband at fs=500

LINE_BG_HALFWIDTH = 6.0                  # local background window, Hz
LINE_PEAK_HALFWIDTH = 1.0                # excluded from background, Hz
LINE_DETECT_DB = 6.0                     # unchanged from the original test
LINE_MAX_REPORT = 8

MODULATION_HALFWIDTH = 4.0               # narrow bandpass around the line

FROZEN_NAME = '02_fhr_detector.py'
SEARCH_DEPTH = 5

# From 03_evaluate.py's record-level rule. CONFIRM this matches that file.
RECORD_OK_FRACTION = 0.50


def _trapz(y, x):
    f = getattr(np, 'trapezoid', None) or np.trapz
    return float(f(y, x))


# ------------------------------------------------------- frozen import
def load_frozen(depth=SEARCH_DEPTH):
    """Locate and import the frozen detector BY PATH. Never copy it."""
    here = os.path.dirname(os.path.abspath(__file__))
    tried = []
    for _ in range(depth + 1):
        cand = os.path.join(here, FROZEN_NAME)
        tried.append(cand)
        if os.path.isfile(cand):
            spec = importlib.util.spec_from_file_location('frozen_fhr', cand)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)   # NOTE: module-level side effect --
                                           # creates plots/fpcgdb/. Harmless.
            print(f"frozen detector : {cand}")
            return mod
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    print(f"could not find {FROZEN_NAME}. Looked in:")
    for t in tried:
        print("   ", t)
    sys.exit(1)


# ------------------------------------------------------------------ io
def load_capture(path):
    t_us, adc = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('t_us'):
                continue
            a, b = line.split(',')
            t_us.append(int(a))
            adc.append(int(b))
    t = np.array(t_us, dtype=np.int64)
    x = np.array(adc, dtype=np.float64)
    fs = (len(t) - 1) / ((t[-1] - t[0]) / 1e6)
    return t, x, fs


# ------------------------------------------------------------- test 1
def detect_lines(xc, fs):
    """Every narrow in-band peak standing above its own local background."""
    nper = min(2500, len(xc))
    f, p = welch(xc, fs=fs, nperseg=nper)

    inb = (f >= BAND_LOW) & (f <= BAND_HIGH)
    e_in = _trapz(p[inb], f[inb])

    idx = np.where(inb)[0]
    lines = []
    for i in idx:
        if i == 0 or i == len(p) - 1:
            continue
        if not (p[i] > p[i - 1] and p[i] >= p[i + 1]):
            continue                                   # local maximum only
        f0 = f[i]
        near = (f >= f0 - LINE_PEAK_HALFWIDTH) & (f <= f0 + LINE_PEAK_HALFWIDTH)
        bg = ((f >= f0 - LINE_BG_HALFWIDTH) & (f <= f0 + LINE_BG_HALFWIDTH)
              & ~near)
        if not bg.any():
            continue
        base = float(np.median(p[bg]))
        db = 10.0 * np.log10(p[i] / (base + 1e-30))
        if db >= LINE_DETECT_DB:
            lines.append({'hz': float(f0), 'db': float(db),
                          'pct_inband': 100.0 * _trapz(p[near], f[near]) / e_in})

    # collapse neighbours belonging to one physical line
    lines.sort(key=lambda d: -d['db'])
    kept = []
    for L in lines:
        if all(abs(L['hz'] - k['hz']) > 2 * LINE_PEAK_HALFWIDTH for k in kept):
            kept.append(L)
    return kept[:LINE_MAX_REPORT], f, p


def modulation_depth(xc, fs, f0):
    """How amplitude-modulated is this line?

    A steady tone is largely flattened by the envelope stage. A modulated one
    is not, and modulation is what the autocorrelation would lock onto.
    """
    lo = max(1.0, f0 - MODULATION_HALFWIDTH)
    hi = min(0.95 * fs / 2, f0 + MODULATION_HALFWIDTH)
    b, a = butter(4, [lo / (fs / 2), hi / (fs / 2)], btype='band')
    narrow = filtfilt(b, a, xc)
    env = np.abs(hilbert(narrow))
    return float(np.mean(env)), float(np.std(env) / (np.mean(env) + 1e-12))


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--glob', default=CAP_GLOB)
    args = ap.parse_args()
    os.makedirs(PLOT_DIR, exist_ok=True)

    paths = sorted(globmod.glob(args.glob))
    if not paths:
        print(f"no captures matching {args.glob}")
        sys.exit(1)

    frozen = load_frozen()
    thr = getattr(frozen, 'CONFIDENCE_THRESHOLD', 0.45)
    print(f"CONFIDENCE_THRESHOLD = {thr}   (imported, not hardcoded)")
    print(f"record OK if >= {RECORD_OK_FRACTION:.0%} of windows clear it\n")

    all_lines, summary = {}, []

    for path in paths:
        name = os.path.basename(path)
        print("=" * 70)
        print(name)
        t, x, fs = load_capture(path)
        xc = x - np.mean(x)
        print(f"  {len(x)} samples, fs={fs:.4f} Hz, "
              f"{(t[-1]-t[0])/1e6:.1f} s, std={np.std(x):.3f} counts")

        # ---- TEST 1 ----
        lines, f, p = detect_lines(xc, fs)
        print(f"\n  --- narrow lines in {BAND_LOW:g}-{BAND_HIGH:g} Hz "
              f"(>{LINE_DETECT_DB:g} dB over local background) ---")
        if not lines:
            print("    none")
        for L in lines:
            print(f"    {L['hz']:7.2f} Hz  {L['db']:+6.1f} dB  "
                  f"{L['pct_inband']:5.2f} % of in-band power")
        all_lines[name] = lines

        if lines:
            f0 = lines[0]['hz']
            amp, cv = modulation_depth(xc, fs, f0)
            print(f"\n    strongest line {f0:.2f} Hz: envelope mean {amp:.2f} "
                  f"counts, CV {cv:.3f}")
            print("    (low CV = steady tone, largely flattened by the 50 ms")
            print("     smoother; high CV = modulated, which survives it)")

        # ---- TEST 2, pre-registered ----
        print(f"\n  --- frozen detector on quiet-room noise [PRE-REGISTERED] ---")
        print("    prediction: mean confidence << threshold, verdict LOW")
        res = frozen.detect_fhr(x, fs)
        conf = np.asarray(res['confidence'])
        bpms = np.asarray(res['bpms'])
        if len(conf) == 0:
            print("    no windows produced")
            continue

        rel = float(np.mean(conf >= thr))
        verdict = 'OK' if rel >= RECORD_OK_FRACTION else 'LOW'
        print(f"    windows          : {len(conf)}")
        print(f"    mean confidence  : {np.mean(conf):.3f}  "
              f"(max {np.max(conf):.3f})")
        print(f"    windows >= {thr}  : {rel:.1%}")
        print(f"    BPM mean/std     : {np.mean(bpms):.1f} / {np.std(bpms):.1f}")
        print(f"    VERDICT          : {verdict}")
        if verdict == 'OK':
            print("\n    *** PREDICTION FAILED. The detector reports a reliable")
            print("        heart rate from a signal containing no heart. This")
            print("        is a critical finding for the FHR channel and must")
            print("        be characterised before any acoustic capture. ***")
        else:
            print("    -> as predicted. The noise floor does not fool the")
            print("       detector, so an OK verdict later carries meaning.")

        summary.append({'name': name, 'fs': fs, 'std': float(np.std(x)),
                        'n_lines': len(lines),
                        'top_hz': lines[0]['hz'] if lines else None,
                        'top_db': lines[0]['db'] if lines else None,
                        'mean_conf': float(np.mean(conf)),
                        'reliable_frac': rel, 'verdict': verdict})
        print()

    # ---- cross-run consistency ----
    print("=" * 70)
    print("CROSS-RUN SUMMARY\n")
    print(f"  {'capture':<45} {'std':>6} {'top line':>10} {'dB':>6} "
          f"{'conf':>6} {'rel':>6}  verdict")
    for s in summary:
        top = f"{s['top_hz']:.2f}" if s['top_hz'] else '-'
        db = f"{s['top_db']:+.1f}" if s['top_db'] else '-'
        print(f"  {s['name'][:45]:<45} {s['std']:6.2f} {top:>10} {db:>6} "
              f"{s['mean_conf']:6.3f} {s['reliable_frac']:6.1%}  {s['verdict']}")

    tops = [s['top_hz'] for s in summary if s['top_hz']]
    if len(tops) > 1:
        print(f"\n  top-line frequency across runs: "
              f"{np.mean(tops):.2f} +/- {np.std(tops):.2f} Hz")
        if np.std(tops) < 0.5:
            print("  -> stable across runs: a persistent environmental or")
            print("     electrical source, not a transient event.")
        else:
            print("  -> varies across runs: not a single fixed source.")

    # ---- overlay plot ----
    fig, ax = plt.subplots(figsize=(12, 5))
    for path in paths:
        _, x, fs = load_capture(path)
        xc = x - np.mean(x)
        f, p = welch(xc, fs=fs, nperseg=min(2500, len(xc)))
        ax.semilogy(f, p, linewidth=0.8, alpha=0.75,
                    label=os.path.basename(path)[-19:-4])
    ax.axvspan(BAND_LOW, BAND_HIGH, alpha=0.08, color='green')
    for h in (50, 100, 150):
        ax.axvline(h, color='r', alpha=0.35, linestyle='--')
    ax.set_xlabel('Hz'); ax.set_ylabel('PSD')
    ax.set_title('Noise-floor PSD, all runs (green = passband, red = nominal mains)')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    plt.tight_layout()
    out = os.path.join(PLOT_DIR, 'quick_check_lines_overlay.png')
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"\noverlay -> {out}")


if __name__ == '__main__':
    main()
