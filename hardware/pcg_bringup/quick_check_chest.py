"""
quick_check_chest.py -- is there heart sound in the battery chest captures?

EXPLORATORY, NOT A FROZEN-PIPELINE RESULT. An adult heart runs 60-100 BPM,
outside the frozen detector's 110-180 BPM lag window, so the frozen PIPELINE
cannot be used and its verdicts do not apply. This script imports the frozen
module's own primitives by path -- bandpass_filter, shannon_energy_envelope,
estimate_fhr_autocorr -- and calls the last one with min_bpm=40, max_bpm=120,
which are EXISTING PARAMETERS of that function, not modifications to it.
Every other parameter stays at its frozen default. Any claim from this script
is a qualitative acquisition demonstration.

PRE-REGISTERED BEFORE FIRST RUN
  - Trim: first 5 s excluded (visible placement transient; fixed in advance).
  - Consistency requires BOTH:
      * median BPM of high-confidence windows within +/-10 of the pulse
        counted at capture time (parsed from the _pNN tag), and
      * >= 50% of windows clearing the imported CONFIDENCE_THRESHOLD.
    A confident rate at some OTHER frequency FAILS -- the pure-tone false
    positive (notes 7.2) established that confident-and-wrong is possible.
  - Known environmental lines, identified in EARLIER non-chest captures and
    therefore not physiological: 47.8, 69.2, 131.6, ~168, ~231 Hz. Periodicity
    attributable to these does not count as heart sound.

TWO METHODS, DELIBERATELY DIFFERENT
  1. Autocorrelation of the Shannon envelope (frozen primitives, adult band).
  2. Envelope peak-picking with inter-peak intervals: no periodicity
     assumption at all. Peaks = local maxima with prominence > 3*MAD and
     spacing >= 0.4 s; intervals outside 0.5-1.5 s (40-120 BPM) discarded and
     the discard fraction reported as a regularity measure.
  Agreement of both with each other AND with the counted pulse is the claim.

USAGE
    python quick_check_chest.py
    python quick_check_chest.py --glob "captures/pcg_char_chest_*batt*.csv"
"""
import os
import re
import sys
import glob as globmod
import argparse
import importlib.util

import numpy as np
from scipy.signal import find_peaks

import matplotlib
matplotlib.use('Agg')          # MUST precede the frozen module's pyplot import
import matplotlib.pyplot as plt

CAP_GLOB = os.path.join('captures', 'pcg_char_chest_*batt*.csv')
PLOT_DIR = 'plots'

TRIM_S = 5.0                   # pre-registered
ADULT_MIN_BPM = 40
ADULT_MAX_BPM = 120
PULSE_TOL_BPM = 10.0           # pre-registered
RELIABLE_FRAC = 0.50

PEAK_MIN_GAP_S = 0.4           # 150 BPM ceiling for peak picking
INTERVAL_LO_S = 0.5            # 120 BPM
INTERVAL_HI_S = 1.5            # 40 BPM

FROZEN_NAME = '02_fhr_detector.py'
SEARCH_DEPTH = 5

ENV_LINES_HZ = [47.8, 69.2, 131.6, 168.0, 231.0]   # known, non-physiological


def load_frozen():
    here = os.path.dirname(os.path.abspath(__file__))
    for _ in range(SEARCH_DEPTH + 1):
        cand = os.path.join(here, FROZEN_NAME)
        if os.path.isfile(cand):
            spec = importlib.util.spec_from_file_location('frozen_fhr', cand)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            print(f"frozen primitives : {cand}")
            return mod
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    print(f"could not find {FROZEN_NAME}")
    sys.exit(1)


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
    return x, fs


def counted_pulse(path):
    m = re.search(r'_p(\d+)_', os.path.basename(path))
    return int(m.group(1)) if m else None


def analyse(frozen, path, thr):
    name = os.path.basename(path)
    print("=" * 70)
    print(name)
    pulse = counted_pulse(path)
    if pulse is None:
        print("  [!] no _pNN_ in the filename; consistency cannot be judged")
    else:
        print(f"  counted pulse   : {pulse} BPM  "
              f"(consistency band {pulse - PULSE_TOL_BPM:.0f}-"
              f"{pulse + PULSE_TOL_BPM:.0f})")

    x, fs = load_capture(path)
    n_trim = int(TRIM_S * fs)
    x = x[n_trim:]
    print(f"  {len(x)} samples after {TRIM_S:g} s trim, fs={fs:.4f} Hz, "
          f"std={np.std(x):.1f} counts")

    # ---- frozen primitives, adult parameters --------------------------
    filt = frozen.bandpass_filter(x, fs)
    env = frozen.shannon_energy_envelope(filt, fs)
    times, bpms, conf = frozen.estimate_fhr_autocorr(
        env, fs, min_bpm=ADULT_MIN_BPM, max_bpm=ADULT_MAX_BPM)

    print(f"\n  --- method 1: envelope autocorrelation "
          f"({ADULT_MIN_BPM}-{ADULT_MAX_BPM} BPM window) [EXPLORATORY] ---")
    ac_bpm = None
    if len(conf) == 0:
        print("    no windows")
    else:
        high = conf >= thr
        frac = float(np.mean(high))
        ac_bpm = float(np.median(bpms[high])) if high.any() \
            else float(np.median(bpms))
        print(f"    windows          : {len(conf)}")
        print(f"    mean confidence  : {np.mean(conf):.3f} "
              f"(max {np.max(conf):.3f})")
        print(f"    windows >= {thr}  : {frac:.1%}")
        print(f"    median BPM       : {ac_bpm:.1f}"
              + ("  (high-conf windows)" if high.any()
                 else "  [FALLBACK median-of-all]"))
        if high.any():
            print(f"    BPM std (high)   : {np.std(bpms[high]):.1f}")

    # ---- method 2: peak intervals, no periodicity assumption ----------
    print(f"\n  --- method 2: envelope peak intervals ---")
    mad = float(np.median(np.abs(env - np.median(env))))
    peaks, _ = find_peaks(env, distance=int(PEAK_MIN_GAP_S * fs),
                          prominence=3.0 * mad)
    pk_bpm = None
    if len(peaks) < 3:
        print(f"    only {len(peaks)} peaks found; no interval statistics")
    else:
        iv = np.diff(peaks) / fs
        ok = (iv >= INTERVAL_LO_S) & (iv <= INTERVAL_HI_S)
        print(f"    peaks            : {len(peaks)} over "
              f"{len(env)/fs:.0f} s  ({len(peaks)/(len(env)/fs):.2f}/s)")
        print(f"    intervals in 0.5-1.5 s: {np.sum(ok)}/{len(iv)} "
              f"({100*np.mean(ok):.0f}%)")
        if ok.any():
            pk_bpm = float(60.0 / np.median(iv[ok]))
            cv = float(np.std(iv[ok]) / np.mean(iv[ok]))
            print(f"    median interval  : {np.median(iv[ok]):.3f} s "
                  f"-> {pk_bpm:.1f} BPM")
            print(f"    interval CV      : {cv:.2f}  "
                  f"(healthy resting sinus rhythm is roughly 0.02-0.1)")

    # ---- verdict against the pre-registered criteria ------------------
    print(f"\n  --- verdict ---")
    if pulse is None:
        print("    no counted pulse; reporting numbers only")
    else:
        for label, val in (("autocorr", ac_bpm), ("peak-interval", pk_bpm)):
            if val is None:
                print(f"    {label:<14}: no estimate")
                continue
            d = abs(val - pulse)
            print(f"    {label:<14}: {val:6.1f} BPM  "
                  f"(|delta| {d:.1f})  "
                  f"{'CONSISTENT' if d <= PULSE_TOL_BPM else 'NOT consistent'}")
        if ac_bpm is not None and pk_bpm is not None:
            print(f"    methods agree to {abs(ac_bpm - pk_bpm):.1f} BPM")

    # ---- plot ---------------------------------------------------------
    stem = name.replace('.csv', '')
    fig, ax = plt.subplots(3, 1, figsize=(12, 9))
    tt = np.arange(len(env)) / fs
    w = tt <= 20.0
    ax[0].plot(tt[w], env[w], linewidth=0.7)
    pw = peaks[peaks < np.sum(w)]
    ax[0].plot(pw / fs, env[pw], 'rx', markersize=6)
    ax[0].set_title(f'{stem} -- envelope, first 20 s after trim '
                    f'(x = detected peaks)')
    ax[0].set_xlabel('Time (s)'); ax[0].grid(alpha=0.3)

    if len(conf):
        hi = conf >= thr
        if hi.any():
            ax[1].plot(times[hi], bpms[hi], 'o-', markersize=3,
                       label=f'conf >= {thr}')
        if (~hi).any():
            ax[1].plot(times[~hi], bpms[~hi], 'x', alpha=0.5,
                       label=f'conf < {thr}')
        if pulse:
            ax[1].axhspan(pulse - PULSE_TOL_BPM, pulse + PULSE_TOL_BPM,
                          alpha=0.15, color='green',
                          label=f'counted pulse {pulse} +/- {PULSE_TOL_BPM:g}')
        ax[1].set_ylim(ADULT_MIN_BPM - 5, ADULT_MAX_BPM + 5)
        ax[1].set_title('autocorr BPM track (EXPLORATORY: adult parameters '
                        'on frozen primitives)')
        ax[1].set_xlabel('Time (s)'); ax[1].set_ylabel('BPM')
        ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    if len(peaks) >= 3:
        iv = np.diff(peaks) / fs
        ax[2].hist(iv, bins=40, range=(0, 2))
        ax[2].axvspan(INTERVAL_LO_S, INTERVAL_HI_S, alpha=0.12, color='green')
        if pulse:
            ax[2].axvline(60.0 / pulse, color='r', linestyle='--',
                          label=f'counted pulse ({60.0/pulse:.2f} s)')
            ax[2].legend(fontsize=8)
        ax[2].set_title('inter-peak intervals (green = 40-120 BPM)')
        ax[2].set_xlabel('interval (s)'); ax[2].grid(alpha=0.3)

    plt.tight_layout()
    out = os.path.join(PLOT_DIR, f'quick_check_chest_{stem}.png')
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"\n  plot -> {out}\n")


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
    print(f"threshold         : {thr} (imported)")
    print(f"pre-registered    : trim {TRIM_S:g} s | consistent if median "
          f"high-conf BPM within +/-{PULSE_TOL_BPM:g} of counted pulse AND "
          f">= {RELIABLE_FRAC:.0%} windows clear threshold")
    print(f"environmental Hz  : {ENV_LINES_HZ} (from earlier non-chest "
          f"captures; not physiological)\n")

    for p in paths:
        analyse(frozen, p, thr)


if __name__ == '__main__':
    main()
