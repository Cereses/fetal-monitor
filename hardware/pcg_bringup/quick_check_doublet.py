"""
quick_check_doublet.py -- are the envelope peaks S1/S2 doublets?

EXPLORATORY. Adult heart sounds, outside the frozen detector's 110-180 BPM
window. Frozen primitives are imported by path and called with existing
parameters; nothing frozen is modified. No frozen-pipeline verdict applies.

THE QUESTION
quick_check_chest.py found evenly spaced envelope peaks at ~0.68 s in the
bare-capsule contact capture (CV 0.05, 100% of windows reliable), with a
visible SMALLER SECONDARY BUMP after each large peak. That is the shape of an
S1/S2 doublet -- but the peak finder used distance=0.4 s, and systole is
roughly 0.3 s, so S2 was SUPPRESSED BY CONSTRUCTION and never counted.

Lowering distance to 0.15 s admits S2. If the events are cardiac, the interval
series becomes ALTERNATING: S1->S2 short (systole), S2->S1 long (diastole).
Systole being shorter than diastole is close to diagnostic -- very little else
that is contact-dependent produces an alternating doublet at a plausible rate.

PRE-REGISTERED BEFORE THE RUN

  Detection: distance 0.15 s. Prominence = 10% of the MEDIAN PROMINENCE of the
  coarse (0.4 s) peaks -- data-derived but scale-free, so the same rule works
  on the ~50-count contact captures and the ~4-count air control without a
  hand-tuned constant.

  D1  peak count 1.6x-2.4x the coarse count
  D2  lag-1 autocorrelation of the interval series <= -0.30
      A strictly alternating short-long sequence has lag-1 correlation near
      -1; a uniform sequence has it near 0. This is the load-bearing test
      because it needs NO threshold on interval length.
  D3  intervals bimodal (Otsu split), short/(short+long) in 0.25-0.45
      Systole occupies roughly a third of the cycle at resting rates.
  D4  cycle period (short+long) agrees with the coarse-peak rate within 5 BPM

  CARDIAC requires D1-D4. Any failure = not established.

  FALSIFIER: the air control runs the identical analysis. If it ALSO shows
  alternation, this test measures nothing and the result is discarded.

USAGE
    python quick_check_doublet.py
    python quick_check_doublet.py --glob "captures/pcg_char_chest_*batt*.csv"
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

TRIM_S = 5.0                   # as in quick_check_chest.py
COARSE_GAP_S = 0.40            # the setting that suppressed S2
FINE_GAP_S = 0.15              # admits S2
PROM_FRACTION = 0.10           # of median coarse prominence

D1_LO, D1_HI = 1.6, 2.4
D2_LAG1_MAX = -0.30
D3_SHORT_FRAC_LO, D3_SHORT_FRAC_HI = 0.25, 0.45
D4_BPM_TOL = 5.0

FROZEN_NAME = '02_fhr_detector.py'
SEARCH_DEPTH = 5


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
    return x, (len(t) - 1) / ((t[-1] - t[0]) / 1e6)


def otsu_split(v):
    """Threshold maximising between-class variance. Returns (thr, lo, hi)."""
    s = np.sort(v)
    best_thr, best_var = None, -1.0
    for i in range(1, len(s)):
        a, b = s[:i], s[i:]
        if len(a) < 2 or len(b) < 2:
            continue
        w = len(a) * len(b) / len(s) ** 2
        var = w * (a.mean() - b.mean()) ** 2
        if var > best_var:
            best_var, best_thr = var, 0.5 * (s[i - 1] + s[i])
    if best_thr is None:
        return None, None, None
    return best_thr, v[v <= best_thr], v[v > best_thr]


def analyse(frozen, path):
    name = os.path.basename(path)
    print("=" * 72)
    print(name)
    m = re.search(r'_p(\d+)_', name)
    pulse = int(m.group(1)) if m else None

    x, fs = load_capture(path)
    x = x[int(TRIM_S * fs):]
    print(f"  {len(x)} samples after {TRIM_S:g} s trim, fs={fs:.4f} Hz, "
          f"std={np.std(x):.1f} counts")

    env = frozen.shannon_energy_envelope(
        frozen.bandpass_filter(x, fs), fs)

    # ---- coarse pass: S2 suppressed by construction --------------------
    mad = float(np.median(np.abs(env - np.median(env))))
    c_pk, c_props = find_peaks(env, distance=int(COARSE_GAP_S * fs),
                               prominence=3.0 * mad)
    if len(c_pk) < 5:
        print(f"  only {len(c_pk)} coarse peaks; cannot proceed")
        return
    c_prom = float(np.median(c_props['prominences']))
    c_iv = np.diff(c_pk) / fs
    ok = (c_iv >= 0.5) & (c_iv <= 1.5)
    c_bpm = float(60.0 / np.median(c_iv[ok])) if ok.any() else float('nan')
    print(f"\n  coarse ({COARSE_GAP_S:g} s gap): {len(c_pk)} peaks, "
          f"median prominence {c_prom:.4f}, {c_bpm:.1f} BPM")

    # ---- fine pass: S2 admitted ----------------------------------------
    prom = PROM_FRACTION * c_prom
    f_pk, _ = find_peaks(env, distance=int(FINE_GAP_S * fs), prominence=prom)
    ratio = len(f_pk) / len(c_pk)
    print(f"  fine   ({FINE_GAP_S:g} s gap): {len(f_pk)} peaks "
          f"(prominence {prom:.4f})")

    if len(f_pk) < 6:
        print("  too few fine peaks; cannot test alternation")
        return
    iv = np.diff(f_pk) / fs

    # ---- D1 -------------------------------------------------------------
    d1 = D1_LO <= ratio <= D1_HI
    print(f"\n  D1 count ratio   : {ratio:.2f}x  "
          f"(need {D1_LO}-{D1_HI})   {'PASS' if d1 else 'FAIL'}")

    # ---- D2: the load-bearing test --------------------------------------
    lag1 = float(np.corrcoef(iv[:-1], iv[1:])[0, 1]) if len(iv) > 2 else 0.0
    d2 = lag1 <= D2_LAG1_MAX
    print(f"  D2 lag-1 autocorr: {lag1:+.3f}  "
          f"(need <= {D2_LAG1_MAX})   {'PASS' if d2 else 'FAIL'}")
    print(f"     (near -1 = strict alternation; near 0 = uniform spacing)")

    # ---- D3 -------------------------------------------------------------
    thr, lo, hi = otsu_split(iv)
    d3 = False
    short = long_ = cycle = None
    if thr is not None and len(lo) >= 2 and len(hi) >= 2:
        short, long_ = float(np.median(lo)), float(np.median(hi))
        cycle = short + long_
        frac = short / cycle
        d3 = D3_SHORT_FRAC_LO <= frac <= D3_SHORT_FRAC_HI
        print(f"  D3 bimodal split : short {short:.3f} s (n={len(lo)}), "
              f"long {long_:.3f} s (n={len(hi)})")
        print(f"     short fraction: {frac:.3f}  "
              f"(need {D3_SHORT_FRAC_LO}-{D3_SHORT_FRAC_HI})   "
              f"{'PASS' if d3 else 'FAIL'}")
    else:
        print("  D3 bimodal split : could not split   FAIL")

    # ---- D4 -------------------------------------------------------------
    d4 = False
    if cycle:
        bpm_cycle = 60.0 / cycle
        d = abs(bpm_cycle - c_bpm)
        d4 = d <= D4_BPM_TOL
        print(f"  D4 cycle rate    : {bpm_cycle:.1f} BPM vs coarse "
              f"{c_bpm:.1f} (delta {d:.1f}, need <= {D4_BPM_TOL:g})   "
              f"{'PASS' if d4 else 'FAIL'}")

    # ---- verdict --------------------------------------------------------
    cardiac = d1 and d2 and d3 and d4
    print(f"\n  VERDICT: {'DOUBLET PATTERN CONSISTENT WITH S1/S2' if cardiac else 'NOT ESTABLISHED'}")
    if cardiac and pulse:
        print(f"  counted pulse {pulse} BPM; cycle rate "
              f"{60.0/cycle:.1f} BPM (delta {abs(60.0/cycle - pulse):.1f})")

    # ---- plot -----------------------------------------------------------
    fig, ax = plt.subplots(3, 1, figsize=(12, 9))
    tt = np.arange(len(env)) / fs
    w = tt <= 8.0
    n_w = int(np.sum(w))
    ax[0].plot(tt[w], env[w], linewidth=0.8)
    cw = c_pk[c_pk < n_w]
    fw = f_pk[f_pk < n_w]
    sec = np.setdiff1d(fw, cw)
    ax[0].plot(cw / fs, env[cw], 'rx', markersize=9, label='coarse (S1?)')
    ax[0].plot(sec / fs, env[sec], 'g+', markersize=9,
               label='additional at fine gap (S2?)')
    ax[0].set_title(f'{name} -- envelope, first 8 s')
    ax[0].set_xlabel('Time (s)'); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    ax[1].plot(iv, 'o-', markersize=3, linewidth=0.8)
    if thr is not None:
        ax[1].axhline(thr, color='r', linestyle='--',
                      label=f'Otsu split {thr:.3f} s')
        ax[1].legend(fontsize=8)
    ax[1].set_title(f'interval sequence (lag-1 autocorr {lag1:+.3f}; '
                    f'alternation = sawtooth)')
    ax[1].set_xlabel('interval index'); ax[1].set_ylabel('s')
    ax[1].grid(alpha=0.3)

    ax[2].hist(iv, bins=50, range=(0, 1.5))
    if thr is not None:
        ax[2].axvline(thr, color='r', linestyle='--')
    ax[2].set_title('interval distribution (two clusters = doublet)')
    ax[2].set_xlabel('interval (s)'); ax[2].grid(alpha=0.3)

    plt.tight_layout()
    out = os.path.join(PLOT_DIR,
                       f"quick_check_doublet_{name.replace('.csv','')}.png")
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"  plot -> {out}\n")


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
    print(f"pre-registered: D1 count {D1_LO}-{D1_HI}x | D2 lag-1 <= "
          f"{D2_LAG1_MAX} | D3 short fraction {D3_SHORT_FRAC_LO}-"
          f"{D3_SHORT_FRAC_HI} | D4 within {D4_BPM_TOL:g} BPM")
    print("the AIR capture is the falsifier: if it also alternates, "
          "this test measures nothing\n")

    for p in paths:
        analyse(frozen, p)


if __name__ == '__main__':
    main()
