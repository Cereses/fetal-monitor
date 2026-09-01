"""
quick_check_stimulus.py -- measure the reference segment before anything is
built on top of it.

Answers two questions that were guessed at rather than measured:

  1. QUANTISATION. 01_stimulus.py reported a range of -21400.0 to 19400.0.
     Both endpoints being exact multiples of 100 suggests wfdb's p_signal is
     scaling raw counts by the header's gain field, i.e. the TRUE raw ADC
     swing is 100x smaller than it appears. If that holds, the reference data
     is itself coarsely quantised, which changes how demanding a target the
     device's ADC actually has to hit. Inferred from two round numbers is not
     the same as measured, so this measures it.

  2. BAND ENERGY. Decides the playback hardware. Fetal S1 energy sits low, and
     small speakers have essentially no output below ~150 Hz. If most of this
     segment's in-band energy is under 100 Hz, a phone or laptop speaker
     cannot reproduce it and sealed earbuds coupled to the capsule are the
     only viable source. Measuring this now avoids attributing a speaker's
     rolloff to the microphone later.

Bands are reported against the FROZEN pipeline's own passband, which at
fs=333 Hz evaluates to 25 Hz .. min(200, 0.95*166.5) = 158.175 Hz. Energy
outside that is discarded by the detector and is reported separately so it
cannot inflate the in-band figures.

USAGE
    python quick_check_stimulus.py
    python quick_check_stimulus.py --npy stimulus/<other>_native333.npy
"""
import os
import glob
import argparse

import numpy as np
from scipy.signal import welch

STIM_DIR = 'stimulus'
FS = 333.0

BAND_LOW = 25.0
BAND_HIGH = min(200.0, 0.95 * (FS / 2))     # 158.175 Hz -- matches the pipeline

SUB_BANDS = [(0.0, 25.0), (25.0, 50.0), (50.0, 100.0),
             (100.0, BAND_HIGH), (BAND_HIGH, FS / 2)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--npy', default=None)
    args = ap.parse_args()

    path = args.npy
    if path is None:
        hits = sorted(glob.glob(os.path.join(STIM_DIR, '*_native333.npy')))
        if not hits:
            print(f"no *_native333.npy in {STIM_DIR}/ -- run 01_stimulus.py first")
            return
        path = hits[0]
        if len(hits) > 1:
            print(f"note: {len(hits)} candidates, using {os.path.basename(path)}")

    x = np.load(path)
    print(f"file     : {path}")
    print(f"samples  : {len(x)}  ({len(x)/FS:.3f} s at {FS:g} Hz)\n")

    # ---------------------------------------------------------- quantisation
    print("=== quantisation ===")
    u = np.unique(x)
    print(f"  distinct values : {len(u)} out of {len(x)} samples")

    if len(u) < 2:
        print("  [!] signal is constant")
        return

    steps = np.diff(u)
    gcd_step = steps.min()
    for s in np.unique(np.round(steps, 9)):
        gcd_step = np.gcd(int(round(gcd_step * 1e6)), int(round(s * 1e6))) / 1e6
    print(f"  smallest gap    : {steps.min():.6f}")
    print(f"  common divisor  : {gcd_step:.6f}")

    if gcd_step > 1.5:
        implied = x / gcd_step
        print(f"  -> values are quantised in steps of {gcd_step:g}")
        print(f"  -> IMPLIED RAW COUNTS: {implied.min():.0f} to {implied.max():.0f}"
              f"  (span {implied.max()-implied.min():.0f})")
        span = implied.max() - implied.min()
        print(f"  -> usable resolution in this segment: "
              f"~{np.log2(span+1):.1f} bits")
        print("  -> CONFIRMS the gain-scaling hypothesis.")
    else:
        print("  -> step size is ~1; values are already raw counts.")
        print("  -> the round min/max were coincidence. Hypothesis REJECTED.")
    print()

    # ---------------------------------------------------------- amplitude
    print("=== amplitude ===")
    xc = x - np.mean(x)
    rms = float(np.sqrt(np.mean(xc ** 2)))
    peak = float(np.max(np.abs(xc)))
    print(f"  mean            : {np.mean(x):.2f}")
    print(f"  rms (DC-removed): {rms:.2f}")
    print(f"  peak            : {peak:.2f}")
    print(f"  crest factor    : {peak/rms:.2f}  ({20*np.log10(peak/rms):.1f} dB)")
    print("  (PCG is bursty: peak-normalising for the WAV sets level from the")
    print("   loudest transient, so average playback level sits well below")
    print("   full scale. A high crest factor means volume must go up, which")
    print("   also raises room noise in the quiet intervals.)")
    print()

    # ---------------------------------------------------------- spectrum
    print("=== band energy ===")
    nper = min(4096, len(xc))
    f, p = welch(xc, fs=FS, nperseg=nper)
    total = float(np.trapezoid(p, f)) if hasattr(np, 'trapezoid') \
        else float(np.trapz(p, f))

    for lo, hi in SUB_BANDS:
        m = (f >= lo) & (f < hi)
        e = float(np.trapezoid(p[m], f[m])) if hasattr(np, 'trapezoid') \
            else float(np.trapz(p[m], f[m]))
        tag = ''
        if hi <= BAND_LOW:
            tag = '  (discarded by pipeline highpass)'
        elif lo >= BAND_HIGH:
            tag = '  (discarded by pipeline lowpass)'
        print(f"  {lo:6.1f} - {hi:6.1f} Hz : {100*e/total:5.1f} %{tag}")

    inb = (f >= BAND_LOW) & (f < BAND_HIGH)
    e_in = float(np.trapezoid(p[inb], f[inb])) if hasattr(np, 'trapezoid') \
        else float(np.trapz(p[inb], f[inb]))
    print(f"\n  in-band total ({BAND_LOW:g}-{BAND_HIGH:.1f} Hz): "
          f"{100*e_in/total:.1f} % of all energy")

    if e_in > 0:
        centroid = float(np.sum(f[inb] * p[inb]) / np.sum(p[inb]))
        csum = np.cumsum(p[inb])
        csum /= csum[-1]
        med = float(f[inb][np.searchsorted(csum, 0.5)])
        print(f"  in-band centroid : {centroid:.1f} Hz")
        print(f"  in-band median   : {med:.1f} Hz  "
              f"(half the in-band energy is below this)")

        below100 = (f >= BAND_LOW) & (f < 100.0)
        e_lo = float(np.trapezoid(p[below100], f[below100])) \
            if hasattr(np, 'trapezoid') else float(np.trapz(p[below100], f[below100]))
        frac = 100 * e_lo / e_in
        print(f"\n  fraction of IN-BAND energy below 100 Hz: {frac:.1f} %")
        if frac > 60:
            print("  -> a phone or laptop speaker will not reproduce most of this.")
            print("     Use sealed earbuds coupled directly to the capsule.")
        elif frac > 35:
            print("  -> marginal for a small speaker. Earbuds preferred;")
            print("     the tone sweep will settle it.")
        else:
            print("  -> energy sits high enough that a small speaker is viable.")


if __name__ == '__main__':
    main()
