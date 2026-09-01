"""
quick_check_sweep.py -- what does the playback chain actually reproduce?

Measures the frequency response of speaker -> air -> capsule -> ADC from a
tone-sweep capture, so that a run C result can be attributed to the right
stage instead of guessed at.

WHAT THIS IS NOT
It is NOT a go/no-go gate on the acoustic loopback. estimate_fhr_autocorr runs
on the Shannon energy ENVELOPE: it detects the TIMING of S1/S2 bursts, not
their spectral content. A heartbeat bandpass-filtered by a mediocre speaker
still has bursts at the same instants, so as long as some of each burst's
energy survives, the envelope keeps its structure and the autocorrelation can
still lock on. A coloured response therefore does not necessarily kill run C.
This measurement tells us what to EXPECT, and lets us explain a failure.

THREE THINGS MEASURED, EACH FOR A REASON

  1. NARROWBAND POWER AT EACH FUNDAMENTAL, not broadband RMS.
     If the subwoofer cannot reproduce 25 Hz it will distort instead, putting
     energy at 50/75/100 Hz. A broadband measurement would then score 25 Hz as
     "present" while the capsule is actually hearing harmonics of it. Power is
     integrated in a narrow window at f0 only, and harmonic content is
     reported separately so distortion is visible rather than hidden.

  2. AMPLITUDE ACROSS THIRDS OF EACH BURST.
     The 75 Hz burst appears in the raw trace to grow over roughly its first
     1.5 s. The stimulus fades are 50 ms, so that is not the stimulus.
     Candidates are a room mode building, slow amplifier dynamics, or a
     bass-boost circuit. Measured, not eyeballed.

  3. ENERGY-WEIGHTED TRANSMISSION OF THE ACTUAL REFERENCE SEGMENT.
     The response at six discrete frequencies is interpolated onto the p21
     segment's own spectrum, and the fraction of its in-band energy that
     survives is computed. That converts an abstract response curve into the
     one number that matters: how much of the signal run C depends on would
     actually get through.

BURST IDENTIFICATION
Bursts are found by short-time RMS, not by assuming a start time. The layout
(6 bursts, 3 s each, 1 s gaps) comes from stimulus_manifest.json. If exactly
six are found they are assigned in frequency order; otherwise the script
reports what it found and stops rather than guessing an alignment.

USAGE
    python quick_check_sweep.py
    python quick_check_sweep.py --csv captures/pcg_char_sweep_sub_*.csv
"""
import os
import sys
import glob as globmod
import json
import argparse

import numpy as np
from scipy.signal import welch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ config
CAP_GLOB = os.path.join('captures', 'pcg_char_sweep_*.csv')
STIM_DIR = 'stimulus'
PLOT_DIR = 'plots'
RES_DIR = 'results'

TONES_HZ = [25.0, 50.0, 75.0, 100.0, 150.0, 200.0]

RMS_WIN_S = 0.10               # short-time RMS window for burst detection
BURST_MIN_S = 1.5              # a real burst is 3 s; reject shorter blips
BURST_THRESH_MULT = 4.0        # above this x the quiet floor counts as burst
EDGE_TRIM_S = 0.25             # drop each burst's edges: fades + settling

NARROW_HALFWIDTH = 2.0         # +/- Hz integrated as "the fundamental"
N_HARMONICS = 4

# p21's own passband at its native rate: min(200, 0.95*166.5)
P21_BAND = (25.0, 158.175)


def _trapz(y, x):
    f = getattr(np, 'trapezoid', None) or np.trapz
    return float(f(y, x))


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


# --------------------------------------------------------- burst finding
def find_bursts(xc, fs):
    """Locate bursts by short-time RMS. Returns list of (i0, i1)."""
    w = max(1, int(RMS_WIN_S * fs))
    n = len(xc) // w
    rms = np.array([np.sqrt(np.mean(xc[i*w:(i+1)*w] ** 2)) for i in range(n)])

    # Quiet floor from the lower quartile of frames: the capture is mostly
    # silence plus six bursts, so the low quartile is reliably background.
    floor = float(np.percentile(rms, 25))
    thr = floor * BURST_THRESH_MULT
    active = rms > thr

    segs, start = [], None
    for i, a in enumerate(active):
        if a and start is None:
            start = i
        elif not a and start is not None:
            segs.append((start, i))
            start = None
    if start is not None:
        segs.append((start, len(active)))

    out = []
    for s, e in segs:
        if (e - s) * RMS_WIN_S >= BURST_MIN_S:
            out.append((s * w, min(e * w, len(xc))))
    return out, floor, thr, rms, w


def band_power(f, p, f0, hw=NARROW_HALFWIDTH):
    m = (f >= f0 - hw) & (f <= f0 + hw)
    if not m.any():
        return 0.0
    return _trapz(p[m], f[m])


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default=None)
    args = ap.parse_args()
    for d in (PLOT_DIR, RES_DIR):
        os.makedirs(d, exist_ok=True)

    path = args.csv
    if path is None:
        hits = sorted(globmod.glob(CAP_GLOB))
        if not hits:
            print(f"no sweep captures matching {CAP_GLOB}")
            sys.exit(1)
        path = hits[-1]
    elif '*' in path:
        hits = sorted(globmod.glob(path))
        if not hits:
            print(f"no match for {path}")
            sys.exit(1)
        path = hits[-1]

    t, x, fs = load_capture(path)
    xc = x - np.mean(x)
    print(f"capture : {path}")
    print(f"          {len(x)} samples, fs={fs:.4f} Hz, "
          f"{(t[-1]-t[0])/1e6:.2f} s\n")

    bursts, floor, thr, rms, w = find_bursts(xc, fs)
    print("=== burst detection ===")
    print(f"  quiet floor      : {floor:.2f} counts rms")
    print(f"  threshold ({BURST_THRESH_MULT:g}x)  : {thr:.2f}")
    print(f"  bursts found     : {len(bursts)}")
    for k, (i0, i1) in enumerate(bursts):
        print(f"    {k+1}: {i0/fs:6.2f} - {i1/fs:6.2f} s  "
              f"({(i1-i0)/fs:.2f} s)")

    if len(bursts) != len(TONES_HZ):
        print(f"\n  [STOP] expected {len(TONES_HZ)} bursts, found {len(bursts)}.")
        print("         Not guessing an alignment. Check the capture, or")
        print("         adjust BURST_THRESH_MULT if bursts merged or split.")
        sys.exit(1)

    starts = [b[0] / fs for b in bursts]
    spacing = np.diff(starts)
    print(f"\n  spacing          : {np.mean(spacing):.2f} +/- "
          f"{np.std(spacing):.2f} s   (manifest: 4.00 s)")
    print(f"  playback began   : ~{starts[0] - 1.0:.2f} s into the capture")
    if abs(np.mean(spacing) - 4.0) > 0.3:
        print("  [!] spacing does not match the manifest; alignment suspect")

    # ------------------------------------------------- per-burst analysis
    print("\n=== per-burst response ===")
    rows = []
    for k, ((i0, i1), f0) in enumerate(zip(bursts, TONES_HZ)):
        trim = int(EDGE_TRIM_S * fs)
        a, b = i0 + trim, i1 - trim
        if b - a < int(0.5 * fs):
            a, b = i0, i1
        seg = xc[a:b]

        nper = min(1024, len(seg))
        f, p = welch(seg, fs=fs, nperseg=nper)

        p_fund = band_power(f, p, f0)
        harms, p_harm_tot = [], 0.0
        for h in range(2, N_HARMONICS + 1):
            fh = f0 * h
            if fh > 0.95 * fs / 2:
                break
            ph = band_power(f, p, fh)
            harms.append((fh, ph))
            p_harm_tot += ph

        thd = 100.0 * p_harm_tot / (p_fund + 1e-30)
        rms_all = float(np.sqrt(np.mean(seg ** 2)))

        third = len(seg) // 3
        r1 = float(np.sqrt(np.mean(seg[:third] ** 2)))
        r3 = float(np.sqrt(np.mean(seg[-third:] ** 2)))
        ramp_db = 20 * np.log10((r3 + 1e-12) / (r1 + 1e-12))

        rows.append({'hz': f0, 'p_fund': p_fund, 'rms': rms_all,
                     'thd_pct': thd, 'ramp_db': ramp_db,
                     'harmonics': harms,
                     't0': i0 / fs, 't1': i1 / fs})

    ref = max(r['p_fund'] for r in rows)
    print(f"  {'Hz':>6} {'rel dB':>8} {'rms':>8} {'harm %':>8} "
          f"{'ramp dB':>8}   dominant harmonic")
    for r in rows:
        r['rel_db'] = 10 * np.log10((r['p_fund'] + 1e-30) / ref)
        if r['harmonics']:
            fh, ph = max(r['harmonics'], key=lambda z: z[1])
            dom = f"{fh:.0f} Hz at {10*np.log10((ph+1e-30)/(r['p_fund']+1e-30)):+.1f} dB"
        else:
            dom = '-'
        print(f"  {r['hz']:>6.0f} {r['rel_db']:>8.1f} {r['rms']:>8.1f} "
              f"{r['thd_pct']:>8.1f} {r['ramp_db']:>+8.1f}   {dom}")

    print("\n  rel dB  : narrowband power at the fundamental, vs the strongest burst")
    print("  harm %  : harmonic power as % of fundamental. HIGH means the driver")
    print("            is distorting, so 'response' at that tone is partly fake.")
    print("  ramp dB : last third vs first third. Large means the level was")
    print("            still changing, so the burst is not a steady measurement.")

    # -------------------------------------- weighted transmission of p21
    print("\n=== what this means for the p21 segment ===")
    npys = sorted(globmod.glob(os.path.join(STIM_DIR, '*_native333.npy')))
    if not npys:
        print("  no reference segment found; skipping")
    else:
        seg = np.load(npys[0])
        fs21 = 333.0
        s = seg - np.mean(seg)
        f21, p21 = welch(s, fs=fs21, nperseg=min(2048, len(s)))
        m = (f21 >= P21_BAND[0]) & (f21 <= P21_BAND[1])
        fb, pb = f21[m], p21[m]

        hz = np.array([r['hz'] for r in rows])
        db = np.array([r['rel_db'] for r in rows])
        # interpolate in log-frequency, clamp outside the measured range
        db_i = np.interp(np.log10(fb), np.log10(hz), db,
                         left=db[0], right=db[-1])
        gain = 10 ** (db_i / 10.0)

        trans = _trapz(pb * gain, fb) / _trapz(pb, fb)
        print(f"  reference       : {os.path.basename(npys[0])}")
        print(f"  in-band energy surviving the measured response: "
              f"{100*trans:.1f} %  ({10*np.log10(trans+1e-30):+.1f} dB)")

        for lo, hi in [(25, 50), (50, 100), (100, 158.175)]:
            mm = (fb >= lo) & (fb < hi)
            if not mm.any():
                continue
            share = 100 * _trapz(pb[mm], fb[mm]) / _trapz(pb, fb)
            surv = 100 * _trapz(pb[mm]*gain[mm], fb[mm]) / _trapz(pb*gain, fb)
            print(f"    {lo:>5.0f}-{hi:<6.0f} Hz : {share:5.1f} % of the original"
                  f"  ->  {surv:5.1f} % of what survives")
        print("\n  The band shares SHIFT: whatever the speaker favours becomes")
        print("  over-represented in what the microphone actually receives.")
        print("  Run C is therefore not a test of the original signal, but of")
        print("  the original signal as filtered by this playback chain -- and")
        print("  the pipeline needs BURST TIMING, which survives filtering, not")
        print("  spectral fidelity, which does not.")

    # ------------------------------------------------------------- plot
    fig, ax = plt.subplots(3, 1, figsize=(12, 10))
    ts = np.arange(len(xc)) / fs
    ax[0].plot(ts, xc, linewidth=0.3)
    for r, tone in zip(rows, TONES_HZ):
        ax[0].axvspan(r['t0'], r['t1'], alpha=0.12, color='green')
        ax[0].text((r['t0']+r['t1'])/2, ax[0].get_ylim()[1]*0.85,
                   f"{tone:.0f}", ha='center', fontsize=9)
    ax[0].set_title('sweep capture with detected bursts')
    ax[0].set_xlabel('Time (s)'); ax[0].grid(alpha=0.3)

    tr = np.arange(len(rms)) * RMS_WIN_S
    ax[1].semilogy(tr, rms + 1e-9, linewidth=1.0)
    ax[1].axhline(thr, color='r', linestyle='--', label='burst threshold')
    ax[1].axhline(floor, color='g', linestyle=':', label='quiet floor')
    ax[1].set_title('short-time RMS (burst detection)')
    ax[1].set_xlabel('Time (s)'); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    ax[2].plot([r['hz'] for r in rows], [r['rel_db'] for r in rows],
               'o-', markersize=7)
    ax[2].axhline(0, color='k', alpha=0.3)
    ax[2].axvspan(25, 50, alpha=0.12, color='orange',
                  label='64% of p21 in-band energy')
    ax[2].set_xlabel('Hz'); ax[2].set_ylabel('relative dB')
    ax[2].set_title('measured playback response')
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)

    plt.tight_layout()
    out = os.path.join(PLOT_DIR, 'quick_check_sweep.png')
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"\nplot -> {out}")

    with open(os.path.join(RES_DIR, 'quick_check_sweep.json'), 'w') as fh:
        json.dump({'capture': path, 'fs': fs, 'quiet_floor_rms': floor,
                   'bursts': [{k: v for k, v in r.items()
                               if k != 'harmonics'} for r in rows]},
                  fh, indent=2)


if __name__ == '__main__':
    main()
