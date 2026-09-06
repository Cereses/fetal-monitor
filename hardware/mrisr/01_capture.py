"""
01_capture.py -- multi-rate ISR characterisation.

Captures the interleaved PCG/UC stream from mrisr_capture.ino, splits it into
two files, and runs the two measurements THIS PHASE EXISTS FOR.

WHY THIS PHASE EXISTS
A fetal heart rate deceleration is interpreted by its PHASE relative to the
contraction: early decelerations (nadir at the contraction peak) are benign,
late decelerations (nadir after the peak) indicate uteroplacental
insufficiency. Same shape, opposite meaning, timing is the only difference.
Independent clocks make that unrecoverable. The shared time base is what makes
the FHR/UC pairing -- cardiotocography -- mean anything.

THE TWO OPEN QUESTIONS, both against baselines already measured

  1. SENSOR_VN ERRATUM UNDER CHANNEL ALTERNATION.
     notes.md 6 cleared GPIO39 using SINGLE-CHANNEL CONTINUOUS sampling and
     stated explicitly that this does NOT clear it for the multi-rate ISR,
     because the erratum is most reported when the ADC alternates channels.
     This firmware alternates on every tick. This is that test.
     Discriminator: a glitch injects ISOLATED OUTLIERS, which inflates p2p far
     more than std. The RATIO is the metric, not either alone.
     Baseline: std 3.82 / 4.11 / 4.00, p2p/std 14.5-14.9.

  2. CROSS-CHANNEL SETTLING.
     The ESP32 has ONE SAR converter behind a multiplexer. A read too soon
     after a channel switch can return a value contaminated by the previous
     channel. PCG sits near 1909 and UC near 1444, so contamination would pull
     each toward the other and raise both noise floors.
     Baselines: PCG std 3.82-4.11 (notes 4), UC std 3.92 (toco phase).

  Neither is assumed. Both are compared against numbers already in hand.

UC IS WRITTEN IN 04_detect.py's EXISTING FORMAT
captures/contractions_<stamp>_uc4hz.csv, columns t_s,uc, float-valued. The
firmware sends the SUM of 125 samples; this script divides by 125.0. That
reproduces the toco phase's float means bit-for-bit, which matters: 04_detect
records that float values make uc_detector's flat_run_mask INERT on device
data. Integer values would silently REACTIVATE it -- same frozen detector,
different behaviour, no error. No new consumer, no format negotiation.

BLOCKING PROMPTS ARE THE KNOWN HAZARD
The toco phase lost 55.47 s of data when input() blocked while the ESP32 kept
streaming and the OS buffer overflowed. This script has no prompt during
capture, and calls reset_input_buffer() before the read begins.

USAGE
    python 01_capture.py
    python 01_capture.py --seconds 90 --tag baseline
"""
import os
import re
import sys
import time
import json
import argparse
from datetime import datetime

import numpy as np
from scipy.signal import welch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import serial

# ------------------------------------------------------------------ config
PORT = 'COM8'
BAUD = 115200
SECONDS = 90.0

EXPECT_PCG_PIN = 39
EXPECT_UC_PIN = 35
EXPECT_FS = 500
EXPECT_DECIM = 125
NOMINAL_FS = 500.0
FS_UC = 4.0

CAP_DIR, RES_DIR, PLOT_DIR = 'captures', 'results', 'plots'

# ---- baselines, all previously measured, none invented ----
PCG_STD_BASELINE = (3.82, 4.11)        # notes.md 4, three quiet-room runs
PCG_RATIO_BASELINE = (14.5, 14.9)      # notes.md 6, G34 and G39
UC_STD_BASELINE = 3.92                 # toco phase, 1 kOhm divider
PCG_BIAS_BASELINE = (1907.0, 1911.0)

STD_TOLERANCE = 1.5                    # pre-registered: >1.5x baseline = contaminated
RATIO_TOLERANCE = 1.3                  # pre-registered: >1.3x baseline = glitching

HEADER_TIMEOUT_S = 25.0


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
        if line.startswith('# mrisr_capture'):
            header = line
            break
    if header is None:
        ser.close()
        print("[FAIL] no header. Is the Arduino Serial Monitor holding the port?")
        sys.exit(1)
    print(f"  {header}")

    m_pcg = re.search(r'pcg=(\d+)@(\d+)', header)
    m_uc = re.search(r'uc=(\d+)@(\d+)', header)
    m_dec = re.search(r'decim=(\d+)', header)
    if not (m_pcg and m_uc and m_dec):
        ser.close()
        print("[FAIL] header present but unparseable")
        sys.exit(1)
    pcg_pin, fs = int(m_pcg.group(1)), int(m_pcg.group(2))
    uc_pin = int(m_uc.group(1))
    decim = int(m_dec.group(1))
    if (pcg_pin != EXPECT_PCG_PIN or uc_pin != EXPECT_UC_PIN
            or fs != EXPECT_FS or decim != EXPECT_DECIM):
        ser.close()
        print(f"[FAIL] firmware reports pcg={pcg_pin} uc={uc_pin} fs={fs} "
              f"decim={decim}; expected {EXPECT_PCG_PIN}/{EXPECT_UC_PIN}/"
              f"{EXPECT_FS}/{EXPECT_DECIM}")
        sys.exit(1)
    print(f"  [ok] pcg={pcg_pin} uc={uc_pin} fs={fs} decim={decim}")
    ser.readline()                     # consume the format comment
    return ser, header


def capture(ser, seconds):
    """Parse both line formats. PCG starts with a digit, UC with 'U'."""
    pt, pv, ut, uv, bad = [], [], [], [], 0
    ser.reset_input_buffer()
    t0, nxt = time.time(), 15.0
    while time.time() - t0 < seconds:
        raw = ser.readline()
        if not raw:
            continue
        line = raw.decode('utf-8', errors='replace').strip()
        if not line or line.startswith('#'):
            continue
        try:
            if line[0] == 'U':
                _, a, b = line.split(',')
                ut.append(int(a))
                uv.append(int(b))
            else:
                a, b = line.split(',')
                pt.append(int(a))
                pv.append(int(b))
        except ValueError:
            bad += 1
            continue
        el = time.time() - t0
        if el >= nxt:
            print(f"    {el:5.0f} s   pcg {len(pt):6d}   uc {len(ut):4d}")
            nxt += 15.0
    return (np.array(pt, dtype=np.int64), np.array(pv, dtype=np.int64),
            np.array(ut, dtype=np.int64), np.array(uv, dtype=np.int64), bad)


def channel_stats(x, label):
    x = x.astype(np.float64)
    d = np.abs(x - np.mean(x))
    return {'label': label, 'n': int(len(x)), 'mean': float(np.mean(x)),
            'std': float(np.std(x)), 'min': float(np.min(x)),
            'max': float(np.max(x)), 'p2p': float(np.max(x) - np.min(x)),
            'ratio': float((np.max(x) - np.min(x)) / (np.std(x) + 1e-12))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', default=PORT)
    ap.add_argument('--seconds', type=float, default=SECONDS)
    ap.add_argument('--tag', default='baseline')
    args = ap.parse_args()
    for d in (CAP_DIR, RES_DIR, PLOT_DIR):
        os.makedirs(d, exist_ok=True)

    ser, header = open_and_sync(args.port, BAUD)
    print(f"\n  capturing {args.seconds:.0f} s -- room quiet, hands off the "
          f"bench, do not disturb the spring clip")
    try:
        pt, pv, ut, uv, bad = capture(ser, args.seconds)
    finally:
        ser.close()

    if len(pt) < 1000 or len(ut) < 10:
        print(f"[FAIL] pcg {len(pt)} samples, uc {len(ut)} samples")
        sys.exit(1)

    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    stem = f'mrisr_{args.tag}_{stamp}'

    pcg_path = os.path.join(CAP_DIR, stem + '_pcg500.csv')
    with open(pcg_path, 'w') as fh:
        fh.write(header + '\n')
        fh.write(f'# {stamp} tag={args.tag}\n')
        fh.write('t_us,adc\n')
        for a, b in zip(pt, pv):
            fh.write(f'{a},{b}\n')

    # UC in 04_detect.py's exact format: t_s,uc float-valued
    uc_path = os.path.join(CAP_DIR, f'contractions_{stamp}_uc4hz.csv')
    uc_mean = uv.astype(np.float64) / float(EXPECT_DECIM)
    uc_t_s = ut.astype(np.float64) / 1e6
    with open(uc_path, 'w') as fh:
        fh.write('t_s,uc\n')
        for a, b in zip(uc_t_s, uc_mean):
            fh.write(f'{a:.6f},{b:.4f}\n')

    print(f"\npcg -> {pcg_path}  ({len(pt)} samples)")
    print(f"uc  -> {uc_path}  ({len(ut)} samples, 04_detect.py format)")

    # ---------------------------------------------------------- integrity
    print("\n=== integrity ===")
    print(f"  malformed lines : {bad}")
    dp = np.diff(pt)
    if np.any(dp <= 0):
        print(f"  [FAIL] {int(np.sum(dp <= 0))} non-increasing PCG timestamps."
              f" Raw files kept; analysis refused.")
        sys.exit(1)
    if np.any(np.diff(ut) <= 0):
        print("  [FAIL] non-increasing UC timestamps.")
        sys.exit(1)
    print("  [ok] both channels strictly increasing")

    fs_p = (len(pt) - 1) / ((pt[-1] - pt[0]) / 1e6)
    fs_u = (len(ut) - 1) / ((ut[-1] - ut[0]) / 1e6)
    print(f"  pcg achieved fs : {fs_p:.4f} Hz  (nominal {NOMINAL_FS:g})")
    print(f"  uc  achieved fs : {fs_u:.4f} Hz  (nominal {FS_UC:g})")
    print(f"  dt mean/std     : {np.mean(dp):.2f} / {np.std(dp):.2f} us")
    print(f"  gaps >2x        : {int(np.sum(dp > 2 * 1e6 / NOMINAL_FS))}")
    print(f"  ratio pcg:uc    : {len(pt)/len(ut):.2f}  "
          f"(expected {EXPECT_DECIM})")

    # The gate 04_detect.py applies. Ours should trim nothing.
    gap = float(np.max(np.diff(uc_t_s)))
    print(f"  largest uc gap  : {gap:.3f} s  "
          f"({'[ok] under 5 s, 04_detect trims nothing' if gap < 5.0 else '[!] 04_detect WILL TRIM'})")

    # ------------------------------------------------ the two measurements
    p = channel_stats(pv, 'PCG')
    u = channel_stats(uv.astype(np.float64) / EXPECT_DECIM, 'UC')

    print("\n" + "=" * 68)
    print("MEASUREMENT 1 -- SENSOR_VN erratum under channel alternation")
    print("  notes.md 6 cleared GPIO39 SINGLE-CHANNEL ONLY and deferred this.")
    print("  Glitches inflate p2p far more than std, so the RATIO decides.")
    print("=" * 68)
    print(f"  {'':<18} {'measured':>10} {'baseline':>16}")
    print(f"  {'PCG std':<18} {p['std']:>10.2f} "
          f"{f'{PCG_STD_BASELINE[0]:.2f}-{PCG_STD_BASELINE[1]:.2f}':>16}")
    print(f"  {'PCG p2p':<18} {p['p2p']:>10.0f} {'120-131':>16}")
    print(f"  {'PCG p2p/std':<18} {p['ratio']:>10.2f} "
          f"{f'{PCG_RATIO_BASELINE[0]:.1f}-{PCG_RATIO_BASELINE[1]:.1f}':>16}")

    ratio_hi = PCG_RATIO_BASELINE[1] * RATIO_TOLERANCE
    erratum = p['ratio'] > ratio_hi
    print(f"\n  pre-registered: erratum ACTIVE if ratio > {ratio_hi:.1f} "
          f"({RATIO_TOLERANCE}x baseline)")
    print(f"  VERDICT: {'ERRATUM ACTIVE' if erratum else 'not detected'}")
    if erratum:
        print("  -> isolated outliers present that single-channel sampling did")
        print("     not produce. Options: move PCG to GPIO36, add a settling")
        print("     delay between reads, or restructure the read order.")

    print("\n" + "=" * 68)
    print("MEASUREMENT 2 -- cross-channel settling")
    print("  One SAR converter behind a mux. A read too soon after switching")
    print("  carries contamination from the previous channel. PCG ~1909 and")
    print("  UC ~1444, so contamination pulls each toward the other.")
    print("=" * 68)
    print(f"  {'':<18} {'measured':>10} {'baseline':>16}")
    print(f"  {'PCG bias':<18} {p['mean']:>10.1f} "
          f"{f'{PCG_BIAS_BASELINE[0]:.0f}-{PCG_BIAS_BASELINE[1]:.0f}':>16}")
    print(f"  {'PCG std':<18} {p['std']:>10.2f} "
          f"{f'{PCG_STD_BASELINE[0]:.2f}-{PCG_STD_BASELINE[1]:.2f}':>16}")
    print(f"  {'UC mean':<18} {u['mean']:>10.1f} {'(preload dependent)':>16}")
    print(f"  {'UC std':<18} {u['std']:>10.2f} "
          f"{f'{UC_STD_BASELINE:.2f}':>16}")

    pcg_hi = PCG_STD_BASELINE[1] * STD_TOLERANCE
    uc_hi = UC_STD_BASELINE * STD_TOLERANCE
    pcg_bad = p['std'] > pcg_hi
    uc_bad = u['std'] > uc_hi
    print(f"\n  pre-registered: contaminated if std > {STD_TOLERANCE}x baseline")
    print(f"    PCG {p['std']:.2f} vs {pcg_hi:.2f}   "
          f"{'CONTAMINATED' if pcg_bad else 'ok'}")
    print(f"    UC  {u['std']:.2f} vs {uc_hi:.2f}   "
          f"{'CONTAMINATED' if uc_bad else 'ok'}")
    if not (pcg_bad or uc_bad):
        print("  -> both floors at their single-channel values. analogRead()")
        print("     settles adequately between channels at this rate; no")
        print("     inter-read delay is needed, and that is now measured.")

    print(f"\n  PCG bias {'ok' if PCG_BIAS_BASELINE[0] <= p['mean'] <= PCG_BIAS_BASELINE[1] else 'OUTSIDE baseline range'}")
    print(f"  UC drift over capture: "
          f"{np.mean(u_last := uc_mean[-len(uc_mean)//4:]) - np.mean(uc_mean[:len(uc_mean)//4]):+.1f} counts "
          f"(last quarter vs first)")

    # ------------------------------------------------------------- plot
    fig, ax = plt.subplots(4, 1, figsize=(12, 12))
    tp = (pt - pt[0]) / 1e6
    tu = (ut - ut[0]) / 1e6

    ax[0].plot(tp, pv, linewidth=0.3)
    ax[0].set_title(f'{stem} -- PCG (GPIO39, 500 Hz)')
    ax[0].set_xlabel('Time (s)'); ax[0].set_ylabel('ADC'); ax[0].grid(alpha=0.3)

    ax[1].plot(tu, uc_mean, 'o-', markersize=2, linewidth=0.8)
    ax[1].set_title(f'UC (GPIO35, {fs_u:.3f} Hz, mean of {EXPECT_DECIM})')
    ax[1].set_xlabel('Time (s)'); ax[1].set_ylabel('counts'); ax[1].grid(alpha=0.3)

    ax[2].hist(dp, bins=60)
    ax[2].set_title(f'PCG sample interval: mean {np.mean(dp):.2f} us, '
                    f'std {np.std(dp):.2f} us')
    ax[2].set_xlabel('dt (us)'); ax[2].set_yscale('log'); ax[2].grid(alpha=0.3)

    f, ps = welch(pv - np.mean(pv), fs=fs_p, nperseg=min(2500, len(pv)))
    ax[3].semilogy(f, ps, linewidth=0.8)
    ax[3].axvspan(25, 200, alpha=0.08, color='green')
    for h in (50, 100, 150):
        ax[3].axvline(h, color='r', alpha=0.35, linestyle='--')
    ax[3].set_title('PCG PSD (green = pipeline passband, red = mains)')
    ax[3].set_xlabel('Hz'); ax[3].grid(alpha=0.3)

    plt.tight_layout()
    ppath = os.path.join(PLOT_DIR, f'01_{stem}.png')
    plt.savefig(ppath, dpi=120)
    plt.close()

    res = {'stem': stem, 'header': header, 'tag': args.tag,
           'malformed': bad, 'pcg_file': pcg_path, 'uc_file': uc_path,
           'fs_pcg': fs_p, 'fs_uc': fs_u,
           'dt_mean_us': float(np.mean(dp)), 'dt_std_us': float(np.std(dp)),
           'n_gaps_2x': int(np.sum(dp > 2 * 1e6 / NOMINAL_FS)),
           'largest_uc_gap_s': gap,
           'pcg': p, 'uc': u,
           'erratum_active': bool(erratum),
           'pcg_contaminated': bool(pcg_bad), 'uc_contaminated': bool(uc_bad),
           'baselines': {'pcg_std': PCG_STD_BASELINE,
                         'pcg_ratio': PCG_RATIO_BASELINE,
                         'uc_std': UC_STD_BASELINE,
                         'pcg_bias': PCG_BIAS_BASELINE},
           'preregistered': {'std_tolerance': STD_TOLERANCE,
                             'ratio_tolerance': RATIO_TOLERANCE}}
    rpath = os.path.join(RES_DIR, stem + '.json')
    with open(rpath, 'w') as fh:
        json.dump(res, fh, indent=2)
    print(f"\nresults -> {rpath}")
    print(f"plot    -> {ppath}")


if __name__ == '__main__':
    main()
