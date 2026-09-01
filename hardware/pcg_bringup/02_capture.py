"""
02_capture.py -- MAX4466 channel characterisation.

Captures 90 s of raw PCG-channel data from pcg_capture.ino, validates it, and
measures the three things nothing so far could establish.

WHY THIS RUN EXISTS

1. ACHIEVED SAMPLE RATE, not nominal. The frozen pipeline computes
   BPM = 60 * fs / peak_lag, so fs scales every result linearly. A 0.2% rate
   error is 0.3 BPM at 140 -- in range, plausible, and wrong. SAMPLE_HZ in the
   firmware is a request; this measures what the hardware actually delivered
   over tens of thousands of samples.

2. NOISE SPECTRUM. pcg_smoke established the noise was ~3.5 counts but had no
   frequency information at all. Ghana mains is 50 Hz, which lies INSIDE the
   pipeline's 25-200 Hz passband -- the bandpass cannot remove it. A strong
   periodic line near the fetal S1 band is a direct false-positive risk for an
   autocorrelation detector, which by construction rewards periodicity.

3. A BASELINE. The RC anti-alias filter, the acoustic chamber and the loopback
   are all judged as deltas from these numbers.

PRE-REGISTERED BEFORE THE RUN
  - Mains harmonic "detected" if peak exceeds local background by >6 dB;
    "dominant" if >20 dB. Fixed now so the result cannot be reclassified after
    it is seen.
  - Settling is MEASURED, not assumed: statistics are reported for the full
    90 s and for the last 80 s separately. If they agree, settling is a
    non-issue and that is stated with evidence rather than asserted.
  - Nominal fs is 500 Hz. Achieved fs within 0.1% is "meets spec".

HARD GATES -- the script refuses to continue
  - Firmware header must match the expected pin and rate. Two pin/firmware
    mismatches were caught by hand today; this makes the check automatic.
  - Timestamps must be strictly increasing. Non-monotonic timestamps mean
    serial corruption or a dropped byte, which is exactly the failure mode
    that forced the drop from 230400 to 115200 baud.

The raw CSV is written before any analysis, so evidence survives even if the
analysis rejects the capture.

USAGE
    python 02_capture.py
    python 02_capture.py --port COM8 --seconds 90 --tag baseline_noRC
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
matplotlib.use('Agg')          # headless: set BEFORE pyplot is imported
import matplotlib.pyplot as plt

import serial

# ------------------------------------------------------------------ config
PORT = 'COM8'
BAUD = 115200
SECONDS = 90.0

EXPECT_PIN = 39
EXPECT_FS = 500
NOMINAL_FS = 500.0

CAP_DIR = 'captures'
RES_DIR = 'results'
PLOT_DIR = 'plots'

# Pipeline passband at fs=500: min(200, 0.95*250) = 200 Hz
BAND_LOW = 25.0
BAND_HIGH = 200.0

MAINS_HZ = [50.0, 100.0, 150.0]
MAINS_HALFWIDTH = 1.0          # +/- Hz counted as "the harmonic"
MAINS_BG_HALFWIDTH = 6.0       # local background window, harmonic excluded
MAINS_DETECT_DB = 6.0          # pre-registered
MAINS_DOMINANT_DB = 20.0       # pre-registered

SETTLE_SPLIT_S = 10.0          # full run vs last (SECONDS - this)

HEADER_TIMEOUT_S = 25.0

# toco_bringup measured the bare ADC noise floor at this, same rate, same
# attenuation. Carried here as the comparison point for whether the MAX4466
# contributes noise above the converter.
TOCO_ADC_NOISE_STD = 3.92


def _trapz(y, x):
    f = getattr(np, 'trapezoid', None) or np.trapz
    return float(f(y, x))


# ------------------------------------------------------------------ serial
def open_and_sync(port, baud):
    """Open the port, resynchronise on the firmware header, validate it.

    Opening the port toggles DTR/RTS, which resets the ESP32 on a CP2102 board,
    so a fresh header should follow. reset_input_buffer() clears the bootloader
    chatter -- without it the toco phase recorded a ~55 s stale-data prefix
    that looked like real signal.

    The header is then used as a SYNC MARKER: everything before it is
    discarded, so a partial line or leftover boot text cannot enter the data.
    """
    print(f"opening {port} at {baud}...")
    ser = serial.Serial(port, baud, timeout=1)
    time.sleep(0.2)
    ser.reset_input_buffer()

    print("waiting for firmware header...")
    t_start = time.time()
    prompted = False
    header = None

    while time.time() - t_start < HEADER_TIMEOUT_S:
        raw = ser.readline()
        if not raw:
            if not prompted and time.time() - t_start > 6.0:
                print("  no header yet -- press EN/RST on the ESP32")
                prompted = True
            continue
        line = raw.decode('utf-8', errors='replace').strip()
        if line.startswith('# pcg_capture'):
            header = line
            break

    if header is None:
        ser.close()
        print(f"\n[FAIL] no header within {HEADER_TIMEOUT_S:.0f} s.")
        print("       Is the Arduino Serial Monitor still holding the port?")
        sys.exit(1)

    print(f"  header: {header}")

    m_pin = re.search(r'pin=(\d+)', header)
    m_fs = re.search(r'fs=(\d+)', header)
    if not m_pin or not m_fs:
        ser.close()
        print("[FAIL] header present but unparseable.")
        sys.exit(1)

    pin, fs = int(m_pin.group(1)), int(m_fs.group(1))
    if pin != EXPECT_PIN or fs != EXPECT_FS:
        ser.close()
        print(f"\n[FAIL] firmware reports pin={pin} fs={fs}; "
              f"expected pin={EXPECT_PIN} fs={EXPECT_FS}.")
        print("       The firmware and the wiring disagree. Nothing below the")
        print("       header can be trusted until they match.")
        sys.exit(1)
    print(f"  [ok] pin={pin} fs={fs} as expected")

    # consume the column header
    col = ser.readline().decode('utf-8', errors='replace').strip()
    if col != 't_us,adc':
        print(f"  [!] expected 't_us,adc', got '{col}' -- continuing, but the")
        print("      line format may have changed")
    return ser, header


def capture(ser, seconds):
    """Read for `seconds`. Returns (t_us, adc, n_malformed, examples)."""
    print(f"\ncapturing {seconds:.0f} s -- keep the room quiet, do not touch "
          f"the mic or the bench")
    t_us, adc, malformed, examples = [], [], 0, []
    t_start = time.time()
    next_report = 10.0

    while time.time() - t_start < seconds:
        raw = ser.readline()
        if not raw:
            continue
        line = raw.decode('utf-8', errors='replace').strip()
        if not line or line.startswith('#'):
            continue
        try:
            a, b = line.split(',')
            t_us.append(int(a))
            adc.append(int(b))
        except ValueError:
            malformed += 1
            if len(examples) < 5:
                examples.append(line)
            continue

        el = time.time() - t_start
        if el >= next_report:
            print(f"  {el:5.0f} s  {len(t_us):7d} samples")
            next_report += 10.0

    return (np.array(t_us, dtype=np.int64),
            np.array(adc, dtype=np.int64), malformed, examples)


# ------------------------------------------------------------------ checks
def assert_monotonic(t_us):
    """Refuse to proceed on non-increasing timestamps."""
    d = np.diff(t_us)
    bad = np.where(d <= 0)[0]
    print("\n=== monotonicity gate ===")
    if len(bad) == 0:
        print(f"  [ok] all {len(d)} intervals strictly positive")
        return
    print(f"  [FAIL] {len(bad)} non-increasing interval(s)")
    for i in bad[:5]:
        print(f"    index {i}: {t_us[i]} -> {t_us[i+1]}  (d={d[i]})")
    print("\n  Non-increasing timestamps mean serial corruption or a dropped")
    print("  byte -- the failure that forced 230400 -> 115200. The raw CSV has")
    print("  been written; the analysis is refused.")
    sys.exit(1)


# ------------------------------------------------------------------ timing
def analyse_timing(t_us, label):
    d = np.diff(t_us).astype(np.float64)
    dur = (t_us[-1] - t_us[0]) / 1e6
    fs = (len(t_us) - 1) / dur
    err_pct = 100.0 * (fs - NOMINAL_FS) / NOMINAL_FS

    out = {
        'label': label, 'n': int(len(t_us)), 'duration_s': dur,
        'fs_achieved': fs, 'fs_error_pct': err_pct,
        'dt_mean_us': float(np.mean(d)), 'dt_std_us': float(np.std(d)),
        'dt_min_us': float(np.min(d)), 'dt_max_us': float(np.max(d)),
        'dt_p99_us': float(np.percentile(d, 99)),
        'n_gaps_2x': int(np.sum(d > 2 * 1e6 / NOMINAL_FS)),
    }

    print(f"\n=== timing [{label}] ===")
    print(f"  samples        : {out['n']}")
    print(f"  duration       : {dur:.3f} s")
    print(f"  achieved fs    : {fs:.4f} Hz   (nominal {NOMINAL_FS:g})")
    print(f"  rate error     : {err_pct:+.4f} %")
    print(f"  dt mean/std    : {out['dt_mean_us']:.2f} / "
          f"{out['dt_std_us']:.2f} us")
    print(f"  dt min/max     : {out['dt_min_us']:.0f} / "
          f"{out['dt_max_us']:.0f} us")
    print(f"  intervals >2x  : {out['n_gaps_2x']}  (dropped-sample gaps)")

    if abs(err_pct) < 0.1:
        print(f"  [ok] within the pre-registered 0.1% band")
    else:
        print(f"  [!] outside 0.1%. USE {fs:.4f} Hz as fs downstream, not 500.")
    return out


# ------------------------------------------------------------------ signal
def analyse_signal(adc, fs, label):
    x = adc.astype(np.float64)
    xc = x - np.mean(x)
    u = np.unique(adc)

    out = {
        'label': label,
        'mean': float(np.mean(x)), 'std': float(np.std(x)),
        'min': int(np.min(adc)), 'max': int(np.max(adc)),
        'p2p': int(np.max(adc) - np.min(adc)),
        'n_distinct': int(len(u)),
    }

    print(f"\n=== signal [{label}] ===")
    print(f"  mean           : {out['mean']:.2f} counts")
    print(f"  std            : {out['std']:.3f} counts")
    print(f"  min/max/p2p    : {out['min']} / {out['max']} / {out['p2p']}")
    print(f"  distinct values: {out['n_distinct']}")
    print(f"  ADC-only floor : {TOCO_ADC_NOISE_STD:.2f} (toco, same rate/att)")
    if out['std'] <= TOCO_ADC_NOISE_STD * 1.15:
        print("  -> at or below the bare-ADC floor: the amplifier contributes")
        print("     little above the converter, so added gain should buy SNR")
        print("     nearly one-for-one until amplifier noise appears.")
    else:
        print("  -> above the bare-ADC floor: something is contributing noise")
        print("     beyond the converter. The spectrum below should show what.")
    return out, xc


def analyse_spectrum(xc, fs, label):
    nper = min(2500, len(xc))          # ~0.2 Hz resolution at 500 Hz
    f, p = welch(xc, fs=fs, nperseg=nper)

    total = _trapz(p, f)
    inb = (f >= BAND_LOW) & (f <= BAND_HIGH)
    e_in = _trapz(p[inb], f[inb])

    print(f"\n=== spectrum [{label}] ===")
    print(f"  resolution     : {f[1]-f[0]:.3f} Hz  (nperseg={nper})")
    print(f"  in-band {BAND_LOW:g}-{BAND_HIGH:g} Hz: "
          f"{100*e_in/total:.1f} % of total power")

    print(f"\n  mains harmonics (pre-registered: >{MAINS_DETECT_DB:g} dB "
          f"detected, >{MAINS_DOMINANT_DB:g} dB dominant)")
    mains = []
    for f0 in MAINS_HZ:
        pk = (f >= f0 - MAINS_HALFWIDTH) & (f <= f0 + MAINS_HALFWIDTH)
        bg = ((f >= f0 - MAINS_BG_HALFWIDTH) & (f <= f0 + MAINS_BG_HALFWIDTH)
              & ~pk)
        if not pk.any() or not bg.any():
            continue
        peak = float(np.max(p[pk]))
        base = float(np.median(p[bg]))
        db = 10.0 * np.log10(peak / (base + 1e-30))
        e_h = _trapz(p[pk], f[pk])
        frac = 100.0 * e_h / e_in if e_in > 0 else 0.0

        if db > MAINS_DOMINANT_DB:
            verdict = 'DOMINANT'
        elif db > MAINS_DETECT_DB:
            verdict = 'detected'
        else:
            verdict = 'not detected'
        print(f"    {f0:6.1f} Hz : {db:+6.1f} dB over background, "
              f"{frac:5.2f} % of in-band power   {verdict}")
        mains.append({'hz': f0, 'db_over_bg': db,
                      'pct_of_inband': frac, 'verdict': verdict})

    tot_pct = sum(m['pct_of_inband'] for m in mains)
    print(f"    all harmonics: {tot_pct:.2f} % of in-band power")

    # Strongest in-band line, whatever it is -- an autocorrelation detector
    # rewards periodicity, so any narrow line in band is worth naming.
    ib = np.where(inb)[0]
    k = ib[int(np.argmax(p[ib]))]
    print(f"\n  strongest in-band line: {f[k]:.2f} Hz")

    return {'label': label, 'inband_pct_of_total': 100 * e_in / total,
            'mains': mains, 'mains_total_pct_inband': tot_pct,
            'peak_inband_hz': float(f[k])}, f, p


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', default=PORT)
    ap.add_argument('--seconds', type=float, default=SECONDS)
    ap.add_argument('--tag', default='baseline')
    args = ap.parse_args()

    for d in (CAP_DIR, RES_DIR, PLOT_DIR):
        os.makedirs(d, exist_ok=True)

    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    stem = f'pcg_char_{args.tag}_{stamp}'

    ser, header = open_and_sync(args.port, BAUD)
    try:
        t_us, adc, malformed, examples = capture(ser, args.seconds)
    finally:
        ser.close()

    if len(t_us) < 100:
        print(f"\n[FAIL] only {len(t_us)} samples captured.")
        sys.exit(1)

    # Write raw BEFORE analysis: evidence survives a rejected capture.
    cap_path = os.path.join(CAP_DIR, stem + '.csv')
    with open(cap_path, 'w') as fh:
        fh.write(header + '\n')
        fh.write(f'# captured {stamp} port={args.port} tag={args.tag}\n')
        fh.write('t_us,adc\n')
        for a, b in zip(t_us, adc):
            fh.write(f'{a},{b}\n')
    print(f"\nraw -> {cap_path}  ({len(t_us)} samples)")

    print("\n=== line integrity ===")
    if malformed:
        print(f"  [!] {malformed} malformed line(s) -- a dropped byte is the")
        print("      230400-baud failure mode reappearing at 115200.")
        for e in examples:
            print(f"      {e!r}")
    else:
        print("  [ok] no malformed lines")

    assert_monotonic(t_us)

    # Full run, then the tail -- settling measured rather than assumed.
    tim_full = analyse_timing(t_us, 'full')
    sig_full, xc_full = analyse_signal(adc, tim_full['fs_achieved'], 'full')
    spec_full, f, p = analyse_spectrum(xc_full, tim_full['fs_achieved'], 'full')

    t0 = t_us[0] + int(SETTLE_SPLIT_S * 1e6)
    m = t_us >= t0
    tim_tail = sig_tail = spec_tail = None
    if np.sum(m) > 1000:
        tail_lbl = f'last {args.seconds - SETTLE_SPLIT_S:.0f}s'
        tim_tail = analyse_timing(t_us[m], tail_lbl)
        sig_tail, xc_tail = analyse_signal(adc[m], tim_tail['fs_achieved'],
                                           tail_lbl)
        spec_tail, _, _ = analyse_spectrum(xc_tail, tim_tail['fs_achieved'],
                                           tail_lbl)

        print("\n=== settling ===")
        ds = abs(sig_tail['std'] - sig_full['std'])
        df = abs(tim_tail['fs_achieved'] - tim_full['fs_achieved'])
        print(f"  std   full {sig_full['std']:.3f} vs tail "
              f"{sig_tail['std']:.3f}   (delta {ds:.3f})")
        print(f"  fs    full {tim_full['fs_achieved']:.4f} vs tail "
              f"{tim_tail['fs_achieved']:.4f}   (delta {df:.4f})")
        if ds < 0.15 * sig_full['std'] and df < 0.05:
            print("  -> first 10 s do not differ materially. No settling")
            print("     exclusion is needed, and that is now measured.")
        else:
            print("  -> first 10 s DO differ. A lead-in exclusion should be")
            print("     pre-registered before any acoustic capture.")

    res = {'capture': stem, 'header': header, 'port': args.port,
           'tag': args.tag, 'requested_s': args.seconds,
           'malformed_lines': malformed,
           'timing_full': tim_full, 'signal_full': sig_full,
           'spectrum_full': spec_full,
           'timing_tail': tim_tail, 'signal_tail': sig_tail,
           'spectrum_tail': spec_tail,
           'preregistered': {
               'mains_detect_db': MAINS_DETECT_DB,
               'mains_dominant_db': MAINS_DOMINANT_DB,
               'fs_tolerance_pct': 0.1,
               'settle_split_s': SETTLE_SPLIT_S}}
    rpath = os.path.join(RES_DIR, stem + '.json')
    with open(rpath, 'w') as fh:
        json.dump(res, fh, indent=2)
    print(f"\nresults -> {rpath}")

    # ---- plot: for spotting the gross, not for measuring ----
    ts = (t_us - t_us[0]) / 1e6
    d = np.diff(t_us).astype(np.float64)
    fig, ax = plt.subplots(4, 1, figsize=(12, 12))

    ax[0].plot(ts, adc, linewidth=0.3)
    ax[0].set_title(f'{stem} -- raw, full capture')
    ax[0].set_xlabel('Time (s)'); ax[0].set_ylabel('ADC counts')
    ax[0].grid(alpha=0.3)

    w = ts <= 1.0
    ax[1].plot(ts[w], adc[w], '.-', linewidth=0.6, markersize=2)
    ax[1].set_title('first 1 s (50 Hz would show ~10-sample periodicity)')
    ax[1].set_xlabel('Time (s)'); ax[1].set_ylabel('ADC counts')
    ax[1].grid(alpha=0.3)

    ax[2].hist(d, bins=60)
    ax[2].set_title(f"sample interval: mean {np.mean(d):.2f} us, "
                    f"std {np.std(d):.2f} us")
    ax[2].set_xlabel('dt (us)'); ax[2].set_yscale('log')
    ax[2].grid(alpha=0.3)

    ax[3].semilogy(f, p, linewidth=0.8)
    for f0 in MAINS_HZ:
        ax[3].axvline(f0, color='r', alpha=0.4, linestyle='--')
    ax[3].axvspan(BAND_LOW, BAND_HIGH, alpha=0.08, color='green')
    ax[3].set_title('Welch PSD (green = pipeline passband, red = mains)')
    ax[3].set_xlabel('Hz'); ax[3].set_ylabel('PSD')
    ax[3].grid(alpha=0.3)

    plt.tight_layout()
    ppath = os.path.join(PLOT_DIR, f'02_char_{args.tag}_{stamp}.png')
    plt.savefig(ppath, dpi=120)
    plt.close()
    print(f"plot    -> {ppath}")


if __name__ == '__main__':
    main()
