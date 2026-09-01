"""
quick_check_level.py -- find the playback level and gain for run C.

Captures a short stretch of p21 playing at the microphone, judges it against
fixed criteria, and prints one verdict plus a suggested adjustment. Run it,
change one thing, run it again.

WHY THIS IS NOT 02_capture.py
02_capture's mains detection and settling split are meaningless on a
music-like, non-stationary signal, and a tuning loop wants a single verdict
rather than a page of statistics. Every attempt is appended to a log CSV, so
the tuning trajectory becomes evidence instead of something remembered.

CRITERIA, FIXED BEFORE THE FIRST RUN

  PEAK EXCURSION 400-900 counts from the DC bias.
    p21's crest factor is 6.27 (15.9 dB), so peaks sit far above RMS and
    headroom must be left for them. Above 900 risks clipping on transients;
    below 400 throws away SNR against a 4-count noise floor for no reason.

  CLIPPING detected two ways.
    ADC rails are the obvious one, but the MAX4466's own output stage clips
    LONG BEFORE the ADC reaches 0 or 4095, and that shows up as FLAT-TOPPING:
    many samples resting at the exact same extreme value. In clean noisy data
    the exact maximum should occur once or twice, not hundreds of times.

  CREST FACTOR tracked across attempts. This is the real clipping detector.
    If the volume goes up 6 dB and RMS follows by 6 dB but peak only rises 3,
    something in the chain is compressing. A FALLING crest factor means
    clipping no matter what the absolute peak value says. A single capture
    cannot show this; the log across attempts can.

  SNR against the measured quiet-room floor of 4.0 counts std
    (baselines measured 3.82 / 4.11 / 4.00). Below 20 dB is worth naming.

ADJUSTMENT ORDER -- volume first, then trimpot.
The sweep already showed the subwoofer distorting badly at 25 Hz (39.3%
harmonic content), and raising volume makes that worse. The trimpot raises
signal against the ADC floor WITHOUT touching speaker distortion, so it is the
better lever. But it also invalidates the noise-floor baseline it was measured
at, so if the trimpot moves, the quiet-room floor must be re-measured.

USAGE
    python quick_check_level.py --note "volume 40%, trimpot factory"
    python quick_check_level.py --seconds 30 --note "volume 60%, trimpot +1 turn"
"""
import os
import re
import sys
import csv
import time
import argparse
from datetime import datetime

import numpy as np
import serial

# ------------------------------------------------------------------ config
PORT = 'COM8'
BAUD = 115200
SECONDS = 20.0

EXPECT_PIN = 39
EXPECT_FS = 500

CAP_DIR = 'captures'
RES_DIR = 'results'
LOG_CSV = os.path.join(RES_DIR, 'level_finding_log.csv')

BASELINE_NOISE_STD = 4.0       # measured: 3.82 / 4.11 / 4.00 counts

PEAK_MIN = 400                 # counts from bias
PEAK_MAX = 900
PEAK_TARGET = 650              # middle of the band

RAIL_LOW, RAIL_HIGH = 50, 4045
FLATTOP_FRAC = 0.001           # >0.1% of samples at the exact extreme

SNR_WARN_DB = 20.0
HEADER_TIMEOUT_S = 25.0


def open_and_sync(port, baud):
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
    t0 = time.time()
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
    return (np.array(t_us, dtype=np.int64),
            np.array(adc, dtype=np.int64), bad)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', default=PORT)
    ap.add_argument('--seconds', type=float, default=SECONDS)
    ap.add_argument('--note', default='', help='volume and trimpot setting')
    args = ap.parse_args()

    for d in (CAP_DIR, RES_DIR):
        os.makedirs(d, exist_ok=True)

    print("Start p21 playing at the microphone, THEN let this run.\n")
    ser, header = open_and_sync(args.port, BAUD)
    print(f"  capturing {args.seconds:.0f} s...")
    try:
        t_us, adc, bad = capture(ser, args.seconds)
    finally:
        ser.close()

    if len(adc) < 100:
        print(f"[FAIL] only {len(adc)} samples")
        sys.exit(1)

    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    stem = f'pcg_level_{stamp}'
    cap = os.path.join(CAP_DIR, stem + '.csv')
    with open(cap, 'w') as fh:
        fh.write(header + '\n')
        fh.write(f'# level-finding {stamp} note={args.note}\n')
        fh.write('t_us,adc\n')
        for a, b in zip(t_us, adc):
            fh.write(f'{a},{b}\n')

    x = adc.astype(np.float64)
    fs = (len(t_us) - 1) / ((t_us[-1] - t_us[0]) / 1e6)
    bias = float(np.mean(x))
    xc = x - bias
    rms = float(np.sqrt(np.mean(xc ** 2)))
    peak = float(np.max(np.abs(xc)))
    crest = peak / (rms + 1e-12)
    snr_db = 20 * np.log10(rms / BASELINE_NOISE_STD)

    vmax, vmin = int(np.max(adc)), int(np.min(adc))
    n_at_max = int(np.sum(adc == vmax))
    n_at_min = int(np.sum(adc == vmin))
    ft_thresh = max(3, int(FLATTOP_FRAC * len(adc)))
    flattop = (n_at_max > ft_thresh) or (n_at_min > ft_thresh)
    railed = (vmin < RAIL_LOW) or (vmax > RAIL_HIGH)

    print(f"\n  note        : {args.note or '(none)'}")
    print(f"  samples     : {len(adc)}  fs={fs:.3f} Hz  "
          f"malformed={bad}")
    print(f"  bias        : {bias:.1f} counts")
    print(f"  rms         : {rms:.1f} counts")
    print(f"  peak excur. : {peak:.0f} counts   "
          f"(target {PEAK_MIN}-{PEAK_MAX})")
    print(f"  crest factor: {crest:.2f}   "
          f"(p21 source is 6.27 -- watch this ACROSS runs)")
    print(f"  min/max     : {vmin} / {vmax}")
    print(f"  at extremes : {n_at_min} at min, {n_at_max} at max  "
          f"(flat-top if > {ft_thresh})")
    print(f"  SNR         : {snr_db:.1f} dB over the {BASELINE_NOISE_STD:g}-count floor")

    if railed:
        verdict = 'CLIPPING (ADC rail)'
    elif flattop:
        verdict = 'CLIPPING (flat-topped)'
    elif peak > PEAK_MAX:
        verdict = 'TOO HIGH'
    elif peak < PEAK_MIN:
        verdict = 'TOO LOW'
    else:
        verdict = 'GOOD'

    print(f"\n  VERDICT     : {verdict}")
    if verdict == 'GOOD':
        print("  -> level is set. Do not touch the volume or the trimpot again")
        print("     before run C, and record both settings in the notes.")
    elif verdict.startswith('CLIPPING'):
        print("  -> come DOWN. Clipped peaks are destroyed transients, and")
        print("     transient TIMING is the only thing the pipeline uses.")
    else:
        adj = 20 * np.log10(PEAK_TARGET / (peak + 1e-12))
        print(f"  -> adjust by about {adj:+.1f} dB "
              f"({'up' if adj > 0 else 'down'}).")
        print("     Volume first, up to where the sub starts distorting;")
        print("     then the trimpot, which does not affect the speaker.")

    if snr_db < SNR_WARN_DB:
        print(f"\n  [!] SNR {snr_db:.1f} dB is below {SNR_WARN_DB:g} dB. If this")
        print("      cannot be improved, record it -- a weak run C would then")
        print("      be explained by level, not by the acquisition chain.")

    row = {'timestamp': stamp, 'note': args.note, 'n': len(adc),
           'fs': round(fs, 3), 'bias': round(bias, 1), 'rms': round(rms, 1),
           'peak_excursion': round(peak, 0), 'crest_factor': round(crest, 2),
           'snr_db': round(snr_db, 1), 'vmin': vmin, 'vmax': vmax,
           'n_at_min': n_at_min, 'n_at_max': n_at_max,
           'verdict': verdict, 'capture': os.path.basename(cap)}
    new = not os.path.exists(LOG_CSV)
    with open(LOG_CSV, 'a', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(row))
        if new:
            wr.writeheader()
        wr.writerow(row)

    print(f"\n  capture -> {cap}")
    print(f"  log     -> {LOG_CSV}")

    if not new:
        with open(LOG_CSV, newline='') as fh:
            rows = list(csv.DictReader(fh))
        if len(rows) > 1:
            print("\n  attempts so far:")
            print(f"    {'rms':>8} {'peak':>8} {'crest':>7} {'snr dB':>8}  note")
            for r in rows[-6:]:
                print(f"    {float(r['rms']):>8.1f} "
                      f"{float(r['peak_excursion']):>8.0f} "
                      f"{float(r['crest_factor']):>7.2f} "
                      f"{float(r['snr_db']):>8.1f}  {r['note']}")
            cf = [float(r['crest_factor']) for r in rows]
            pk = [float(r['peak_excursion']) for r in rows]
            if len(cf) >= 2 and pk[-1] > pk[-2] and cf[-1] < cf[-2] * 0.85:
                print("\n    [!] crest factor FELL while level rose. Something is")
                print("        compressing -- speaker, amplifier, or MAX4466.")
                print("        Come down regardless of the peak value.")


if __name__ == '__main__':
    main()
