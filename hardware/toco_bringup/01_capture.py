"""
01_capture.py
FSR402 tocodynamometer bring-up: structured 90 s characterisation capture.

Reads the toco_capture.ino serial stream (t_us,adc @ 500 Hz), writes a
timestamped CSV, logs the cue protocol, and saves a diagnostic plot.
"""

import os
import sys
import time
from datetime import datetime

import matplotlib
matplotlib.use("Agg")            # save-only, no window (avoids the 3.14 crash)
import matplotlib.pyplot as plt
import numpy as np
import serial

# ---- frozen constants -------------------------------------------------
PORT       = "COM8"
BAUD       = 115200
FS_NOMINAL = 500.0               # Hz, must match the sketch
DURATION_S = 90.0

CUES = [                         # (t_s, label)
    ( 0.0, "REST     baseline, hands off"),
    (20.0, "PRESS    light, hold steady"),
    (35.0, "RELEASE"),
    (45.0, "PRESS    light, hold steady"),
    (60.0, "RELEASE"),
    (70.0, "PRESS    light, hold steady"),
    (85.0, "RELEASE  hands off until end"),
]

CAP_DIR  = "captures"
PLOT_DIR = "plots"


def wait_for_header(ser, timeout_s=5.0):
    """Opening the port resets the ESP32, so the sketch reprints its header.
    Consuming it here discards boot noise and starts us on clean data."""
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        line = ser.readline().decode("utf-8", errors="replace").strip()
        if line.startswith("#"):
            print(f"  firmware: {line}")
        if line == "t_us,adc":
            return True
    print("  WARNING: header not seen; proceeding anyway")
    return False


def capture(ser):
    rows, cue_log = [], []
    next_cue = 0
    t_start = time.time()

    while True:
        elapsed = time.time() - t_start
        if elapsed >= DURATION_S:
            break

        if next_cue < len(CUES) and elapsed >= CUES[next_cue][0]:
            t_cue, label = CUES[next_cue]
            print(f"  [{elapsed:5.1f}s] {label}", flush=True)
            cue_log.append((elapsed, label.split()[0]))
            next_cue += 1

        line = ser.readline().decode("utf-8", errors="replace").strip()
        if not line or line.startswith("#") or line.startswith("t_us"):
            continue
        try:
            t_us, adc = line.split(",")
            rows.append((int(t_us), int(adc)))
        except ValueError:
            continue                      # partial line, drop it

    return np.array(rows, dtype=np.int64), cue_log


def report(t_us, adc):
    t_s = (t_us - t_us[0]) / 1e6
    dt_ms = np.diff(t_us) / 1000.0

    print(f"\n  samples        {len(adc)}")
    print(f"  duration       {t_s[-1]:.2f} s")
    print(f"  actual rate    {len(adc) / t_s[-1]:.2f} Hz  (nominal {FS_NOMINAL:.0f})")
    print(f"  interval       mean {dt_ms.mean():.3f} ms, std {dt_ms.std():.3f} ms")
    print(f"  dropped slots  {int((dt_ms > 3.0).sum())}  (intervals > 3 ms)")
    print(f"  adc range      {adc.min()} .. {adc.max()}")
    print(f"  at zero        {100.0 * (adc == 0).mean():.2f} %")
    print(f"  saturated      {100.0 * (adc == 4095).mean():.2f} %")
    return t_s


def plot(t_s, adc, cue_log, stem):
    # 125-sample box-car -> 4 Hz, the rate the CTU-CHB pipeline expects
    n = (len(adc) // 125) * 125
    dec = adc[:n].reshape(-1, 125).mean(axis=1)
    t_dec = t_s[:n:125] + (125 / FS_NOMINAL) / 2.0

    fig, ax = plt.subplots(2, 1, figsize=(13, 7), sharex=True)

    ax[0].plot(t_s, adc, lw=0.4)
    ax[0].set_ylabel("ADC (raw, 500 Hz)")
    ax[0].set_title(stem)

    ax[1].plot(t_dec, dec, lw=1.2, color="tab:red")
    ax[1].set_ylabel("ADC (averaged, 4 Hz)")
    ax[1].set_xlabel("time (s)")

    for a in ax:
        for t_cue, label in cue_log:
            a.axvline(t_cue, color="0.6", ls="--", lw=0.7)
        a.grid(alpha=0.3)

    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f"{stem}.png")
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


def main():
    os.makedirs(CAP_DIR, exist_ok=True)
    os.makedirs(PLOT_DIR, exist_ok=True)
    stem = "toco_" + datetime.now().strftime("%Y%m%d_%H%M%S")

    print(f"Opening {PORT} at {BAUD}...")
    try:
        ser = serial.Serial(PORT, BAUD, timeout=1.0)
    except serial.SerialException as e:
        print(f"  FAILED: {e}")
        print("  Is the Arduino Serial Monitor still open? Close it and retry.")
        sys.exit(1)

    with ser:
        time.sleep(2.0)                   # let the board finish resetting
        ser.reset_input_buffer()
        wait_for_header(ser)
        print(f"\nCapturing {DURATION_S:.0f} s. Follow the cues:\n")
        data, cue_log = capture(ser)

    if len(data) < 100:
        print("\n  FAILED: almost no data. Check the port and the sketch.")
        sys.exit(1)

    t_us, adc = data[:, 0], data[:, 1]
    t_s = report(t_us, adc)

    csv_path = os.path.join(CAP_DIR, f"{stem}.csv")
    np.savetxt(csv_path, data, fmt="%d", delimiter=",",
               header="t_us,adc", comments="")

    cue_path = os.path.join(CAP_DIR, f"{stem}_cues.csv")
    with open(cue_path, "w") as f:
        f.write("t_s,cue\n")
        for t_cue, label in cue_log:
            f.write(f"{t_cue:.3f},{label}\n")

    png_path = plot(t_s, adc, cue_log, stem)
    print(f"\n  wrote {csv_path}")
    print(f"  wrote {cue_path}")
    print(f"  wrote {png_path}")


if __name__ == "__main__":
    main()