"""
02_preload.py
FSR402 bring-up: static preload characterisation.

Answers two questions the finger-press capture could not:
  1. Where does the resting level sit under a realistic belt preload?
     (determines whether R_fixed gives a usable operating window)
  2. How much does the sensor creep under sustained static load?
     (creep lives in the same 30-60 s band as a real contraction)

Records one run per object. Nothing may touch the sensor during a run.
"""

import os
import sys
import time
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import serial

# ---- frozen constants -------------------------------------------------
PORT       = "COM8"
BAUD       = 115200
FS         = 500.0               # Hz, must match toco_capture.ino
R_FIXED    = 1000.0              # ohm, the pull-down actually fitted
VREF       = 3.3
FULL_SCALE = 4095.0

EMPTY_S    = 10.0                # unloaded, before placement
SETTLE_S   = 10.0                # discarded after onset: placement transient
LOADED_S   = 80.0                # loaded, creep measured after SETTLE_S

ONSET_COUNTS = 300               # load is present above this; unloaded reads 0
                                 # (ADC dead zone) so any real load clears it

# a usable resting window: clear of the ADC dead zone, with room above
WINDOW_LO, WINDOW_HI = 600, 2000

CAP_DIR, PLOT_DIR, RES_DIR = "captures", "plots", "results"


def adc_to_ohms(counts):
    """Invert the divider. counts -> FSR resistance in ohms."""
    counts = np.asarray(counts, dtype=float)
    v = counts / FULL_SCALE * VREF
    with np.errstate(divide="ignore", invalid="ignore"):
        r = R_FIXED * (VREF - v) / v
    return np.where(v > 0, r, np.inf)


def find_onset(t, x):
    """First sample where load is present. Returns t of onset, or 0.0."""
    hit = np.flatnonzero(x > ONSET_COUNTS)
    return float(t[hit[0]]) if len(hit) else 0.0


def open_port():
    try:
        ser = serial.Serial(PORT, BAUD, timeout=1.0)
    except serial.SerialException as e:
        print(f"  FAILED to open {PORT}: {e}")
        print("  Is the Arduino Serial Monitor still open?")
        sys.exit(1)
    time.sleep(2.0)                      # board resets on port open
    ser.reset_input_buffer()
    deadline = time.time() + 5.0
    while time.time() < deadline:        # consume the sketch header
        line = ser.readline().decode("utf-8", errors="replace").strip()
        if line.startswith("#"):
            print(f"  firmware: {line}")
        if line == "t_us,adc":
            return ser
    print("  WARNING: header not seen; proceeding anyway")
    return ser


def read_for(ser, seconds, banner=None):
    # The board streams continuously, but this script blocks on input() between
    # reads. The OS serial buffer fills during that wait and then overflows,
    # so without this flush a capture begins with stale backlog followed by a
    # gap of lost samples. Discard everything queued before timing starts.
    ser.reset_input_buffer()

    if banner:
        print(f"  {banner}", flush=True)
    rows = []
    t_start = time.time()
    last_tick = -1
    while True:
        elapsed = time.time() - t_start
        if elapsed >= seconds:
            break
        tick = int(elapsed) // 10
        if tick != last_tick:
            print(f"    {int(elapsed):3d}/{int(seconds)} s", flush=True)
            last_tick = tick
        line = ser.readline().decode("utf-8", errors="replace").strip()
        if not line or line.startswith("#") or line.startswith("t_us"):
            continue
        try:
            t_us, adc = line.split(",")
            rows.append((int(t_us), int(adc)))
        except ValueError:
            continue
    return np.array(rows, dtype=np.int64)


def check_continuity(label, data):
    """Report dropped samples. A gap means the host fell behind the board."""
    if len(data) < 2:
        return
    dt_ms = np.diff(data[:, 0]) / 1000.0
    nominal = 1000.0 / FS
    gaps = np.flatnonzero(dt_ms > 10 * nominal)
    span = (data[-1, 0] - data[0, 0]) / 1e6
    if len(gaps):
        print(f"    WARNING [{label}]: {len(gaps)} gap(s), "
              f"largest {dt_ms[gaps].max() / 1000.0:.2f} s")
    print(f"    [{label}] {len(data)} samples over {span:.1f} s "
          f"= {len(data) / max(span, 1e-9):.1f} Hz")


def analyse(label, empty, loaded):
    """Returns a dict of the numbers that decide the resistor choice."""
    e_adc = empty[:, 1].astype(float)
    t = (loaded[:, 0] - loaded[0, 0]) / 1e6
    x = loaded[:, 1].astype(float)

    # find load onset from the data rather than assuming it at t=0
    t_onset = find_onset(t, x)
    keep = t >= t_onset + SETTLE_S
    if keep.sum() < 100:                 # onset too late, or never loaded
        print(f"    WARNING [{label}]: only {int(keep.sum())} samples after "
              f"onset+{SETTLE_S:.0f}s; falling back to whole record")
        keep = np.ones_like(t, dtype=bool)
        t_onset = 0.0
    tc, xc = t[keep], x[keep]

    settled = xc[tc >= tc[-1] - 10.0]    # final 10 s
    onset   = xc[tc <= tc[0] + 5.0]      # first 5 s after settle window opens

    fit   = np.polyfit(tc, xc, 1)
    slope = fit[0]                       # counts/s
    resid = xc - np.polyval(fit, tc)

    r_settled = adc_to_ohms(settled.mean())
    headroom  = FULL_SCALE - settled.mean()

    return {
        "label":        label,
        "t_onset":      t_onset,
        "empty_mean":   e_adc.mean(),
        "empty_std":    e_adc.std(),
        "onset_mean":   onset.mean(),
        "settled_mean": settled.mean(),
        "settled_std":  settled.std(),
        "creep_cps":    slope,
        "creep_total":  slope * (tc[-1] - tc[0]),
        "creep_pct":    100.0 * slope * (tc[-1] - tc[0]) / max(settled.mean(), 1),
        "noise_std":    resid.std(),
        "r_ohms":       r_settled,
        "headroom":     headroom,
        "in_window":    WINDOW_LO <= settled.mean() <= WINDOW_HI,
    }


def print_row(a):
    print(f"\n  --- {a['label']} ---")
    print(f"  load onset at   {a['t_onset']:7.2f} s into the loaded record")
    print(f"  unloaded        {a['empty_mean']:7.1f}  (std {a['empty_std']:.2f})")
    print(f"  at settle+0s    {a['onset_mean']:7.1f}")
    print(f"  settled (last 10s) {a['settled_mean']:7.1f}  (std {a['settled_std']:.2f})")
    print(f"  implied R_FSR   {a['r_ohms']:7.0f} ohm")
    print(f"  creep           {a['creep_cps']:+7.3f} counts/s"
          f"  -> {a['creep_total']:+.1f} counts ({a['creep_pct']:+.1f} %)")
    print(f"  noise std       {a['noise_std']:7.2f} counts (detrended)")
    print(f"  headroom to FS  {a['headroom']:7.0f} counts")
    print(f"  resting window  {'OK' if a['in_window'] else 'OUT OF RANGE'}"
          f"  (target {WINDOW_LO}-{WINDOW_HI})")
    if a["empty_std"] > 50:
        print( "  NOTE: unloaded segment is not flat. The load was probably "
               "applied before the unloaded recording finished.")


def plot(runs, stem):
    fig, ax = plt.subplots(2, 1, figsize=(13, 8))

    for label, empty, loaded in runs:
        t = (loaded[:, 0] - loaded[0, 0]) / 1e6
        x = loaded[:, 1].astype(float)
        t_analysis = find_onset(t, x) + SETTLE_S

        line, = ax[0].plot(t, x, lw=0.6, label=label)

        n = (len(x) // 125) * 125         # 500 Hz -> 4 Hz box-car
        dec = x[:n].reshape(-1, 125).mean(axis=1)
        t_d = t[:n:125] + (125 / FS) / 2.0
        ax[1].plot(t_d, dec, lw=1.3, label=label, color=line.get_color())

        # mark where this run's analysis window actually opens, not a fixed
        # offset: onset is detected per run and differs between them
        for a in ax:
            a.axvline(t_analysis, color=line.get_color(), ls="--", lw=0.9)

    for a in ax:
        a.axhspan(WINDOW_LO, WINDOW_HI, color="tab:green", alpha=0.08)
        a.grid(alpha=0.3)
        a.legend(fontsize=8)

    ax[0].set_ylabel("ADC (raw, 500 Hz)")
    ax[0].set_title(f"{stem}   shaded = target resting window,"
                    f"  dashed = start of analysis window (onset + "
                    f"{SETTLE_S:.0f} s)")
    ax[1].set_ylabel("ADC (averaged, 4 Hz)")
    ax[1].set_xlabel("time since recording started (s)")

    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f"{stem}.png")
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


def main():
    for d in (CAP_DIR, PLOT_DIR, RES_DIR):
        os.makedirs(d, exist_ok=True)
    stem = "preload_" + datetime.now().strftime("%Y%m%d_%H%M%S")

    print("\nFSR402 static preload characterisation")
    print(f"R_fixed = {R_FIXED:.0f} ohm")
    print("Nothing may touch the sensor once an object is placed.\n")

    ser = open_port()
    runs, results = [], []

    with ser:
        while True:
            label = input("\nObject label (blank to finish): ").strip()
            if not label:
                break

            print(f"\n[{label}]  clear the sensor completely.")
            input("  Press Enter when the pad is bare and untouched...")
            empty = read_for(ser, EMPTY_S, f"recording {EMPTY_S:.0f} s unloaded")

            print(f"\n  Place '{label}' on the pad now, then take your hands off.")
            input("  Press Enter the moment it is placed and released...")
            loaded = read_for(ser, SETTLE_S + LOADED_S,
                              f"recording {SETTLE_S + LOADED_S:.0f} s loaded")

            if len(empty) < 100 or len(loaded) < 100:
                print("  Not enough data for this run, skipping.")
                continue

            check_continuity("unloaded", empty)
            check_continuity("loaded", loaded)

            a = analyse(label, empty, loaded)
            print_row(a)
            runs.append((label, empty, loaded))
            results.append(a)

            safe = "".join(c if c.isalnum() else "_" for c in label)
            np.savetxt(os.path.join(CAP_DIR, f"{stem}_{safe}.csv"),
                       loaded, fmt="%d", delimiter=",",
                       header="t_us,adc", comments="")

    if not results:
        print("\nNo runs recorded.")
        return

    res_path = os.path.join(RES_DIR, f"{stem}.csv")
    keys = list(results[0].keys())
    with open(res_path, "w") as f:
        f.write(",".join(keys) + "\n")
        for a in results:
            f.write(",".join(str(a[k]) for k in keys) + "\n")

    png = plot(runs, stem)
    print(f"\n  wrote {res_path}")
    print(f"  wrote {png}")
    print(f"  wrote {len(runs)} capture file(s) to {CAP_DIR}/")


if __name__ == "__main__":
    main()
