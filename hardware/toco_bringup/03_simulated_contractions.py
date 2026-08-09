"""
03_simulated_contractions.py
FSR402 bring-up: contraction-shaped stimulus capture.

Spring clip supplies static preload; hand pressure traces slow ramp-hold-release
cycles approximating uterine contraction morphology. Writes both the raw
500 Hz stream and a 4 Hz decimated channel for the frozen CTU-CHB detector.

PROTOCOL LENGTH MATTERS. uc_detector.BASELINE_WIN_S is 600 s, so a record
shorter than 10 minutes makes percentile_filter degenerate to a global
percentile and the rolling baseline cannot roll. The 5-minute run on
2026-08-08 tracked only 14% of a 410-count hysteresis ratchet for exactly this
reason. LONG is the default so the baseline mechanism is genuinely exercised.
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
PORT, BAUD = "COM8", 115200
FS         = 500.0
DECIM      = 125                 # 500 Hz -> 4 Hz, matches CTU-CHB
R_FIXED, VREF, FULL = 1000.0, 3.3, 4095.0

PRELOAD_LO, PRELOAD_HI = 600, 1500      # resting window with room above
PEAK_TARGET            = 3000           # aim for ramps to reach roughly here

# --- protocol ----------------------------------------------------------
# SHORT (5 min, superseded):  30, 3, 15, 20, 15, 40
# LONG  (20 min, default):    60, 7, 20, 25, 20, 100
BASELINE_S = 60.0
N_CYCLES   = 7
RAMP_S, HOLD_S, FALL_S, REST_S = 20.0, 25.0, 20.0, 100.0

GAP_WARN_S = 0.5                 # a timestamp step larger than this is a drop

CAP_DIR, PLOT_DIR = "captures", "plots"


def open_port():
    try:
        ser = serial.Serial(PORT, BAUD, timeout=1.0)
    except serial.SerialException as e:
        print(f"  FAILED to open {PORT}: {e}")
        print("  Is the Arduino Serial Monitor still open?")
        sys.exit(1)
    time.sleep(2.0)
    ser.reset_input_buffer()
    deadline = time.time() + 5.0
    while time.time() < deadline:
        line = ser.readline().decode("utf-8", errors="replace").strip()
        if line.startswith("#"):
            print(f"  firmware: {line}")
        if line == "t_us,adc":
            return ser
    print("  WARNING: header not seen; proceeding anyway")
    return ser


def read_for(ser, seconds, phase_fn=None):
    """Read for `seconds`. phase_fn(elapsed) may print a cue; called ~4 Hz."""
    # The board streams continuously while this script blocks on input(). The
    # OS serial buffer fills during that wait and then overflows, so without
    # this flush a capture opens with stale backlog followed by lost samples.
    # Observed on 2026-08-08 as a 55.47 s gap and an 86.35 s cue offset.
    ser.reset_input_buffer()

    rows, t_start, last = [], time.time(), -1.0
    while True:
        elapsed = time.time() - t_start
        if elapsed >= seconds:
            break
        if phase_fn and elapsed - last >= 0.25:
            phase_fn(elapsed)
            last = elapsed
        line = ser.readline().decode("utf-8", errors="replace").strip()
        if not line or line[0] in "#t":
            continue
        try:
            a, b = line.split(",")
            rows.append((int(a), int(b)))
        except ValueError:
            continue
    return np.array(rows, dtype=np.int64)


def check_continuity(data):
    """Report dropped samples. Over 20 minutes a single OS hiccup can open a
    gap, and a gap invalidates the decimation across it."""
    if len(data) < 2:
        return
    dt = np.diff(data[:, 0]) / 1e6
    span = (data[-1, 0] - data[0, 0]) / 1e6
    gaps = np.flatnonzero(dt > GAP_WARN_S)
    print(f"  continuity   {len(data)} samples over {span:.1f} s "
          f"= {len(data)/max(span,1e-9):.2f} Hz")
    if len(gaps):
        print(f"  WARNING: {len(gaps)} gap(s) over {GAP_WARN_S} s, "
              f"largest {dt[gaps].max():.2f} s")
        for i in gaps[:5]:
            t0 = (data[i, 0] - data[0, 0]) / 1e6
            print(f"    gap at t={t0:.1f} s, {dt[i]:.2f} s lost")
    else:
        print(f"  continuity   no gaps over {GAP_WARN_S} s")


def adjust_preload(ser):
    """Reposition the clip until the resting level leaves room for ramps."""
    print("\n--- preload adjustment ---")
    print(f"  target resting level {PRELOAD_LO}-{PRELOAD_HI} counts")
    print("  slide the clip partly off the pad to reduce it\n")
    while True:
        input("  Hands off, then press Enter to sample 3 s...")
        d = read_for(ser, 3.0)
        if len(d) < 100:
            print("  no data; check the port")
            continue
        x = d[:, 1].astype(float)
        v = x.mean() / FULL * VREF
        r = R_FIXED * (VREF - v) / v if v > 0 else float("inf")
        print(f"    mean {x.mean():7.1f}   std {x.std():5.2f}"
              f"   range {x.min():.0f}-{x.max():.0f}   R_FSR ~ {r:.0f} ohm")
        if x.mean() < PRELOAD_LO:
            print("    TOO LOW  near the ADC dead zone: tighten the clip\n")
        elif x.mean() > PRELOAD_HI:
            print(f"    TOO HIGH  only {FULL - x.mean():.0f} counts of headroom:"
                  f" slide the clip further off\n")
        else:
            print(f"    IN RANGE  {FULL - x.mean():.0f} counts of headroom\n")
        if input("  Accept and start the protocol? [y/N] ").strip().lower() == "y":
            return x.mean()


def build_schedule():
    """(t_start, t_end, phase) covering the whole protocol."""
    sched, t = [], 0.0
    sched.append((t, t + BASELINE_S, "REST")); t += BASELINE_S
    for i in range(N_CYCLES):
        sched.append((t, t + RAMP_S, f"RAMP {i+1}"));  t += RAMP_S
        sched.append((t, t + HOLD_S, f"HOLD {i+1}"));  t += HOLD_S
        sched.append((t, t + FALL_S, f"FALL {i+1}"));  t += FALL_S
        sched.append((t, t + REST_S, f"REST {i+1}"));  t += REST_S
    return sched, t


def main():
    for d in (CAP_DIR, PLOT_DIR):
        os.makedirs(d, exist_ok=True)
    stem = "contractions_" + datetime.now().strftime("%Y%m%d_%H%M%S")

    ser = open_port()
    with ser:
        rest_level = adjust_preload(ser)
        sched, total = build_schedule()

        print(f"\n--- protocol: {total:.0f} s ({total/60:.1f} min), "
              f"{N_CYCLES} cycles ---")
        print(f"  RAMP: increase pressure smoothly toward ~{PEAK_TARGET} counts")
        print( "  HOLD: keep it steady    FALL: release smoothly, not abruptly")
        print( "  REST: hands completely off the sensor")
        print(f"  expected rate: {N_CYCLES/(total/60)*10:.2f} per 10 min\n")
        input("  Press Enter to begin...")

        state = {"i": -1}

        def cue(elapsed):
            for k, (a, b, name) in enumerate(sched):
                if a <= elapsed < b:
                    if k != state["i"]:
                        rem = (total - elapsed) / 60.0
                        print(f"  [{elapsed:7.1f}s] {name:9s} "
                              f"({rem:4.1f} min left)", flush=True)
                        state["i"] = k
                    return

        data = read_for(ser, total, cue)

    if len(data) < 1000:
        print("\n  FAILED: almost no data.")
        sys.exit(1)

    t_us, adc = data[:, 0], data[:, 1]
    t = (t_us - t_us[0]) / 1e6
    x = adc.astype(float)

    n = (len(x) // DECIM) * DECIM
    uc = x[:n].reshape(-1, DECIM).mean(axis=1)
    t_uc = t[:n:DECIM] + (DECIM / FS) / 2.0

    print(f"\n=== capture ===")
    check_continuity(data)
    print(f"  duration     {t[-1]:.1f} s ({t[-1]/60:.2f} min)")
    print(f"  raw range    {x.min():.0f} .. {x.max():.0f}")
    print(f"  saturated    {100*(x >= FULL).mean():.2f} %")
    print(f"  at zero      {100*(x == 0).mean():.2f} %")
    print(f"  preload ref  {rest_level:.0f}")
    print(f"  4 Hz samples {len(uc)}   range {uc.min():.0f} .. {uc.max():.0f}")

    # Baseline ratchet: compare each REST phase to the first.
    print(f"\n=== baseline by rest phase ===")
    rest_means = []
    for a, b, name in sched:
        if not name.startswith("REST"):
            continue
        m = (t >= a + 10) & (t < b - 2)      # skip recovery transient
        if m.sum() > 100:
            rest_means.append((name, float(x[m].mean())))
    if rest_means:
        ref = rest_means[0][1]
        for name, v in rest_means:
            print(f"  {name:9s} {v:7.1f}   {v-ref:+7.1f} from first")

    raw_path = os.path.join(CAP_DIR, f"{stem}_raw500.csv")
    np.savetxt(raw_path, data, fmt="%d", delimiter=",",
               header="t_us,adc", comments="")

    uc_path = os.path.join(CAP_DIR, f"{stem}_uc4hz.csv")
    np.savetxt(uc_path, np.column_stack([t_uc, uc]), fmt="%.4f",
               delimiter=",", header="t_s,uc", comments="")

    cue_path = os.path.join(CAP_DIR, f"{stem}_protocol.csv")
    with open(cue_path, "w") as f:
        f.write("t_start_s,t_end_s,phase\n")
        for a, b, name in sched:
            f.write(f"{a:.3f},{b:.3f},{name}\n")

    fig, ax = plt.subplots(2, 1, figsize=(15, 8), sharex=True)
    ax[0].plot(t, x, lw=0.3)
    ax[0].set_ylabel("ADC (raw, 500 Hz)")
    ax[0].set_title(stem)
    ax[1].plot(t_uc, uc, lw=1.0, color="tab:red")
    ax[1].set_ylabel("UC (averaged, 4 Hz)")
    ax[1].set_xlabel("time (s)")
    for a_ in ax:
        for s, e, name in sched:
            if name.startswith(("RAMP", "HOLD", "FALL")):
                a_.axvspan(s, e, color="tab:orange", alpha=0.08)
        a_.axhline(rest_level, color="0.4", ls=":", lw=0.8)
        a_.grid(alpha=0.3)
    fig.tight_layout()
    png = os.path.join(PLOT_DIR, f"{stem}.png")
    fig.savefig(png, dpi=110)
    plt.close(fig)

    size_mb = os.path.getsize(raw_path) / 1e6
    print(f"\n  wrote {raw_path}  ({size_mb:.1f} MB)")
    print(f"  wrote {uc_path}")
    print(f"  wrote {cue_path}")
    print(f"  wrote {png}")
    if size_mb > 5:
        print(f"\n  NOTE: the raw file is {size_mb:.0f} MB. Commit the 4 Hz "
              f"file and the plot;")
        print( "  keep the raw one out of git unless it is needed as evidence.")


if __name__ == "__main__":
    main()
