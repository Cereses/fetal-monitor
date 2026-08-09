"""
04_detect.py — run the FROZEN CTU-CHB contraction detector against FSR402
device data, and score it against the capture protocol.

This imports uc_detector unmodified. Nothing about the algorithm or its
constants is redefined here: if the detector is edited, this script inherits
the edit. Reimplementing it would silently make it a different algorithm and
void the comparison, which is the whole point of the exercise.

WHAT THIS DOES AND DOES NOT ESTABLISH
  DOES     shows whether device-produced signals have morphology the frozen
           detector recognises, and where it breaks on device-specific
           artifacts (baseline ratchet, short records, float-valued samples).
  DOES NOT produce a performance figure comparable to the CTU-CHB result. The
           protocol is a scripted hand stimulus, not expert-annotated labour.
           The 77.32% F1 belongs to CTU-CHB and cannot be claimed here.

Usage:
    python 04_detect.py                     # newest capture in captures/
    python 04_detect.py <stem>              # e.g. contractions_20260808_151139
"""

import os
import sys
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def _locate_detector():
    """Find uc_detector.py without copying it.

    It lives in the software pipeline, not in hardware/toco_bringup/. Importing
    the real file (rather than a duplicate) is what guarantees this script runs
    the frozen algorithm: a copy would silently fork the moment either is edited.
    Set UC_DETECTOR_DIR to override if the layout changes.
    """
    tried = []
    env = os.environ.get("UC_DETECTOR_DIR")
    if env:
        tried.append(env)
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(5):                       # walk up toward the repo root
        tried += [d, os.path.join(d, "src"), os.path.join(d, "scripts")]
        d = os.path.dirname(d)
    for c in tried:
        if os.path.isfile(os.path.join(c, "uc_detector.py")):
            sys.path.insert(0, c)
            return c
    print("Could not find uc_detector.py. Looked in:")
    for c in tried:
        print(f"  {c}")
    print("\nSet the path explicitly, e.g.")
    print('  $env:UC_DETECTOR_DIR = "C:\\Users\\Jaydenn\\Documents\\fetal-monitor"')
    sys.exit(1)


_DETECTOR_DIR = _locate_detector()
import uc_detector as UD

# ---- frozen constants -------------------------------------------------
FS_UC       = 4.0       # decimated device rate; matches CTU-CHB UC channel
MIN_GAP_S   = 5.0       # a time step larger than this is a capture dropout
OVERLAP_MIN = 0.30      # scoring rule: a detection matches a reference event
                        # if their intersection covers at least this fraction
                        # of the reference. Stated here because it is NOT
                        # necessarily the rule score_uc.py uses on CTU-CHB;
                        # with 3 events the rule matters less than declaring it.

CAP_DIR, PLOT_DIR, RES_DIR = "captures", "plots", "results"


# ----------------------------------------------------------------------- #
def newest_stem():
    hits = sorted(glob.glob(os.path.join(CAP_DIR, "contractions_*_uc4hz.csv")))
    if not hits:
        print(f"No captures matching {CAP_DIR}/contractions_*_uc4hz.csv")
        sys.exit(1)
    return os.path.basename(hits[-1])[:-len("_uc4hz.csv")]


def load_uc(stem):
    path = os.path.join(CAP_DIR, f"{stem}_uc4hz.csv")
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return d[:, 0], d[:, 1]


def find_trim_point(t):
    """Return the time after the largest capture gap.

    03_simulated_contractions.py (pre-flush) began recording with stale serial
    backlog, then lost samples when the OS buffer overflowed. The protocol
    proper starts after that discontinuity. Detect it from the timestamps
    rather than hardcoding an offset, so a clean capture trims to nothing.
    """
    if t.size < 2:
        return t[0] if t.size else 0.0
    dt = np.diff(t)
    i = int(np.argmax(dt))
    if dt[i] > MIN_GAP_S:
        return float(t[i + 1]), float(dt[i])
    return float(t[0]), 0.0


def load_protocol(stem, offset):
    """Reference contraction events: one per RAMP->HOLD->FALL group, shifted
    into device time by `offset`. REST phases are not events."""
    path = os.path.join(CAP_DIR, f"{stem}_protocol.csv")
    phases = []
    with open(path) as f:
        next(f)
        for line in f:
            a, b, name = line.strip().split(",")
            phases.append((float(a), float(b), name))

    events, cur = [], None
    for a, b, name in phases:
        kind = name.split()[0]
        if kind == "RAMP":
            cur = [a, b, name.split()[-1]]
        elif kind in ("HOLD", "FALL") and cur is not None:
            cur[1] = b
            if kind == "FALL":
                events.append((cur[0] + offset, cur[1] + offset, cur[2]))
                cur = None
    return events, phases


def score(det_spans_s, ref_events):
    """Greedy one-to-one match by overlap fraction of the reference event."""
    used, rows = set(), []
    for rs, re_, tag in ref_events:
        best, best_ov = None, 0.0
        for k, (ds, de) in enumerate(det_spans_s):
            if k in used:
                continue
            ov = max(0.0, min(re_, de) - max(rs, ds))
            if ov > best_ov:
                best, best_ov = k, ov
        frac = best_ov / (re_ - rs)
        if best is not None and frac >= OVERLAP_MIN:
            used.add(best)
            ds, de = det_spans_s[best]
            union = max(re_, de) - min(rs, ds)
            rows.append(dict(tag=tag, matched=True, ref=(rs, re_),
                             det=(ds, de), cover=frac,
                             iou=best_ov / union if union > 0 else 0.0))
        else:
            rows.append(dict(tag=tag, matched=False, ref=(rs, re_),
                             det=None, cover=frac, iou=0.0))
    fp = [det_spans_s[k] for k in range(len(det_spans_s)) if k not in used]
    return rows, fp


# ----------------------------------------------------------------------- #
def main():
    for d in (PLOT_DIR, RES_DIR):
        os.makedirs(d, exist_ok=True)
    stem = sys.argv[1] if len(sys.argv) > 1 else newest_stem()
    print(f"capture: {stem}\n")
    print(f"detector: {os.path.join(_DETECTOR_DIR, 'uc_detector.py')}\n")

    t_all, uc_all = load_uc(stem)
    t_trim, gap = find_trim_point(t_all)
    keep = t_all >= t_trim
    t, uc = t_all[keep], uc_all[keep]

    print("=== input ===")
    print(f"  raw file        {len(t_all)} samples, {t_all[0]:.1f}-{t_all[-1]:.1f} s")
    if gap > 0:
        print(f"  capture gap     {gap:.2f} s -> trimming to t >= {t_trim:.2f} s")
    print(f"  analysed        {len(t)} samples, {t[-1]-t[0]:.1f} s "
          f"({(t[-1]-t[0])/60:.1f} min)")
    print(f"  range           {uc.min():.0f} .. {uc.max():.0f} counts")
    print(f"  effective rate  {len(t)/(t[-1]-t[0]):.3f} Hz (expected {FS_UC})")

    # --- the detector, unmodified ------------------------------------------
    spans, diag = UD.detect_contractions(uc, fs=FS_UC)
    conf, frac, rate = UD.confidence(spans, diag["valid"], fs=FS_UC)

    print("\n=== detector state ===")
    print(f"  amplitude ref   {diag['amp']:.1f} counts (p{UD.AMP_PCTL} of detrended)")
    print(f"  threshold       {diag['thr']:.1f} counts "
          f"({UD.THRESH_FRAC:.2f} x amp)")
    print(f"  analysable      {100*frac:.1f} %")
    print(f"  rate            {rate:.2f} per 10 min")
    print(f"  confidence      {conf}")

    # Device data is float-valued (means of 125 raw samples), so exact-equality
    # flat runs essentially never occur. On CTU-CHB the UC channel is integer
    # quantised, which is what makes flat_run_mask fire. Report it rather than
    # let a silently inert dropout mask look like a clean result.
    dead_frac = 1.0 - frac
    print(f"  dropout masked  {100*dead_frac:.2f} %"
          + ("   (mask effectively inert on float data)"
             if dead_frac < 0.001 else ""))

    # Baseline window vs record length: the rolling baseline cannot roll if the
    # record is shorter than its window.
    win_min = UD.BASELINE_WIN_S / 60.0
    rec_min = (t[-1] - t[0]) / 60.0
    if rec_min < win_min:
        print(f"\n  WARNING: BASELINE_WIN_S is {win_min:.1f} min but the record "
              f"is {rec_min:.1f} min.")
        print( "  The rolling baseline degenerates to a global percentile, so "
               "any slow")
        print( "  drift (FSR hysteresis) is NOT tracked out. Not a valid test "
               "of the")
        print( "  baseline mechanism; it tests thresholding only.")

    det_s = [(t[0] + s / FS_UC, t[0] + e / FS_UC) for s, e in spans]
    print(f"\n=== detections ({len(det_s)}) ===")
    for i, (a, b) in enumerate(det_s, 1):
        print(f"  {i}  {a:7.1f} - {b:7.1f} s   dur {b-a:5.1f} s")

    ref, phases = load_protocol(stem, t_trim)
    rows, fp = score(det_s, ref)

    print(f"\n=== scoring vs protocol (overlap >= {OVERLAP_MIN:.2f}) ===")
    for r in rows:
        rs, re_ = r["ref"]
        if r["matched"]:
            ds, de = r["det"]
            print(f"  cycle {r['tag']}  ref {rs:6.1f}-{re_:6.1f}  "
                  f"det {ds:6.1f}-{de:6.1f}  cover {100*r['cover']:5.1f}%  "
                  f"IoU {r['iou']:.2f}  MATCH")
        else:
            print(f"  cycle {r['tag']}  ref {rs:6.1f}-{re_:6.1f}  "
                  f"best cover {100*r['cover']:5.1f}%  MISS")
    for a, b in fp:
        print(f"  unmatched detection {a:6.1f}-{b:6.1f} s  FALSE POSITIVE")

    tp = sum(r["matched"] for r in rows)
    fn, nfp = len(rows) - tp, len(fp)
    se = tp / len(rows) if rows else 0.0
    pp = tp / (tp + nfp) if (tp + nfp) else 0.0
    f1 = 2 * se * pp / (se + pp) if (se + pp) else 0.0
    print(f"\n  TP {tp}  FN {fn}  FP {nfp}   Se {100*se:.1f}%  "
          f"PPV {100*pp:.1f}%  F1 {100*f1:.1f}%")
    print( "  (3 events on scripted stimulus: descriptive only, NOT a "
           "performance claim)")

    with open(os.path.join(RES_DIR, f"{stem}_detect.csv"), "w") as f:
        f.write("cycle,ref_start_s,ref_end_s,det_start_s,det_end_s,"
                "cover,iou,matched\n")
        for r in rows:
            ds, de = r["det"] if r["det"] else ("", "")
            f.write(f"{r['tag']},{r['ref'][0]:.2f},{r['ref'][1]:.2f},"
                    f"{ds if ds=='' else f'{ds:.2f}'},"
                    f"{de if de=='' else f'{de:.2f}'},"
                    f"{r['cover']:.4f},{r['iou']:.4f},{int(r['matched'])}\n")

    # --- plot ---------------------------------------------------------------
    fig, ax = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    for i, (rs, re_, tag) in enumerate(ref):
        ax[0].axvspan(rs, re_, color="tab:orange", alpha=0.25,
                      label="protocol event" if i == 0 else None)
    ax[0].plot(t, uc, lw=1.2, color="tab:blue", label="UC (4 Hz)")
    ax[0].plot(t, diag["baseline"], lw=1.2, color="tab:green", label="baseline")
    if np.isfinite(diag["thr"]):
        ax[0].plot(t, diag["baseline"] + diag["thr"], lw=1.0, ls="--",
                   color="tab:red", label="threshold")
    ymax = uc.max()
    for i, (a, b) in enumerate(det_s):
        ax[0].plot([a, b], [ymax * 1.03] * 2, lw=5, color="tab:purple",
                   solid_capstyle="butt", label="detected" if i == 0 else None)
    ax[0].set_ylabel("UC (counts)")
    ax[0].set_title(f"{stem} — frozen uc_detector on FSR402 device data — {conf}")
    ax[0].legend(fontsize=8, ncol=5, loc="upper left")

    ax[1].plot(t, diag["detr"], lw=1.2, color="tab:blue", label="detrended")
    if np.isfinite(diag["thr"]):
        ax[1].axhline(diag["thr"], ls="--", color="tab:red", lw=1.0,
                      label=f"thr {diag['thr']:.0f}")
        ax[1].axhline(UD.PEAK_FRAC * diag["amp"], ls=":", color="tab:brown",
                      lw=1.0, label=f"peak gate {UD.PEAK_FRAC*diag['amp']:.0f}")
    ax[1].set_ylabel("detrended (counts)")
    ax[1].set_xlabel("time (s, device clock)")
    ax[1].legend(fontsize=8, ncol=3)

    for a_ in ax:
        a_.grid(alpha=0.3)
    fig.tight_layout()
    png = os.path.join(PLOT_DIR, f"{stem}_detect.png")
    fig.savefig(png, dpi=115)
    plt.close(fig)

    print(f"\n  wrote {os.path.join(RES_DIR, f'{stem}_detect.csv')}")
    print(f"  wrote {png}")


if __name__ == "__main__":
    main()
