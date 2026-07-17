"""
02_qrs_detect.py — Pan-Tompkins QRS detection, implemented from scratch.

Five stages, each solving one problem:
  1. Bandpass 5-15 Hz          : isolate the QRS frequency band
  2. Derivative                : a QRS is defined by steepness, not height
  3. Squaring                  : rectify + nonlinearly boost large values
  4. Moving-window integration : one smooth blob per QRS (~150 ms wide)
  5. Adaptive threshold        : learn signal/noise levels, decide which blobs are beats

Deliberate deviation from the 1985 paper: every stage here is ZERO-PHASE /
CENTERED (filtfilt, centered derivative, centered integration). The original
used causal real-time filters and had to hand-correct their group delay. We are
doing offline analysis on a file, so we can afford zero-phase filtering, and
detected peaks land where the beats actually are. That matters because E1 scores
POSITIONS against the .atr annotations.

Run mode: single record + diagnostic plot. Batch scoring lives in 03.
"""
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt, find_peaks

import wfdb

DATASET = 'mitdb'
DATA_DIR = os.path.join('data', DATASET)
PLOT_DIR = 'plots'

RECORD = '207'        # record to inspect
PLOT_START_S = 300.0    # window of the record to draw (seconds)
PLOT_DUR_S = 5.0

# --- algorithm constants (seconds; converted to samples using each record's fs) ---
BP_LOW_HZ = 5.0
BP_HIGH_HZ = 15.0
BP_ORDER = 3
MWI_WIDTH_S = 0.150       # widest a QRS complex realistically gets
REFRACTORY_S = 0.200      # physiological floor: no two beats closer than this
REFINE_WIN_S = 0.075      # search radius when snapping to the true R-peak
TWAVE_WIN_S = 0.360       # peaks closer than this may be T-waves, not beats


# ---------------------------------------------------------------- utilities
def find_channel(rec, *targets):
    """Index of first matching lead by name; falls back to 0.

    Record 114 stores its leads reversed (['V5', 'MLII']), so index 0 is NOT
    reliably MLII. Matching by name is what catches that.
    """
    cleaned = [n.strip() for n in rec.sig_name]
    for t in targets:
        if t in cleaned:
            return cleaned.index(t)
    return 0


def clean_nonfinite(x):
    """Linearly interpolate over NaN/inf before filtering.

    Hard-won lesson: a single non-finite sample makes scipy's filtfilt return
    all-NaN, silently poisoning the entire downstream pipeline.
    """
    x = np.asarray(x, dtype=float)
    bad = ~np.isfinite(x)
    if bad.any():
        idx = np.arange(x.size)
        x = x.copy()
        x[bad] = np.interp(idx[bad], idx[~bad], x[~bad])
    return x


# ---------------------------------------------------------------- stages
def stage1_bandpass(x, fs):
    """5-15 Hz Butterworth bandpass, zero-phase.

    QRS energy concentrates around 10-15 Hz. This drops baseline wander
    (<1 Hz, from breathing), the slow P/T waves, and muscle/mains noise
    (>20 Hz). SOS form for numerical stability; Nyquist clamp for safety.
    """
    nyq = 0.5 * fs
    high = min(BP_HIGH_HZ, 0.95 * nyq)
    sos = butter(BP_ORDER, [BP_LOW_HZ / nyq, high / nyq], btype='band', output='sos')
    return sosfiltfilt(sos, x)


def stage2_derivative(x, fs):
    """5-point centered derivative: measures slope.

    The point of the whole algorithm in one line. A QRS is not necessarily the
    TALLEST thing in the trace, but it is reliably the STEEPEST. P and T waves
    are slow humps, so differentiating flattens them while the QRS survives.
    Kernel is centered, so no delay is introduced.
    """
    kernel = np.array([-1, -2, 0, 2, 1], dtype=float) * (fs / 8.0)
    return np.convolve(x, kernel[::-1], mode='same')


def stage3_square(x):
    """Point-by-point squaring.

    Two jobs: (a) rectify, since the R deflection can be negative-going in some
    leads, and (b) squaring is nonlinear, so it amplifies large values far more
    than small ones, widening the gap between QRS and residual noise.
    """
    return x ** 2


def stage4_integrate(x, fs):
    """Centered moving-window mean, ~150 ms wide.

    The squared signal is a burst of spikes per QRS. Integrating merges each
    burst into ONE smooth blob, so we can count beats instead of counting
    spikes. Window is centered ('same'), so the blob sits over the beat rather
    than lagging it.
    """
    n = max(1, int(round(MWI_WIDTH_S * fs)))
    return np.convolve(x, np.ones(n) / n, mode='same')


def stage5_adaptive_threshold(integ, band, fs):
    """Learn running signal/noise estimates and decide which blobs are beats.

    Why adaptive: over 30 minutes amplitude drifts with posture, electrode
    contact and breathing. Any fixed threshold either misses quiet stretches or
    fires on noisy ones. Pan-Tompkins keeps two running estimates:

        SPKI : running peak level of things believed to be BEATS
        NPKI : running peak level of things believed to be NOISE
        THR  = NPKI + 0.25 * (SPKI - NPKI)      i.e. sit 25% up the gap

    Each accepted/rejected peak nudges its estimate by 12.5% toward the new
    value, so the threshold tracks the signal slowly and does not lurch.

    Also implemented:
      - refractory period (200 ms), a physiological floor on beat spacing
      - T-wave discrimination: a peak arriving <360 ms after a beat is only a
        beat if it is STEEPER than half the previous beat's steepness.
        A tall-but-lazy T-wave gets rejected on slope.
      - searchback: if no beat for 1.66x the recent average RR, re-scan that
        gap at half threshold to recover a missed beat.

    Returns integrated-signal peak indices judged to be QRS complexes.
    """
    refractory = int(round(REFRACTORY_S * fs))
    twave_win = int(round(TWAVE_WIN_S * fs))

    # Candidate blobs. `distance` enforces the refractory period up front.
    cand, _ = find_peaks(integ, distance=refractory)
    if cand.size == 0:
        return np.array([], dtype=int)

    # Initialise estimates from the first 2 s, the paper's learning phase.
    init = integ[:min(integ.size, int(2 * fs))]
    spki = float(np.max(init)) * 0.25 if init.size else 0.0
    npki = float(np.mean(init)) * 0.5 if init.size else 0.0

    # Slope proxy for T-wave discrimination, measured on the bandpassed signal.
    slope = np.abs(np.diff(band, prepend=band[0]))

    def local_slope(i):
        a = max(0, i - int(0.05 * fs))
        b = min(slope.size, i + int(0.05 * fs))
        return float(np.max(slope[a:b])) if b > a else 0.0

    qrs = []
    rr = []           # recent accepted RR intervals (samples)
    last_slope = 0.0

    for i in cand:
        peak = float(integ[i])
        thr_i1 = npki + 0.25 * (spki - npki)
        thr_i2 = 0.5 * thr_i1

        is_qrs = False
        if peak > thr_i1:
            is_qrs = True
            # T-wave check: too soon after the last beat?
            if qrs and (i - qrs[-1]) < twave_win:
                s = local_slope(i)
                if s < 0.5 * last_slope:
                    is_qrs = False          # lazy slope -> it is a T-wave
                    npki = 0.125 * peak + 0.875 * npki

        if is_qrs:
            # --- searchback: did we skip a beat before this one? ---
            if len(rr) >= 2:
                rr_avg = float(np.mean(rr[-8:]))
                if qrs and (i - qrs[-1]) > 1.66 * rr_avg:
                    lo, hi = qrs[-1] + refractory, i - refractory
                    if hi > lo:
                        gap = cand[(cand > lo) & (cand < hi)]
                        gap = gap[integ[gap] > thr_i2]
                        if gap.size:
                            m = gap[np.argmax(integ[gap])]
                            qrs.append(int(m))
                            rr.append(m - qrs[-2]) if len(qrs) >= 2 else None
                            # recovered beats update SPKI more cautiously
                            spki = 0.25 * float(integ[m]) + 0.75 * spki

            if qrs:
                rr.append(i - qrs[-1])
            qrs.append(int(i))
            last_slope = local_slope(i)
            spki = 0.125 * peak + 0.875 * spki
        else:
            npki = 0.125 * peak + 0.875 * npki

    return np.array(sorted(set(qrs)), dtype=int)


def refine_to_r_peak(band, idxs, fs):
    """Snap each accepted blob to the true R-peak in the bandpassed signal.

    The integration blob is smooth, so its maximum is near the beat but not
    exactly on it. We search a small window around it and take the largest
    absolute deflection. Skipping this costs tens of ms of position error on
    every single beat, which directly hurts E1 scoring.
    """
    r = int(round(REFINE_WIN_S * fs))
    out = []
    for i in idxs:
        a, b = max(0, i - r), min(band.size, i + r + 1)
        if b > a:
            out.append(a + int(np.argmax(np.abs(band[a:b]))))
    return np.array(sorted(set(out)), dtype=int)


# ---------------------------------------------------------------- pipeline
def detect_qrs(sig, fs):
    """Full Pan-Tompkins. Returns (r_peaks, stage_dict) for plotting."""
    x = clean_nonfinite(sig)
    band = stage1_bandpass(x, fs)
    deriv = stage2_derivative(band, fs)
    sq = stage3_square(deriv)
    integ = stage4_integrate(sq, fs)
    blobs = stage5_adaptive_threshold(integ, band, fs)
    peaks = refine_to_r_peak(band, blobs, fs)
    return peaks, {'raw': x, 'band': band, 'deriv': deriv, 'sq': sq, 'integ': integ}


# ---------------------------------------------------------------- diagnostics
def plot_stages(stages, peaks, fs, name, ann_samples=None):
    os.makedirs(PLOT_DIR, exist_ok=True)
    a = int(PLOT_START_S * fs)
    b = min(stages['raw'].size, a + int(PLOT_DUR_S * fs))
    t = np.arange(a, b) / fs

    panels = [
        ('raw', 'Raw MLII (mV)'),
        ('band', 'Stage 1: bandpass 5-15 Hz'),
        ('deriv', 'Stage 2: derivative (slope)'),
        ('sq', 'Stage 3: squared'),
        ('integ', 'Stage 4: moving-window integration (150 ms)'),
    ]
    fig, axes = plt.subplots(len(panels) + 1, 1, figsize=(13, 13), sharex=True)

    for ax, (key, title) in zip(axes, panels):
        ax.plot(t, stages[key][a:b], lw=0.8)
        ax.set_title(title, loc='left', fontsize=9)
        ax.grid(alpha=0.3)

    # Final panel: raw trace with detections vs reference annotations.
    ax = axes[-1]
    ax.plot(t, stages['raw'][a:b], lw=0.8, color='0.4')
    sel = peaks[(peaks >= a) & (peaks < b)]
    ax.plot(sel / fs, stages['raw'][sel], 'v', ms=7, color='tab:red',
            label=f'detected ({sel.size})')
    if ann_samples is not None:
        asel = ann_samples[(ann_samples >= a) & (ann_samples < b)]
        ax.plot(asel / fs, stages['raw'][asel], 'o', ms=9, mfc='none',
                color='tab:green', label=f'reference .atr ({asel.size})')
    ax.set_title('Stage 5: detections vs reference', loc='left', fontsize=9)
    ax.set_xlabel('time (s)')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle(f'Pan-Tompkins stages — {DATASET} record {name} @ {fs} Hz')
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f'{DATASET}_{name}_pantompkins_stages.png')
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  plot -> {out}")


# ---------------------------------------------------------------- main
# Non-beat annotation symbols. These mark rhythm changes, signal quality and
# episode boundaries, NOT heartbeats. Concrete proof they must go: record 100's
# first annotation sits at sample 18, only 59 samples (0.16 s) before the next.
# That is ~375 BPM, physically impossible. It is a '+' rhythm marker. Left in,
# it would count as a phantom reference beat and corrupt the score.
NON_BEAT = set('+~|"[]!x')


def main():
    path = os.path.join(DATA_DIR, RECORD)
    rec = wfdb.rdrecord(path)
    fs = rec.fs
    ch = find_channel(rec, 'MLII', 'II')
    sig = rec.p_signal[:, ch]

    print(f"=== {DATASET} record {RECORD} ===")
    print(f"  fs       : {fs} Hz")
    print(f"  channels : {rec.sig_name}  -> using index {ch} ({rec.sig_name[ch].strip()})")
    print(f"  duration : {sig.size / fs / 60:.1f} min")

    peaks, stages = detect_qrs(sig, fs)

    ann = wfdb.rdann(path, 'atr')
    keep = np.array([s not in NON_BEAT for s in ann.symbol])
    ref = ann.sample[keep]
    print(f"\n  annotations total : {len(ann.sample)}")
    print(f"  non-beat filtered : {int((~keep).sum())}")
    print(f"  reference beats   : {ref.size}")
    print(f"  detected  beats   : {peaks.size}")
    print(f"  difference        : {peaks.size - ref.size:+d}")

    if peaks.size > 1:
        bpm = 60.0 * fs / np.diff(peaks)
        print(f"\n  median HR (detected)  : {np.median(bpm):.1f} BPM")
    if ref.size > 1:
        rbpm = 60.0 * fs / np.diff(ref)
        print(f"  median HR (reference) : {np.median(rbpm):.1f} BPM")

    plot_stages(stages, peaks, fs, RECORD, ann_samples=ref)
    print("\n  Counts matching is necessary but NOT sufficient — a detector can")
    print("  hit the right count with wrong positions. Sensitivity/PPV in 03.")


if __name__ == '__main__':
    main()
