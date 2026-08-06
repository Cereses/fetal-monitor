"""
uc_detector.py — detect uterine contractions in the CTU-CHB UC channel.

Every design choice below is forced by something visible in
plots/uc_annotation_overlay.png, not chosen a priori:

  ROLLING BASELINE      Record 1025's resting tone falls from ~15 to ~5 across
                        the recording; 1029 wanders 0-25. A fixed offset cannot
                        track that, so the baseline is a rolling low percentile.

  RECORD-RELATIVE       The UC channel is 'nd', uncalibrated, so absolute units
  THRESHOLD             mean nothing across records. Scale comes from the
                        record's own UPPER TAIL (AMP_PCTL), not its typical
                        variation. Record 1029 is why: its unmarked wobbles
                        reach 20-40, which in record 1025 would be genuine
                        contractions. Anchoring to a median-like statistic
                        would drown 1029 in false positives.

  MEDIAN PREFILTER      Records 1003 and 1025 contain single-sample spikes to
                        100 that nobody annotated. A short median filter kills
                        them before they can set the amplitude scale.

  DROPOUT MASK          Record 1003's first 21 minutes are flat and carry zero
                        annotations: the expert did not score dead signal. Any
                        detection there is a false positive we inflicted on
                        ourselves. Masked regions are also interpolated before
                        baseline estimation so a long zero stretch cannot drag
                        the baseline down for the live signal around it.

  GAP MERGE             Quantisation and noise chop a single contraction into
                        several above-threshold fragments. Merge, then gate on
                        duration, in that order.

NOTHING HERE IS TUNED YET. The constants are first guesses. Run this, look at
plots/uc_detection_overlay.png, and change them once with evidence. Do not tune
against the final score.
"""
import os
import csv

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter, uniform_filter1d, percentile_filter
from scipy.signal import find_peaks
import wfdb

DATASET = 'ctu-uhb-ctgdb'
DATA_DIR = os.path.join('data', DATASET)
RESULTS_DIR = 'results'
PLOT_DIR = 'plots'

FS = 4.0

# --- dropout -----------------------------------------------------------------
DROPOUT_MIN_S = 10.0            # matches audit_uc_quality.py; keep them in step

# --- conditioning ------------------------------------------------------------
MEDIAN_S = 5.0                  # spike killer, shorter than any real contraction
SMOOTH_S = 15.0                 # quantisation/noise smoother

# --- baseline ----------------------------------------------------------------
BASELINE_WIN_S = 600.0          # 10 min: spans several contraction cycles, so the
                                # low percentile lands on inter-contraction tone
BASELINE_PCTL = 10

# --- threshold ---------------------------------------------------------------
AMP_PCTL = 98                   # upper-tail reference for "how big is a real
                                # contraction in THIS record". Raised from 95:
                                # in record 1029 the single real contraction
                                # occupies only ~6% of the trace, so p95 landed
                                # on the unmarked wobbles instead of on it.
THRESH_FRAC = 0.30              # detect above this fraction of that reference
PEAK_FRAC = 0.50                # a span must also PEAK above this fraction,
                                # which is what rejects 1029's small wobbles

# --- trough splitting --------------------------------------------------------
# Record 1025 detected 26 against 37 annotated because between minutes 5 and 9
# the inter-contraction troughs never fall below threshold, so four or five
# contractions arrive as one span. Raising the threshold would separate them at
# the cost of the shallow contractions around minutes 12-17, so split on
# internal structure instead.
MIN_PEAK_SEP_S = 40.0           # Set from the reference distribution, not by
                                # eye: annotated start-to-start intervals have
                                # p0.5 = 39 s, so a 40 s floor discards 0.5% of
                                # real pairs. The previous 60 s discarded 2%
                                # corpus-wide and 17% of record 1025, which is
                                # a fast record (median interval 88 s vs the
                                # corpus 135 s, minimum 38 s).
SPLIT_PROM_FRAC = 0.20          # a second peak must be this prominent, relative
                                # to the record amplitude, to justify a split

# --- span gating -------------------------------------------------------------
MERGE_GAP_S = 20.0
MIN_DUR_S = 20.0                # Same rule as MIN_PEAK_SEP_S: annotated
                                # durations have p1 = 18 s and p2 = 23 s, so a
                                # 20 s floor discards 1.3% of real contractions
                                # against 2.6% at the previous 25 s.
MAX_DUR_S = 300.0               # generous: 1029's single annotation runs ~4 min

# --- per-record confidence ---------------------------------------------------
MIN_ANALYSABLE_FRAC = 0.40
PLAUSIBLE_RATE = (1.0, 8.0)     # contractions per 10 analysable minutes


def find_channel(rec, *targets):
    cleaned = [n.strip() for n in rec.sig_name]
    for t in targets:
        if t in cleaned:
            return cleaned.index(t)
    return 0


def flat_run_mask(x, min_n):
    """True where the sample sits inside a run of >= min_n identical values."""
    mask = np.zeros(x.size, bool)
    if x.size == 0:
        return mask
    change = np.flatnonzero(np.diff(x) != 0) + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [x.size]))
    for s, e in zip(starts, ends):
        if e - s >= min_n:
            mask[s:e] = True
    return mask


def spans_from_mask(mask):
    """Contiguous True regions as (start, end) with end exclusive."""
    if not mask.any():
        return []
    d = np.diff(mask.astype(np.int8))
    starts = list(np.flatnonzero(d == 1) + 1)
    ends = list(np.flatnonzero(d == -1) + 1)
    if mask[0]:
        starts.insert(0, 0)
    if mask[-1]:
        ends.append(mask.size)
    return list(zip(starts, ends))


def merge_spans(spans, max_gap):
    if not spans:
        return []
    out = [list(spans[0])]
    for s, e in spans[1:]:
        if s - out[-1][1] <= max_gap:
            out[-1][1] = e
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def split_at_troughs(detr, s, e, amp, fs):
    """Split one above-threshold span wherever it holds several prominent peaks.

    Contractions that follow each other closely never return to baseline in
    between, so thresholding alone fuses them. Peaks separated by at least
    MIN_PEAK_SEP_S and prominent enough to be real are treated as separate
    contractions, with the boundary placed at the minimum between them.
    """
    seg = detr[s:e]
    if seg.size < 3:
        return [(s, e)]
    peaks, _ = find_peaks(seg,
                          distance=max(1, int(MIN_PEAK_SEP_S * fs)),
                          prominence=SPLIT_PROM_FRAC * amp)
    if peaks.size <= 1:
        return [(s, e)]
    bounds = [s]
    for a, b in zip(peaks[:-1], peaks[1:]):
        bounds.append(s + a + int(np.argmin(seg[a:b])))
    bounds.append(e)
    return [(u, v) for u, v in zip(bounds[:-1], bounds[1:]) if v > u]


def detect_contractions(uc, fs=FS):
    """Detect contractions. Returns (spans, diag) where spans is a list of
    (start_sample, end_sample) and diag holds the intermediate signals for
    plotting and debugging."""
    n = uc.size
    finite = np.isfinite(uc)
    x = np.where(finite, uc, np.nan)

    # --- dropout mask -------------------------------------------------------
    dead = ~finite.copy()
    if finite.any():
        vals = np.where(finite, uc, np.inf)   # sentinel breaks runs at NaN gaps
        dead |= flat_run_mask(vals, int(DROPOUT_MIN_S * fs))
    valid = ~dead

    if valid.sum() < int(60 * fs):            # under a minute of live signal
        return [], dict(valid=valid, smooth=np.full(n, np.nan),
                        baseline=np.full(n, np.nan), thr=np.nan, amp=np.nan)

    # --- interpolate across dead regions before filtering -------------------
    # Filtering across a zeroed stretch would pull the baseline down for the
    # live signal beside it. Interpolation keeps the baseline honest; the mask
    # still suppresses detection inside those regions later.
    idx = np.arange(n)
    filled = np.interp(idx, idx[valid], uc[valid])

    # --- condition ----------------------------------------------------------
    med = median_filter(filled, size=max(3, int(MEDIAN_S * fs)), mode='nearest')
    smooth = uniform_filter1d(med, size=max(3, int(SMOOTH_S * fs)), mode='nearest')

    # --- baseline and detrend ----------------------------------------------
    baseline = percentile_filter(
        smooth, percentile=BASELINE_PCTL,
        size=max(3, int(BASELINE_WIN_S * fs)), mode='nearest')
    detr = np.clip(smooth - baseline, 0, None)

    # --- record-relative scale ---------------------------------------------
    amp = float(np.percentile(detr[valid], AMP_PCTL))
    if amp <= 0:
        return [], dict(valid=valid, smooth=smooth, baseline=baseline,
                        thr=np.nan, amp=amp)
    thr = THRESH_FRAC * amp

    # --- spans --------------------------------------------------------------
    above = (detr > thr) & valid
    spans = merge_spans(spans_from_mask(above), int(MERGE_GAP_S * fs))

    # Split fused contractions BEFORE gating, so that duration limits apply to
    # individual contractions rather than to accidental concatenations.
    split = []
    for s, e in spans:
        split.extend(split_at_troughs(detr, s, e, amp, fs))

    keep = []
    for s, e in split:
        dur = (e - s) / fs
        if dur < MIN_DUR_S or dur > MAX_DUR_S:
            continue
        if detr[s:e].max() < PEAK_FRAC * amp:   # rejects low wobbles
            continue
        keep.append((s, e))

    return keep, dict(valid=valid, smooth=smooth, baseline=baseline,
                      thr=thr, amp=amp, detr=detr)


def confidence(spans, valid, fs=FS):
    frac = float(valid.mean())
    minutes = valid.sum() / fs / 60.0
    rate = len(spans) / minutes * 10 if minutes > 0 else 0.0
    ok = (frac >= MIN_ANALYSABLE_FRAC
          and PLAUSIBLE_RATE[0] <= rate <= PLAUSIBLE_RATE[1])
    return ('HIGH' if ok else 'LOW'), frac, rate


# --------------------------------------------------------------------------- #
def load_contractions():
    path = os.path.join(RESULTS_DIR, 'ctu_contractions.csv')
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            out.setdefault(r['record'], []).append(
                (int(r['start_sample']), int(r['end_sample'])))
    return out


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    ann = load_contractions()
    available = sorted(f[:-4] for f in os.listdir(DATA_DIR) if f.endswith('.hea'))

    # Same four records as the annotation overlay, so the two figures compare
    # directly: most annotated, typical, sparsest, unannotated.
    have = sorted((len(ann.get(r, [])), r) for r in available if ann.get(r))
    none = [r for r in available if not ann.get(r)]
    picks = [have[-1][1], have[len(have) // 2][1], have[0][1]]
    if none:
        picks.append(none[0])

    fig, axes = plt.subplots(len(picks), 1, figsize=(15, 3.4 * len(picks)))
    if len(picks) == 1:
        axes = [axes]

    print(f"{'record':<10}{'ref':<7}{'det':<7}{'conf':<7}{'analysable':<13}"
          f"{'rate/10min':<12}{'amp':<8}{'thr'}")
    for ax, name in zip(axes, picks):
        rec = wfdb.rdrecord(os.path.join(DATA_DIR, name))
        uc = rec.p_signal[:, find_channel(rec, 'UC')].astype(float)
        spans, d = detect_contractions(uc)
        conf, frac, rate = confidence(spans, d['valid'])
        ref = ann.get(name, [])
        print(f"{name:<10}{len(ref):<7}{len(spans):<7}{conf:<7}"
              f"{100*frac:>5.1f}%{'':<7}{rate:>6.2f}{'':<6}"
              f"{d['amp']:>6.1f}  {d['thr']:>6.1f}")

        t = np.arange(uc.size) / FS / 60.0
        ax.axhspan(0, 0, color='none')
        for s, e in spans_from_mask(~d['valid']):
            ax.axvspan(s / FS / 60, e / FS / 60, color='0.88', zorder=0)
        for i, (s, e) in enumerate(ref):
            ax.axvspan(s / FS / 60, e / FS / 60, color='tab:orange', alpha=0.30,
                       zorder=1, label='annotated' if i == 0 else None)
        ax.plot(t, uc, lw=0.6, color='tab:blue', alpha=0.55, zorder=2)
        ax.plot(t, d['baseline'], lw=1.0, color='tab:green', zorder=3,
                label='baseline')
        if np.isfinite(d['thr']):
            ax.plot(t, d['baseline'] + d['thr'], lw=0.9, ls='--',
                    color='tab:red', zorder=3, label='threshold')
        ymax = np.nanmax(uc) if np.isfinite(uc).any() else 1
        for i, (s, e) in enumerate(spans):
            ax.plot([s / FS / 60, e / FS / 60], [ymax * 1.04] * 2, lw=4,
                    color='tab:purple', solid_capstyle='butt', zorder=4,
                    label='detected' if i == 0 else None)
        ax.set_xlim(t[0], t[-1])
        ax.set_ylabel('UC (nd)')
        ax.set_title(f"record {name} — {len(ref)} annotated, "
                     f"{len(spans)} detected — {conf}", fontsize=10)
        ax.legend(loc='upper right', fontsize=8, ncol=4)

    axes[-1].set_xlabel('time (minutes)')
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, 'uc_detection_overlay.png')
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out}")


if __name__ == '__main__':
    main()
