"""
02_mhr_detector.py — Maternal heart rate from a PPG signal.

Unlike the FHR detector (heart *sounds*, low SNR, autocorrelation), maternal
PPG is a high-SNR direct measurement: large, well-separated systolic peaks.
So the pipeline is:

    bandpass 0.5-5 Hz  ->  auto-orient peaks up  ->  per-window:
    z-normalize -> find_peaks -> inter-beat intervals -> BPM + confidence

Confidence is based on the regularity (coefficient of variation) of the
inter-beat intervals within a window: a clean pulse train is evenly spaced,
noise-driven false peaks are not. Two-tier confidence matches the FHR
detector (per-window >= CONF_THRESH, record reliable if >= RELIABLE_FRAC
of windows clear it), so the evaluation scaffolding carries over.

All spatial quantities are expressed via fs, so the same code works at
125 Hz (BIDMC) and 25 Hz (aromring).
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import wfdb
from scipy.signal import butter, sosfiltfilt, find_peaks
from scipy.stats import skew

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
DATASET = 'bidmc'
DATA_DIR = os.path.join('data', DATASET)
PLOTS_DIR = os.path.join('plots', DATASET)
os.makedirs(PLOTS_DIR, exist_ok=True)

PPG_CHANNEL = 'PLETH'        # channel name to look for in BIDMC signal records

LOW_HZ, HIGH_HZ = 0.5, 5.0   # cardiac band: 0.5-5 Hz ~= 30-300 BPM
FILTER_ORDER = 2             # low normalized cutoff -> keep order modest + use SOS

MIN_BPM, MAX_BPM = 40, 200   # plausible HR bounds
WINDOW_S, HOP_S = 8.0, 4.0   # ~12 beats per 8 s window at 90 BPM; 50% overlap
MIN_IBI = 3                  # need >=3 intervals (>=4 peaks) to judge a window
CV_REF = 0.25                # IBI coeff-of-variation that maps to zero confidence
PROMINENCE = 0.35            # peak prominence in std units (after z-normalization)

DEDOUBLE_MIN_PEAKS = 6       # need enough peaks to judge doubling reliably
DEDOUBLE_MIN_BPM = 150       # only suspect doubling above this apparent rate
DEDOUBLE_RATIO = 0.60        # dropped/kept prominence below this = dicrotic doubling

CONF_THRESH = 0.45           # per-window high-confidence cutoff (matches FHR)
RELIABLE_FRAC = 0.50         # >=50% high-conf windows -> record reliable (matches FHR)

NORMAL_BPM = (60, 100)       # typical adult resting band (shaded in plot only)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def find_channel(rec, target):
    """Index of `target` channel, tolerating trailing commas/whitespace
    in BIDMC channel names (e.g. 'PLETH,' instead of 'PLETH')."""
    cleaned = [name.strip().strip(',').strip() for name in rec.sig_name]
    return cleaned.index(target) if target in cleaned else None


def bandpass_filter(x, fs, low_hz=LOW_HZ, high_hz=HIGH_HZ, order=FILTER_ORDER):
    """Butterworth bandpass in SOS form (stable at low normalized cutoffs).
    Cleans non-finite samples first (sosfiltfilt would propagate them) and
    clamps the high cutoff below Nyquist so the same call works at any fs."""
    x = np.asarray(x, dtype=np.float64)

    bad = ~np.isfinite(x)
    if bad.any():
        n_bad = int(bad.sum())
        print(f"  [bandpass] cleaning {n_bad} non-finite sample(s) "
              f"({100 * n_bad / x.size:.3f}%) before filtering")
        good_idx = np.where(~bad)[0]
        if len(good_idx) < 2:
            raise ValueError("Signal has too few finite samples to interpolate")
        x = x.copy()
        x[bad] = np.interp(np.where(bad)[0], good_idx, x[good_idx])

    nyq = 0.5 * fs
    high_hz = min(high_hz, 0.95 * nyq)
    low_hz = min(low_hz, high_hz - 0.1)
    print(f"  [bandpass] fs={fs} Hz, using {low_hz:.2f}-{high_hz:.2f} Hz")
    sos = butter(order, [low_hz / nyq, high_hz / nyq], btype='band', output='sos')
    return sosfiltfilt(sos, x)


def orient_peaks_up(x):
    """PPG systolic peaks should point up. A signal with sharp peaks is
    skewed toward them; if skew is negative the pulse is inverted, so flip."""
    return -x if skew(x) < 0 else x


def _window_bounds(n, fs):
    """(start, end) sample index pairs. If the record is shorter than one
    window (e.g. the 4 s aromring clip), return a single full-length window."""
    win = int(round(WINDOW_S * fs))
    hop = int(round(HOP_S * fs))
    if n <= win:
        return [(0, n)]
    bounds = []
    start = 0
    while start + win <= n:
        bounds.append((start, start + win))
        start += hop
    # tail window if a meaningful remainder is left over
    if bounds and bounds[-1][1] < n and (n - bounds[-1][1]) >= 0.5 * win:
        bounds.append((n - win, n))
    return bounds


def _dedouble_peaks(wn, peaks, prominences, fs):
    """Correct harmonic doubling robustly. When the apparent rate is
    implausibly high, re-detect peaks with the spacing widened to ~1.5x the
    current median interval (forces skipping the dicrotic hump between true
    beats). Commit only if the dropped peaks were markedly smaller than the
    kept ones -- that confirms dicrotic doubling rather than a genuine fast
    rhythm (where all peaks are similar height).

    `wn` is the z-normalized window; re-detecting on it (rather than parity-
    subsetting the existing peaks) is robust to missed humps and uneven hump
    heights, which broke the earlier even/odd approach."""
    if len(peaks) < DEDOUBLE_MIN_PEAKS:
        return peaks

    med_ibi = np.median(np.diff(peaks))            # samples
    bpm = 60.0 / (med_ibi / fs)
    if bpm < DEDOUBLE_MIN_BPM:                     # plausible rate -> trust it
        return peaks

    # Candidate true-beat set: re-detect with spacing near the halved interval.
    new_dist = max(1, int(round(1.5 * med_ibi)))
    kept, _ = find_peaks(wn, distance=new_dist, prominence=PROMINENCE)
    if len(kept) < MIN_IBI + 1:
        return peaks

    # kept is a subset of peaks (same prominence filter, only stricter spacing).
    prom_by_peak = dict(zip(peaks.tolist(), np.asarray(prominences, float).tolist()))
    kept_set = set(kept.tolist())
    dropped = [p for p in peaks.tolist() if p not in kept_set]
    if not dropped:
        return peaks

    kept_prom = np.median([prom_by_peak[p] for p in kept.tolist() if p in prom_by_peak])
    drop_prom = np.median([prom_by_peak[p] for p in dropped])
    if kept_prom <= 0:
        return peaks

    # Dropped peaks clearly smaller -> dicrotic doubling -> accept de-doubled set.
    if drop_prom / kept_prom <= DEDOUBLE_RATIO:
        return kept
    return peaks


# ----------------------------------------------------------------------
# Core detector
# ----------------------------------------------------------------------
def detect_mhr(raw, fs, auto_orient=True, invert=False):
    """Estimate per-window maternal heart rate from a raw PPG signal."""
    raw = np.asarray(raw, dtype=np.float64)
    filt = bandpass_filter(raw, fs)
    if auto_orient:
        filt = orient_peaks_up(filt)
    if invert:
        filt = -filt

    min_dist = max(1, int(round((60.0 / MAX_BPM) * fs)))

    bpms, confs, times, all_peaks = [], [], [], []
    for a, b in _window_bounds(len(filt), fs):
        w = filt[a:b]
        center_t = 0.5 * (a + b) / fs
        sd = w.std()
        if sd < 1e-12:
            bpms.append(np.nan); confs.append(0.0); times.append(center_t)
            continue

        wn = (w - w.mean()) / sd
        peaks, props = find_peaks(wn, distance=min_dist, prominence=PROMINENCE)

        if len(peaks) < MIN_IBI + 1:
            all_peaks.extend((peaks + a).tolist())
            bpms.append(np.nan); confs.append(0.0); times.append(center_t)
            continue

        # Correct harmonic doubling before measuring the rate.
        peaks = _dedouble_peaks(wn, peaks, props['prominences'], fs)
        all_peaks.extend((peaks + a).tolist())

        if len(peaks) < MIN_IBI + 1:          # de-doubling may leave too few
            bpms.append(np.nan); confs.append(0.0); times.append(center_t)
            continue

        ibi = np.diff(peaks) / fs
        bpm = 60.0 / np.median(ibi)
        cv = np.std(ibi) / np.mean(ibi)
        conf = float(np.clip(1.0 - cv / CV_REF, 0.0, 1.0))
        if not (MIN_BPM <= bpm <= MAX_BPM):
            conf = 0.0

        bpms.append(bpm); confs.append(conf); times.append(center_t)

    bpms = np.asarray(bpms)
    confs = np.asarray(confs)
    times = np.asarray(times)
    high_frac = float((confs >= CONF_THRESH).mean()) if confs.size else 0.0

    return {
        'fs': fs,
        'filtered': filt,
        'bpms': bpms,
        'confidence': confs,
        'window_times': times,
        'peaks': np.array(sorted(set(all_peaks)), dtype=int),
        'high_conf_fraction': high_frac,
        'reliable': high_frac >= RELIABLE_FRAC,
    }


# ----------------------------------------------------------------------
# Diagnostic plot
# ----------------------------------------------------------------------
def plot_diagnostic(raw, result, fs, name, out, zoom_s=20):
    filt = result['filtered']
    peaks = result['peaks']
    times = result['window_times']
    bpms = result['bpms']
    conf = result['confidence']
    high = conf >= CONF_THRESH

    nz = min(int(zoom_s * fs), len(raw))
    tz = np.arange(nz) / fs

    fig, ax = plt.subplots(4, 1, figsize=(12, 11))

    ax[0].plot(tz, raw[:nz], lw=0.8)
    ax[0].set_title(f'{name} - raw PPG (first {zoom_s} s)')
    ax[0].set_xlabel('Time (s)'); ax[0].set_ylabel('PPG')

    ax[1].plot(tz, filt[:nz], lw=0.8)
    pk = peaks[peaks < nz]
    ax[1].plot(pk / fs, filt[pk], 'r.', ms=7, label='detected peaks')
    ax[1].set_title(f'Bandpass {LOW_HZ}-{HIGH_HZ} Hz + detected peaks (first {zoom_s} s)')
    ax[1].set_xlabel('Time (s)'); ax[1].set_ylabel('Amplitude')
    ax[1].legend(loc='upper right')

    ax[2].axhspan(*NORMAL_BPM, color='green', alpha=0.10, label='Typical resting')
    ax[2].plot(times[high], bpms[high], 'o-', color='tab:blue', ms=4,
               label=f'High conf (>={CONF_THRESH})')
    ax[2].plot(times[~high], bpms[~high], 'x', color='tab:orange', ms=6,
               label=f'Low conf (<{CONF_THRESH})')
    ax[2].set_ylim(40, 180)
    ax[2].set_title(f'Estimated MHR over time  '
                    f'(reliable: {100 * result["high_conf_fraction"]:.0f}% of windows)')
    ax[2].set_xlabel('Time (s)'); ax[2].set_ylabel('BPM')
    ax[2].legend(loc='upper right', fontsize=8)

    ax[3].plot(times, conf, '.-', color='tab:purple', lw=1.0)
    ax[3].axhline(CONF_THRESH, ls='--', color='gray', label=f'threshold {CONF_THRESH}')
    ax[3].set_ylim(0, 1)
    ax[3].set_title('Per-window confidence (inter-beat-interval regularity)')
    ax[3].set_xlabel('Time (s)'); ax[3].set_ylabel('Confidence')
    ax[3].legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    plt.savefig(out, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------
# Main: run on one BIDMC record
# ----------------------------------------------------------------------
def main():
    headers = sorted(f for f in os.listdir(DATA_DIR) if f.endswith('.hea'))
    sig_headers = [h for h in headers if not h[:-4].endswith('n')]   # skip numerics
    if not sig_headers:
        print(f"No signal records in {DATA_DIR}/. Run the BIDMC download first.")
        return

    # Pick a record by index - change [0] to inspect a different one.
    record_name = os.path.splitext(sig_headers[23])[0]
    rec = wfdb.rdrecord(os.path.join(DATA_DIR, record_name))

    idx = find_channel(rec, PPG_CHANNEL)
    if idx is None:
        print(f"No '{PPG_CHANNEL}' channel in {record_name}; channels: {rec.sig_name}")
        return

    raw = rec.p_signal[:, idx]
    fs = rec.fs
    result = detect_mhr(raw, fs)

    valid = result['bpms'][np.isfinite(result['bpms'])]
    print(f"Record: {record_name}  (fs={fs} Hz, {len(raw)/fs:.1f} s)")
    print(f"Windows analyzed     : {len(result['bpms'])}")
    if valid.size:
        print(f"Median estimated MHR : {np.median(valid):.1f} BPM")
        print(f"Std of MHR estimates : {np.std(valid):.1f} BPM")
    print(f"Mean confidence      : {np.mean(result['confidence']):.2f}")
    print(f"High-conf windows    : {100*result['high_conf_fraction']:.0f}%  "
          f"-> {'OK' if result['reliable'] else 'LOW'}")

    out = os.path.join(PLOTS_DIR, f'02_mhr_{record_name}.png')
    plot_diagnostic(raw, result, fs, record_name, out)
    print(f"\nDiagnostic plot saved to: {out}")


if __name__ == '__main__':
    main()
