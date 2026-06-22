"""02_fhr_detector.py
Full FHR detection pipeline:
  raw PCG -> bandpass filter -> Shannon-energy envelope -> autocorrelation -> BPM
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
import wfdb

DATASET = 'fpcgdb'
DATA_DIR = os.path.join('data', DATASET)
PLOTS_DIR = os.path.join('plots', DATASET)
os.makedirs(PLOTS_DIR, exist_ok=True)

CONFIDENCE_THRESHOLD = 0.45 #Below this, treat output as unreliable

# ----- Signal processing primitives -------------------------------------------

def bandpass_filter(x, fs, low_hz=25.0, high_hz=200.0, order=4):
    """Butterworth bandpass. Default 25-200 Hz isolates the fetal heart-sound band.
    High cutoff is clamped below Nyquist so the filter is valid across sampling rates."""
    # Clean NaN/Inf values via linear interpolation — filtfilt propagates them globally
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
    low_hz = min(low_hz, high_hz - 1.0)   # safety: keep low < high
    print(f"  [bandpass] fs={fs} Hz, using {low_hz:.1f}-{high_hz:.1f} Hz")
    b, a = butter(order, [low_hz / nyq, high_hz / nyq], btype='band')
    return filtfilt(b, a, x)


def shannon_energy_envelope(x, fs, smooth_ms=50):
    """Smoothed normalized Shannon energy. Highlights heart-sound bursts (S1/S2)
    and suppresses low-amplitude noise."""
    x = x / (np.max(np.abs(x)) + 1e-12)
    energy = -(x ** 2) * np.log(x ** 2 + 1e-12)
    # Smooth with a moving-average window
    win = max(1, int(smooth_ms / 1000.0 * fs))
    kernel = np.ones(win) / win
    return np.convolve(energy, kernel, mode='same')


def estimate_fhr_autocorr(envelope, fs, win_sec=4.0, hop_sec=1.0,
                          min_bpm=110, max_bpm=180):
    """Sliding-window autocorrelation. For each window, find the autocorrelation
    peak in the lag range corresponding to plausible fetal BPMs."""
    win = int(win_sec * fs)
    hop = int(hop_sec * fs)
    min_lag = int(60.0 / max_bpm * fs)
    max_lag = int(60.0 / min_bpm * fs)

    times, bpms, confidences = [], [], []
    for start in range(0, len(envelope) - win, hop):
        seg = envelope[start:start + win]
        seg = seg - np.mean(seg)
        if np.std(seg) < 1e-9:
            continue
        ac = np.correlate(seg, seg, mode='full')
        ac = ac[len(ac) // 2:]
        ac = ac / (ac[0] + 1e-12)  # normalize so ac[0] = 1
        if max_lag >= len(ac):
            continue
        peak_lag = int(np.argmax(ac[min_lag:max_lag])) + min_lag
        bpm = 60.0 * fs / peak_lag
        times.append((start + win / 2) / fs)
        bpms.append(bpm)
        confidences.append(float(ac[peak_lag]))  # 0..1, higher = more periodic

    return np.array(times), np.array(bpms), np.array(confidences)


# ----- End-to-end pipeline ----------------------------------------------------

def detect_fhr(raw, fs):
    filtered = bandpass_filter(raw, fs)
    envelope = shannon_energy_envelope(filtered, fs)
    times, bpms, conf = estimate_fhr_autocorr(envelope, fs)
    return {
        'filtered': filtered,
        'envelope': envelope,
        'times': times,
        'bpms': bpms,
        'confidence': conf,
    }


# ----- Diagnostic plot --------------------------------------------------------

def plot_diagnostic(raw, result, fs, record_name, save_path,
                    conf_threshold=CONFIDENCE_THRESHOLD):
    t = np.arange(len(raw)) / fs
    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=False)

    axes[0].plot(t, raw, linewidth=0.6)
    axes[0].set_title(f'{record_name} — raw PCG')
    axes[0].set_ylabel('Amplitude')
    axes[0].grid(alpha=0.3)

    axes[1].plot(t, result['filtered'], linewidth=0.6)
    axes[1].set_title('Bandpass 25–200 Hz')
    axes[1].set_ylabel('Amplitude')
    axes[1].grid(alpha=0.3)

    axes[2].plot(t, result['envelope'], linewidth=0.7)
    axes[2].set_title('Shannon-energy envelope')
    axes[2].set_ylabel('Energy')
    axes[2].set_xlabel('Time (s)')
    axes[2].grid(alpha=0.3)

    # ----- BPM panel with confidence-based styling --------------------------
    ax_bpm = axes[3]
    times = np.asarray(result['times'])
    bpms = np.asarray(result['bpms'])
    conf = np.asarray(result['confidence'])

    high = conf >= conf_threshold
    low = ~high

    # High-confidence: solid blue line + dots
    if high.any():
        ax_bpm.plot(times[high], bpms[high], 'o-', color='steelblue',
                    markersize=4, linewidth=1.2,
                    label=f'High confidence (≥{conf_threshold:.2f})')
    # Low-confidence: faded orange X markers, no connecting line
    if low.any():
        ax_bpm.plot(times[low], bpms[low], 'x', color='darkorange',
                    markersize=5, alpha=0.6,
                    label=f'Low confidence (<{conf_threshold:.2f})')

    ax_bpm.axhspan(110, 160, alpha=0.1, color='green', label='Normal range')
    ax_bpm.set_ylim(80, 200)
    ax_bpm.set_xlim(0, t[-1])
    ax_bpm.set_xlabel('Time (s)')
    ax_bpm.set_ylabel('BPM')

    # Title shows what fraction of the recording was trustworthy
    high_frac = float(np.mean(high)) if len(conf) else 0.0
    ax_bpm.set_title(
        f'Estimated FHR over time  '
        f'(reliable: {high_frac:.0%} of windows)'
    )
    ax_bpm.grid(alpha=0.3)
    ax_bpm.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    plt.savefig(save_path, dpi=120)
    plt.close()

# ----- Main -------------------------------------------------------------------

def main():
    headers = sorted(f for f in os.listdir(DATA_DIR) if f.endswith('.hea'))
    if not headers:
        print(f"No records in {DATA_DIR}/. Run download_data.py first.")
        return

    # Pick a record by index — change [0] to inspect a different one
    record_name = os.path.splitext(headers[20])[0]

    record = wfdb.rdrecord(os.path.join(DATA_DIR, record_name))
    raw = record.p_signal[:, 0]
    fs = record.fs

    result = detect_fhr(raw, fs)

    print(f"Record: {record_name}  (fs={fs} Hz, {len(raw)/fs:.1f} s)")
    print(f"Windows analyzed     : {len(result['bpms'])}")
    if len(result['bpms']):
        print(f"Mean estimated FHR   : {np.mean(result['bpms']):.1f} BPM")
        print(f"Std of FHR estimates : {np.std(result['bpms']):.1f} BPM")
        print(f"Mean confidence      : {np.mean(result['confidence']):.2f}")

    out = os.path.join(PLOTS_DIR, f'02_fhr_{record_name}.png')
    plot_diagnostic(raw, result, fs, record_name, out)
    print(f"\nDiagnostic plot saved to: {out}")


if __name__ == '__main__':
    main()
