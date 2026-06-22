"""01_explore.py
Load one fetal PCG record, print metadata, plot the raw waveform.
This is just to confirm the data loads and to see what you're working with.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import wfdb

DATASET = 'fpcgdb'
DATA_DIR = os.path.join('data', DATASET)
PLOTS_DIR = os.path.join('plots', DATASET)
os.makedirs(PLOTS_DIR, exist_ok=True)


def list_records(data_dir):
    """Return sorted list of record base names (no extension) in data_dir."""
    bases = set()
    for f in os.listdir(data_dir):
        if f.endswith('.hea'):
            bases.add(os.path.splitext(f)[0])
    return sorted(bases)


def load_record(data_dir, record_name):
    """Load a WFDB record and any annotation files we can find."""
    path = os.path.join(data_dir, record_name)
    record = wfdb.rdrecord(path)

    annotations = {}
    for ext in ('qrs', 'atr', 'fqrs'):
        try:
            ann = wfdb.rdann(path, ext)
            annotations[ext] = ann
        except (FileNotFoundError, Exception):
            pass

    return record, annotations


def plot_record(record, annotations, save_path, seconds=10):
    """Plot the first `seconds` of signal, mark annotations if present."""
    fs = record.fs
    n = min(int(seconds * fs), len(record.p_signal))
    t = np.arange(n) / fs
    sig = record.p_signal[:n, 0]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(t, sig, linewidth=0.7)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude')
    ax.set_title(f'{record.record_name}  |  fs = {fs} Hz  |  '
                 f'duration = {len(record.p_signal)/fs:.1f} s')
    ax.grid(True, alpha=0.3)

    # Overlay annotations if any fall in our window
    for ext, ann in annotations.items():
        beat_times = ann.sample / fs
        beat_times = beat_times[beat_times < seconds]
        if len(beat_times):
            ax.vlines(beat_times, sig.min(), sig.max(),
                      colors='r', alpha=0.4, label=f'.{ext} ({len(ann.sample)} total)')
    if annotations:
        ax.legend()

    plt.tight_layout()
    plt.savefig(save_path, dpi=120)
    plt.close()


def main():
    records = list_records(DATA_DIR)
    print(f"Found {len(records)} records in {DATA_DIR}")
    print("First 5:", records[:5])
    if not records:
        print("\nNo records found. Run `python download_data.py` first.")
        return

    # Inspect a record by index — change [0] to look at a different one
    name = records[16]
    record, annotations = load_record(DATA_DIR, name)

    print(f"\nRecord: {name}")
    print(f"  Sampling rate : {record.fs} Hz")
    print(f"  Duration      : {len(record.p_signal)/record.fs:.2f} s")
    print(f"  Channels      : {record.sig_name}")
    print(f"  Units         : {record.units}")
    print(f"  Annotations   : {list(annotations.keys()) or 'none'}")

    if 'comments' in dir(record) and record.comments:
        print(f"  Header notes  : {record.comments}")

    out = os.path.join(PLOTS_DIR, f'01_explore_{name}.png')
    plot_record(record, annotations, out)
    print(f"\nPlot saved to: {out}")


if __name__ == '__main__':
    main()
