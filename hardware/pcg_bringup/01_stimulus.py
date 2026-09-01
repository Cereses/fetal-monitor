"""
01_stimulus.py -- build the acoustic loopback stimuli.

Produces three things:

  1. stimulus/<record>_<start>s_<dur>s.wav
       An fpcgdb segment resampled to the playback rate, for playing at the
       microphone.

  2. stimulus/<record>_<start>s_<dur>s_native333.npy
       THE SAME SEGMENT at fpcgdb's native rate. This is the run-A reference.
       The published p21 figures (0.76 confidence, 97% reliable) are for the
       full ~20 minute record, NOT for a 60 s excerpt, so the frozen pipeline
       must be re-run on this exact array to get the number run C is compared
       against. Comparing a segment capture to a whole-record figure would be
       comparing two different things.

  3. stimulus/tone_sweep.wav
       Equal-amplitude bursts at 25/50/75/100/150/200 Hz. Measures the low-end
       response of the WHOLE playback-plus-acquisition chain, so that any loss
       seen in the fpcgdb loopback can be attributed to a stage instead of
       guessed at.

WHY THE RESAMPLING RATIO IS ASSERTED, NOT TRUSTED
fpcgdb is 333 Hz; sound hardware is not. Resampling is unavoidable, and a wrong
ratio changes the PLAYBACK SPEED, which changes the HEART RATE, silently.
Upsampling by 100x and labelling the file 44100 Hz plays it 1.32x fast: a
140 BPM record becomes 185 BPM. Every sample is valid, the file opens, the
waveform looks correct, and the answer is wrong.

44100 / 333 = 4900 / 37 exactly, so resample_poly(x, 4900, 37) is exact.
This script computes the ratio from the ACTUAL header rate via gcd rather than
hardcoding it, and refuses to write unless input and output durations agree to
under a millisecond.

NO FILTERING IS APPLIED TO THE STIMULUS
DC removal and peak normalisation only. Bandpassing the WAV before playback
would pre-clean the signal and contaminate the comparison; the frozen
pipeline's own bandpass must be the only filter in the path.

USAGE
    python 01_stimulus.py
    python 01_stimulus.py --record fetal_PCG_p02_GW_31 --start 0 --dur 120
"""
import os
import sys
import json
import math
import argparse

import numpy as np
from scipy.signal import resample_poly
from scipy.io import wavfile

import matplotlib
matplotlib.use('Agg')          # headless: Python 3.14 windowing bug on this setup
import matplotlib.pyplot as plt

import wfdb

# ------------------------------------------------------------------ config
DATASET = 'fpcgdb'
SEARCH_DEPTH = 5


def resolve_data_dir(dataset, depth=SEARCH_DEPTH):
    """Find data/<dataset>/ in this directory or an ancestor.

    hardware/pcg_bringup/ sits below the repo root where data/ lives. Walking
    up rather than hardcoding '../../data' means the script survives being run
    from a different depth or a different working directory. Outputs stay
    local to this directory so bring-up artefacts never mix with the frozen
    offline pipeline's results/ and plots/.

    Set FPCGDB_DIR to override.
    """
    env = os.environ.get('FPCGDB_DIR')
    if env:
        if os.path.isdir(env):
            return env
        print(f"FPCGDB_DIR is set to '{env}' but that is not a directory")
        sys.exit(1)

    tried = []
    here = os.path.dirname(os.path.abspath(__file__))
    for _ in range(depth + 1):
        cand = os.path.join(here, 'data', dataset)
        tried.append(cand)
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent

    print(f"could not find data/{dataset}/. Looked in:")
    for c in tried:
        print("   ", c)
    print("\nSet the path explicitly, e.g.")
    print(r'  $env:FPCGDB_DIR = "C:\Users\Jaydenn\Documents\fetal-monitor\data\fpcgdb"')
    sys.exit(1)


DATA_DIR = resolve_data_dir(DATASET)

STIM_DIR = 'stimulus'
PLOT_DIR = 'plots'

# Pre-registered primary stimulus. p21 is the BEST record in the frozen fpcgdb
# evaluation (0.76 confidence, 97% reliable). Chosen deliberately: if the
# acquisition chain cannot reproduce the cleanest record in the set, it
# certainly cannot reproduce a median one, so this gives an unambiguous
# pass/fail rather than an equivocal one.
DEFAULT_RECORD = 'fetal_PCG_p21_GW_39'
DEFAULT_START_S = 60.0         # skip the opening; settling and handling noise
DEFAULT_DUR_S = 60.0           # 60 s at 4 s window / 1 s hop = ~57 windows

PLAYBACK_FS = 44100            # standard rate; every sound card handles it
PEAK = 0.9                     # normalisation target, leaves headroom

TONES_HZ = [25, 50, 75, 100, 150, 200]
TONE_S = 3.0                   # per burst
TONE_GAP_S = 1.0               # silence between bursts, for segmentation
TONE_FADE_S = 0.05             # raised-cosine edges; a hard edge is broadband
                               # and would smear energy across the whole band,
                               # destroying the response measurement


# ------------------------------------------------------------------ helpers
def exact_ratio(fs_in, fs_out):
    """Reduce fs_out/fs_in to lowest terms. Returns (up, down)."""
    g = math.gcd(int(fs_in), int(fs_out))
    return int(fs_out) // g, int(fs_in) // g


def to_int16(x):
    """DC-remove, peak-normalise, convert to int16. No filtering."""
    x = np.asarray(x, dtype=np.float64)
    x = x - np.mean(x)
    peak = np.max(np.abs(x))
    if peak < 1e-12:
        raise ValueError("signal is flat; nothing to write")
    x = x / peak * PEAK
    return (x * 32767.0).astype(np.int16)


def fade(x, fs, fade_s=TONE_FADE_S):
    n = int(fade_s * fs)
    if n * 2 >= len(x):
        return x
    w = 0.5 * (1 - np.cos(np.pi * np.arange(n) / n))   # raised cosine
    x = x.copy()
    x[:n] *= w
    x[-n:] *= w[::-1]
    return x


# ------------------------------------------------------------------ fpcgdb
def load_segment(record, start_s, dur_s):
    """Read one channel of an fpcgdb record. Returns (seg, fs, meta)."""
    path = os.path.join(DATA_DIR, record)
    if not os.path.exists(path + '.hea'):
        avail = sorted(f[:-4] for f in os.listdir(DATA_DIR) if f.endswith('.hea'))
        print(f"record '{record}' not found in {DATA_DIR}/")
        print("available:")
        for a in avail:
            print("   ", a)
        sys.exit(1)

    rec = wfdb.rdrecord(path)
    fs = float(rec.fs)
    sig = np.asarray(rec.p_signal[:, 0], dtype=np.float64)

    # fpcgdb headers declare units of "mV" but the values are raw ADC counts.
    # Recorded here for the audit trail; it does not affect anything, since
    # every stage downstream is amplitude-relative.
    units = rec.units[0] if rec.units else '?'

    i0 = int(round(start_s * fs))
    i1 = i0 + int(round(dur_s * fs))
    if i1 > len(sig):
        print(f"requested {start_s}+{dur_s}s but record is only "
              f"{len(sig)/fs:.1f}s ({len(sig)} samples at {fs} Hz)")
        sys.exit(1)

    seg = sig[i0:i1]
    n_nonfinite = int(np.sum(~np.isfinite(seg)))

    meta = {
        'record': record, 'fs_native': fs, 'declared_units': units,
        'record_len_s': len(sig) / fs, 'start_s': start_s, 'dur_s': dur_s,
        'i0': i0, 'i1': i1, 'n_samples': len(seg),
        'n_nonfinite': n_nonfinite,
        'sig_name': rec.sig_name[0] if rec.sig_name else '?',
        'n_channels': rec.p_signal.shape[1],
    }
    return seg, fs, meta


def write_fpcg_wav(seg, fs_native, meta, out_wav):
    up, down = exact_ratio(fs_native, PLAYBACK_FS)
    print(f"  resample ratio : up={up} down={down}  "
          f"({fs_native} x {up}/{down} = {fs_native * up / down:.6f} Hz)")

    if abs(fs_native * up / down - PLAYBACK_FS) > 1e-9:
        print("  [FAIL] ratio does not land exactly on the playback rate.")
        sys.exit(1)

    res = resample_poly(seg, up, down)

    # --- the assertion that matters -------------------------------------
    dur_in = len(seg) / fs_native
    dur_out = len(res) / PLAYBACK_FS
    print(f"  duration in    : {dur_in:.6f} s ({len(seg)} samples)")
    print(f"  duration out   : {dur_out:.6f} s ({len(res)} samples)")
    print(f"  difference     : {abs(dur_out - dur_in)*1000:.4f} ms")
    if abs(dur_out - dur_in) > 1e-3:
        print("  [FAIL] duration changed by more than 1 ms. Playback speed, and")
        print("         therefore reported BPM, would be wrong. Not writing.")
        sys.exit(1)
    print("  [ok] duration preserved -> playback speed is correct")

    wavfile.write(out_wav, PLAYBACK_FS, to_int16(res))
    meta.update({'playback_fs': PLAYBACK_FS, 'up': up, 'down': down,
                 'dur_in_s': dur_in, 'dur_out_s': dur_out})
    return res


# ------------------------------------------------------------------ tones
def write_tone_sweep(out_wav):
    parts, layout, t = [], [], 0.0
    gap = np.zeros(int(TONE_GAP_S * PLAYBACK_FS))
    parts.append(gap)
    t += TONE_GAP_S
    for f0 in TONES_HZ:
        n = int(TONE_S * PLAYBACK_FS)
        tt = np.arange(n) / PLAYBACK_FS
        parts.append(fade(np.sin(2 * np.pi * f0 * tt), PLAYBACK_FS))
        layout.append({'hz': f0, 'start_s': round(t, 4),
                       'end_s': round(t + TONE_S, 4)})
        t += TONE_S
        parts.append(gap)
        t += TONE_GAP_S

    x = np.concatenate(parts)
    # Equal amplitude per burst BEFORE normalisation, so relative response is
    # the only thing the measurement can be reading.
    x = x / (np.max(np.abs(x)) + 1e-12) * PEAK
    wavfile.write(out_wav, PLAYBACK_FS, (x * 32767).astype(np.int16))
    return layout, len(x) / PLAYBACK_FS


# ------------------------------------------------------------------ main
print(f"data dir : {DATA_DIR}\n")
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--record', default=DEFAULT_RECORD)
    ap.add_argument('--start', type=float, default=DEFAULT_START_S)
    ap.add_argument('--dur', type=float, default=DEFAULT_DUR_S)
    args = ap.parse_args()

    os.makedirs(STIM_DIR, exist_ok=True)
    os.makedirs(PLOT_DIR, exist_ok=True)

    stem = f"{args.record}_{int(args.start)}s_{int(args.dur)}s"

    print("=== fpcgdb segment ===")
    seg, fs, meta = load_segment(args.record, args.start, args.dur)
    print(f"  record         : {meta['record']}")
    print(f"  channels       : {meta['n_channels']}  (using '{meta['sig_name']}')")
    print(f"  native fs      : {fs} Hz")
    print(f"  record length  : {meta['record_len_s']:.1f} s")
    print(f"  declared units : {meta['declared_units']}  "
          f"(headers say mV; values are raw counts)")
    print(f"  segment        : {args.start:.1f}-{args.start+args.dur:.1f} s, "
          f"{meta['n_samples']} samples")
    print(f"  non-finite     : {meta['n_nonfinite']}")
    print(f"  range          : {np.nanmin(seg):.1f} to {np.nanmax(seg):.1f}")

    native_npy = os.path.join(STIM_DIR, stem + '_native333.npy')
    np.save(native_npy, seg)
    print(f"  run-A reference -> {native_npy}")

    out_wav = os.path.join(STIM_DIR, stem + '.wav')
    res = write_fpcg_wav(seg, fs, meta, out_wav)
    print(f"  wav            -> {out_wav}")

    print("\n=== tone sweep ===")
    tone_wav = os.path.join(STIM_DIR, 'tone_sweep.wav')
    layout, tone_dur = write_tone_sweep(tone_wav)
    for L in layout:
        print(f"  {L['hz']:>3} Hz  {L['start_s']:>6.2f} - {L['end_s']:>6.2f} s")
    print(f"  total {tone_dur:.2f} s -> {tone_wav}")

    manifest = {'fpcg': meta, 'tones': {'layout': layout,
                                        'tone_s': TONE_S, 'gap_s': TONE_GAP_S,
                                        'fade_s': TONE_FADE_S,
                                        'total_s': tone_dur}}
    mpath = os.path.join(STIM_DIR, 'stimulus_manifest.json')
    with open(mpath, 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(f"\nmanifest -> {mpath}")

    # Verification plot. The printed numbers above are the evidence; this is
    # only for spotting something grossly wrong, like a segment of silence.
    fig, ax = plt.subplots(2, 1, figsize=(12, 6))
    ax[0].plot(np.arange(len(seg)) / fs, seg, linewidth=0.5)
    ax[0].set_title(f"{stem} -- native {fs:g} Hz (run-A reference)")
    ax[0].set_xlabel('Time (s)')
    ax[0].grid(alpha=0.3)
    ax[1].plot(np.arange(len(res)) / PLAYBACK_FS, res, linewidth=0.3)
    ax[1].set_title(f"resampled to {PLAYBACK_FS} Hz (what gets played)")
    ax[1].set_xlabel('Time (s)')
    ax[1].grid(alpha=0.3)
    plt.tight_layout()
    p = os.path.join(PLOT_DIR, f'01_stimulus_{stem}.png')
    plt.savefig(p, dpi=120)
    plt.close()
    print(f"plot -> {p}")


if __name__ == '__main__':
    main()
