"""
02_detect_device.py — run the FROZEN Pan-Tompkins detector on AD8232 device data.

THE QUESTION THIS ANSWERS
E1 scored Se 99.52% / PPV 99.35% on MIT-BIH: 360 Hz, clinical amplifiers,
cardiologist annotations, 46 records. None of that is this prototype. This
script asks the only question that has been outstanding since the smoke test:
does the same algorithm, unmodified, find beats in a signal from OUR hardware?

Everything that could differ, differs:
    MIT-BIH                     this capture
    360 Hz                      250 Hz
    clinical front end          AD8232 breakout, 3-lead
    12-bit, well scaled         12-bit, R waves SATURATE (1.11% of samples)
    .atr expert annotations     31 hand-reviewed plateau centres
    30 min x 46 records         30 s x 1 record
A pass here does not make the device clinical. It establishes that the frozen
pipeline transfers to the hardware, which is the specific claim the writeup
needs and currently cannot make.

NOTHING IS COPIED. Both the detector and the scorer are imported BY PATH:
    detect_qrs                     <- 02_qrs_detect.py
    match_beats, TOLERANCE_S       <- 03_qrs_evaluate.py
Copying either would mean this script could silently score a different
algorithm than the one that produced the MIT-BIH headline. The toco phase hit
exactly this (04_detect.py could not see uc_detector.py) and the fix was the
same: resolve the path, import the original, never duplicate. If the import
fails, that is a REAL failure and the script stops rather than falling back to
a local reimplementation.

INTEGRITY GATE
The manifest carries SHA-256 of both the capture and the frozen annotations.
Both are re-checked here before anything runs. This is what makes the audit
trail mean something: it proves the detector scored the SAME bytes that were
annotated, and that the reference was not edited after the freeze commit. A
mismatch stops the run.

TIMING ERROR, AND THE CORRECTION
match_beats() works in integer samples and returns signed errors as
    error = detection - reference_integer
The frozen reference carries two columns (see 01_annotate.py):
    sample        floor(plateau centre)  -> what match_beats() receives
    sample_exact  true float centre      -> what the timing error should be
                                            measured against
Ten of 31 beats sit on even-width plateaus, so their true centre is half a
sample (2 ms) later than the integer. Correcting for that:
    corrected_error = raw_error - (sample_exact - sample)
When every reference beat matches, the errors array is in reference order and
the correction is applied EXACTLY, per beat. If any reference is missed the
pairing is no longer positional, and rather than reimplementing the matching
(which would defeat the purpose of importing it) the script falls back to the
scalar mean offset and SAYS SO in the output. An approximate correction that
announces itself beats an exact one that required forking the scorer.

PRE-REGISTERED PREDICTION (frozen in the manifest before this script existed)
    signed timing error small, within about +/-1 sample (4 ms)
Recorded honestly as the WEAKER of two predictions considered. The stronger
one — that clipped plateaus bias refined peaks early, because np.argmax returns
the first index of a maximum — was withdrawn before the freeze: refine_to_r_peak
operates on the BANDPASSED signal, and a 3-sample plateau (12 ms) is far
shorter than the 5-15 Hz passband period (67-200 ms), so the flat top is
smoothed into a single rounded peak. The prediction is checked below either
way, and a miss is reported as a miss.

WHAT A FAILURE WOULD MEAN
Se/PPV below the MIT-BIH figures is not automatically a defect. Possible
causes, in the order worth checking:
    saturation      R waves are clipped, so the peak the detector refines to
                    is truncated. Affects morphology far more than slope.
    50 Hz mains     3.9% of signal power, but the 5-15 Hz bandpass puts it deep
                    in the stopband, so this should NOT be the cause.
    fs = 250 Hz     every detector constant is a duration converted via fs
                    (audited before this run), so rate alone should not matter.
    n = 31 beats    ONE missed beat costs 3.2 percentage points of sensitivity.
                    MIT-BIH's 46 records carry ~110,000 beats. These numbers
                    are not comparably precise and must not be quoted as if
                    they were.
Diagnose structurally before touching any constant. Nothing in 02_qrs_detect.py
or 03_qrs_evaluate.py is to be modified in response to this run: they are
frozen, and tuning them here would invalidate the MIT-BIH result they produced.

USAGE
    python 02_detect_device.py
    python 02_detect_device.py --capture raw_ecg_capture_20260729_180936.csv
    python 02_detect_device.py --pipeline-dir ../../ecg
"""
import os
import csv
import json
import argparse
import importlib.util

import numpy as np
import matplotlib
matplotlib.use('Agg')          # headless: Python 3.14 windowing bug on this setup
import matplotlib.pyplot as plt

import importlib.machinery
import hashlib

CAPTURE_DIR = 'captures'
PLOT_DIR = 'plots'
RESULTS_DIR = 'results'

DEFAULT_CAPTURE = 'raw_ecg_capture_20260729_180936.csv'

DETECTOR_FILE = '02_qrs_detect.py'
EVALUATOR_FILE = '03_qrs_evaluate.py'

PLOT_START_S = 6.0             # window drawn in the stage figure
PLOT_DUR_S = 6.0

# MIT-BIH headline, for context only. NOT a target and NOT a pass mark.
MITDB_SE = 99.52
MITDB_PPV = 99.35

PREDICTED_ABS_SAMPLES = 1.0    # pre-registered: signed error within ~1 sample


# ---------------------------------------------------------------- imports
def resolve_pipeline_file(name, explicit_dir=None):
    """Find a frozen pipeline script without copying it.

    Walks upward from this file, checking each ancestor and a few likely
    subdirectory names. Same lesson as the toco phase: the detector must be
    IMPORTED from its original location, so a change there is reflected here
    rather than silently diverging.
    """
    if explicit_dir:
        p = os.path.join(explicit_dir, name)
        if os.path.exists(p):
            return os.path.abspath(p)
        raise FileNotFoundError(f"{name} not found in --pipeline-dir {explicit_dir}")

    here = os.path.dirname(os.path.abspath(__file__)) or os.getcwd()
    subdirs = ('', 'ecg', 'pipeline', 'pipelines', 'src', 'scripts',
               os.path.join('pipelines', 'ecg'), os.path.join('src', 'ecg'))
    d = here
    for _ in range(6):
        for sub in subdirs:
            p = os.path.join(d, sub, name)
            if os.path.exists(p):
                return os.path.abspath(p)
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    raise FileNotFoundError(
        f"could not locate {name} by walking up from {here}.\n"
        f"  Pass --pipeline-dir pointing at the directory holding it.\n"
        f"  Do NOT copy the file here: this script must score the SAME code\n"
        f"  that produced the MIT-BIH result.")


def load_module(path, name):
    """Import a script whose filename starts with a digit.

    `import 02_qrs_detect` is a syntax error, so load by path. Both files guard
    main() behind __name__, so importing executes nothing.
    """
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except ImportError as e:
        raise ImportError(
            f"failed importing {os.path.basename(path)}: {e}\n"
            f"  The pipeline scripts import wfdb at module level. Activate the\n"
            f"  project venv, or: pip install wfdb") from e
    return mod


# ---------------------------------------------------------------- io
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def load_capture(path):
    ms, adc = [], []
    with open(path, newline='') as fh:
        rdr = csv.reader(fh)
        next(rdr)
        for row in rdr:
            if row:
                ms.append(int(row[0]))
                adc.append(int(row[1]))
    return np.array(ms, dtype=np.int64), np.array(adc, dtype=float)


def load_reference(path):
    sample, exact = [], []
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            sample.append(int(row['sample']))
            exact.append(float(row.get('sample_exact') or row['sample']))
    return (np.array(sample, dtype=np.int64), np.array(exact, dtype=float))


def integrity_gate(cap_path, ref_path, man_path):
    """Refuse to score unless the bytes match what was frozen."""
    print("=== Integrity gate ===")
    with open(man_path) as fh:
        man = json.load(fh)

    ok = True
    for label, path, key in (('capture   ', cap_path, 'capture_sha256'),
                             ('annotation', ref_path, 'annotation_sha256')):
        actual = sha256_file(path)
        expect = man.get(key, '')
        if actual == expect:
            print(f"  [ok]   {label} sha256 matches manifest ({actual[:16]}...)")
        else:
            ok = False
            print(f"  [FAIL] {label} sha256 MISMATCH")
            print(f"         manifest: {expect[:32]}")
            print(f"         actual  : {actual[:32]}")

    if not ok:
        print("\n  The data or the reference changed after the freeze.")
        print("  Re-freeze and re-commit, or restore the committed files.")
        print("  Refusing to score.")
        raise SystemExit(1)

    print(f"  [ok]   frozen reference: {man['n_reference_beats']} beats, "
          f"{man['n_manual_entries']} manual")
    pred = man.get('pre_registration', {}).get('prediction', '(none recorded)')
    print(f"  pre-registered prediction:\n    {pred}\n")
    return man


# ---------------------------------------------------------------- scoring
def score(ev, ref_int, ref_exact, det, fs):
    """Score with the IMPORTED matcher. Returns a result dict."""
    tol = int(round(ev.TOLERANCE_S * fs))
    tp, fp, fn, err, claimed = ev.match_beats(ref_int, det, tol)

    se = tp / (tp + fn) if (tp + fn) else 0.0
    ppv = tp / (tp + fp) if (tp + fp) else 0.0
    der = (fp + fn) / ref_int.size if ref_int.size else 0.0

    offsets = ref_exact - ref_int
    if tp == ref_int.size and err.size == ref_int.size:
        # every reference matched, so err is positionally aligned with ref
        corrected = err - offsets
        mode = 'exact (per beat)'
    else:
        corrected = err - float(offsets.mean())
        mode = f'scalar (mean offset {offsets.mean():.4f} samples) — APPROXIMATE'

    return {
        'tol_samples': tol, 'tp': tp, 'fp': fp, 'fn': fn,
        'se': se, 'ppv': ppv, 'der': der,
        'err_raw': err, 'err_corrected': corrected,
        'correction_mode': mode, 'claimed': claimed,
    }


def report(r, ref_int, det, fs):
    n = ref_int.size
    print("=== Detection vs frozen reference ===")
    print(f"  reference beats  : {n}")
    print(f"  detected         : {det.size}")
    print(f"  tolerance        : +/-{r['tol_samples']} samples "
          f"({1000*r['tol_samples']/fs:.0f} ms)")
    print(f"  TP / FP / FN     : {r['tp']} / {r['fp']} / {r['fn']}")
    print(f"\n  sensitivity      : {100*r['se']:6.2f}%   "
          f"(MIT-BIH headline {MITDB_SE}%)")
    print(f"  PPV              : {100*r['ppv']:6.2f}%   "
          f"(MIT-BIH headline {MITDB_PPV}%)")
    print(f"  DER              : {100*r['der']:6.2f}%")
    print(f"\n  One beat is worth {100.0/n:.2f} percentage points of sensitivity")
    print(f"  here, against ~110,000 beats on MIT-BIH. These figures are NOT")
    print(f"  comparably precise and must not be quoted side by side as if")
    print(f"  they were.")

    e = r['err_corrected']
    if e.size:
        ms_ = e / fs * 1000.0
        print(f"\n=== Timing ===")
        print(f"  correction       : {r['correction_mode']}")
        print(f"  signed mean      : {e.mean():+.3f} samples ({ms_.mean():+.2f} ms)")
        print(f"  signed median    : {np.median(e):+.3f} samples "
              f"({np.median(ms_):+.2f} ms)")
        print(f"  MAE              : {np.abs(e).mean():.3f} samples "
              f"({np.abs(ms_).mean():.2f} ms)")
        print(f"  worst            : {e[np.argmax(np.abs(e))]:+.1f} samples "
              f"({ms_[np.argmax(np.abs(ms_))]:+.1f} ms)")
        raw_ms = r['err_raw'] / fs * 1000.0
        print(f"  (uncorrected signed mean: {raw_ms.mean():+.2f} ms — the "
              f"difference is the plateau-centre offset)")

        print(f"\n  Pre-registered prediction: signed error within about "
              f"+/-{PREDICTED_ABS_SAMPLES:.0f} sample "
              f"({1000*PREDICTED_ABS_SAMPLES/fs:.0f} ms)")
        if abs(e.mean()) <= PREDICTED_ABS_SAMPLES:
            print(f"  HELD: |{e.mean():+.3f}| <= {PREDICTED_ABS_SAMPLES:.0f} samples")
        else:
            print(f"  MISSED: |{e.mean():+.3f}| > {PREDICTED_ABS_SAMPLES:.0f} samples.")
            print("  Report the miss. A prediction that only counts when it")
            print("  succeeds is not evidence of anything.")

    if r['fn'] or r['fp']:
        print(f"\n=== Where the errors are ===")
        det_arr = np.asarray(det)
        fps = det_arr[~r['claimed']]
        if fps.size:
            print(f"  false positives at: "
                  f"{', '.join(f'{s/fs:.3f}s' for s in fps[:10])}"
                  f"{' ...' if fps.size > 10 else ''}")
        if r['fn']:
            tol = r['tol_samples']
            missed = [s for s in ref_int
                      if not np.any(np.abs(det_arr - s) <= tol)]
            print(f"  missed beats at   : "
                  f"{', '.join(f'{s/fs:.3f}s' for s in missed[:10])}"
                  f"{' ...' if len(missed) > 10 else ''}")
        print("  Go and LOOK at these in the plot before changing anything.")


def hr_report(ref_exact, det, fs):
    print(f"\n=== Heart rate ===")
    for label, arr in (('reference', ref_exact), ('detected ', np.asarray(det, float))):
        if arr.size > 1:
            rr = np.diff(arr) / fs
            print(f"  {label}: median {60/np.median(rr):5.1f} BPM, "
                  f"range {60/rr.max():.1f}-{60/rr.min():.1f} BPM")
    print("  The range is respiratory sinus arrhythmia, not instability:")
    print("  rate rises on inhalation and falls on exhalation.")


# ---------------------------------------------------------------- plots
def plot_stages(stages, det, ref_int, fs, stem):
    os.makedirs(PLOT_DIR, exist_ok=True)
    a = int(PLOT_START_S * fs)
    b = min(stages['raw'].size, a + int(PLOT_DUR_S * fs))
    t = np.arange(a, b) / fs

    panels = [('raw', 'Raw ADC (AD8232 -> ESP32, 250 Hz)'),
              ('band', 'Stage 1: bandpass 5-15 Hz'),
              ('deriv', 'Stage 2: derivative (slope)'),
              ('sq', 'Stage 3: squared'),
              ('integ', 'Stage 4: moving-window integration (150 ms)')]
    fig, axes = plt.subplots(len(panels) + 1, 1, figsize=(14, 14), sharex=True)
    for ax, (key, title) in zip(axes, panels):
        ax.plot(t, stages[key][a:b], lw=0.8)
        ax.set_title(title, loc='left', fontsize=9)
        ax.grid(alpha=0.3)

    ax = axes[-1]
    ax.plot(t, stages['raw'][a:b], lw=0.8, color='0.4')
    d = np.asarray(det)
    dsel = d[(d >= a) & (d < b)]
    rsel = ref_int[(ref_int >= a) & (ref_int < b)]
    ax.plot(dsel / fs, stages['raw'][dsel], 'v', ms=8, color='tab:red',
            label=f'detected ({dsel.size})')
    ax.plot(rsel / fs, stages['raw'][rsel], 'o', ms=11, mfc='none',
            color='tab:green', label=f'frozen reference ({rsel.size})')
    ax.set_title('Stage 5: detections vs frozen reference', loc='left', fontsize=9)
    ax.set_xlabel('seconds')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle(f'Pan-Tompkins on device data — {stem} @ {fs:.0f} Hz')
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f'{stem}_device_stages.png')
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  plot -> {out}")


def plot_overview(adc, det, ref_int, err_corrected, fs, stem):
    fig, axes = plt.subplots(2, 1, figsize=(15, 7),
                             gridspec_kw={'height_ratios': [2, 1]})
    t = np.arange(adc.size) / fs
    ax = axes[0]
    ax.plot(t, adc, lw=0.5, color='0.4')
    d = np.asarray(det)
    ax.plot(d / fs, adc[d], 'v', ms=6, color='tab:red',
            label=f'detected ({d.size})')
    ax.plot(ref_int / fs, adc[ref_int], 'o', ms=9, mfc='none',
            color='tab:green', label=f'reference ({ref_int.size})')
    ax.set_ylabel('raw ADC')
    ax.set_title(f'{stem} — full record')
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    if err_corrected.size:
        ms_ = err_corrected / fs * 1000.0
        ax.hist(ms_, bins=15, color='tab:blue', alpha=0.75, edgecolor='white')
        ax.axvline(0, color='0.3', lw=1)
        ax.axvline(ms_.mean(), color='tab:red', ls='--', lw=1.2,
                   label=f'mean {ms_.mean():+.2f} ms')
        ax.axvline(1000.0 / fs, color='tab:orange', ls=':', lw=1,
                   label=f'+/-1 sample ({1000.0/fs:.0f} ms)')
        ax.axvline(-1000.0 / fs, color='tab:orange', ls=':', lw=1)
        ax.legend(fontsize=8)
    ax.set_xlabel('signed timing error, corrected (ms)   detection - reference')
    ax.set_ylabel('beats')
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f'{stem}_device_scoring.png')
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  plot -> {out}")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--capture', default=DEFAULT_CAPTURE)
    ap.add_argument('--pipeline-dir', default=None,
                    help='directory holding the frozen 02_/03_ pipeline scripts')
    args = ap.parse_args()

    stem = os.path.splitext(os.path.basename(args.capture))[0]
    cap_path = os.path.join(CAPTURE_DIR, args.capture)
    ref_path = os.path.join(RESULTS_DIR, f'{stem}_annotations_frozen.csv')
    man_path = os.path.join(RESULTS_DIR, f'{stem}_annotations_manifest.json')

    for p in (cap_path, ref_path, man_path):
        if not os.path.exists(p):
            print(f"Missing {p}")
            if 'annotations' in p:
                print("Run 01_annotate.py (and --freeze) first.")
            raise SystemExit(1)

    print(f"=== 02_detect_device.py — {stem} ===\n")

    det_path = resolve_pipeline_file(DETECTOR_FILE, args.pipeline_dir)
    ev_path = resolve_pipeline_file(EVALUATOR_FILE, args.pipeline_dir)
    print("=== Frozen code (imported, not copied) ===")
    print(f"  detector : {det_path}")
    print(f"  scorer   : {ev_path}")
    print(f"  sha256   : {sha256_file(det_path)[:16]}... / "
          f"{sha256_file(ev_path)[:16]}...\n")

    qrs = load_module(det_path, 'qrs_detect')
    ev = load_module(ev_path, 'qrs_evaluate')

    man = integrity_gate(cap_path, ref_path, man_path)

    ms, adc = load_capture(cap_path)
    ref_int, ref_exact = load_reference(ref_path)
    fs = 1000.0 / float(np.median(np.diff(ms)))
    print(f"  capture  : {adc.size} samples, {fs:.1f} Hz, "
          f"{(ms[-1]-ms[0])/1000:.1f} s\n")

    det, stages = qrs.detect_qrs(adc, fs)

    r = score(ev, ref_int, ref_exact, det, fs)
    report(r, ref_int, det, fs)
    hr_report(ref_exact, det, fs)

    print()
    plot_stages(stages, det, ref_int, fs, stem)
    plot_overview(adc, det, ref_int, r['err_corrected'], fs, stem)

    out = os.path.join(RESULTS_DIR, f'{stem}_device_eval.json')
    ms_err = r['err_corrected'] / fs * 1000.0
    with open(out, 'w') as fh:
        json.dump({
            'capture': os.path.basename(cap_path),
            'capture_sha256': man['capture_sha256'],
            'annotation_sha256': man['annotation_sha256'],
            'detector_file': det_path,
            'detector_sha256': sha256_file(det_path),
            'evaluator_sha256': sha256_file(ev_path),
            'fs_hz': fs,
            'tolerance_ms': ev.TOLERANCE_S * 1000.0,
            'reference_beats': int(ref_int.size),
            'detected': int(np.asarray(det).size),
            'TP': r['tp'], 'FP': r['fp'], 'FN': r['fn'],
            'sensitivity': r['se'], 'PPV': r['ppv'], 'DER': r['der'],
            'timing_correction_mode': r['correction_mode'],
            'timing_signed_mean_ms': float(ms_err.mean()) if ms_err.size else None,
            'timing_mae_ms': float(np.abs(ms_err).mean()) if ms_err.size else None,
            'prediction_held': (bool(abs(r['err_corrected'].mean())
                                     <= PREDICTED_ABS_SAMPLES)
                                if r['err_corrected'].size else None),
        }, fh, indent=2)
    print(f"  results -> {out}")
    print("\n  Frozen code was IMPORTED. If any number here is disappointing,")
    print("  diagnose it — do not edit 02_qrs_detect.py or 03_qrs_evaluate.py.")
    print("  Those files produced the MIT-BIH result and changing them now")
    print("  would invalidate it.")


if __name__ == '__main__':
    main()
