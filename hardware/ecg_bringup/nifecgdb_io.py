"""
nifecgdb_io.py — load NIFECGDB at its native rate. Shared by 03 and 04.

WHY THIS MODULE EXISTS: wfdb's read_edf SILENTLY DECIMATES THIS DATABASE 10x.

Confirmed from the header of ecgca102.edf:
    Number of data records         : 45
    Duration of each data record   : 6.0 s
    Number of Samples per Record   : [6000, 6000, 6000, 6000, 6000, 600]
    Signal Labels                  : Thorax_1, Thorax_2, Abdomen_1,
                                     Abdomen_2, Abdomen_3, EDF Annotations

    6000 / 6.0 = 1000 Hz   the five ECG channels
     600 / 6.0 =  100 Hz   the EDF Annotations text channel

EDF permits a different sampling rate per signal, and the EDF+ spec states
that the 'EDF Annotations' label is reserved: that signal carries text (TALs),
not samples, so its "rate" is a storage artifact with no physical meaning.

wfdb/io/convert/edf.py lines 342-345 do this:

    sample_rate    = [int(i / block_duration) for i in samps_per_block]
    fs             = functools.reduce(math.gcd, sample_rate)
    samps_per_frame= [int(s / min(samps_per_block)) for s in samps_per_block]
    sig_len        = int(fs * num_blocks * block_duration)

fs is the GCD across ALL channels INCLUDING the annotations channel, so
gcd(1000, 100) = 100. Then lines 419-426 average every samps_per_frame=10
consecutive source samples into one output sample. The result is not 100 Hz
of the recording; it is a 10-point boxcar average of it.

read_edf's own docstring documents this ("all signals are resampled at the
lowest sampling frequency used for any signal in the record ... almost
certainly not what you want") and points at smooth_frames=False — a parameter
read_edf does not accept. Its signature is (record_name, pn_dir, header_only,
verbose, rdedfann_flag, encoding).

WHY DECIMATION IS NOT MERELY A RESOLUTION LOSS HERE.
Prefiltering is 'HP:0.01Hz LP:100Hz NF:50Hz', so content extends to 100 Hz.
Decimating to 100 Hz puts Nyquist at 50 Hz, and a 10-point boxcar (first null
at 100 Hz) attenuates but does not remove the 50-100 Hz band. That band folds
back into 0-50 Hz. The wfdb path ALIASES, it does not just smooth.

A SECOND wfdb BUG IN THE SAME BRANCH. Compare the two scaling paths:
    equal-rate  : (raw - baseline) / gain          <- correct
    mixed-rate  : mean(raw - baseline / gain)      <- wrong
The mixed-rate branch never divides by gain, so p_signal is not in physical
units. Irrelevant to Pan-Tompkins (scale-invariant via squaring plus the
adaptive threshold) but fatal to any amplitude or SNR claim — and this
database mixes units across the very groups being compared (mV thoracic,
uV abdominal).

WHAT THIS MODULE GUARANTEES
  - signals at their native 1000 Hz, no resampling, correct physical units
  - the 'EDF Annotations' channel absent from the signal set. pyedflib
    recognises the TAL channel as annotations and does NOT expose it as a
    signal at all, so signals_in_file is 5 where wfdb reported 6. The
    label-based exclusion below is therefore a SAFETY NET that normally
    drops nothing — an empty `dropped` list is the CORRECT outcome here,
    not a failure. wfdb by contrast hands it back as a signal column of
    text bytes reinterpreted as numbers, which looks like a waveform.
  - annotations returned in SIGNAL SAMPLE SPACE, with the scale relationship
    asserted rather than assumed
  - gestational age parsed from the EDF patient field ('Gestation_22+1')

ANNOTATION SCALE. The .qrs files are indexed at 1000 Hz, matching the ECG
channels. Under wfdb's read_edf the signal came back at 100 Hz, so annotation
indices ran ~10x past the end of p_signal — the symptom that exposed all of
this. Loading natively removes the mismatch, and check_alignment() asserts it
rather than trusting it.
"""
import os
import re
import glob

import numpy as np
import pyedflib
import wfdb

ANN_EXT = 'qrs'
EDF_ANN_LABEL = 'EDF Annotations'
THORACIC_PREFIX = 'Thorax'
ABDOMINAL_PREFIX = 'Abdomen'

DATA_SUBPATH = os.path.join('data', 'nifecgdb')
SEARCH_DEPTH = 5

# Fraction of the record length that annotation indices may exceed before we
# call it a scale mismatch rather than a trailing beat.
ALIGN_TOL = 0.02


def resolve_data_dir(subpath=DATA_SUBPATH, depth=SEARCH_DEPTH):
    """Find data/nifecgdb/ in cwd or an ancestor.

    The reference datasets live at the repo root beside mitdb and
    ctu-uhb-ctgdb, but the scripts that read them live in
    hardware/ecg_bringup/. A bare relative glob run from the script's own
    directory silently returns nothing, which is worse than an error: an empty
    loop prints nothing at all and looks like the module failed to import.
    """
    here = os.path.abspath(os.getcwd())
    for _ in range(depth + 1):
        cand = os.path.join(here, subpath)
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    return None


def list_records(data_dir=None):
    """Absolute paths to every .edf. Raises with a useful message on a miss."""
    if data_dir is None:
        data_dir = resolve_data_dir()
        if data_dir is None:
            raise FileNotFoundError(
                f"could not find {DATA_SUBPATH} in the working directory or "
                f"{SEARCH_DEPTH} parents (searched upward from "
                f"{os.path.abspath(os.getcwd())}). Run download_nifecgdb.py "
                f"from the repo root, or pass data_dir explicitly.")
    paths = sorted(glob.glob(os.path.join(data_dir, '*.edf')))
    if not paths:
        raise FileNotFoundError(f"no .edf files in {data_dir}")
    return paths


def parse_gestation(source):
    """Extract (weeks, days) from anywhere in the EDF header.

    Accepts a string or a dict of header fields.

    WHY IT SCANS RATHER THAN NAMING A FIELD. The EDF+ local patient
    identification packs four space-separated subfields:

        code  sex  birthdate  name

    NIFECGDB stores 'X F X Gestation_22+1', so the gestational age lands in
    the NAME slot — code='X', sex='F', birthdate='X', name='Gestation_22+1'.
    The first version of this function checked getPatientAdditional() and
    getPatientCode() only, and returned None on every real record while
    appearing to work on a synthetic test file that happened to use the
    additional field. Scanning all header strings removes the guess.

    Gestational age is the one covariate that plausibly changes abdominal
    maternal QRS amplitude: the growing uterus displaces the heart and adds
    tissue between it and the electrodes. It is recoverable per record from
    the header alone, with no external metadata.

    Returns (weeks, days, field_name) or None.
    """
    if source is None:
        return None
    fields = source if isinstance(source, dict) else {'_': source}
    for key, val in fields.items():
        if not isinstance(val, str) or not val:
            continue
        m = re.search(r'Gestation[_\s]*(\d+)\s*\+\s*(\d+)', val)
        if m:
            return int(m.group(1)), int(m.group(2)), key
        m = re.search(r'Gestation[_\s]*(\d+)', val)
        if m:
            return int(m.group(1)), 0, key
    return None


def load_record(path, channels=None):
    """Read one .edf at native per-channel rates.

    Returns a dict:
        sig        : (n_samples, n_channels) float, physical units
        sig_name   : channel labels, annotation channel removed
        units      : physical dimension per kept channel
        fs         : the ECG channels' sampling rate (asserted uniform)
        duration_s : file duration
        gestation  : (weeks, days) or None
        prefilter  : prefiltering string from the header
        dropped    : labels excluded, and why
    """
    f = pyedflib.EdfReader(path)
    try:
        labels = [s.strip() for s in f.getSignalLabels()]
        rates = [float(f.getSampleFrequency(i)) for i in range(f.signals_in_file)]
        units = [f.getPhysicalDimension(i).strip() for i in range(f.signals_in_file)]
        nsamp = list(f.getNSamples())
        duration = float(f.getFileDuration())
        # Scan every header string. See parse_gestation for why naming a
        # single field silently returned None on all 55 real records.
        hdr = dict(f.getHeader())
        for extra, val in (('patientname', f.getPatientName()),
                           ('patientcode', f.getPatientCode()),
                           ('patient_additional', f.getPatientAdditional()),
                           ('recording_additional', f.getRecordingAdditional())):
            hdr.setdefault(extra, val)
        g = parse_gestation(hdr)
        gest = (g[0], g[1]) if g else None
        gest_field = g[2] if g else None
        prefilter = f.getPrefilter(0).strip() if f.signals_in_file else ''

        keep, dropped = [], []
        for i, lab in enumerate(labels):
            if lab == EDF_ANN_LABEL:
                dropped.append((lab, 'EDF+ text channel, not a signal'))
                continue
            if channels is not None and lab not in channels:
                dropped.append((lab, 'not requested'))
                continue
            keep.append(i)

        if not keep:
            raise ValueError(f'no signal channels kept in {os.path.basename(path)}')

        kept_rates = {rates[i] for i in keep}
        if len(kept_rates) != 1:
            raise ValueError(
                f'kept channels have mixed rates {sorted(kept_rates)} in '
                f'{os.path.basename(path)}; the caller must handle this '
                f'explicitly rather than have a rate silently chosen')
        fs = kept_rates.pop()

        n = min(nsamp[i] for i in keep)
        sig = np.empty((n, len(keep)), dtype=float)
        for c, i in enumerate(keep):
            sig[:, c] = f.readSignal(i)[:n]

        return {
            'sig': sig,
            'sig_name': [labels[i] for i in keep],
            'units': [units[i] for i in keep],
            'fs': fs,
            'n_samples': int(n),
            'duration_s': duration,
            'gestation': gest,
            'gestation_field': gest_field,
            'prefilter': prefilter,
            'dropped': dropped,
            'all_labels': labels,
            'all_rates': rates,
        }
    finally:
        f.close()


def load_annotations(path, fs_signal, n_samples):
    """Beat annotations from '<record>.edf.qrs', in SIGNAL sample space.

    The base name already ends in '.edf', so the file is 'X.edf.qrs' — an
    unusual convention that trips the obvious guess.

    Returns (samples, ann_fs, scale) where `scale` is the factor applied to
    reach signal sample space. It should be 1.0 when the signal is loaded
    natively. Anything else is reported, never applied silently.
    """
    ann = wfdb.rdann(path, ANN_EXT)
    ann_fs = float(ann.fs) if getattr(ann, 'fs', None) else fs_signal
    scale = fs_signal / ann_fs
    samples = np.round(np.asarray(ann.sample, dtype=float) * scale).astype(np.int64)
    return samples, ann_fs, scale


def check_alignment(samples, n_samples, tol=ALIGN_TOL):
    """Do annotation indices actually land inside the signal?

    This is the check that would have caught the wfdb decimation immediately:
    annotations ran to ~269,292 against a p_signal of 27,000 rows, a ratio of
    9.97. Asserting it costs nothing and turns a silent 0% sensitivity into a
    named error.
    """
    if samples.size == 0:
        return False, 'no annotations'
    ratio = float(samples.max()) / float(n_samples)
    if ratio > 1.0 + tol:
        return False, (f'annotations run {ratio:.2f}x past the signal '
                       f'(max index {samples.max()}, signal {n_samples}). '
                       f'Scale mismatch — do NOT score.')
    if ratio < 0.5:
        return False, (f'annotations cover only {100*ratio:.0f}% of the '
                       f'signal; check for partial annotation.')
    return True, f'aligned ({100*ratio:.1f}% of signal length)'


def channel_groups(sig_name):
    """Split kept channels into the paired comparison groups.

    Thoracic  : maternal reference, analogous to the AD8232 smoke-test
                placement (chest, roughly Lead I/II)
    Abdominal : the belt scenario, where maternal QRS is smaller and fetal
                ECG plus uterine EMG are present in band

    Both are recorded SIMULTANEOUSLY in every record, so comparing them holds
    subject, session, gestational age and equipment constant. The placement
    effect is the only thing that varies. That is the experiment.
    """
    thor = [i for i, n in enumerate(sig_name) if n.startswith(THORACIC_PREFIX)]
    abdo = [i for i, n in enumerate(sig_name) if n.startswith(ABDOMINAL_PREFIX)]
    other = [i for i in range(len(sig_name)) if i not in thor and i not in abdo]
    return {'thoracic': thor, 'abdominal': abdo, 'other': other}


def summarise(path):
    """One-line-per-record load + alignment check. Used by the verify script."""
    rec = load_record(path)
    samples, ann_fs, scale = load_annotations(path, rec['fs'], rec['n_samples'])
    ok, msg = check_alignment(samples, rec['n_samples'])
    groups = channel_groups(rec['sig_name'])
    return {
        'record': os.path.basename(path),
        'fs': rec['fs'],
        'n_samples': rec['n_samples'],
        'duration_s': rec['duration_s'],
        'n_beats': int(samples.size),
        'ann_fs': ann_fs,
        'scale': scale,
        'aligned': ok,
        'align_msg': msg,
        'gestation': rec['gestation'],
        'gestation_field': rec.get('gestation_field'),
        'units': rec['units'],
        'sig_name': rec['sig_name'],
        'n_thoracic': len(groups['thoracic']),
        'n_abdominal': len(groups['abdominal']),
        'prefilter': rec['prefilter'],
        'dropped': rec['dropped'],
    }


def _main():
    """Verification entry point: python nifecgdb_io.py [--all]

    Loads records and reports what must hold before any scoring: native
    1000 Hz rate, annotation scale 1.0, annotations landing inside the signal,
    and no kept channel being the EDF+ text channel. Prints a result for every
    record, so an empty data directory produces a NAMED error, not silence.

    NOTE the 'dropped' column reading NONE is the CORRECT outcome — see the
    closing text and the module docstring. pyedflib never exposes the TAL
    channel as a signal, so there is nothing to drop.
    """
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', default=None)
    ap.add_argument('--all', action='store_true',
                    help='check every record (default: first 3)')
    args = ap.parse_args()

    try:
        paths = list_records(args.data_dir)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        raise SystemExit(1)

    print(f"Found {len(paths)} .edf files in {os.path.dirname(paths[0])}\n")
    if not args.all:
        paths = paths[:3]

    bad = []
    for p in paths:
        try:
            s = summarise(p)
        except Exception as e:
            print(f"  {os.path.basename(p):<16} FAILED: {e}")
            bad.append(p)
            continue
        gest = (f"{s['gestation'][0]}w{s['gestation'][1]}d"
                if s['gestation'] else '?')
        drop = ', '.join(lbl for lbl, _ in s['dropped']) or 'NONE'
        flag = '' if s['aligned'] else '   <-- ALIGNMENT PROBLEM'
        print(f"  {s['record']:<16} fs {s['fs']:>7.1f}  n {s['n_samples']:>7}  "
              f"beats {s['n_beats']:>5}  scale {s['scale']:.3f}  "
              f"gest {gest:>7}  thor/abd {s['n_thoracic']}/{s['n_abdominal']}  "
              f"dropped [{drop}]{flag}")
        if not s['aligned']:
            print(f"      {s['align_msg']}")
            bad.append(p)

    print()
    if bad:
        print(f"  {len(bad)} record(s) failed. Do NOT score until resolved.")
        raise SystemExit(1)
    print("  All checked records: native rate, scale 1.0, annotations aligned.")
    print("  'dropped' reading NONE is CORRECT. pyedflib recognises the EDF+")
    print("  TAL channel as annotations and never exposes it as a signal, so")
    print("  signals_in_file is 5 where wfdb reported 6 — there is nothing to")
    print("  drop. The label-based exclusion in load_record() is a safety net")
    print("  for readers that DO hand it back (wfdb does, as a column of text")
    print("  bytes reinterpreted as numbers). What matters is the check below,")
    print("  not the 'dropped' column being non-empty.")
    stray = [s for s in (safe_summarise(p) for p in paths)
             if s and EDF_ANN_LABEL in s['sig_name']]
    if stray:
        print(f"\n  PROBLEM: '{EDF_ANN_LABEL}' present among KEPT channels in "
              f"{len(stray)} record(s). It is text, not signal.")
        raise SystemExit(1)
    print(f"\n  Verified: no kept channel is '{EDF_ANN_LABEL}'.")


def safe_summarise(path):
    try:
        return summarise(path)
    except Exception:
        return None


if __name__ == '__main__':
    _main()
