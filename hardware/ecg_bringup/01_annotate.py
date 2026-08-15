"""
01_annotate.py — build the frozen R-peak reference for AD8232 device data.

WHY THIS SCRIPT EXISTS
Every ECG number so far (E1: Se 99.52%, PPV 99.35%) was scored against MIT-BIH
.atr files: expert-produced, external, free. Device data has no annotations.
Without a reference, running detect_qrs() on a capture yields "31 beats, median
61 BPM" — a PLAUSIBILITY CHECK, not a validation. It cannot separate a detector
that found every beat from one that missed three and invented three, because
both produce the same count. This script builds the arbiter.

TWO-PHASE BY DESIGN, AND THE ORDER IS THE POINT
    pass 1 (default)  : propose candidates + render review plots
    -- human reviews, edits the draft CSV --
    pass 2 (--freeze) : validate, hash, write the frozen reference

The reference must be fixed BEFORE detect_qrs() ever runs on this capture. The
freeze step records a SHA-256 of both the capture and the annotation set, so
"annotated blind to detector output" is VERIFIABLE rather than asserted. Same
discipline as pre-registering the flutter-audit exclusions and the mitdb 102/104
lead-mismatch decision: the convention is fixed before the number is seen.

WHY CANDIDATES ARE MACHINE-PROPOSED, AND WHY THAT IS NOT CIRCULAR
Clicking 31 peaks by hand is fine; clicking 31 peaks by hand on a headless
matplotlib (Agg, forced by the Python 3.14 windowing bug on this setup) is not
possible. So candidates are proposed automatically and REVIEWED by a human.

That is exactly how MIT-BIH itself was built: a simple slope detector proposed
beats, cardiologists corrected them. The human remains the arbiter; the machine
only saves the clicking.

The proposal method is deliberately chosen to share NOTHING with the detector
under test:
    proposal  : contiguous runs of ADC == 4095 (the saturation plateau),
                centre of each run
    detector  : 5-15 Hz bandpass -> derivative -> square -> 150 ms integration
                -> adaptive SPKI/NPKI threshold -> refine to local max
No filter, no threshold, no parameter is shared. The proposal exploits a
property of THIS hardware (the AD8232 output saturates the ESP32 ADC on every
R wave) that MIT-BIH does not have. If Pan-Tompkins is broken, the proposal is
unaffected, which is the whole requirement for an independent reference.

Verified on capture 20260729_180936 before writing this script:
    31 plateaus / 30.0 s, median RR 0.988 s (60.7 BPM), range 55.4-79.8 BPM
    0 gaps above 1.5x median RR, 0 spacings below the 200 ms refractory floor
Every beat saturates, so the proposal has full coverage on this capture. That
is a PROPERTY OF THIS RECORDING, not a general assumption — the review plots
exist precisely so a missed beat is visible, and --freeze re-runs the gap and
refractory checks on the edited set.

PRE-REGISTRATION (fixed before any detector run)
    capture    : raw_ecg_capture_20260729_180936.csv, full 30 s, no exclusions
    reference  : plateau centres, human-reviewed, frozen here
    scoring    : +/-150 ms tolerance, greedy one-to-one matching, both IMPORTED
                 from 03_qrs_evaluate.py, never reimplemented
    prediction : signed timing error (detection - reference) should be SMALL,
                 within about +/-1 sample.
                 NOTE this is a CORRECTION to an earlier prediction of a
                 systematic early bias. That argument was: plateaus are flat,
                 np.argmax returns the FIRST index of a maximum, so refined
                 peaks land early. The mechanism is weak, because
                 refine_to_r_peak() operates on the BANDPASSED signal, and a
                 3-sample plateau (12 ms at 250 Hz) is far shorter than the
                 5-15 Hz passband period (67-200 ms). The bandpass smooths the
                 flat top into a single rounded peak near the plateau centre.
                 Recording the weaker prediction rather than quietly keeping
                 the stronger one: a prediction is only evidence if it was
                 written down before the measurement.

VALIDATION GATE
Runs before anything else, and refuses to proceed on failure. Three checks,
each from a defect actually found in this project rather than imagined:
    monotonic ms  : capture 180817 carries row 6874 as "7510,1776" where
                    27510 was meant — a single dropped leading byte at
                    230400 baud. The ADC value was intact; only the timestamp
                    lost a character. Silent, and plausible-looking.
    no -1 rows    : the firmware writes val = lo ? -1 : analogRead(). The -1
                    lead-off sentinel is FINITE, so clean_nonfinite() (which
                    only interpolates NaN/inf) passes it straight through as a
                    ~1800-count negative excursion — comparable to an R wave,
                    opposite sign, two huge slope events, guaranteed false
                    positives with no warning raised anywhere.
    uniform dt    : the sketch schedules on micros() but timestamps on
                    millis(), so sub-millisecond jitter is invisible by
                    construction. This check therefore proves only that no
                    interval deviates by >=1 ms. It does NOT establish that
                    jitter is zero, and the notes should not claim it does.

USAGE
    python 01_annotate.py                     # propose + plot
    python 01_annotate.py --freeze            # validate + hash + freeze
    python 01_annotate.py --capture NAME.csv  # a different capture

EDITING THE DRAFT
results/<capture>_annotations_draft.csv, columns: sample, ms, adc, source, keep
    drop a false candidate : set keep to 0 (keep the row — the audit trail
                             should show what was rejected and why)
    add a missed beat      : append a row with ms filled in and sample blank;
                             read the time off the review plot. --freeze snaps
                             it to the local raw maximum within +/-50 ms and
                             REPORTS the snap distance for every manual entry,
                             so the adjustment is auditable rather than hidden.
"""
import os
import sys
import csv
import json
import hashlib
import argparse
import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')          # headless: Python 3.14 windowing bug on this setup
import matplotlib.pyplot as plt

CAPTURE_DIR = 'captures'
PLOT_DIR = 'plots'
RESULTS_DIR = 'results'

DEFAULT_CAPTURE = 'raw_ecg_capture_20260729_180936.csv'

FS_NOMINAL = 250.0             # firmware: PERIOD_US = 1000000 / 250
ADC_MAX = 4095                 # 12-bit, analogReadResolution(12)
LEAD_OFF_SENTINEL = -1

PANEL_S = 6.0                  # seconds per review panel
SNAP_WIN_S = 0.050             # manual-entry snap radius
REFRACTORY_S = 0.200           # physiological floor, matches the detector's
GAP_FACTOR = 1.5               # flag RR intervals above this x median


# ---------------------------------------------------------------- io
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def load_capture(path):
    """Read ms,adc. Returns (ms int array, adc int array)."""
    ms, adc = [], []
    with open(path, newline='') as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        if [c.strip().lower() for c in header] != ['ms', 'adc']:
            raise ValueError(f"unexpected header {header}, expected ['ms','adc']")
        for i, row in enumerate(rdr, start=2):
            if not row:
                continue
            if len(row) != 2:
                raise ValueError(f"line {i}: expected 2 fields, got {len(row)}: {row}")
            ms.append(int(row[0]))
            adc.append(int(row[1]))
    return np.array(ms, dtype=np.int64), np.array(adc, dtype=np.int64)


# ---------------------------------------------------------------- validation
def validate(ms, adc, path):
    """Hard gate. Returns a report dict; raises SystemExit on any failure.

    Loud failure is the entire point. Each of these defects produced data that
    looked completely normal until it was specifically checked for.
    """
    print("=== Validation gate ===")
    failures = []

    # --- monotonic timestamps (capture 180817's dropped byte) ---
    dt = np.diff(ms)
    bad = np.where(dt <= 0)[0]
    if bad.size:
        failures.append(f"{bad.size} non-increasing timestamp(s)")
        print(f"  [FAIL] timestamps not strictly increasing at {bad.size} row(s):")
        for b in bad[:5]:
            print(f"         row {b+2}: ms {ms[b]} -> {ms[b+1]} (delta {dt[b]:+d})")
        print("         a single dropped serial byte turns 27510 into 7510.")
        print("         drop the row, or re-capture at 115200 baud.")
    else:
        print(f"  [ok]   timestamps strictly increasing ({ms.size} rows)")

    # --- lead-off sentinel ---
    n_lo = int((adc == LEAD_OFF_SENTINEL).sum())
    if n_lo:
        failures.append(f"{n_lo} lead-off sentinel row(s)")
        first = np.where(adc == LEAD_OFF_SENTINEL)[0][:5]
        print(f"  [FAIL] {n_lo} row(s) carry the -1 lead-off sentinel")
        print(f"         first at ms {[int(ms[i]) for i in first]}")
        print("         -1 is FINITE: clean_nonfinite() will not catch it, and it")
        print("         reads as a ~1800-count negative spike into the detector.")
        print("         electrodes lost contact — re-capture, do not patch.")
    else:
        print(f"  [ok]   no lead-off sentinels")

    # --- sampling interval ---
    if dt.size:
        uniq = np.unique(dt)
        expected = int(round(1000.0 / FS_NOMINAL))
        if uniq.size == 1 and uniq[0] == expected:
            print(f"  [ok]   sample interval uniform at {expected} ms "
                  f"({FS_NOMINAL:.0f} Hz nominal)")
        else:
            good = dt[dt > 0]
            print(f"  [warn] non-uniform interval: values {uniq[:8].tolist()}"
                  f"{' ...' if uniq.size > 8 else ''}")
            if good.size:
                print(f"         median {np.median(good):.0f} ms, "
                      f"achieved {1000.0/np.median(good):.2f} Hz")
        print("         NOTE: timestamps come from millis() (1 ms resolution)")
        print("         while scheduling uses micros(). Sub-ms jitter is invisible")
        print("         here. This proves no interval deviates by >=1 ms; it does")
        print("         NOT establish that jitter is zero.")

    # --- clipping census (characterisation, not a failure) ---
    n_clip = int((adc >= ADC_MAX).sum())
    print(f"  [info] {n_clip} sample(s) at ADC ceiling {ADC_MAX} "
          f"({100.0*n_clip/adc.size:.2f}%)")
    if n_clip == 0:
        failures.append("no clipped samples: plateau proposal has nothing to find")
        print("         [FAIL] the plateau proposal needs saturation to work.")
        print("         this capture does not saturate — annotate by another means.")

    if failures:
        print(f"\n  VALIDATION FAILED: {'; '.join(failures)}")
        print("  Refusing to proceed. Fix the capture rather than the checks.")
        raise SystemExit(1)

    print("  Validation passed.\n")
    return {'n_rows': int(ms.size), 'n_clipped': n_clip,
            'sha256': sha256_file(path)}


# ---------------------------------------------------------------- proposal
def propose_from_plateaus(adc):
    """Centre of each contiguous run of ADC == ADC_MAX.

    Shares no filter, threshold or parameter with detect_qrs(). Centre rather
    than first/last index because the plateau is a symmetric artifact of
    saturation, so its midpoint is the least biased estimate of where the true
    (unsaturated) peak would have been.

    AMENDMENT (made before any detector run — see module docstring).
    Returns BOTH an exact float centre and a floor-based integer index.

    The first version used int(round(mean)). For an even-width plateau the
    mean lands on x.5, and Python's round() is BANKER'S ROUNDING: it goes to
    the nearest EVEN integer. So round(3000.5)=3000 but round(3001.5)=3002.
    The reference position therefore depended on the PARITY OF THE SAMPLE
    INDEX, which is an accident of when the capture happened to start.

    On capture 180936 that affects 10 of 31 plateaus (all width 2). Five
    rounded up, five down, so the net effect on the mean signed error was
    exactly 0.000 ms — coincidence, not design, and it would not cancel on a
    different capture.

    Irrelevant to Se/PPV: the EC57 tolerance is 150 ms = 37.5 samples at
    250 Hz, and this is half of one sample (2 ms). But we pre-registered a
    prediction about the SIGNED TIMING ERROR at roughly +/-1 sample (4 ms),
    and parity-dependent noise of +/-0.5 sample is the same order as the
    effect being measured. Measuring a 4 ms effect through 2 ms of avoidable
    noise is a bad position to argue from.

    So:
        sample_exact : true float centre (797.5), used for timing-error
                       reporting only
        sample       : floor(centre), used for matching only. Deterministic
                       and stated, rather than inherited from a rounding mode
                       most readers would not expect. Its residual offset from
                       the exact centre is 0.0 or 0.5 samples, known per beat
                       and recorded, so it can be corrected exactly.
    03_qrs_evaluate.py is NOT modified: match_beats() still receives integers.
    """
    idx = np.where(adc >= ADC_MAX)[0]
    if idx.size == 0:
        return (np.array([], dtype=float), np.array([], dtype=int),
                np.array([], dtype=int))
    runs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
    exact = np.array([float(r.mean()) for r in runs], dtype=float)
    centres = np.floor(exact).astype(int)
    widths = np.array([len(r) for r in runs], dtype=int)
    return exact, centres, widths


def sanity_report(samples, fs, label):
    """Physiological plausibility. Flags, never edits — the human decides."""
    if samples.size < 2:
        print(f"  {label}: fewer than 2 beats, nothing to check")
        return
    rr = np.diff(samples) / fs
    med = float(np.median(rr))
    print(f"  {label}: {samples.size} beats")
    print(f"    RR median {med:.3f} s ({60/med:.1f} BPM), "
          f"range {rr.min():.3f}-{rr.max():.3f} s "
          f"({60/rr.max():.1f}-{60/rr.min():.1f} BPM)")
    gaps = np.where(rr > GAP_FACTOR * med)[0]
    if gaps.size:
        print(f"    [check] {gaps.size} interval(s) above {GAP_FACTOR}x median "
              f"— a beat may be MISSING here:")
        for g in gaps:
            print(f"            {samples[g]/fs:7.3f} s -> {samples[g+1]/fs:7.3f} s "
                  f"({rr[g]:.3f} s)")
    else:
        print(f"    [ok]    no interval above {GAP_FACTOR}x median")
    short = np.where(rr < REFRACTORY_S)[0]
    if short.size:
        print(f"    [check] {short.size} interval(s) below the {REFRACTORY_S*1000:.0f} ms "
              f"refractory floor — likely a DOUBLE mark:")
        for s in short:
            print(f"            {samples[s]/fs:7.3f} s -> {samples[s+1]/fs:7.3f} s "
                  f"({rr[s]:.3f} s)")
    else:
        print(f"    [ok]    no interval below the refractory floor")


# ---------------------------------------------------------------- plots
def plot_review(ms, adc, samples, fs, stem):
    """Overview + per-panel detail. These are what the human actually reviews.

    Every candidate is labelled with its index so a specific mark can be named
    in the draft CSV. Panels are short enough that a missing mark is obvious
    rather than lost in a 30-second squeeze.
    """
    os.makedirs(PLOT_DIR, exist_ok=True)
    t = (ms - ms[0]) / 1000.0
    dur = t[-1]

    # --- overview ---
    fig, ax = plt.subplots(figsize=(15, 4))
    ax.plot(t, adc, lw=0.5, color='0.35')
    if samples.size:
        ax.plot(t[samples], adc[samples], 'v', ms=6, color='tab:red',
                label=f'candidates ({samples.size})')
        ax.legend(loc='lower right', fontsize=8)
    ax.axhline(ADC_MAX, ls='--', lw=0.7, color='tab:orange',
               label=f'ADC ceiling {ADC_MAX}')
    ax.set_xlabel('seconds')
    ax.set_ylabel('raw ADC')
    ax.set_title(f'{stem} — full record, {samples.size} candidates')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f'{stem}_annot_overview.png')
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  plot -> {out}")

    # --- panels ---
    n_panels = int(np.ceil(dur / PANEL_S))
    fig, axes = plt.subplots(n_panels, 1, figsize=(15, 2.4 * n_panels))
    axes = np.atleast_1d(axes)
    for k, ax in enumerate(axes):
        a, b = k * PANEL_S, min((k + 1) * PANEL_S, dur)
        m = (t >= a) & (t <= b)
        ax.plot(t[m], adc[m], lw=0.8, color='0.35')
        sel = samples[(t[samples] >= a) & (t[samples] <= b)] if samples.size else np.array([], int)
        for s in sel:
            j = int(np.where(samples == s)[0][0])
            ax.plot(t[s], adc[s], 'v', ms=7, color='tab:red')
            ax.annotate(f'{j}', (t[s], adc[s]), textcoords='offset points',
                        xytext=(0, 9), ha='center', fontsize=7, color='tab:red')
        ax.axhline(ADC_MAX, ls='--', lw=0.6, color='tab:orange')
        ax.set_xlim(a, b)
        # Headroom above the ceiling so the index labels sit clear of the
        # title. Without it the marks at the panel edges are unreadable, and
        # an unreadable review plot defeats the point of the review.
        lo_y = float(np.min(adc[m])) if m.any() else 0.0
        ax.set_ylim(lo_y - 100, ADC_MAX + 450)
        ax.set_ylabel('ADC')
        ax.grid(alpha=0.3)
        ax.set_title(f'{a:.0f}-{b:.0f} s   ({sel.size} candidates)',
                     loc='left', fontsize=9)
    axes[-1].set_xlabel('seconds')
    fig.suptitle(f'{stem} — review panels: is every beat marked exactly once?')
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f'{stem}_annot_panels.png')
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  plot -> {out}")


# ---------------------------------------------------------------- draft io
def write_draft(path, ms, adc, exact, samples, widths):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['sample', 'sample_exact', 'ms', 'adc', 'plateau_width',
                    'source', 'keep'])
        for e, s, pw in zip(exact, samples, widths):
            w.writerow([int(s), f'{e:.1f}', int(ms[s]), int(adc[s]), int(pw),
                        'auto', 1])
    print(f"  draft -> {path}")
    print("    'sample' feeds match_beats(); 'sample_exact' is for timing error only.")
    print("    Leave sample_exact alone when editing — added rows may omit it.")


def read_draft(path, ms, adc, fs):
    """Read the edited draft. Snaps manual entries, reports every snap.

    A manual entry snaps to an ACTUAL sample index, so its exact centre IS that
    integer — there is no plateau to take a midpoint of. Rows omitting
    sample_exact therefore fall back to the integer, which is right for manual
    entries and harmless elsewhere.
    """
    kept, manual_log = [], []
    snap = int(round(SNAP_WIN_S * fs))

    with open(path, newline='') as fh:
        for i, row in enumerate(csv.DictReader(fh), start=2):
            if str(row.get('keep', '1')).strip() in ('0', 'false', 'False', ''):
                continue
            src = (row.get('source') or 'auto').strip()
            raw_sample = (row.get('sample') or '').strip()
            raw_exact = (row.get('sample_exact') or '').strip()

            if raw_sample:
                s = int(raw_sample)
                e = float(raw_exact) if raw_exact else float(s)
                # sample must be floor(sample_exact) or the correction applied
                # downstream is wrong. Catch a hand-edit that breaks the pair.
                if not (0.0 <= e - s < 1.0):
                    raise ValueError(
                        f"draft line {i}: sample_exact {e} not in [{s}, {s+1}) "
                        f"— 'sample' must be floor('sample_exact')")
            else:
                if not (row.get('ms') or '').strip():
                    raise ValueError(f"draft line {i}: needs 'sample' or 'ms'")
                target = int(row['ms'])
                s0 = int(np.argmin(np.abs(ms - target)))
                a, b = max(0, s0 - snap), min(adc.size, s0 + snap + 1)
                s = a + int(np.argmax(adc[a:b]))
                e = float(s)
                manual_log.append({'line': i, 'requested_ms': target,
                                   'snapped_sample': s, 'snapped_ms': int(ms[s]),
                                   'snap_ms': int(ms[s]) - target})
                src = 'manual'
            kept.append((s, e, src))

    kept.sort()
    samples = np.array([k[0] for k in kept], dtype=int)
    exact = np.array([k[1] for k in kept], dtype=float)
    sources = [k[2] for k in kept]

    if samples.size != np.unique(samples).size:
        dup = [int(x) for x in samples if list(samples).count(x) > 1]
        raise ValueError(f"duplicate sample indices in draft: {sorted(set(dup))}")

    if manual_log:
        print(f"\n  Manual entries snapped to local raw maximum "
              f"(+/-{SNAP_WIN_S*1000:.0f} ms):")
        for m in manual_log:
            print(f"    line {m['line']}: requested {m['requested_ms']} ms -> "
                  f"sample {m['snapped_sample']} ({m['snapped_ms']} ms), "
                  f"moved {m['snap_ms']:+d} ms")
    return samples, exact, sources, manual_log


# ---------------------------------------------------------------- passes
def do_propose(cap_path, stem):
    ms, adc = load_capture(cap_path)
    meta = validate(ms, adc, cap_path)

    fs = 1000.0 / float(np.median(np.diff(ms)))
    exact, samples, widths = propose_from_plateaus(adc)

    print("=== Candidate proposal (saturation plateaus) ===")
    print(f"  plateaus       : {samples.size}")
    if widths.size:
        print(f"  plateau width  : {widths.min()}-{widths.max()} samples "
              f"(median {np.median(widths):.1f})")
        n_even = int((widths % 2 == 0).sum())
        print(f"  even-width     : {n_even}/{widths.size} — centre falls between "
              f"two samples")
        print(f"                   ({0.5/fs*1000:.0f} ms), recorded exactly in "
              f"'sample_exact'")
    sanity_report(exact, fs, 'candidates')

    plot_review(ms, adc, samples, fs, stem)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    draft = os.path.join(RESULTS_DIR, f'{stem}_annotations_draft.csv')
    if os.path.exists(draft):
        print(f"\n  {draft} already exists — NOT overwriting.")
        print("  Delete it if you want to regenerate from scratch.")
    else:
        write_draft(draft, ms, adc, exact, samples, widths)

    print("\n  NEXT: open the panel plot and the draft side by side.")
    print("  Confirm every beat carries exactly one mark. Set keep=0 to reject")
    print("  a candidate; append a row with ms filled in to add a missed beat.")
    print("  Then:  python 01_annotate.py --freeze")


def do_freeze(cap_path, stem):
    ms, adc = load_capture(cap_path)
    meta = validate(ms, adc, cap_path)
    fs = 1000.0 / float(np.median(np.diff(ms)))

    draft = os.path.join(RESULTS_DIR, f'{stem}_annotations_draft.csv')
    if not os.path.exists(draft):
        print(f"Missing {draft}. Run without --freeze first.")
        raise SystemExit(1)

    samples, exact, sources, manual_log = read_draft(draft, ms, adc, fs)

    print("\n=== Frozen reference ===")
    sanity_report(exact, fs, 'reference')
    n_manual = sum(1 for s in sources if s == 'manual')
    print(f"  auto-proposed kept : {samples.size - n_manual}")
    print(f"  manually added     : {n_manual}")

    # Offset between the integer index used for matching and the true centre.
    # 0.0 for odd-width plateaus and manual entries, 0.5 for even-width.
    offsets = exact - samples
    mean_off = float(offsets.mean()) if offsets.size else 0.0
    print(f"\n  centre offset (sample_exact - sample):")
    print(f"    beats with a 0.5-sample offset : {int((offsets > 0).sum())}"
          f"/{offsets.size}")
    print(f"    mean offset                    : {mean_off:.4f} samples "
          f"({mean_off/fs*1000:.3f} ms)")
    print("    02 subtracts this from the raw signed error to report the")
    print("    timing bias against the true plateau centres.")

    frozen = os.path.join(RESULTS_DIR, f'{stem}_annotations_frozen.csv')
    with open(frozen, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['sample', 'sample_exact', 'ms', 'adc', 'source'])
        for s, e, src in zip(samples, exact, sources):
            w.writerow([int(s), f'{e:.1f}', int(ms[s]), int(adc[s]), src])

    ann_hash = sha256_file(frozen)
    manifest = {
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'capture_file': os.path.basename(cap_path),
        'capture_sha256': meta['sha256'],
        'annotation_file': os.path.basename(frozen),
        'annotation_sha256': ann_hash,
        'n_reference_beats': int(samples.size),
        'n_manual_entries': n_manual,
        'manual_snaps': manual_log,
        'fs_hz': fs,
        'sample_period_ms': 1000.0 / fs,
        'proposal_method': 'centre of contiguous ADC==4095 saturation runs',
        'proposal_independence': ('shares no filter, threshold or parameter '
                                  'with 02_qrs_detect.py'),
        'centre_convention': {
            'sample': 'floor(plateau centre) — used for match_beats() only',
            'sample_exact': 'true float centre — used for timing error only',
            'n_offset_beats': int((offsets > 0).sum()),
            'mean_offset_samples': mean_off,
            'mean_offset_ms': mean_off / fs * 1000.0,
        },
        'amendments': [{
            'when': 'before any detector run on this capture',
            'change': ('centre changed from int(round(mean)) to a float '
                       'sample_exact plus floor-based integer sample'),
            'reason': ("int(round()) is banker's rounding, so an even-width "
                       'plateau resolved to the nearest EVEN index. The '
                       'reference position depended on the parity of the '
                       'sample index, an accident of capture start time. '
                       'Affected 10/31 plateaus on capture 180936: 5 up, '
                       '5 down, net 0.000 ms by coincidence rather than '
                       'design. Irrelevant to Se/PPV (tolerance 37.5 samples) '
                       'but the same order as the pre-registered +/-1 sample '
                       'timing prediction.'),
            'frozen_code_modified': 'none — 03_qrs_evaluate.py untouched',
        }],
        'pre_registration': {
            'capture': 'raw_ecg_capture_20260729_180936.csv, full 30 s, no exclusions',
            'scoring': ('+/-150 ms tolerance, greedy one-to-one matching, '
                        'imported from 03_qrs_evaluate.py'),
            'prediction': ('signed timing error small, within about +/-1 sample; '
                           'bandpass smooths the 12 ms plateau so the earlier '
                           'predicted early bias should be weak or absent'),
        },
        'frozen_before_detector_run': True,
    }
    man_path = os.path.join(RESULTS_DIR, f'{stem}_annotations_manifest.json')
    with open(man_path, 'w') as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\n  reference -> {frozen}")
    print(f"  manifest  -> {man_path}")
    print(f"  capture    sha256: {meta['sha256'][:16]}...")
    print(f"  annotation sha256: {ann_hash[:16]}...")
    print("\n  Commit BOTH files now, before running the detector. The commit")
    print("  timestamp is what makes 'frozen before detection' checkable by")
    print("  someone who was not in the room.")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--capture', default=DEFAULT_CAPTURE)
    ap.add_argument('--freeze', action='store_true',
                    help='validate the edited draft and write the frozen reference')
    args = ap.parse_args()

    cap_path = os.path.join(CAPTURE_DIR, args.capture)
    if not os.path.exists(cap_path):
        print(f"Missing {cap_path}")
        print(f"Place the capture in {CAPTURE_DIR}/ and retry.")
        raise SystemExit(1)

    stem = os.path.splitext(os.path.basename(cap_path))[0]
    print(f"=== 01_annotate.py — {stem} ===")
    print(f"  mode: {'FREEZE' if args.freeze else 'PROPOSE'}\n")

    if args.freeze:
        do_freeze(cap_path, stem)
    else:
        do_propose(cap_path, stem)


if __name__ == '__main__':
    main()
