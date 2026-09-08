"""
02_ctg.py -- score one accepted capture and draw the cardiotocograph. Phase 6.

Governed by notes/prereg.md. Run ONCE per accepted capture (prereg 4.4).

WHAT THIS PRODUCES
The figure the project exists to produce: fetal heart rate and contraction
timing on ONE shared device time axis from ONE concurrent capture, each channel
scored by the detector already validated for it.

WHAT IT DOES NOT SHOW, AND THE FIGURE SAYS SO
The FHR channel is a recorded fetal PCG played through a speaker. The UC channel
is a hand-pressed force sensor. There is no physiological relationship between
them. The figure demonstrates concurrent acquisition on a shared time base. It
does not demonstrate deceleration detection, and no early/late classification
can be drawn from it. That sentence is drawn INTO the figure, not left to the
surrounding text.

FROZEN CODE IS IMPORTED BY PATH, NEVER COPIED
  02_fhr_detector.py   detect_fhr(raw, fs)          whole record, one call
  uc_detector.py       detect_contractions, confidence
Resolved paths are printed and written to the results JSON.

UC SCORING
The authoritative UC scorer for device data is hardware/toco_bringup/04_detect.py.
If its output CSV is supplied with --uc-detect-csv, C3 is read from it. If not,
this script reproduces its DOCUMENTED primary matching rule (overlap of at least
half the shorter span, greedy one-to-one) against the cue windows on device
time, and says so in the output. The detector itself is the frozen module in
both cases; only the matching harness differs.

USAGE
    python 02_ctg.py --stem mrisr_ctg_gold_20260907_143000
    python 02_ctg.py --stem <stem> --uc-detect-csv ../toco_bringup/results/x.csv
"""
import os
import sys
import json
import glob
import argparse
import importlib.util
import subprocess

import numpy as np
from scipy.signal import welch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ config
CAP_DIR, RES_DIR, PLOT_DIR = 'captures', 'results', 'plots'
FS_PCG = 500.0
FS_UC = 4.0

# ---- frozen pipeline files, relative to the repo root (two levels up) ----
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
FHR_PATH = os.path.join(REPO_ROOT, '02_fhr_detector.py')
UC_PATH = os.path.join(REPO_ROOT, 'uc_detector.py')

# ---- pre-registered criteria, prereg 6. Not tunable here. ----
C1_MIN_RELIABLE_FRAC = 0.50          # 03_evaluate.py record rule
C2_REFERENCE_BPM = 133.93            # run C
C2_TOLERANCE_BPM = 3.0               # run C's own tolerance
C3_EXPECTED_CONTRACTIONS = 7
NULL_RELIABLE_FRAC = 0.0             # 0 of 261 windows, prereg 5

# ---- diagnostics, prereg 9. Reported, never deciding. ----
BAND = (25.0, 200.0)
LINE_HALF_WIDTH_HZ = 1.0
CAPTION = ("FHR: recorded fetal PCG played through a speaker. UC: hand-pressed "
           "force sensor.\nNo physiological relationship exists between the "
           "two traces. This figure demonstrates concurrent\nacquisition on a "
           "shared time base only. It does not demonstrate deceleration "
           "detection.")


def load_frozen(name, path):
    if not os.path.exists(path):
        print(f'[FAIL] frozen file not found: {path}')
        sys.exit(1)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print(f'  frozen {name:<14} -> {spec.origin}')
    return mod


def git_head():
    try:
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None


def trapz(y, x):
    f = getattr(np, 'trapezoid', None) or np.trapz
    return float(f(y, x))


# ------------------------------------------------------------------ loaders
def load_pcg(path):
    d = np.loadtxt(path, delimiter=',', skiprows=3)
    return d[:, 0].astype(np.int64), d[:, 1].astype(np.float64)


def load_uc(path):
    d = np.loadtxt(path, delimiter=',', skiprows=1)
    return d[:, 0].astype(np.float64), d[:, 1].astype(np.float64)


def load_cues(path):
    cues = []
    with open(path) as fh:
        hdr = fh.readline().strip().split(',')
        for line in fh:
            v = line.strip().split(',')
            if len(v) != len(hdr):
                continue
            c = dict(zip(hdr, v))
            cues.append(dict(scheduled_s=float(c['scheduled_s']),
                             host_s=float(c['host_s']),
                             device_s=float(c['device_s']),
                             phase=c['phase'], cycle=int(c['cycle'])))
    return cues


def cue_windows(cues):
    """Per cycle: RAMP start .. FALL end, HOLD start .. HOLD end, REST span.
    All on DEVICE time, taken from the device timestamp recorded at each cue."""
    by = {}
    for c in cues:
        by.setdefault(c['cycle'], {})[c['phase']] = c['device_s']
    out = []
    for k in sorted(by):
        if k == 0:
            continue
        ph = by[k]
        if not all(p in ph for p in ('RAMP', 'HOLD', 'FALL', 'REST')):
            continue
        out.append(dict(cycle=k,
                        press=(ph['RAMP'], ph['REST']),      # whole press
                        hold=(ph['HOLD'], ph['FALL']),
                        rest_start=ph['REST']))
    return out


# ------------------------------------------------------------------ matching
def overlap(a, b):
    return max(0.0, min(a[1], b[1]) - max(a[0], b[0]))


def match_spans(refs, dets):
    """Documented primary rule from 04_detect.py: overlap >= half the shorter
    span, greedy one-to-one. Returns list of (ref_idx, det_idx, iou)."""
    pairs = []
    for i, r in enumerate(refs):
        for j, d in enumerate(dets):
            ov = overlap(r, d)
            shorter = min(r[1] - r[0], d[1] - d[0])
            if shorter > 0 and ov >= 0.5 * shorter:
                union = (r[1] - r[0]) + (d[1] - d[0]) - ov
                pairs.append((ov, i, j, ov / union if union > 0 else 0.0))
    pairs.sort(reverse=True)
    used_r, used_d, out = set(), set(), []
    for ov, i, j, iou in pairs:
        if i in used_r or j in used_d:
            continue
        used_r.add(i); used_d.add(j)
        out.append((i, j, iou))
    return sorted(out)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stem', required=True,
                    help='e.g. mrisr_ctg_gold_20260907_143000')
    ap.add_argument('--uc-detect-csv', default=None,
                    help='output of toco_bringup/04_detect.py for this capture; '
                         'if given, C3 is read from it')
    args = ap.parse_args()
    os.makedirs(RES_DIR, exist_ok=True)
    os.makedirs(PLOT_DIR, exist_ok=True)

    print(f'script  : {os.path.abspath(__file__)}')
    head = git_head()
    print(f'git HEAD: {head or "(unavailable)"}')

    # ---- locate the three capture files from the capture results JSON ----
    cap_json = os.path.join(RES_DIR, args.stem + '.json')
    if not os.path.exists(cap_json):
        print(f'[FAIL] capture results not found: {cap_json}')
        sys.exit(1)
    with open(cap_json) as fh:
        cap = json.load(fh)
    if not cap.get('c5', {}).get('pass', False):
        print('[FAIL] this capture did NOT pass C5. It is void (prereg 8). '
              'Refusing to score it.')
        sys.exit(1)
    if cap.get('detectors_run', False):
        print('[!] capture JSON says detectors were already run. prereg 4.4: '
              'one capture, one analysis. Proceeding, but record why.')
    pcg_path, uc_path, cue_path = (cap['files']['pcg'], cap['files']['uc'],
                                   cap['files']['cues'])
    validated = bool(cap.get('validated_protocol', False))
    profile = cap.get('profile', '?')

    print(f'\ncapture : {args.stem}  profile={profile}  '
          f'validated_protocol={validated}')
    print('frozen modules:')
    fhr = load_frozen('fhr_detector', FHR_PATH)
    ucd = load_frozen('uc_detector', UC_PATH)

    pt_us, pv = load_pcg(pcg_path)
    ut_s, uv = load_uc(uc_path)
    cues = load_cues(cue_path)
    t0_us = pt_us[0]
    tp = (pt_us - t0_us) / 1e6                  # shared device axis, seconds
    tu = ut_s - t0_us / 1e6
    fs_p = (len(pt_us) - 1) / (tp[-1] - tp[0])
    windows = cue_windows(cues)
    # cue device times were recorded as absolute device us; re-zero them
    for w in windows:
        w['press'] = tuple(x - t0_us / 1e6 for x in w['press'])
        w['hold'] = tuple(x - t0_us / 1e6 for x in w['hold'])
        w['rest_start'] -= t0_us / 1e6
    print(f'  pcg {len(pv)} samples, {tp[-1]:.1f} s at {fs_p:.4f} Hz')
    print(f'  uc  {len(uv)} samples, {len(windows)} complete press cycles')

    # ================================================== FHR, one call
    print('\n' + '=' * 70)
    print('FHR -- frozen 02_fhr_detector.detect_fhr, ONE call on the whole record')
    print('=' * 70)
    r = fhr.detect_fhr(pv, fs_p)
    ft = np.asarray(r['times']); fb = np.asarray(r['bpms'])
    fc = np.asarray(r['confidence'])
    thr = float(fhr.CONFIDENCE_THRESHOLD)
    high = fc >= thr
    frac_high = float(np.mean(high)) if len(fc) else 0.0
    med_high = float(np.median(fb[high])) if high.any() else float('nan')
    med_all = float(np.median(fb)) if len(fb) else float('nan')
    print(f'  windows            : {len(fb)}')
    print(f'  reliable (>= {thr}) : {100*frac_high:.1f}%   '
          f'(null 0.0% of 261, run C 100%)')
    print(f'  median BPM, high   : {med_high:.2f}')
    print(f'  median BPM, all    : {med_all:.2f}')
    print(f'  median conf        : {np.median(fc):.3f}   max {fc.max():.3f}')

    c1 = frac_high >= C1_MIN_RELIABLE_FRAC
    c2_delta = abs(med_high - C2_REFERENCE_BPM)
    c2 = c2_delta <= C2_TOLERANCE_BPM
    print(f'\n  C1 verdict OK (>= {int(100*C1_MIN_RELIABLE_FRAC)}% reliable) : '
          f'{"PASS" if c1 else "FAIL"}')
    print(f'  C2 |median - {C2_REFERENCE_BPM}| <= {C2_TOLERANCE_BPM}       : '
          f'{c2_delta:.2f} BPM  {"PASS" if c2 else "FAIL"}')

    # ================================================== UC, frozen detector
    print('\n' + '=' * 70)
    print('UC -- frozen uc_detector.detect_contractions')
    print('=' * 70)
    n_uc = len(uv)
    win_ratio = (n_uc / FS_UC) / ucd.BASELINE_WIN_S
    print(f'  record / BASELINE_WIN_S : {win_ratio:.2f}x')
    spans, d = ucd.detect_contractions(uv, fs=FS_UC)
    conf, frac, rate = ucd.confidence(spans, d['valid'], fs=FS_UC)
    baseline = d['baseline']
    n_distinct = int(len(np.unique(np.round(baseline, 6)))) \
        if np.isfinite(baseline).any() else 0
    det = [(tu[s], tu[min(e, n_uc - 1)]) for s, e in spans]
    print(f'  detected spans          : {len(spans)}')
    print(f'  confidence()            : {conf}  analysable {100*frac:.1f}%  '
          f'rate {rate:.2f}/10min')
    print(f'  baseline distinct values: {n_distinct}  '
          f'{"(DEGENERATE: rolling baseline could not roll)" if n_distinct <= 1 else ""}')

    uc_detector_valid = validated and win_ratio >= 1.0
    if not uc_detector_valid:
        print('  [!] uc_detector output on this record is NOT a validated '
              'configuration (prereg 12). Reported for the plot only; '
              'C3 is not evaluated.')

    # ---- C3 ----
    c3 = None; matched = []; fp = None; iou_list = []
    refs = [w['press'] for w in windows]
    src_c3 = None
    if args.uc_detect_csv and os.path.exists(args.uc_detect_csv):
        src_c3 = os.path.abspath(args.uc_detect_csv)
        rows = np.genfromtxt(args.uc_detect_csv, delimiter=',', names=True)
        rows = np.atleast_1d(rows)
        n_matched = int(np.sum(rows['matched'] == 1))
        n_det_csv = len(rows)
        iou_list = [float(x) for x in rows['iou']]
        fp = n_det_csv - n_matched
        c3 = (n_matched == C3_EXPECTED_CONTRACTIONS and fp == 0
              and conf == 'HIGH')
        print(f'\n  C3 from 04_detect.py CSV: matched {n_matched}/'
              f'{C3_EXPECTED_CONTRACTIONS}, FP {fp}  '
              f'{"PASS" if c3 else "FAIL"}')
    elif uc_detector_valid:
        src_c3 = 'reimplemented primary rule (04_detect.py documented rule)'
        pairs = match_spans(refs, det)
        matched = pairs
        n_matched = len(pairs)
        fp = len(det) - n_matched
        iou_list = [p[2] for p in pairs]
        c3 = (n_matched == C3_EXPECTED_CONTRACTIONS and fp == 0
              and conf == 'HIGH')
        print(f'\n  C3 (documented primary rule, reimplemented here): '
              f'matched {n_matched}/{C3_EXPECTED_CONTRACTIONS}, FP {fp}, '
              f'IoU {" ".join(f"{x:.2f}" for x in iou_list)}  '
              f'{"PASS" if c3 else "FAIL"}')
        print('  -> run toco_bringup/04_detect.py on this UC file and re-run '
              'with --uc-detect-csv for the authoritative figure.')

    # ---- 9.3 cue-to-detection offset ----
    offsets = []
    if matched:
        for i, j, _ in matched:
            offsets.append(det[j][0] - refs[i][0])
        print(f'  cue-to-onset offset     : ' +
              ' '.join(f'{o:+.1f}' for o in offsets) + ' s  (device time)')

    # ================================================== C4, hold vs rest
    print('\n' + '=' * 70)
    print('C4 -- FHR confidence during HOLD vs REST. Reported. No threshold.')
    print('=' * 70)
    in_hold = np.zeros(len(ft), bool); in_rest = np.zeros(len(ft), bool)
    for k, w in enumerate(windows):
        in_hold |= (ft >= w['hold'][0]) & (ft < w['hold'][1])
        rest_end = windows[k + 1]['press'][0] if k + 1 < len(windows) else tp[-1]
        in_rest |= (ft >= w['rest_start']) & (ft < rest_end)
    def summ(m):
        if not m.any():
            return dict(n=0)
        return dict(n=int(m.sum()), median_conf=float(np.median(fc[m])),
                    mean_conf=float(np.mean(fc[m])),
                    p10_conf=float(np.percentile(fc[m], 10)),
                    frac_high=float(np.mean(fc[m] >= thr)),
                    median_bpm=float(np.median(fb[m])))
    hold_s, rest_s = summ(in_hold), summ(in_rest)
    for name, s in (('HOLD', hold_s), ('REST', rest_s)):
        if s['n']:
            print(f'  {name:<5} n={s["n"]:4d}  median conf {s["median_conf"]:.3f}'
                  f'  p10 {s["p10_conf"]:.3f}  reliable {100*s["frac_high"]:.1f}%'
                  f'  median BPM {s["median_bpm"]:.2f}')
    if hold_s['n'] and rest_s['n']:
        dconf = hold_s['median_conf'] - rest_s['median_conf']
        print(f'  HOLD - REST median confidence: {dconf:+.3f}  '
              f'(prediction: no systematic difference)')

    # ================================================== 9.1 dominant line
    print('\n' + '=' * 70)
    print('9.1 -- dominant in-band line, LOCATED, and its share. Reported.')
    print('=' * 70)
    x = pv - np.mean(pv)
    fr, ps = welch(x, fs=fs_p, nperseg=min(2500, len(x)))
    inb = (fr >= BAND[0]) & (fr <= BAND[1])
    e_in = trapz(ps[inb], fr[inb]); e_tot = trapz(ps, fr)
    ib = np.where(inb)[0]
    k = ib[int(np.argmax(ps[ib]))]
    f0 = float(fr[k])
    line = (fr >= f0 - LINE_HALF_WIDTH_HZ) & (fr <= f0 + LINE_HALF_WIDTH_HZ)
    share = 100 * trapz(ps[line], fr[line]) / e_in
    pk50 = (fr >= 49) & (fr <= 51); bg50 = (fr >= 44) & (fr <= 56) & ~pk50
    db50 = 10 * np.log10(ps[pk50].max() / (np.median(ps[bg50]) + 1e-30))
    print(f'  in-band share of total  : {100*e_in/e_tot:.1f}%')
    print(f'  strongest in-band line  : {f0:.2f} Hz')
    print(f'  its share of in-band    : {share:.2f}%   '
          f'(quiet floor mrisr 8.21-8.59%, pcg_bringup 2.58-4.43%, '
          f'occluded 97.91%)')
    print(f'  50 Hz over background   : {db50:+.1f} dB')
    print(f'  prediction (prereg 7)   : line at 47.80 Hz, share below 8.21%  '
          f'-> {"HELD" if (abs(f0-47.80) < 0.5 and share < 8.21) else "DID NOT HOLD"}')

    # ================================================== figure
    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[3, 2.4, 0.5], hspace=0.12)
    ax1 = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)

    ax1.axhspan(110, 160, alpha=0.08, color='green', label='110-160 BPM')
    ax1.axhline(C2_REFERENCE_BPM, color='k', lw=0.8, ls=':',
                label=f'run C reference {C2_REFERENCE_BPM} BPM')
    if (~high).any():
        ax1.plot(ft[~high], fb[~high], 'o', ms=3, mfc='none', mec='darkorange',
                 alpha=0.6, label=f'window confidence < {thr}')
    if high.any():
        ax1.plot(ft[high], fb[high], 'o', ms=3, color='steelblue',
                 label=f'window confidence >= {thr}')
    ax1.set_ylim(100, 190); ax1.set_ylabel('FHR (BPM)')
    ax1.set_title(f'{args.stem}  --  FHR (frozen detector, {len(fb)} windows, '
                  f'{100*frac_high:.0f}% reliable, median {med_high:.2f} BPM)',
                  fontsize=10)
    ax1.legend(loc='upper right', fontsize=8, ncol=2); ax1.grid(alpha=0.3)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2.plot(tu, uv, lw=0.7, color='tab:blue', alpha=0.7, label='UC (4 Hz)')
    if np.isfinite(baseline).any():
        ax2.plot(tu, baseline, lw=1.0, color='tab:green', label='rolling baseline')
        if np.isfinite(d['thr']):
            ax2.plot(tu, baseline + d['thr'], lw=0.9, ls='--', color='tab:red',
                     label='threshold')
    ymax = float(np.nanmax(uv))
    for i, (s, e) in enumerate(det):
        ax2.axvspan(s, e, alpha=0.18, color='tab:purple',
                    label='detected contraction' if i == 0 else None)
    lab = 'UC (frozen detector' + ('' if uc_detector_valid else
                                   ', NOT a validated configuration at this length')
    ax2.set_title(f'{lab}: {len(det)} spans, confidence {conf})', fontsize=10)
    ax2.set_ylabel('UC (counts)'); ax2.legend(loc='upper right', fontsize=8, ncol=4)
    ax2.grid(alpha=0.3); plt.setp(ax2.get_xticklabels(), visible=False)

    for i, w in enumerate(windows):
        ax3.axvspan(w['press'][0], w['press'][1], color='0.75',
                    label='cue window (press)' if i == 0 else None)
        ax3.axvspan(w['hold'][0], w['hold'][1], color='0.45',
                    label='HOLD' if i == 0 else None)
    ax3.set_yticks([]); ax3.set_ylabel('cues', rotation=0, labelpad=18)
    ax3.set_xlabel('Time (s, device clock, shared)')
    ax3.legend(loc='upper right', fontsize=8, ncol=2)
    ax3.set_xlim(0, tp[-1])

    fig.text(0.5, 0.005, CAPTION, ha='center', va='bottom', fontsize=8.5,
             family='monospace')
    fig.subplots_adjust(bottom=0.12, top=0.95)
    ppath = os.path.join(PLOT_DIR, f'ctg_{args.stem}.png')
    fig.savefig(ppath, dpi=130); plt.close(fig)

    # ================================================== results
    res = {
        'stem': args.stem, 'capture_json': os.path.abspath(cap_json),
        'script': os.path.abspath(__file__), 'git_head': head,
        'profile': profile, 'validated_protocol': validated,
        'frozen': {'fhr_detector': FHR_PATH, 'uc_detector': UC_PATH},
        'fhr': {'n_windows': int(len(fb)), 'threshold': thr,
                'frac_high': frac_high, 'median_bpm_high': med_high,
                'median_bpm_all': med_all,
                'median_conf': float(np.median(fc)) if len(fc) else None,
                'max_conf': float(fc.max()) if len(fc) else None,
                'whole_record_single_call': True},
        'uc': {'n_samples': int(n_uc), 'record_over_baseline_win': win_ratio,
               'n_spans': len(det), 'spans_s': det, 'confidence': conf,
               'analysable_frac': float(frac), 'rate_per_10min': float(rate),
               'baseline_distinct_values': n_distinct,
               'validated_configuration': bool(uc_detector_valid),
               'amp': float(d['amp']) if np.isfinite(d['amp']) else None,
               'thr': float(d['thr']) if np.isfinite(d['thr']) else None},
        'cues': {'n_cycles': len(windows),
                 'press_windows_s': [w['press'] for w in windows],
                 'hold_windows_s': [w['hold'] for w in windows]},
        'criteria': {
            'C1': {'pass': bool(c1), 'frac_high': frac_high,
                   'rule': f'>= {C1_MIN_RELIABLE_FRAC} of windows >= {thr}',
                   'null_frac_high': NULL_RELIABLE_FRAC},
            'C2': {'pass': bool(c2), 'delta_bpm': c2_delta,
                   'reference': C2_REFERENCE_BPM, 'tolerance': C2_TOLERANCE_BPM},
            'C3': {'pass': c3, 'source': src_c3, 'fp': fp,
                   'iou': iou_list, 'evaluated': c3 is not None},
            'C4': {'hold': hold_s, 'rest': rest_s, 'threshold': None,
                   'prediction': 'no systematic difference'},
        },
        'diagnostics': {
            'dominant_line_hz': f0, 'dominant_line_share_pct': share,
            'inband_pct_of_total': 100 * e_in / e_tot,
            'db_50hz_over_bg': float(db50),
            'cue_to_onset_offset_s': offsets,
        },
        'plot': ppath,
        'caption': CAPTION.replace('\n', ' '),
    }
    rpath = os.path.join(RES_DIR, f'ctg_{args.stem}.json')
    with open(rpath, 'w') as fh:
        json.dump(res, fh, indent=2)
    print(f'\nresults -> {rpath}')
    print(f'figure  -> {ppath}')
    print('\nOne capture, one analysis (prereg 4.4). If this is re-run, record why.')


if __name__ == '__main__':
    main()
