"""
03_demo_bundle.py -- assemble one demo-data.js for the FetoPulse viewer.

WHY A .js FILE AND NOT .json
The viewer has to run from the filesystem with no network and no server.
Chrome blocks fetch() against file:// URLs, so a page that loaded demo.json
would work on your machine behind a dev server and fail in the exam room. A
<script src="demo-data.js"> tag has no such restriction. The file assigns one
global and nothing else.

WHAT IT BUNDLES
  1. The CTG capture (01_capture.py output) and its scoring (02_ctg.py output).
  2. Optionally a SEPARATE maternal ECG capture, which cannot share the CTG
     capture: mrisr_capture.ino streams two channels, and adding ECG at 250 Hz
     takes the serial link to 102.6% of capacity. Two captures, minutes apart,
     and the viewer says so on screen.

The ECG R-peaks must come from a file produced by the frozen ECG pipeline. This
script does NOT detect beats. If no peak file is given the viewer shows the
trace and states that the detector was not run.

USAGE
    python 03_demo_bundle.py --stem mrisr_ctg_demo_20260907_220000
    python 03_demo_bundle.py --stem <stem> \
        --ecg-csv ../ecg_bringup/captures/ecg_xxx.csv \
        --ecg-peaks ../ecg_bringup/results/ecg_xxx_peaks.csv
"""
import os
import re
import sys
import json
import glob
import argparse
import importlib.util
from datetime import datetime

import numpy as np

RES_DIR = 'results'
OUT = 'viewer'
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
UC_PATH = os.path.join(REPO_ROOT, 'uc_detector.py')

MAX_PCG_POINTS = 6000          # decimate only for drawing; never for analysis


def load_frozen(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def read_csv(path, skip):
    return np.loadtxt(path, delimiter=',', skiprows=skip)


def thin(t, y, n=MAX_PCG_POINTS):
    """Min/max envelope decimation. Keeps visible extremes, unlike stride."""
    t = np.asarray(t, float); y = np.asarray(y, float)
    if t.size <= n:
        return t.tolist(), y.tolist()
    step = int(np.ceil(t.size / (n / 2)))
    ts, ys = [], []
    for i in range(0, t.size, step):
        seg = y[i:i + step]
        if seg.size == 0:
            continue
        tt = t[i:i + step]
        a, b = int(np.argmin(seg)), int(np.argmax(seg))
        for j in sorted((a, b)):
            ts.append(float(tt[j])); ys.append(float(seg[j]))
    return ts, ys


def r3(a):
    return [round(float(x), 3) for x in a]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stem', required=True)
    ap.add_argument('--ecg-csv', default=None)
    ap.add_argument('--ecg-peaks', default=None,
                    help='one column of R-peak sample indices or times, from '
                         'the frozen ECG pipeline. This script detects nothing.')
    ap.add_argument('--ecg-time-unit', default='auto',
                    choices=['auto', 's', 'ms', 'us'],
                    help='unit of column 0 in the ECG capture. 02_detect_device.py '
                         'reads that column as MILLISECONDS.')
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    cap_path = os.path.join(RES_DIR, args.stem + '.json')
    ctg_path = os.path.join(RES_DIR, 'ctg_' + args.stem + '.json')
    if not os.path.exists(cap_path):
        sys.exit(f'[FAIL] no capture results at {cap_path}')
    cap = json.load(open(cap_path))
    ctg = json.load(open(ctg_path)) if os.path.exists(ctg_path) else None
    print(f'capture : {args.stem}  profile={cap.get("profile")}  '
          f'validated={cap.get("validated_protocol")}')
    print(f'scoring : {"found" if ctg else "NOT FOUND (run 02_ctg.py first)"}')

    # ---------------------------------------------------------------- UC
    ucd = load_frozen('uc_detector', UC_PATH)
    uc = read_csv(cap['files']['uc'], 1)
    ut, uv = uc[:, 0], uc[:, 1]
    ut = ut - ut[0]
    spans, d = ucd.detect_contractions(uv, fs=4.0)
    conf, frac, rate = ucd.confidence(spans, d['valid'], fs=4.0)
    n_uc = uv.size
    win_ratio = (n_uc / 4.0) / ucd.BASELINE_WIN_S
    uc_valid = bool(cap.get('validated_protocol')) and win_ratio >= 1.0
    baseline = d['baseline']
    distinct = int(len(np.unique(np.round(baseline, 6)))) \
        if np.isfinite(baseline).any() else 0

    uc_block = {
        't': r3(ut), 'v': r3(uv),
        'baseline': r3(baseline) if np.isfinite(baseline).any() else None,
        'detector_run': uc_valid,
        'spans': [[round(float(ut[s]), 2),
                   round(float(ut[min(e, n_uc - 1)]), 2)] for s, e in spans]
                 if uc_valid else [],
        'record_over_window': round(win_ratio, 3),
        'baseline_window_s': float(ucd.BASELINE_WIN_S),
        'baseline_distinct': distinct,
        'reason': None if uc_valid else (
            f'Contraction detection is not run on this capture. The validated '
            f'detector estimates its baseline over a {ucd.BASELINE_WIN_S:.0f} second '
            f'rolling window, and this record is {n_uc/4.0:.0f} seconds, or '
            f'{win_ratio:.2f} times that window. Below one, the baseline '
            f'collapses to a single constant value and the detector is running '
            f'outside the configuration it was validated in. The trace is shown '
            f'as recorded.')
    }
    print(f'UC      : {n_uc} samples, {win_ratio:.2f}x baseline window, '
          f'detector {"run" if uc_valid else "NOT run"}')

    # --------------------------------------------------------------- FHR
    # The FHR detector is valid at demo length, so it runs here. But a capture
    # taken under the validated protocol is the GOLD capture, and prereg 4.4
    # says it is scored once, by 02_ctg.py. This refuses to be the thing that
    # scores it first.
    if cap.get('validated_protocol') and ctg is None:
        sys.exit('[FAIL] this capture used the validated protocol, so it is a '
                 'gold capture. Score it with 02_ctg.py first (prereg 4.4), '
                 'then bundle.')

    pcg = read_csv(cap['files']['pcg'], 3)
    pt = (pcg[:, 0] - pcg[0, 0]) / 1e6
    pv = pcg[:, 1]
    fs = (pt.size - 1) / (pt[-1] - pt[0])
    fdet = load_frozen('fhr_detector', os.path.join(REPO_ROOT, '02_fhr_detector.py'))
    r = fdet.detect_fhr(pv, fs)
    ft, fb, fc = (np.asarray(r[k]) for k in ('times', 'bpms', 'confidence'))
    thr = float(fdet.CONFIDENCE_THRESHOLD)
    high = fc >= thr
    fhr_block = {
        't': r3(ft), 'bpm': r3(fb), 'conf': [round(float(x), 3) for x in fc],
        'threshold': thr,
        'median_high': float(np.median(fb[high])) if high.any() else None,
        'frac_high': float(np.mean(high)),
        'n_windows': int(fb.size),
        'window_s': 4.0, 'hop_s': 1.0, 'search_lo': 110, 'search_hi': 180,
        'from_02_ctg': ctg is not None,
    }
    if ctg:   # cross-check against the authoritative scoring
        f = ctg['fhr']
        for k, v in (('median_high', f.get('median_bpm_high')),
                     ('frac_high', f.get('frac_high'))):
            if v is not None and fhr_block[k] is not None and abs(fhr_block[k] - v) > 1e-6:
                print(f'  [!] {k} differs from 02_ctg.py: '
                      f'{fhr_block[k]} here vs {v} there. Using 02_ctg.py.')
            if v is not None:
                fhr_block[k] = v
    ptx, pvx = thin(pt, pv)
    pcg_block = {'t': [round(x, 3) for x in ptx], 'v': [round(x, 1) for x in pvx],
                 'decimated': pt.size > MAX_PCG_POINTS, 'n_original': int(pt.size)}
    print(f'FHR     : {fb.size} windows, {100*float(np.mean(high)):.1f}% reliable, '
          f'median ' + (f'{float(np.median(fb[high])):.2f} BPM' if high.any() else 'n/a'))

    # -------------------------------------------------------------- cues
    cues = []
    if os.path.exists(cap['files']['cues']):
        rows = [l.strip().split(',') for l in
                open(cap['files']['cues']).read().splitlines()[1:] if l.strip()]
        by = {}
        for r_ in rows:
            by.setdefault(int(r_[5]), {})[r_[4]] = float(r_[3])
        t0 = float(rows[0][3])
        for k in sorted(by):
            if k == 0 or not all(p in by[k] for p in ('RAMP', 'HOLD', 'FALL', 'REST')):
                continue
            cues.append({'cycle': k,
                         'press': [round(by[k]['RAMP'] - t0, 2),
                                   round(by[k]['REST'] - t0, 2)],
                         'hold': [round(by[k]['HOLD'] - t0, 2),
                                  round(by[k]['FALL'] - t0, 2)]})

    # --------------------------------------------------------------- ECG
    ecg_block = None
    if args.ecg_csv and os.path.exists(args.ecg_csv):
        raw = np.genfromtxt(args.ecg_csv, delimiter=',', skip_header=1)
        raw = raw[np.isfinite(raw).all(axis=1)]
        if raw.shape[1] >= 2:
            t_raw, ev = raw[:, 0], raw[:, 1]
            # 02_detect_device.py's load_capture() reads column 0 as ms.
            # Decide by median sample interval rather than by total span, which
            # is ambiguous: 30,000 could be 30 s in ms or 30 ms in us.
            unit = args.ecg_time_unit
            if unit == 'auto':
                dt = float(np.median(np.diff(t_raw)))
                unit = 'us' if dt > 100 else ('ms' if dt > 0.5 else 's')
            div = {'s': 1.0, 'ms': 1e3, 'us': 1e6}[unit]
            et = (t_raw - t_raw[0]) / div
            fs = (et.size - 1) / (et[-1] - et[0])
            print(f'ECG     : time column read as {unit}, '
                  f'{et[-1]:.1f} s at {fs:.1f} Hz')
            peaks, mhr, src = None, None, None
            if args.ecg_peaks and os.path.exists(args.ecg_peaks):
                hdr = open(args.ecg_peaks).readline().strip().lower()
                cols = [c.strip() for c in hdr.split(',')]
                use = 0
                for want in ('sample_exact', 'sample'):
                    if want in cols:
                        use = cols.index(want); break
                p = np.atleast_1d(np.loadtxt(args.ecg_peaks, delimiter=',',
                                             skiprows=1, usecols=use))
                pt_s = p / fs if p.max() > et[-1] * 2 else p
                peaks = [round(float(x), 3) for x in pt_s if 0 <= x <= et[-1]]
                src = os.path.basename(args.ecg_peaks)
                if len(peaks) > 1:
                    mhr = round(60.0 / float(np.median(np.diff(peaks))), 1)
            etx, evx = thin(et, ev, 4000)
            is_ref = bool(src and 'annotation' in src.lower())
            ecg_block = {
                't': [round(x, 3) for x in etx], 'v': [round(x, 1) for x in evx],
                'fs': round(float(fs), 2), 'duration_s': round(float(et[-1]), 2),
                'peaks': peaks, 'mhr_bpm': mhr,
                'n_beats': len(peaks) if peaks else None,
                'detector_run': peaks is not None and not is_ref,
                'peaks_source': src,
                'peaks_are_reference': is_ref,
                'reason': (None if (peaks and not is_ref) else
                           ('Marks are the frozen manual annotations for this '
                            'capture, not detector output.' if is_ref else
                            'R-peak detection was not run for this capture. '
                            'The trace is shown as recorded.'))}
            print(f'          {len(peaks) if peaks else 0} marks from '
                  f'{src or "no file"}'
                  + (' (REFERENCE annotations, not detections)' if is_ref else '')
                  + (f', rate {mhr} bpm' if mhr else ''))
    if ecg_block is None:
        print('ECG     : not bundled')

    data = {
        'generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
        'capture': {
            'stem': args.stem, 'stamp': cap.get('stamp'),
            'profile': cap.get('profile'),
            'validated_protocol': bool(cap.get('validated_protocol')),
            'profile_note': cap.get('profile_note'),
            'duration_s': round(float(cap['timing']['pcg_span_s']), 2),
            'fs_pcg': round(float(cap['timing']['fs_pcg']), 4),
            'fs_uc': round(float(cap['timing']['fs_uc']), 4),
            'malformed': cap['counts']['malformed'],
            'total_lines': cap['counts']['total_lines'],
            'c5_pass': cap['c5']['pass'],
            'git_head': (cap.get('git_head') or '')[:10],
        },
        'pcg': pcg_block, 'fhr': fhr_block, 'uc': uc_block,
        'cues': cues, 'ecg': ecg_block,
        'scored': ctg is not None,
    }
    path = os.path.join(args.out, 'demo-data.js')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('window.FETOPULSE = ')
        json.dump(data, fh, separators=(',', ':'))
        fh.write(';\n')
    print(f'\nwrote {path}  ({os.path.getsize(path)/1024:.0f} KB)')
    print('Open viewer/index.html in a browser. No server needed.')


if __name__ == '__main__':
    main()
