"""
00_prereg_measurements.py -- persist the measurements cited in notes/prereg.md.

WHY THIS FILE EXISTS
The provenance audit (notes/provenance_audit.md, 2026-09-07) found that the
pre-registration cites measurements that no committed script produces: the
null-run detector output on three mrisr captures, the strongest in-band line
and its share on the same captures, the 50 Hz level, the seam-normalisation
test, the 90 s baseline-degeneracy test, and the press-consistency table. They
were computed in an analysis sandbox and typed into the document. That is the
same failure the PCG notes record as 11.5, and persisting the numbers is the
only fix that does not involve editing a committed pre-registration.

This script is committed AFTER the pre-registration and BEFORE any Phase 6
capture. The git history establishes that ordering. It computes nothing new and
changes nothing; it writes to disk what the pre-registration already states.

INPUTS
Three committed mrisr captures (gitignored raw files; a targeted exception is
required for the three named below), and the two frozen detectors, imported by
path.

USAGE
    python 00_prereg_measurements.py --repo-root ../..
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

FS = 500.0
BAND = (25.0, 200.0)
NULL_CAPTURES = ['mrisr_baseline_20260906_003840',
                 'mrisr_baseline_nofan_20260906_005014',
                 'mrisr_paired_multi_20260906_005816']


def load_frozen(name, path):
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


def spectral(a, fs):
    x = a - a.mean()
    fr, ps = welch(x, fs=fs, nperseg=min(2500, len(x)))
    inb = (fr >= BAND[0]) & (fr <= BAND[1])
    e_in = trapz(ps[inb], fr[inb]); e_tot = trapz(ps, fr)
    ib = np.where(inb)[0]; k = ib[int(np.argmax(ps[ib]))]; f0 = float(fr[k])
    line = (fr >= f0 - 1) & (fr <= f0 + 1)
    share = 100 * trapz(ps[line], fr[line]) / e_in
    pk = (fr >= 49) & (fr <= 51); bg = (fr >= 44) & (fr <= 56) & ~pk
    db = 10 * np.log10(ps[pk].max() / (np.median(ps[bg]) + 1e-30))
    return dict(peak_inband_hz=f0, peak_line_share_pct=float(share),
                inband_pct_of_total=100 * e_in / e_tot,
                db_50hz_over_bg=float(db),
                method='welch nperseg=2500; line = +-1 Hz around the strongest '
                       'bin in 25-200 Hz; 50 Hz level = max in 49-51 Hz over '
                       'median of 44-56 Hz excluding 49-51')


# ---------------------------------------------------------- synthetic helpers
def make_loop(rng, bpm=133.93, loop_s=60.0, sysfrac=0.365):
    n = int(loop_s * FS); x = np.zeros(n); period = 60.0 / bpm; nb = int(0.05 * FS)
    tt = np.arange(nb) / FS
    w = np.exp(-0.5 * ((tt - tt[-1] / 2) / (tt[-1] / 6)) ** 2)
    s1 = np.sin(2 * np.pi * 70 * tt) * w; s2 = 0.6 * np.sin(2 * np.pi * 110 * tt) * w
    k = 0
    while True:
        i1 = int(k * period * FS); i2 = int((k * period + sysfrac * period) * FS)
        if i2 + nb >= n:
            break
        x[i1:i1 + nb] += s1; x[i2:i2 + nb] += s2; k += 1
    return x + 0.05 * rng.standard_normal(n)


def uc_synth(rng, amps, total_s=1215.0, noise=3.92, drift=30.0, fs=4.0):
    n = int(total_s * fs); t = np.arange(n) / fs
    x = 1085 + drift * (t / total_s) + noise * rng.standard_normal(n)
    for k, A in enumerate(amps):
        a = 60.0 + k * 165.0
        r = (t >= a) & (t < a + 20); h = (t >= a + 20) & (t < a + 45)
        f = (t >= a + 45) & (t < a + 65)
        x[r] += A * (t[r] - a) / 20; x[h] += A; x[f] += A * (1 - (t[f] - a - 45) / 20)
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', default=os.path.join(os.path.dirname(__file__), '..', '..'))
    ap.add_argument('--captures-dir', default=None,
                    help='override: directory holding the three mrisr pcg500 CSVs')
    ap.add_argument('--fhr', default=None); ap.add_argument('--uc', default=None)
    args = ap.parse_args()
    root = os.path.abspath(args.repo_root)
    fhr_path = args.fhr or os.path.join(root, '02_fhr_detector.py')
    uc_path = args.uc or os.path.join(root, 'uc_detector.py')
    cap_dir = args.captures_dir or os.path.join(root, 'hardware', 'mrisr', 'captures')
    os.makedirs('results', exist_ok=True)

    print(f'script  : {os.path.abspath(__file__)}')
    head = git_head(); print(f'git HEAD: {head}')
    print('frozen modules:')
    fhr = load_frozen('fhr_detector', fhr_path)
    ucd = load_frozen('uc_detector', uc_path)
    thr = float(fhr.CONFIDENCE_THRESHOLD)
    out = {'script': os.path.abspath(__file__), 'git_head': head,
           'frozen': {'fhr_detector': fhr_path, 'uc_detector': uc_path},
           'purpose': 'persist measurements cited in notes/prereg.md (5, 4.1, 7, 9.1, 12)'}

    # ================================================ 1. null runs (prereg 5, 9.1)
    print('\n1. null runs on committed mrisr captures (prereg 5) and their spectra (9.1)')
    null = []; total_win = 0; total_high = 0
    for stem in NULL_CAPTURES:
        hits = glob.glob(os.path.join(cap_dir, stem + '_pcg500.csv'))
        if not hits:
            print(f'  [MISSING] {stem}_pcg500.csv -- add a .gitignore exception')
            null.append({'capture': stem, 'missing': True}); continue
        d = np.loadtxt(hits[0], delimiter=',', skiprows=3)
        t, a = d[:, 0], d[:, 1]
        fs = (len(t) - 1) / ((t[-1] - t[0]) / 1e6)
        r = fhr.detect_fhr(a.astype(float), fs)
        b, c = np.asarray(r['bpms']), np.asarray(r['confidence'])
        sp = spectral(a, fs)
        row = dict(capture=stem, file=os.path.abspath(hits[0]), n_samples=int(len(t)),
                   fs_achieved=float(fs), n_windows=int(len(b)),
                   median_bpm=float(np.median(b)),
                   bpm_iqr=float(np.percentile(b, 75) - np.percentile(b, 25)),
                   bpm_std=float(np.std(b)),
                   median_conf=float(np.median(c)), mean_conf=float(np.mean(c)),
                   max_conf=float(c.max()),
                   frac_high=float(np.mean(c >= thr)),
                   pcg_std=float(np.std(a)), pcg_p2p=float(a.max() - a.min()),
                   spectrum=sp)
        total_win += len(b); total_high += int(np.sum(c >= thr))
        null.append(row)
        print(f'  {stem:<40} win {len(b):3d}  medBPM {row["median_bpm"]:7.2f}  '
              f'maxConf {row["max_conf"]:.3f}  >=thr {100*row["frac_high"]:.1f}%  '
              f'line {sp["peak_inband_hz"]:.2f} Hz {sp["peak_line_share_pct"]:.2f}%  '
              f'50Hz {sp["db_50hz_over_bg"]:+.1f} dB')
    out['null_runs'] = {'threshold': thr, 'captures': null,
                        'total_windows': total_win, 'total_high': total_high,
                        'frac_high_pooled': (total_high / total_win) if total_win else None}
    print(f'  pooled: {total_high} of {total_win} windows >= {thr}')

    # ================================================ 2. seam test (prereg 4.1)
    print('\n2. seam-normalisation test, synthetic (prereg 4.1)')
    rng = np.random.default_rng(0)
    loop = make_loop(rng); peak = np.max(np.abs(loop))
    seam = []
    for K in (2.0, 20.0, 100.0, 1000.0):
        segs = []
        for i in range(5):
            s = loop.copy()
            if i == 2:
                s[len(s) // 2] += K * peak
            segs.append(s)
        full = np.concatenate(segs)
        rf = fhr.detect_fhr(full, FS)
        tf, bf, cf = (np.asarray(rf[k]) for k in ('times', 'bpms', 'confidence'))
        ts, bs, cs = [], [], []
        for i, s in enumerate(segs):
            rs = fhr.detect_fhr(s, FS)
            ts.append(np.asarray(rs['times']) + i * 60.0)
            bs.append(np.asarray(rs['bpms'])); cs.append(np.asarray(rs['confidence']))
        ts, bs, cs = map(np.concatenate, (ts, bs, cs))
        idx = {round(x, 3): j for j, x in enumerate(ts)}
        pr = [(j, idx[round(x, 3)]) for j, x in enumerate(tf) if round(x, 3) in idx]
        jf = np.array([p[0] for p in pr]); js = np.array([p[1] for p in pr])
        db = np.abs(bf[jf] - bs[js]); dc = np.abs(cf[jf] - cs[js])
        cross = float(np.mean((cf[jf] >= thr) != (cs[js] >= thr)))
        row = dict(transient_x_peak=K, n_matched=len(pr), median_abs_dbpm=float(np.median(db)),
                   max_abs_dbpm=float(db.max()), max_abs_dconf=float(dc.max()),
                   frac_windows_crossing_threshold=cross)
        seam.append(row)
        print(f'  K={K:6.0f}  max|dBPM| {row["max_abs_dbpm"]:.4f}  '
              f'max|dconf| {row["max_abs_dconf"]:.3f}  crossing {100*cross:.1f}%')
    out['seam_test'] = {'design': '5 x 60 s synthetic S1/S2 doublet at 133.93 BPM, one '
                                  'transient mid-record, whole-record vs 60 s segmented '
                                  'calls matched on window centre', 'rows': seam}

    # ================================================ 3. 90 s degeneracy (prereg 12)
    print('\n3. baseline degeneracy at 90 s vs 1215 s, synthetic (prereg 12)')
    deg = []
    for label, total, cyc in (('gold_1215s', 1215.0, 7), ('demo_90s', 90.0, 1)):
        rng = np.random.default_rng(1)
        x = uc_synth(rng, [400] * cyc, total_s=total)
        spans, d = ucd.detect_contractions(x)
        conf, frac, rate = ucd.confidence(spans, d['valid'])
        b = d['baseline']
        row = dict(record=label, n_samples=int(x.size),
                   filter_size=int(ucd.BASELINE_WIN_S * ucd.FS),
                   baseline_distinct_values=int(len(np.unique(np.round(b, 6)))),
                   detected=len(spans), expected=cyc, confidence=conf,
                   rate_per_10min=float(rate))
        deg.append(row)
        print(f'  {label:<12} distinct baseline values {row["baseline_distinct_values"]:5d}  '
              f'detected {len(spans)}/{cyc}  confidence {conf}')
    out['baseline_degeneracy'] = deg

    # ================================================ 4. press consistency (notes)
    print('\n4. press-consistency tolerance, synthetic')
    rng = np.random.default_rng(7)
    cases = [('all identical', [400] * 7),
             ('+/-10% random', list(400 * (1 + 0.10 * rng.uniform(-1, 1, 7)))),
             ('+/-25% random', list(400 * (1 + 0.25 * rng.uniform(-1, 1, 7)))),
             ('+/-50% random', list(400 * (1 + 0.50 * rng.uniform(-1, 1, 7)))),
             ('one press 2x', [400, 400, 800, 400, 400, 400, 400]),
             ('one press 3x', [400, 400, 1200, 400, 400, 400, 400]),
             ('one press 0.5x', [400, 400, 200, 400, 400, 400, 400]),
             ('fade 500->250', list(np.linspace(500, 250, 7))),
             ('climb 250->500', list(np.linspace(250, 500, 7)))]
    for A in (10, 20, 40, 80, 150, 400, 1000):
        cases.append((f'absolute {A} counts', [A] * 7))
    press = []
    for name, amps in cases:
        spans, d = ucd.detect_contractions(uc_synth(rng, amps))
        conf, _, _ = ucd.confidence(spans, d['valid'])
        press.append(dict(condition=name, amplitudes=[float(a) for a in amps],
                          detected=len(spans), confidence=conf, amp=float(d['amp'])))
        print(f'  {name:<24} detected {len(spans)}  {conf}')
    out['press_consistency'] = press

    path = os.path.join('results', 'prereg_measurements.json')
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f'\nresults -> {path}')


if __name__ == '__main__':
    main()
