"""
01_capture.py -- concurrent CTG capture. Phase 6.

Governed by notes/prereg.md, committed before this script ran for the first
time. Read it before changing anything here.

WHAT THIS SCRIPT DOES
Captures the interleaved PCG/UC stream from mrisr_capture.ino for the full
1,215 s toco protocol with BOTH channels live -- speaker looping the fetal PCG
WAV at the microphone, FSR402 pressed by hand on the scripted cue schedule --
and writes the two channels in the exact formats the frozen detectors already
consume.

WHAT THIS SCRIPT DOES NOT DO
It runs NO detector. Not the FHR detector, not the UC detector, not a partial
one, not "just to check". Scoring happens in 02_ctg.py, once, after the capture
is accepted. Splitting them is deliberate: prereg 8 says a capture may be voided
only for an operational reason, never because of a detector result, and that
rule is easier to keep when the capture tool cannot produce a detector result.

WHY THE PROTOCOL IS UNCHANGED
The 20-minute protocol (60 s baseline, then 7 cycles of 20 s ramp / 25 s hold /
20 s fall / 100 s rest) is reused with its timings untouched. The single-channel
run of this exact protocol returned 7/7 with 0 false positives. Repeating it
with the PCG channel live makes this a REPLICATION with a known prior answer, so
any deviation is attributable to the second channel rather than to a protocol
invented for this phase. 1,215 s is also 2.02x uc_detector's BASELINE_WIN_S of
600 s, which is the condition under which its rolling baseline can actually
roll.

CUES ARE NAVIGATION AIDS, NOT GROUND TRUTH
The toco phase lost 55.47 s of data to a blocking input() while the ESP32 kept
streaming, and measured an 86.35 s host-clock offset from the same cause. So:

  - The preflight checklist blocks BEFORE the serial port is opened.
  - Nothing blocks once capture begins.
  - reset_input_buffer() is called immediately before the read loop.

Each cue records BOTH its host elapsed time and the device timestamp of the most
recent PCG sample received at that instant. The difference is host-to-device
buffer lag, measured rather than assumed. Contraction onset is still determined
from the signal in 02_ctg.py, never from cue timing.

USAGE
    python 01_capture.py --profile gold      # 1,215 s, the real thing
    python 01_capture.py --profile dry       #   105 s, setup check
    python 01_capture.py --profile demo      #    90 s, live demonstration
"""
import os
import re
import sys
import time
import json
import argparse
import subprocess
from datetime import datetime

import numpy as np
from scipy.signal import welch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import serial

# ------------------------------------------------------------------ config
PORT = 'COM8'
BAUD = 115200

EXPECT_PCG_PIN = 39
EXPECT_UC_PIN = 35
EXPECT_FS = 500
EXPECT_DECIM = 125
NOMINAL_FS_PCG = 500.0
NOMINAL_FS_UC = 4.0

CAP_DIR, RES_DIR, PLOT_DIR = 'captures', 'results', 'plots'
HEADER_TIMEOUT_S = 25.0

# ---- C5 integrity gates, prereg 6. Objective and signal-independent. ----
FS_TOLERANCE_FRAC = 0.001         # 0.1% of nominal, both channels
MAX_UC_GAP_S = 5.0                # above this, 04_detect.py trims

# ---- prereg 7 predictions, printed against measurement, never enforced ----
PREDICTED_MALFORMED_RATE = 1 / 20000.0

# ------------------------------------------------------------------ profiles
# Only 'gold' is the validated protocol. The other two exist to check the rig
# and to drive the demonstration, and both are labelled as NOT the validated
# protocol in the results JSON so no downstream reader can mistake them for it.
PROFILES = {
    'gold': dict(baseline_s=60.0, cycles=7,
                 ramp_s=20.0, hold_s=25.0, fall_s=20.0, rest_s=100.0,
                 validated=True,
                 note='the frozen toco protocol, timings unchanged'),
    'dry':  dict(baseline_s=20.0, cycles=1,
                 ramp_s=20.0, hold_s=25.0, fall_s=20.0, rest_s=20.0,
                 validated=False,
                 note='rig check only. NOT the validated protocol.'),
    'demo': dict(baseline_s=10.0, cycles=1,
                 ramp_s=15.0, hold_s=20.0, fall_s=15.0, rest_s=30.0,
                 validated=False,
                 note='live demonstration. NOT the validated protocol. '
                      'uc_detector MUST NOT be run on this: 90 s is 0.15x '
                      'BASELINE_WIN_S and the rolling baseline degenerates '
                      'to a constant. See prereg 12.'),
}


def build_schedule(p):
    """Cue list as (t_s, phase, cycle). Total duration is exact by construction."""
    cues, t = [], 0.0
    cues.append((t, 'BASELINE', 0))
    t += p['baseline_s']
    for k in range(1, p['cycles'] + 1):
        for phase in ('RAMP', 'HOLD', 'FALL', 'REST'):
            cues.append((t, phase, k))
            t += p[phase.lower() + '_s']
    return cues, t


PHASE_HELP = {
    'BASELINE': 'hands off. preload only.',
    'RAMP':     'press, increasing steadily',
    'HOLD':     'hold firm and steady',
    'FALL':     'release steadily',
    'REST':     'hands off completely',
}


# ------------------------------------------------------------------ preflight
def preflight(profile_name, p, total_s):
    """Blocking checklist. Runs BEFORE the port is opened, deliberately."""
    print('=' * 70)
    print(f'CTG CONCURRENT CAPTURE -- profile: {profile_name}')
    print(f'  {p["note"]}')
    print(f'  duration {total_s:.0f} s ({total_s/60:.1f} min), '
          f'{p["cycles"]} contraction cycle(s)')
    if p['cycles'] > 0 and total_s > 0:
        print(f'  rate {p["cycles"] / (total_s/60) * 10:.2f} per 10 min')
    print('=' * 70)
    print('\nPreflight. Every item is a documented failure mode.\n')
    items = [
        'Laptop UNPLUGGED from mains (USB power raises MAX4466 mains '
        'content by ~40 dB)',
        'Speaker looping the fetal PCG WAV, volume 40%, ALREADY PLAYING',
        'Microphone at the run C distance, bare capsule, trimpot untouched',
        'FSR402 under spring-clip preload (unloaded it reads exactly 0)',
        'FSR on a SEPARATE SURFACE from the microphone (C4 confound)',
        'Arduino Serial Monitor CLOSED (it holds the port)',
        'Nothing else will need this laptop for the next '
        f'{total_s/60:.0f} minutes',
    ]
    for i, s in enumerate(items, 1):
        print(f'  {i}. {s}')
    print('\nOnce capture starts, DO NOT type in this window. A blocking '
          'prompt\ncost the toco phase 55.47 s of data.\n')
    ans = input('All checked? Type GO to begin, anything else to abort: ')
    if ans.strip() != 'GO':
        print('aborted, nothing captured')
        sys.exit(0)
    print()


# ------------------------------------------------------------------ serial
def open_and_sync(port, baud):
    print(f'opening {port}...')
    ser = serial.Serial(port, baud, timeout=1)
    time.sleep(0.2)
    ser.reset_input_buffer()

    t0, header, prompted = time.time(), None, False
    while time.time() - t0 < HEADER_TIMEOUT_S:
        raw = ser.readline()
        if not raw:
            if not prompted and time.time() - t0 > 6.0:
                print('  no header yet -- press EN/RST on the ESP32')
                prompted = True
            continue
        line = raw.decode('utf-8', errors='replace').strip()
        if line.startswith('# mrisr_capture'):
            header = line
            break
    if header is None:
        ser.close()
        print('[FAIL] no header. Is the Arduino Serial Monitor holding the port?')
        sys.exit(1)
    print(f'  {header}')

    m_pcg = re.search(r'pcg=(\d+)@(\d+)', header)
    m_uc = re.search(r'uc=(\d+)@(\d+)', header)
    m_dec = re.search(r'decim=(\d+)', header)
    if not (m_pcg and m_uc and m_dec):
        ser.close()
        print('[FAIL] header present but unparseable')
        sys.exit(1)
    pcg_pin, fs = int(m_pcg.group(1)), int(m_pcg.group(2))
    uc_pin, decim = int(m_uc.group(1)), int(m_dec.group(1))
    if (pcg_pin != EXPECT_PCG_PIN or uc_pin != EXPECT_UC_PIN
            or fs != EXPECT_FS or decim != EXPECT_DECIM):
        ser.close()
        print(f'[FAIL] firmware reports pcg={pcg_pin} uc={uc_pin} fs={fs} '
              f'decim={decim}; expected {EXPECT_PCG_PIN}/{EXPECT_UC_PIN}/'
              f'{EXPECT_FS}/{EXPECT_DECIM}')
        print('       the firmware must NOT be modified for this phase')
        sys.exit(1)
    print(f'  [ok] pcg={pcg_pin} uc={uc_pin} fs={fs} decim={decim}')
    ser.readline()                      # consume the format comment
    return ser, header


def capture(ser, seconds, cues):
    """Read both line formats and fire cues, without ever blocking.

    PCG lines start with a digit, UC lines with 'U'. Each fired cue records the
    host elapsed time AND the device timestamp of the most recent PCG sample,
    so buffer lag is measured rather than assumed.
    """
    pt, pv, ut, uv, bad = [], [], [], [], 0
    fired = []
    ci = 0
    ser.reset_input_buffer()
    t0 = time.time()
    nxt_tick = 30.0

    while True:
        el = time.time() - t0
        if el >= seconds:
            break

        # ---- cues first, so a slow readline cannot delay one ----
        while ci < len(cues) and el >= cues[ci][0]:
            ct, phase, cyc = cues[ci]
            dev_us = pt[-1] if pt else 0
            fired.append(dict(scheduled_s=ct, host_s=el,
                              device_us=int(dev_us), phase=phase, cycle=cyc))
            tag = f'cycle {cyc}' if cyc else ''
            print(f'  [{el:6.1f}s] {phase:<8} {tag:<8} '
                  f'{PHASE_HELP[phase]}', flush=True)
            ci += 1

        raw = ser.readline()
        if not raw:
            continue
        line = raw.decode('utf-8', errors='replace').strip()
        if not line or line.startswith('#'):
            continue
        try:
            if line[0] == 'U':
                _, a, b = line.split(',')
                ut.append(int(a))
                uv.append(int(b))
            else:
                a, b = line.split(',')
                pt.append(int(a))
                pv.append(int(b))
        except ValueError:
            bad += 1
            continue

        if el >= nxt_tick:
            uc_now = uv[-1] / float(EXPECT_DECIM) if uv else 0.0
            flag = '  <-- UC AT ZERO' if uc_now == 0.0 else ''
            print(f'      .. {el:5.0f}/{seconds:.0f} s   '
                  f'pcg {len(pt):7d}   uc {len(ut):5d}   '
                  f'uc now {uc_now:7.1f}{flag}', flush=True)
            nxt_tick += 30.0

    return (np.array(pt, dtype=np.int64), np.array(pv, dtype=np.int64),
            np.array(ut, dtype=np.int64), np.array(uv, dtype=np.int64),
            bad, fired, time.time() - t0)


def git_head():
    try:
        return subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', default=PORT)
    ap.add_argument('--profile', default='gold', choices=sorted(PROFILES))
    ap.add_argument('--tag', default=None,
                    help='defaults to the profile name')
    args = ap.parse_args()

    p = PROFILES[args.profile]
    tag = args.tag or args.profile
    cues, total_s = build_schedule(p)

    for d in (CAP_DIR, RES_DIR, PLOT_DIR):
        os.makedirs(d, exist_ok=True)

    print(f'script  : {os.path.abspath(__file__)}')
    head = git_head()
    print(f'git HEAD: {head or "(not a git repo / git unavailable)"}\n')

    preflight(args.profile, p, total_s)

    ser, header = open_and_sync(args.port, BAUD)
    print(f'\ncapturing {total_s:.0f} s. Follow the cues.\n')
    aborted = False
    try:
        pt, pv, ut, uv, bad, fired, wall_s = capture(ser, total_s, cues)
    except KeyboardInterrupt:
        aborted = True
        pt, pv, ut, uv, bad, fired, wall_s = capture.partial
        print('\n[ABORTED] files will still be written and marked aborted.')
    finally:
        ser.close()

    if len(pt) < 1000 or len(ut) < 10:
        print(f'[FAIL] pcg {len(pt)} samples, uc {len(ut)} samples. '
              f'Nothing usable; treat as VOID (operational).')
        sys.exit(1)

    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    stem = f'mrisr_ctg_{tag}_{stamp}'

    # ---------------------------------------------------------- write files
    pcg_path = os.path.join(CAP_DIR, stem + '_pcg500.csv')
    with open(pcg_path, 'w') as fh:
        fh.write(header + '\n')
        fh.write(f'# {stamp} tag={tag} profile={args.profile} '
                 f'validated_protocol={p["validated"]}\n')
        fh.write('t_us,adc\n')
        for a, b in zip(pt, pv):
            fh.write(f'{a},{b}\n')

    # UC in 04_detect.py's exact format: t_s,uc float-valued. The firmware
    # sends the SUM of 125 samples; dividing by 125.0 reproduces the toco
    # phase's float means bit-for-bit. Integer values would silently
    # REACTIVATE uc_detector's flat_run_mask -- same frozen detector,
    # different behaviour, no error.
    uc_path = os.path.join(CAP_DIR, f'contractions_{stamp}_uc4hz.csv')
    uc_mean = uv.astype(np.float64) / float(EXPECT_DECIM)
    uc_t_s = ut.astype(np.float64) / 1e6
    with open(uc_path, 'w') as fh:
        fh.write('t_s,uc\n')
        for a, b in zip(uc_t_s, uc_mean):
            fh.write(f'{a:.6f},{b:.4f}\n')

    cue_path = os.path.join(CAP_DIR, f'cues_{stamp}.csv')
    with open(cue_path, 'w') as fh:
        fh.write('scheduled_s,host_s,device_us,device_s,phase,cycle\n')
        for c in fired:
            fh.write(f'{c["scheduled_s"]:.3f},{c["host_s"]:.3f},'
                     f'{c["device_us"]},{c["device_us"]/1e6:.3f},'
                     f'{c["phase"]},{c["cycle"]}\n')

    print(f'\npcg  -> {pcg_path}  ({len(pt)} samples, '
          f'{os.path.getsize(pcg_path)/1e6:.1f} MB)')
    print(f'uc   -> {uc_path}  ({len(ut)} samples, 04_detect.py format)')
    print(f'cues -> {cue_path}  ({len(fired)} cues)')

    # ------------------------------------------------------- C5 gates
    print('\n' + '=' * 70)
    print('C5 -- CAPTURE INTEGRITY. Objective, signal-independent, prereg 6.')
    print('=' * 70)

    dp = np.diff(pt)
    du = np.diff(ut)
    mono_pcg = bool(np.all(dp > 0))
    mono_uc = bool(np.all(du > 0))

    span_p = (pt[-1] - pt[0]) / 1e6
    span_u = (ut[-1] - ut[0]) / 1e6
    fs_p = (len(pt) - 1) / span_p if span_p > 0 else 0.0
    fs_u = (len(ut) - 1) / span_u if span_u > 0 else 0.0
    fs_p_ok = abs(fs_p - NOMINAL_FS_PCG) <= FS_TOLERANCE_FRAC * NOMINAL_FS_PCG
    fs_u_ok = abs(fs_u - NOMINAL_FS_UC) <= FS_TOLERANCE_FRAC * NOMINAL_FS_UC

    gap = float(np.max(np.diff(uc_t_s))) if len(uc_t_s) > 1 else 0.0
    gap_ok = gap < MAX_UC_GAP_S

    completed = len(fired) == len(cues) and wall_s >= total_s - 1.0

    checks = [
        ('PCG timestamps strictly increasing', mono_pcg,
         f'{int(np.sum(dp <= 0))} violations'),
        ('UC timestamps strictly increasing', mono_uc,
         f'{int(np.sum(du <= 0))} violations'),
        (f'PCG fs within 0.1% of {NOMINAL_FS_PCG:g} Hz', fs_p_ok,
         f'{fs_p:.4f} Hz'),
        (f'UC fs within 0.1% of {NOMINAL_FS_UC:g} Hz', fs_u_ok,
         f'{fs_u:.4f} Hz'),
        (f'largest UC gap under {MAX_UC_GAP_S:g} s', gap_ok,
         f'{gap:.3f} s'),
        ('protocol completed', completed,
         f'{len(fired)}/{len(cues)} cues, {wall_s:.1f}/{total_s:.0f} s'),
    ]
    for label, ok, detail in checks:
        print(f'  [{"ok" if ok else "FAIL"}]  {label:<42} {detail}')

    c5 = all(ok for _, ok, _ in checks)
    print(f'\n  C5: {"PASS" if c5 else "FAIL"}')
    if not c5:
        print('  -> VOID for an operational reason. Record it in notes with')
        print('     the reason, KEEP the files, and repeat. prereg 8.')
    else:
        print('  -> eligible. If this is the FIRST accepted capture of this')
        print('     profile it is the gold capture. Not the best one. prereg 8.')

    # ------------------------------------------------------- diagnostics
    print('\n' + '=' * 70)
    print('DIAGNOSTICS. Reported, never deciding. prereg 9.')
    print('=' * 70)

    n_lines = len(pt) + len(ut) + bad
    rate = bad / n_lines if n_lines else 0.0
    print(f'  malformed lines   : {bad} of {n_lines}  (1 in '
          f'{int(1/rate) if rate else float("inf"):,})')
    print(f'    predicted        : ~{n_lines * PREDICTED_MALFORMED_RATE:.0f} '
          f'at 1 in 20,000')
    print(f'  dt mean/std       : {np.mean(dp):.2f} / {np.std(dp):.2f} us')
    print(f'  gaps > 2x nominal : '
          f'{int(np.sum(dp > 2 * 1e6 / NOMINAL_FS_PCG))}')
    print(f'  pcg:uc line ratio : {len(pt)/len(ut):.2f}  '
          f'(expected {EXPECT_DECIM})')
    print(f'  host vs device    : wall {wall_s:.2f} s, device span '
          f'{span_p:.2f} s, drift {wall_s - span_p:+.2f} s')

    if fired:
        lag = [c['host_s'] - c['device_us'] / 1e6 for c in fired
               if c['device_us'] > 0]
        if lag:
            print(f'  cue buffer lag    : mean {np.mean(lag):+.3f} s, '
                  f'max {np.max(lag):+.3f} s')
            print('    host-clock cues run AHEAD of device time by this much.')
            print('    Cues are navigation aids. Onset comes from the signal.')

    uc_span = float(np.max(uc_mean) - np.min(uc_mean))
    print(f'  UC range          : {np.min(uc_mean):.1f} .. '
          f'{np.max(uc_mean):.1f} counts (span {uc_span:.1f})')
    if np.min(uc_mean) == 0.0:
        print('    [!] UC touched exactly 0. The FSR reads 0 with no preload.')
        print('        Check the spring clip before trusting this capture.')
    print(f'  PCG bias          : {np.mean(pv):.1f} counts, '
          f'std {np.std(pv):.2f}  (run C rms 87.2)')
    print(f'  PCG peak excursion: '
          f'{np.max(np.abs(pv - np.mean(pv))):.0f} counts  '
          f'(level-finding GOOD window 400-900, run C 642)')

    if not p['validated']:
        print(f'\n  [!] profile "{args.profile}" is NOT the validated protocol.')
        print(f'      {p["note"]}')

    # ------------------------------------------------------------- plot
    fig, ax = plt.subplots(4, 1, figsize=(14, 13))
    tp = (pt - pt[0]) / 1e6
    tu = (ut - ut[0]) / 1e6

    ax[0].plot(tp, pv, linewidth=0.2)
    ax[0].set_title(f'{stem} -- PCG (GPIO39, {fs_p:.3f} Hz)')
    ax[0].set_ylabel('ADC')
    ax[0].grid(alpha=0.3)

    ax[1].plot(tu, uc_mean, '-', linewidth=0.9, color='tab:blue')
    for c in fired:
        if c['phase'] in ('RAMP', 'HOLD', 'FALL'):
            t_end = c['host_s'] + p[c['phase'].lower() + '_s']
            ax[1].axvspan(c['host_s'], t_end, alpha=0.10,
                          color='tab:orange', zorder=0)
    ax[1].set_title(f'UC (GPIO35, {fs_u:.4f} Hz, mean of {EXPECT_DECIM}) '
                    f'-- shading = cue windows on the HOST clock, '
                    f'not detections')
    ax[1].set_ylabel('counts')
    ax[1].grid(alpha=0.3)

    ax[2].hist(dp, bins=80)
    ax[2].set_title(f'PCG sample interval: mean {np.mean(dp):.2f} us, '
                    f'std {np.std(dp):.2f} us')
    ax[2].set_xlabel('dt (us)')
    ax[2].set_yscale('log')
    ax[2].grid(alpha=0.3)

    f, ps = welch(pv - np.mean(pv), fs=fs_p, nperseg=min(2500, len(pv)))
    ax[3].semilogy(f, ps, linewidth=0.8)
    ax[3].axvspan(25, 200, alpha=0.08, color='green')
    ax[3].axvline(47.80, color='tab:purple', alpha=0.6, linestyle=':',
                  label='47.80 Hz (measured bench line)')
    for h in (50, 100, 150):
        ax[3].axvline(h, color='r', alpha=0.30, linestyle='--')
    ax[3].set_title('PCG PSD (green = 25-200 Hz pipeline band, '
                    'red = mains harmonics)')
    ax[3].set_xlabel('Hz')
    ax[3].legend(fontsize=8)
    ax[3].grid(alpha=0.3)

    for a in ax[:2]:
        a.set_xlabel('Time (s)')

    plt.tight_layout()
    ppath = os.path.join(PLOT_DIR, f'01_{stem}.png')
    plt.savefig(ppath, dpi=110)
    plt.close()

    # ------------------------------------------------------------- results
    res = {
        'stem': stem, 'stamp': stamp, 'tag': tag,
        'profile': args.profile,
        'validated_protocol': bool(p['validated']),
        'profile_note': p['note'],
        'script': os.path.abspath(__file__),
        'git_head': head,
        'header': header,
        'protocol': dict(p, total_s=total_s, n_cues=len(cues)),
        'files': {'pcg': pcg_path, 'uc': uc_path, 'cues': cue_path,
                  'plot': ppath,
                  'pcg_bytes': os.path.getsize(pcg_path)},
        'counts': {'pcg': int(len(pt)), 'uc': int(len(ut)),
                   'malformed': int(bad), 'total_lines': int(n_lines),
                   'malformed_rate': float(rate)},
        'timing': {'wall_s': float(wall_s),
                   'pcg_span_s': float(span_p), 'uc_span_s': float(span_u),
                   'fs_pcg': float(fs_p), 'fs_uc': float(fs_u),
                   'dt_mean_us': float(np.mean(dp)),
                   'dt_std_us': float(np.std(dp)),
                   'n_gaps_2x': int(np.sum(dp > 2 * 1e6 / NOMINAL_FS_PCG)),
                   'largest_uc_gap_s': float(gap),
                   'host_device_drift_s': float(wall_s - span_p)},
        'signal': {'pcg_mean': float(np.mean(pv)), 'pcg_std': float(np.std(pv)),
                   'pcg_min': float(np.min(pv)), 'pcg_max': float(np.max(pv)),
                   'pcg_p2p': float(np.max(pv) - np.min(pv)),
                   'pcg_peak_excursion': float(np.max(np.abs(pv - np.mean(pv)))),
                   'uc_min': float(np.min(uc_mean)),
                   'uc_max': float(np.max(uc_mean)),
                   'uc_span': uc_span},
        'cues': fired,
        'c5': {'pass': bool(c5),
               'checks': [{'name': n, 'pass': bool(o), 'detail': d}
                          for n, o, d in checks],
               'tolerances': {'fs_tolerance_frac': FS_TOLERANCE_FRAC,
                              'max_uc_gap_s': MAX_UC_GAP_S}},
        'detectors_run': False,
    }
    rpath = os.path.join(RES_DIR, stem + '.json')
    with open(rpath, 'w') as fh:
        json.dump(res, fh, indent=2)

    print(f'\nresults -> {rpath}')
    print(f'plot    -> {ppath}')
    print('\nNo detector has been run. Score with 02_ctg.py, once, '
          'after accepting this capture.')


if __name__ == '__main__':
    main()
