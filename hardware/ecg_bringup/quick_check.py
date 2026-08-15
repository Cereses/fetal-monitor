"""
quick_check.py — verify nifecgdb_io reads the database at its native rate.

Runs from anywhere in the repo. Every failure path prints WHY, because a
script that can finish silently on a missing input is worse than one that
crashes: the previous version globbed a relative path, found nothing, and
exited 0 with no output at all.

Checks, in order:
  1. pyedflib imports and reports a version
  2. data/nifecgdb/ is found (walking up from cwd)
  3. .edf files are present, and .edf.qrs alongside them
  4. each record loads at fs = 1000 Hz, NOT the 100 Hz wfdb reports
  5. annotation scale is 1.0 and indices land inside the signal
  6. no kept channel is the 'EDF Annotations' text channel
     (pyedflib never exposes it, so `dropped` being empty is CORRECT)
  7. gestational age parses from the EDF patient field
"""
import os
import sys
import glob
import argparse
import collections

SUBPATH = os.path.join('data', 'nifecgdb')
DEPTH = 5
N_SHOW = 5
EDF_ANN_LABEL = 'EDF Annotations'


def resolve_upward(subpath=SUBPATH, depth=DEPTH):
    here = os.path.abspath(os.getcwd())
    tried = []
    for _ in range(depth + 1):
        cand = os.path.join(here, subpath)
        tried.append(cand)
        if os.path.isdir(cand):
            return cand, tried
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    return None, tried


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=0,
                    help='records to check; 0 = all (default)')
    args = ap.parse_args()
    print("=== quick_check: nifecgdb native-rate loader ===")
    print(f"  cwd: {os.path.abspath(os.getcwd())}\n")

    # --- 1. dependency ---
    try:
        import pyedflib
        print(f"  [ok]   pyedflib {pyedflib.__version__}")
    except ImportError as e:
        print(f"  [FAIL] pyedflib not importable: {e}")
        print("         pip install pyedflib")
        return 1

    try:
        import nifecgdb_io as nio
        print(f"  [ok]   nifecgdb_io imported from {nio.__file__}")
    except ImportError as e:
        print(f"  [FAIL] nifecgdb_io not importable: {e}")
        print("         it must sit in the same directory as this script,")
        print("         or on PYTHONPATH.")
        return 1

    # --- 2. data ---
    data_dir, tried = resolve_upward()
    if data_dir is None:
        print("  [FAIL] could not find data/nifecgdb/. Looked in:")
        for t in tried:
            print(f"           {t}")
        print("         Run download_nifecgdb.py from the repo root.")
        return 1
    print(f"  [ok]   data dir: {data_dir}")

    # --- 3. files ---
    edfs = sorted(glob.glob(os.path.join(data_dir, '*.edf')))
    qrs = sorted(glob.glob(os.path.join(data_dir, '*.edf.qrs')))
    print(f"  [{'ok' if edfs else 'FAIL'}]   {len(edfs)} .edf file(s), "
          f"{len(qrs)} .edf.qrs file(s)")
    if not edfs:
        listing = os.listdir(data_dir)[:10]
        print(f"         directory contains: {listing}")
        return 1
    if len(qrs) != len(edfs):
        print(f"         MISMATCH: records without annotations cannot be")
        print(f"         scored and must be a pre-registered exclusion.")

    # --- 4-7. load ---
    todo = edfs if args.n <= 0 else edfs[:args.n]
    print(f"\n  Loading {len(todo)} of {len(edfs)} record(s):\n")
    print(f"{'record':>16} {'fs':>7} {'samples':>9} {'beats':>6} {'scale':>6} "
          f"{'align':>6} {'gest':>7} {'thor':>4}{'abd':>4}")
    problems = []
    stats = []
    for p in todo:
        name = os.path.basename(p)
        try:
            s = nio.summarise(p)
        except Exception as e:
            print(f"{name:>16}  LOAD FAILED: {type(e).__name__}: {e}")
            problems.append(f"{name}: {e}")
            continue
        g = s['gestation']
        gs = f"{g[0]}+{g[1]}" if g else '-'
        print(f"{name:>16} {s['fs']:>7.0f} {s['n_samples']:>9} "
              f"{s['n_beats']:>6} {s['scale']:>6.2f} "
              f"{str(s['aligned']):>6} {gs:>7} {s['n_thoracic']:>4}"
              f"{s['n_abdominal']:>4}")
        if s['fs'] != 1000.0:
            problems.append(f"{name}: fs is {s['fs']}, expected 1000")
        if abs(s['scale'] - 1.0) > 1e-9:
            problems.append(f"{name}: annotation scale {s['scale']}, expected 1.0")
        if not s['aligned']:
            problems.append(f"{name}: {s['align_msg']}")
        # NOT a problem if nothing was dropped. pyedflib parses the TAL
        # channel as annotations and never exposes it as a signal, so
        # signals_in_file is 5 where wfdb reported 6. What matters is that no
        # KEPT channel is the text channel.
        if EDF_ANN_LABEL in s['sig_name']:
            problems.append(f"{name}: '{EDF_ANN_LABEL}' present among kept "
                            f"channels — it is text, not signal")
        if s['n_thoracic'] < 1 or s['n_abdominal'] < 1:
            problems.append(f"{name}: channel groups incomplete "
                            f"({s['n_thoracic']} thoracic, "
                            f"{s['n_abdominal']} abdominal)")
        if s['gestation'] is None:
            problems.append(f"{name}: gestational age not found in any "
                            f"header field")
        stats.append(s)

    # --- aggregate: this is what the pre-registration is written against ---
    if stats:
        print(f"\n  === Aggregate over {len(stats)} record(s) ===")
        fss = sorted({x['fs'] for x in stats})
        print(f"    sampling rates      : {fss}")
        scales = sorted({round(x['scale'], 6) for x in stats})
        print(f"    annotation scales   : {scales}")
        print(f"    all aligned         : {all(x['aligned'] for x in stats)}")
        gs = [x['gestation'] for x in stats if x['gestation']]
        if gs:
            wk = [g[0] + g[1] / 7.0 for g in gs]
            print(f"    gestation           : {min(wk):.1f}-{max(wk):.1f} weeks "
                  f"({len(gs)}/{len(stats)} parsed)")
        grp = collections.Counter((x['n_thoracic'], x['n_abdominal'])
                                  for x in stats)
        print(f"    (thoracic, abdominal): {dict(grp)}")
        tot_b = sum(x['n_beats'] for x in stats)
        tot_s = sum(x['duration_s'] for x in stats)
        print(f"    total beats         : {tot_b:,}")
        print(f"    total duration      : {tot_s/60:.1f} min "
              f"({tot_s/3600:.2f} h)")
        print(f"    shortest / longest  : {min(x['duration_s'] for x in stats):.0f} s"
              f" / {max(x['duration_s'] for x in stats):.0f} s")
        labs = collections.Counter(tuple(x['sig_name']) for x in stats)
        print(f"    distinct channel sets: {len(labs)}")
        for lab, c in labs.most_common():
            print(f"      {c:>3} x {list(lab)}")

    # --- detail on one record ---
    s = None
    try:
        s = nio.summarise(edfs[0])
    except Exception:
        pass
    if s:
        print(f"\n  Detail — {s['record']}:")
        print(f"    channels   : {s['sig_name']}")
        print(f"    units      : {s['units']}")
        print(f"    thoracic   : {s['n_thoracic']}   abdominal: {s['n_abdominal']}")
        print(f"    duration   : {s['duration_s']:.1f} s")
        print(f"    prefilter  : {s['prefilter']}")
        print(f"    alignment  : {s['align_msg']}")
        gf = s.get('gestation_field')
        print(f"    gestation  : {s['gestation']} (header field: {gf})")
        print(f"    dropped    : {s['dropped'] or 'none — correct, pyedflib '
                                  'parses the TAL channel as annotations'}")

    print()
    if problems:
        print(f"  {len(problems)} PROBLEM(S):")
        for p in problems:
            print(f"    - {p}")
        return 1

    print("  All checks passed. Signals are at native 1000 Hz, annotations")
    print("  index the signal directly, and the EDF+ text channel is excluded.")
    print("  Next: pre-register the thoracic-vs-abdominal comparison.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
