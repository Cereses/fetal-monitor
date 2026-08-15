"""download_nifecgdb.py — fetch the Non-Invasive Fetal ECG Database.

WHY THIS IS SEPARATE FROM download_data.py
The generic script uses wfdb.dl_database(), which infers each record's files
by reading its WFDB header:

    record = rdheader(rec, ...)         # fetches "<rec>.hea"
    all_files.append(rec + ".hea")

NIFECGDB has NO .hea files. It is EDF+, and each record is a single
'<name>.edf' with beat annotations in '<name>.edf.qrs'. dl_database would
request 'ecgca102.edf.hea', get nothing, and fail. So this uses dl_files(),
which takes an explicit file list rather than inferring one.

download_data.py is left alone: it is correct for every WFDB database in the
project (mitdb, ctu-uhb-ctgdb, fpcgdb, bidmc). Breaking it to accommodate one
EDF outlier would be the wrong trade.

WHAT LANDS
    data/nifecgdb/ecgca###.edf        signals: 2 thoracic + 3-4 abdominal,
                                      1 kHz, 16-bit, 0-100 Hz bandpassed
    data/nifecgdb/ecgca###.edf.qrs    beat annotations (maternal or fetal is
                                      NOT established — 03_explore settles it)

55 records, one subject, gestational weeks 21-40. Roughly 330 MB total.

ANNOTATIONS ARE FETCHED TOLERANTLY. It is not documented that every record
carries a .qrs file, so a missing one is reported and skipped rather than
aborting the run. Records without annotations cannot be scored and become a
pre-registered exclusion.

STATUS: written without network access to PhysioNet from the authoring
sandbox, so the dl_files path is untested against the live server. Failures
are caught and named per file.
"""
import os
import argparse

import wfdb

DATASET = 'nifecgdb'
DL_DIR = os.path.join('data', DATASET)
ANN_EXT = 'qrs'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dl-dir', default=DL_DIR)
    ap.add_argument('--limit', type=int, default=None,
                    help='fetch only the first N records (default: all)')
    args = ap.parse_args()

    os.makedirs(args.dl_dir, exist_ok=True)

    print(f"=== {DATASET} ===")
    recs = wfdb.get_record_list(DATASET)
    print(f"  records listed by PhysioNet : {len(recs)}")
    print(f"  first few                   : {recs[:3]}")

    # Record names in this database already end in '.edf'. If they do not,
    # the assumption below about the annotation filename is wrong too, so say
    # so loudly rather than silently building bad paths.
    if not all(r.endswith('.edf') for r in recs):
        print("\n  NOTE: not every record name ends in '.edf'. The annotation")
        print("  filename convention assumed here ('<rec>.qrs') may be wrong.")
        print("  Inspect the RECORDS file before trusting what lands.")

    if args.limit:
        recs = recs[:args.limit]
        print(f"  limited to                  : {len(recs)}")

    # --- signals: one call, since every record must have its .edf ---
    print(f"\n  downloading {len(recs)} signal file(s) into {args.dl_dir}/ ...")
    try:
        wfdb.dl_files(DATASET, args.dl_dir, recs, keep_subdirs=False)
    except Exception as e:
        print(f"  dl_files failed on the signal files: {e}")
        print("  nothing further will work until this is resolved.")
        raise SystemExit(1)

    # --- annotations: individually, so one missing file does not abort ---
    print(f"  downloading annotation file(s) ...")
    got, missing = [], []
    for r in recs:
        f = f'{r}.{ANN_EXT}'
        try:
            wfdb.dl_files(DATASET, args.dl_dir, [f], keep_subdirs=False)
            got.append(f)
        except Exception:
            missing.append(f)

    files = sorted(os.listdir(args.dl_dir))
    n_edf = sum(1 for f in files if f.endswith('.edf'))
    n_ann = sum(1 for f in files if f.endswith('.' + ANN_EXT))
    size_mb = sum(os.path.getsize(os.path.join(args.dl_dir, f))
                  for f in files) / 1e6

    print(f"\n=== Done ===")
    print(f"  signal files (.edf)      : {n_edf}")
    print(f"  annotation files (.{ANN_EXT})  : {n_ann}")
    print(f"  total size               : {size_mb:.0f} MB")
    print(f"  location                 : {os.path.abspath(args.dl_dir)}")
    if missing:
        print(f"\n  {len(missing)} record(s) have NO annotation file:")
        for f in missing[:10]:
            print(f"    {f}")
        if len(missing) > 10:
            print(f"    ... and {len(missing)-10} more")
        print("  These cannot be scored. List them as a PRE-REGISTERED")
        print("  exclusion before running any detector.")

    print("\n  NEXT: python hardware/ecg_bringup/03_explore_nifecgdb.py")
    print("  That establishes whether the annotations are MATERNAL or FETAL.")
    print("  Nothing downstream is decidable until that number is known.")


if __name__ == '__main__':
    main()
