"""Download a PhysioNet database into data/<DATASET>/.

Set DATASET to switch sources with a one-line edit (the established pattern).
For CTU-CHB the full set is 552 records; SUBSET_N caps it, because 30-50 is
plenty to characterise the contraction pipeline and downloads far faster.
Set SUBSET_N = None to fetch everything.
"""
import os

import wfdb

DATASET = 'ctu-uhb-ctgdb'          # CTU-CHB intrapartum CTG database
DL_DIR = os.path.join('data', DATASET)
SUBSET_N = 200                       # None = all records


def record_list(dataset):
    """Ask PhysioNet which records exist, without downloading signal data.

    dl_database needs record NAMES to fetch a subset. get_record_list returns
    them cheaply (a text index, not the .dat files), so we can slice before
    pulling anything heavy.
    """
    recs = wfdb.get_record_list(dataset)
    return recs


if __name__ == '__main__':
    os.makedirs(DL_DIR, exist_ok=True)

    if SUBSET_N is None:
        print(f"Downloading ALL of {DATASET} into {DL_DIR}/ ...")
        wfdb.dl_database(DATASET, dl_dir=DL_DIR)
    else:
        all_recs = record_list(DATASET)
        subset = all_recs[:SUBSET_N]
        print(f"{DATASET}: {len(all_recs)} records available, "
              f"downloading first {len(subset)} into {DL_DIR}/ ...")
        wfdb.dl_database(DATASET, dl_dir=DL_DIR, records=subset)

    files = sorted(os.listdir(DL_DIR))
    print("Done. Files in:", DL_DIR)
    print("Sample files:", files[:10])
    print(f"Total files: {len(files)}")
