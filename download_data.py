"""Download the simulated fetal PCG database from PhysioNet."""
import os
import wfdb

DATASET = 'mitdb'
DL_DIR = os.path.join('data', DATASET)

if __name__ == '__main__':
    os.makedirs(DL_DIR, exist_ok=True)
    print(f"Downloading {DATASET} into {DL_DIR}/ ...")
    wfdb.dl_database(DATASET, dl_dir=DL_DIR)
    print("Done. Files in:", DL_DIR)
    print("Sample files:", sorted(os.listdir(DL_DIR))[:10])
