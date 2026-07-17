"""
explore_mitdb.py — discover MIT-BIH structure before building the ECG pipeline.
Unlike the PPG/PCG datasets, MIT-BIH ships expert per-beat annotations (.atr):
exact R-peak positions AND beat-type labels. That gives true ground truth for
both QRS detection (positions) and beat classification (labels).
"""
import os
import collections

import wfdb

DATA_DIR = os.path.join('data', 'mitdb')


def find_channel(rec, *targets):
    """Index of first matching lead by name (e.g. 'MLII'); fall back to 0.
    MIT-BIH is usually MLII + V1, but a few records substitute another lead,
    so match by name rather than assuming position."""
    cleaned = [n.strip() for n in rec.sig_name]
    for t in targets:
        if t in cleaned:
            return cleaned.index(t)
    return 0


def list_records(data_dir):
    return sorted(f[:-4] for f in os.listdir(data_dir) if f.endswith('.hea'))


def main():
    recs = list_records(DATA_DIR)
    if not recs:
        print(f"No records in {DATA_DIR}/. Set DATASET='mitdb' in download_data.py and run it.")
        return
    print(f"Records found : {len(recs)}")
    print(f"First 5       : {recs[:5]}")

    # --- detailed look at one record ---
    name = recs[0]
    rec = wfdb.rdrecord(os.path.join(DATA_DIR, name))
    print(f"\n=== Signal record: {name} ===")
    print(f"  fs       : {rec.fs} Hz")
    print(f"  duration : {rec.sig_len / rec.fs / 60:.1f} min")
    print(f"  channels : {rec.sig_name}")
    print(f"  units    : {rec.units}")
    mlii = find_channel(rec, 'MLII', 'II')
    print(f"  MLII (or fallback) at channel index {mlii}")

    # --- annotations for that record ---
    try:
        ann = wfdb.rdann(os.path.join(DATA_DIR, name), 'atr')
        print(f"\n=== Annotations: {name}.atr ===")
        print(f"  total annotations : {len(ann.sample)}")
        print(f"  first 5 positions : {ann.sample[:5].tolist()}")
        print(f"  symbol counts     : {dict(collections.Counter(ann.symbol))}")
    except Exception as e:
        print(f"\n  WARNING: could not read {name}.atr annotations: {e}")
        print("  If missing, the .atr files may need fetching separately (like BIDMC numerics).")
        return

    # --- global beat-type distribution across ALL records ---
    # This surfaces the class-imbalance reality up front (central to E2).
    print(f"\n=== Global symbol distribution (all {len(recs)} records) ===")
    global_counts = collections.Counter()
    for r in recs:
        try:
            global_counts.update(wfdb.rdann(os.path.join(DATA_DIR, r), 'atr').symbol)
        except Exception:
            pass
    total = sum(global_counts.values())
    for sym, cnt in global_counts.most_common():
        print(f"  {sym:>3} : {cnt:>7}  ({100*cnt/total:5.2f}%)")
    print(f"  total: {total}")
    print("\n  Note: not all symbols are beats — e.g. '+' (rhythm change), '~'")
    print("  (signal quality), '|' are non-beat annotations we'll filter out later.")


if __name__ == '__main__':
    main()