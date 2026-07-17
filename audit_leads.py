"""
audit_leads.py — confirm which lead find_channel picks for every MIT-BIH record
before E1. MIT-BIH's first lead is usually MLII, but not always: a few records
lack it entirely and at least one stores the leads reversed. Checking all 48 up
front means Pan-Tompkins never runs silently on the wrong lead. Uses rdheader
(metadata only), so no signal is loaded — it's instant.
"""
import os
import wfdb

DATA_DIR = os.path.join('data', 'mitdb')


def find_channel(rec, *targets):
    """Index of first matching lead by name; falls back to 0. Mirrors explore_mitdb.py."""
    cleaned = [n.strip() for n in rec.sig_name]
    for t in targets:
        if t in cleaned:
            return cleaned.index(t)
    return 0


def list_records(data_dir):
    return sorted(f[:-4] for f in os.listdir(data_dir) if f.endswith('.hea'))


def main():
    recs = list_records(DATA_DIR)
    no_mlii = []
    for name in recs:
        hea = wfdb.rdheader(os.path.join(DATA_DIR, name))
        leads = [n.strip() for n in hea.sig_name]
        idx = find_channel(hea, 'MLII', 'II')
        has = ('MLII' in leads) or ('II' in leads)
        if not has:
            no_mlii.append(name)
        flag = '' if has else '   <-- NO MLII, using index 0 instead'
        print(f"  {name}:  {str(leads):24s} -> idx {idx} ({leads[idx]}){flag}")

    print(f"\nRecords total        : {len(recs)}")
    print(f"Records without MLII : {no_mlii if no_mlii else 'none'}")


if __name__ == '__main__':
    main()