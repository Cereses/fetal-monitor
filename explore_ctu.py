"""
explore_ctu.py — discover CTU-CHB structure before building the FSR402
contraction pipeline. Same discipline as explore_mitdb.py: look before you build.

Unlike MIT-BIH, CTU-CHB ships NO event annotations on PhysioNet. The uterine
contraction ground truth lives in a separate Data in Brief supplement
(Romagnoli et al. 2020, DOI 10.1016/j.dib.2020.105690), downloaded by hand.
So this script characterises the SIGNAL only; annotation parsing comes later.

Questions this needs to answer, because each shapes the pipeline:
  - Are the channels really [FHR, UC] at 4 Hz, as the headers suggested?
  - What is the UC signal's range? It is 'nd' (no dimension): a relative
    tocograph value, exactly like the FSR402 will produce. Not calibrated mmHg.
  - How much of each UC channel is FLAT or ZERO? At least one published group
    dropped the UC channel for quality reasons, so the usable-record count
    matters before we commit to a subset.
"""
import os
import collections

import numpy as np
import wfdb

DATASET = 'ctu-uhb-ctgdb'
DATA_DIR = os.path.join('data', DATASET)

# A UC sample counts as "flat" if it barely moves from its neighbour. At 4 Hz a
# real contraction ramps over tens of seconds, so genuine signal always has
# local variation; long flat runs mean sensor dropout or a lead-off stretch.
FLAT_EPS = 1e-6


def find_channel(rec, *targets):
    """Index of first matching channel by name; fall back to 0.
    Match by NAME because assuming position is the exact mistake MIT-BIH's
    record 114 punished. Do not assume UC is channel 1."""
    cleaned = [n.strip() for n in rec.sig_name]
    for t in targets:
        if t in cleaned:
            return cleaned.index(t)
    return 0


def list_records(data_dir):
    return sorted(f[:-4] for f in os.listdir(data_dir) if f.endswith('.hea'))


def flat_fraction(x):
    """Fraction of samples that are essentially unchanged from the previous one.
    A proxy for dropout: long constant runs are sensor loss, not physiology."""
    if x.size < 2:
        return 1.0
    d = np.abs(np.diff(x))
    return float((d < FLAT_EPS).mean())


def zero_fraction(x):
    """Fraction of exactly-zero samples: a common dropout sentinel in CTG."""
    return float((x == 0).mean()) if x.size else 1.0


def main():
    recs = list_records(DATA_DIR)
    if not recs:
        print(f"No records in {DATA_DIR}/. Set DATASET='{DATASET}' in "
              "download_data.py and run it.")
        return
    print(f"Records found : {len(recs)}")
    print(f"First 5       : {recs[:5]}")

    # --- detailed look at one record ---
    name = recs[0]
    rec = wfdb.rdrecord(os.path.join(DATA_DIR, name))
    print(f"\n=== Signal record: {name} ===")
    print(f"  fs        : {rec.fs} Hz")
    print(f"  duration  : {rec.sig_len / rec.fs / 60:.1f} min")
    print(f"  channels  : {rec.sig_name}")
    print(f"  units     : {rec.units}")
    uc = find_channel(rec, 'UC')
    fhr = find_channel(rec, 'FHR')
    print(f"  FHR at index {fhr}, UC at index {uc}")

    uc_sig = rec.p_signal[:, uc]
    print(f"\n  UC signal (record {name}):")
    print(f"    min / max     : {np.nanmin(uc_sig):.2f} / {np.nanmax(uc_sig):.2f}")
    print(f"    mean / std    : {np.nanmean(uc_sig):.2f} / {np.nanstd(uc_sig):.2f}")
    print(f"    flat fraction : {100*flat_fraction(uc_sig):.1f}%")
    print(f"    zero fraction : {100*zero_fraction(uc_sig):.1f}%")
    print(f"    NaN fraction  : {100*np.isnan(uc_sig).mean():.1f}%")

    # --- channel-layout audit across all records ---
    # Confirm every record is [FHR, UC] at 4 Hz before trusting a fixed index.
    print(f"\n=== Channel-layout audit (all {len(recs)} records) ===")
    layouts = collections.Counter()
    fs_seen = collections.Counter()
    uc_not_1 = []
    for r in recs:
        hea = wfdb.rdheader(os.path.join(DATA_DIR, r))
        names = tuple(n.strip() for n in hea.sig_name)
        layouts[names] += 1
        fs_seen[hea.fs] += 1
        idx = find_channel(hea, 'UC')
        if idx != 1:
            uc_not_1.append((r, names))
    for layout, cnt in layouts.most_common():
        print(f"  {layout} : {cnt}")
    print(f"  sampling rates seen: {dict(fs_seen)}")
    print(f"  records where UC is NOT index 1: "
          f"{uc_not_1 if uc_not_1 else 'none'}")

    # --- UC quality distribution across all records ---
    # This is the decision-driver: how many records carry a usable UC channel?
    print(f"\n=== UC quality across all {len(recs)} records ===")
    stats = []
    for r in recs:
        rr = wfdb.rdrecord(os.path.join(DATA_DIR, r))
        idx = find_channel(rr, 'UC')
        x = rr.p_signal[:, idx]
        stats.append({
            'record': r,
            'flat': flat_fraction(x),
            'zero': zero_fraction(x),
            'nan': float(np.isnan(x).mean()),
            'std': float(np.nanstd(x)),
        })

    flats = np.array([s['flat'] for s in stats])
    zeros = np.array([s['zero'] for s in stats])
    nans = np.array([s['nan'] for s in stats])

    print(f"  flat fraction : median {100*np.median(flats):.1f}%  "
          f"max {100*flats.max():.1f}%")
    print(f"  zero fraction : median {100*np.median(zeros):.1f}%  "
          f"max {100*zeros.max():.1f}%")
    print(f"  NaN  fraction : median {100*np.median(nans):.1f}%  "
          f"max {100*nans.max():.1f}%")

    # A rough usability gate: too much flat/zero/NaN means no scoreable signal.
    bad = [s['record'] for s in stats
           if s['flat'] > 0.5 or s['zero'] > 0.5 or s['nan'] > 0.5]
    print(f"\n  records >50% flat/zero/NaN (likely unusable): "
          f"{len(bad)} -> {bad[:10]}{' ...' if len(bad) > 10 else ''}")
    print(f"  usable records (rough)                       : "
          f"{len(recs) - len(bad)} / {len(recs)}")

    print("\n  Next: obtain the Romagnoli 2020 UC annotations (Data in Brief")
    print("  supplement) so contractions can be scored by sensitivity/")
    print("  specificity. Until then the pipeline is built and checked on the")
    print("  signal's physiological plausibility, as with fpcgdb.")


if __name__ == '__main__':
    main()
