"""
parse_ctu_annotations.py — turn the Romagnoli et al. (2020) annotation .mat
files into two tidy CSVs the contraction pipeline can score against.

Source: Data in Brief 31 (2020) 105690, DOI 10.1016/j.dib.2020.105690,
supplementary file mmc2.zip -> 552 files named annotation_<record>.mat.

WHY THIS SCRIPT EXISTS
----------------------
The .mat files are not a friendly format. Each holds two variables:

  dataloss : 1x3 double
  ann      : 5xN cell array, N = number of samples in the record (4 Hz)

Annotations are stored as MARKER POSITIONS, not per-sample labels: cell (r, c)
is non-empty when an event of type r starts or ends at sample c. The contents
of those cells are MATLAB `string` objects, which serialise as MCOS object
references into the file's subsystem. scipy cannot resolve them, and we do not
need to: the ROW INDEX already encodes the event type, and the marker POSITION
is the ground truth we are after.

ROW -> EVENT MAPPING
--------------------
The paper does not state the row order explicitly, so it was inferred and then
verified against four independent lines of evidence. Documented here because
this inference has to survive a viva:

  row 0  bradycardia    10/552 records, median episode 600 s
  row 1  tachycardia    11/552 records, median episode 660 s
  row 2  acceleration   ~166 events
  row 3  deceleration   ~490 events
  row 4  UTERINE CONTRACTION  <- the one we need

  1. The paper lists exactly five annotated event types, in this order:
     bradycardia, tachycardia, acceleration, deceleration (FHR), then uterine
     contraction (tocogram). Five types, five rows, same order.
  2. Row 4 durations: median 60.0 s, IQR 46.5-76.5 s. Textbook contraction.
  3. Row 4 start-to-start interval: median 141 s = 4.3 per 10 min. Textbook
     active-labour frequency.
  4. Published CTU-CHB event proportions are ~12.2% deceleration vs ~4.4%
     acceleration, a ratio near 2.8:1. Observed rows 3:2 = 490:166 = 2.95:1.

  Rows 0 and 1 carry episodes an order of magnitude longer than rows 2-4,
  consistent with bradycardia and tachycardia being sustained baseline shifts
  rather than transient events.

  REMAINING CHECK (do this once the signals are downloaded): overlay the row-4
  spans on the actual UC waveform for two or three records. If the spans sit on
  the visible tocograph bumps, the mapping is confirmed beyond argument. That
  visual check is worth more than any string label would have been.

PAIRING
-------
Markers come in start/end pairs, so a clean row has an even count. 12 of 552
records have an odd count. In every one of those 12 the unpaired marker falls
97.7-99.6% of the way through the record, 21-90 s from the end: a contraction
that began before the recording stopped at delivery. Those trailing starts are
dropped and the record is kept, with the drop recorded in the manifest.

PRE-REGISTERED EXCLUSIONS
-------------------------
Written to the manifest BEFORE any detector runs, so the exclusion set cannot
be tuned to flatter a result later. Same discipline as the mitdb 102/104 lead
exclusions and the record 207 flutter audit.

  no_annotation   79 records carry zero contraction markers. Note that this is
                  NOT explained by data loss alone: record 1069 has very low
                  loss (1.5 / 0.1 / 5.7) and still has no markers, so the
                  absence reflects the annotators' own inclusion decisions.
                  Excluded because absence of markers cannot be distinguished
                  from absence of contractions, and scoring against it would
                  manufacture false positives.

  Everything else is scorable: 473 records, ~6,750 annotated contractions.

DATALOSS
--------
Three unlabelled columns. Column 0 is close to the figure the paper quotes for
record 1016 (38.97 here vs 38.5 quoted), so it is probably FHR-related, but the
mapping of all three is NOT established. Do not cite these as UC quality until
verified: once explore_ctu.py has run, correlate each column against the
measured flat/zero fraction of the UC channel per record. Whichever column
tracks it is the UC column. Carried through to the manifest untouched so that
check is possible.
"""
import os
import glob

import numpy as np
import scipy.io as sio

# --- config -----------------------------------------------------------------
DATASET = 'ctu-uhb-ctgdb'
ANN_DIR = os.path.join('data', 'ctu-annotations')   # unzip mmc2.zip here
RESULTS_DIR = 'results'

FS = 4.0            # Hz, fixed for CTU-CHB
ROW_UC = 4          # uterine contraction row; see ROW -> EVENT MAPPING above

# Sanity band for reporting only. NOT a filter: annotations are the ground
# truth and we do not get to overrule the gynecologist. Contractions outside
# this band are counted and reported so the number is known, not hidden.
PLAUSIBLE_MIN_S = 20.0
PLAUSIBLE_MAX_S = 120.0


def marker_positions(ann, row):
    """Sample indices where row `row` has a non-empty cell.

    Non-empty is the only test that matters. The cell CONTENTS are MCOS string
    references that scipy renders as opaque structs; their presence is the
    signal, their value is not needed.
    """
    return np.array(
        [c for c in range(ann.shape[1]) if np.asarray(ann[row, c]).size > 0],
        dtype=int,
    )


def record_id(path):
    """annotation_1005.mat -> '1005'. Matches the wfdb record name."""
    return os.path.basename(path)[len('annotation_'):-len('.mat')]


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    files = sorted(glob.glob(os.path.join(ANN_DIR, 'annotation_*.mat')))
    if not files:
        raise SystemExit(
            f"No annotation files in {ANN_DIR}/.\n"
            f"Unzip the Data in Brief supplement (mmc2.zip) into that folder."
        )
    print(f"Found {len(files)} annotation files in {ANN_DIR}/")

    contractions = []   # one row per contraction
    manifest = []       # one row per record

    for path in files:
        rid = record_id(path)
        m = sio.loadmat(path)
        ann = m['ann']
        n_samples = ann.shape[1]
        dl = m['dataloss'].ravel()

        pos = marker_positions(ann, ROW_UC)

        truncated = 0
        if pos.size % 2 == 1:
            # Verify the assumption rather than trusting it blindly: the
            # unpaired marker must be near the end for the truncation story to
            # hold. If it is not, this record needs a human look, so say so.
            frac = pos[-1] / n_samples
            if frac < 0.95:
                print(f"  WARNING {rid}: odd marker count but last marker at "
                      f"{100*frac:.1f}% of record, not near the end. Excluded.")
                manifest.append(dict(
                    record=rid, n_samples=n_samples,
                    duration_min=n_samples / FS / 60,
                    n_contractions=0, n_truncated=0,
                    dataloss_0=dl[0], dataloss_1=dl[1], dataloss_2=dl[2],
                    excluded=True, exclusion_reason='unpaired_marker_midrecord',
                ))
                continue
            pos = pos[:-1]
            truncated = 1

        starts, ends = pos[0::2], pos[1::2]

        if starts.size == 0:
            manifest.append(dict(
                record=rid, n_samples=n_samples,
                duration_min=n_samples / FS / 60,
                n_contractions=0, n_truncated=0,
                dataloss_0=dl[0], dataloss_1=dl[1], dataloss_2=dl[2],
                excluded=True, exclusion_reason='no_annotation',
            ))
            continue

        # Ordering guard. If any end precedes its start, the pairing assumption
        # is wrong for this record and every downstream metric would be garbage.
        if np.any(ends <= starts):
            bad = int((ends <= starts).sum())
            print(f"  WARNING {rid}: {bad} pair(s) with end <= start. Excluded.")
            manifest.append(dict(
                record=rid, n_samples=n_samples,
                duration_min=n_samples / FS / 60,
                n_contractions=0, n_truncated=0,
                dataloss_0=dl[0], dataloss_1=dl[1], dataloss_2=dl[2],
                excluded=True, exclusion_reason='pairing_violation',
            ))
            continue

        for i, (s, e) in enumerate(zip(starts, ends)):
            contractions.append(dict(
                record=rid, index=i,
                start_sample=int(s), end_sample=int(e),
                start_s=s / FS, end_s=e / FS,
                duration_s=(e - s) / FS,
            ))

        manifest.append(dict(
            record=rid, n_samples=n_samples,
            duration_min=n_samples / FS / 60,
            n_contractions=int(starts.size), n_truncated=truncated,
            dataloss_0=dl[0], dataloss_1=dl[1], dataloss_2=dl[2],
            excluded=False, exclusion_reason='',
        ))

    # --- write ---------------------------------------------------------------
    import csv

    c_path = os.path.join(RESULTS_DIR, 'ctu_contractions.csv')
    with open(c_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(contractions[0].keys()))
        w.writeheader()
        w.writerows(contractions)

    m_path = os.path.join(RESULTS_DIR, 'ctu_ann_manifest.csv')
    with open(m_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(manifest[0].keys()))
        w.writeheader()
        w.writerows(manifest)

    # --- report --------------------------------------------------------------
    dur = np.array([c['duration_s'] for c in contractions])
    kept = [r for r in manifest if not r['excluded']]
    excl = [r for r in manifest if r['excluded']]
    n_trunc = sum(r['n_truncated'] for r in kept)

    print(f"\n{'='*62}\nCTU-CHB CONTRACTION ANNOTATIONS\n{'='*62}")
    print(f"  records total          : {len(manifest)}")
    print(f"  records scorable       : {len(kept)}")
    print(f"  records excluded       : {len(excl)}")
    from collections import Counter
    for reason, n in Counter(r['exclusion_reason'] for r in excl).items():
        print(f"      {reason:<28}: {n}")
    print(f"  truncated final contraction dropped in {n_trunc} record(s)")

    print(f"\n  contractions annotated : {len(contractions)}")
    print(f"  per record             : median "
          f"{np.median([r['n_contractions'] for r in kept]):.0f}")

    print(f"\n  duration (s)  p5 {np.percentile(dur,5):6.1f} | "
          f"p25 {np.percentile(dur,25):6.1f} | med {np.median(dur):6.1f} | "
          f"p75 {np.percentile(dur,75):6.1f} | p95 {np.percentile(dur,95):6.1f}")
    out = (dur < PLAUSIBLE_MIN_S) | (dur > PLAUSIBLE_MAX_S)
    print(f"  outside {PLAUSIBLE_MIN_S:.0f}-{PLAUSIBLE_MAX_S:.0f}s "
          f": {out.sum()} ({100*out.mean():.2f}%)  [reported, not filtered]")

    print(f"\n  wrote {c_path}")
    print(f"  wrote {m_path}")


if __name__ == '__main__':
    main()
