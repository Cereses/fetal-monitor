"""
plot_uc_annotations.py — overlay the expert contraction annotations on the raw
UC waveform.

Two jobs:

  1. CLOSE THE ROW-MAPPING QUESTION. The claim that ann row 4 is uterine
     contractions was inferred from event counts, durations and rates, not read
     from documentation. Strong evidence, but circumstantial. If the shaded
     spans land on the visible tocograph bumps, the inference is confirmed by
     direct observation and the figure goes straight into the report.

  2. SHOW WHAT THE DETECTOR HAS TO COPE WITH before a line of it is written:
     how far the resting tone drifts, how large a contraction is relative to
     that drift, and what dropout looks like next to real signal.

Records are chosen from the manifest at runtime, never hardcoded: one heavily
annotated, one typical, one sparse, one with no annotation at all. Hardcoded
record names have bitten this project before.

Run parse_ctu_annotations.py and audit_uc_quality.py first.
"""
import os
import csv

import numpy as np
import matplotlib
matplotlib.use('Agg')            # save only, no window (Python 3.14 windowing crash)
import matplotlib.pyplot as plt
import wfdb

DATASET = 'ctu-uhb-ctgdb'
DATA_DIR = os.path.join('data', DATASET)
RESULTS_DIR = 'results'
PLOT_DIR = 'plots'

FS = 4.0
DROPOUT_MIN_N = int(10.0 * FS)   # same definition as audit_uc_quality.py


def find_channel(rec, *targets):
    cleaned = [n.strip() for n in rec.sig_name]
    for t in targets:
        if t in cleaned:
            return cleaned.index(t)
    return 0


def flat_runs(x):
    if x.size == 0:
        return np.empty(0, int), np.empty(0, int)
    change = np.flatnonzero(np.diff(x) != 0) + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [x.size]))
    return starts, ends - starts


def load_manifest():
    path = os.path.join(RESULTS_DIR, 'ctu_ann_manifest.csv')
    with open(path) as f:
        return {r['record']: int(r['n_contractions']) for r in csv.DictReader(f)}


def load_contractions():
    """record -> list of (start_sample, end_sample)."""
    path = os.path.join(RESULTS_DIR, 'ctu_contractions.csv')
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            out.setdefault(r['record'], []).append(
                (int(r['start_sample']), int(r['end_sample']))
            )
    return out


def pick_records(available, ann_n):
    """One heavy, one typical, one sparse, one unannotated. Chosen, not hardcoded."""
    have = sorted((ann_n.get(r, 0), r) for r in available if ann_n.get(r, 0) > 0)
    none = [r for r in available if ann_n.get(r, 0) == 0]
    picks = []
    if have:
        picks.append((have[-1][1], 'most annotated'))
        picks.append((have[len(have) // 2][1], 'typical'))
        picks.append((have[0][1], 'sparsest'))
    if none:
        picks.append((none[0], 'no annotation'))
    seen, out = set(), []
    for r, label in picks:
        if r not in seen:
            seen.add(r)
            out.append((r, label))
    return out


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    available = sorted(f[:-4] for f in os.listdir(DATA_DIR) if f.endswith('.hea'))
    ann_n = load_manifest()
    spans = load_contractions()

    picks = pick_records(available, ann_n)
    print("Plotting:")
    for r, label in picks:
        print(f"  {r}  ({label}, {ann_n.get(r,0)} contractions)")

    fig, axes = plt.subplots(len(picks), 1, figsize=(15, 3.2 * len(picks)))
    if len(picks) == 1:
        axes = [axes]

    for ax, (name, label) in zip(axes, picks):
        rec = wfdb.rdrecord(os.path.join(DATA_DIR, name))
        uc = rec.p_signal[:, find_channel(rec, 'UC')].astype(float)
        t = np.arange(uc.size) / FS / 60.0        # minutes

        # dropout first, so it sits underneath everything
        v = np.where(np.isfinite(uc), uc, np.nan)
        starts, lens = flat_runs(v[np.isfinite(v)])
        for s, L in zip(starts, lens):
            if L >= DROPOUT_MIN_N:
                ax.axvspan(s / FS / 60, (s + L) / FS / 60,
                           color='0.85', zorder=0)

        for i, (s, e) in enumerate(spans.get(name, [])):
            ax.axvspan(s / FS / 60, e / FS / 60, color='tab:orange',
                       alpha=0.30, zorder=1,
                       label='annotated contraction' if i == 0 else None)

        ax.plot(t, uc, lw=0.8, color='tab:blue', zorder=2)
        ax.set_xlim(t[0], t[-1])
        ax.set_ylabel('UC (nd)')
        ax.set_title(f"record {name} — {label} — "
                     f"{ann_n.get(name,0)} annotated contractions  "
                     f"(grey = flat run >= 10 s)", fontsize=10)
        if spans.get(name):
            ax.legend(loc='upper right', fontsize=8)

    axes[-1].set_xlabel('time (minutes)')
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, 'uc_annotation_overlay.png')
    fig.savefig(out, dpi=130)
    print(f"\nwrote {out}")
    print("\nCheck, in order:")
    print("  1. Do the orange spans sit on visible bumps? That confirms row 4.")
    print("  2. Does the resting tone between contractions drift? That decides")
    print("     whether the detector needs a rolling baseline or a fixed one.")
    print("  3. How large is a contraction against that drift? That sets whether")
    print("     the threshold can be absolute or has to be record-relative.")
    print("  4. Does the unannotated record look genuinely dead, or does it have")
    print("     contractions nobody marked? The second case would matter.")


if __name__ == '__main__':
    main()
