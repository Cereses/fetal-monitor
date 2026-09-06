# §14 — acoustic chamber and surface conduction

Section 14 of the PCG bring-up notes, kept as a separate file alongside
`notes.md` (same arrangement as `ecg_bringup/notes/notes_nifecgdb.md`).
`notes.md` §0, §12 and §13 have been amended to reference it.

Work completed 2026-09-01 and 2026-09-05, after commit `98efe36`.

---

## 14. Acoustic chamber and surface conduction

§0 of the original notes recorded that the phase established nothing about
surface-conducted sound. **That gap is now closed.** This section is the
evidence, and the things it still does not establish.

### 14.1 The chamber

A green PET bottle cap, ~28 mm, with a hole melted through the crown using a
screwdriver tip heated with a lighter. The capsule is pushed through from
outside so its face opens into the cavity; the PCB and all wiring stay outside
the cap, so the wall is crossed at exactly one joint. That joint is sealed with
PVC electrical tape, and the seal was checked by pressing the rim against a
palm and feeling for resistance and release.

The geometry is a **bell**, not a diaphragm: rigid rim against skin, trapped
air, capsule opening into the trapped volume. A bell is the correct choice for
low-frequency sound — a tensioned diaphragm deliberately attenuates the bottom
of the band, which is where fetal S1 energy sits.

Named limitations. The cap has visible flex, and a wall that gives absorbs
displacement instead of converting it to pressure. The seal is tape rather than
a moulded fit. No soldering iron was used for the hole (PET fumes, and the tip
is needed for later work), so the hole is rough and oversize with the tape
doing the sealing.

### 14.2 Body contact plus USB power buries the channel in mains

| | 50 Hz over background | mains share of in-band | std (counts) |
|---|---|---|---|
| bare capsule on chest, USB | +46.8 dB | **97.4%** | 210.4 |
| chamber on chest, USB | +39.8 dB | **94.9%** | 260.9 |
| chamber on chest, **battery** | +5.0 dB, not detected | **1.86%** | 54.7 |

The §7.3 prediction held emphatically. Pressing anything against the chest with
the laptop on charger drives the channel to ~95–97% mains; the first-second
panel of the USB bare capture is a textbook 50 Hz sinusoid. Unplugging the
charger collapsed it by roughly **40 dB**.

This upgrades §7.3 from a characterised interference mode to a characterised
interference mode **with a validated mitigation**. It also made the first two
captures unusable for any chamber comparison — both are saturated, so
bare-versus-cap there would compare two ruined captures. Everything after this
point was captured on battery.

### 14.3 The signal is contact-dependent, not airborne

| | std (counts) |
|---|---|
| quiet-room baselines (§4) | 3.82 / 4.11 / 4.00 |
| **near chest, no contact** | **3.92** |
| bare capsule on chest | 50.1 |
| chamber on chest | 54.7 |

The no-contact control sits at **exactly the bare-ADC noise floor**. Not
slightly above — 3.92 against a 3.92 reference. Airborne room pickup and
breathing sound are ruled out as the source of the ~50-count contact signals.

One weakness, recorded: the control was specified as *the chamber held off the
chest* and executed as *a bare capsule held off the chest* (see §14.7f). It
therefore cannot rule out the chamber acting as a resonant cavity for airborne
sound. It does answer the contact-dependence question, which was the more
important one.

### 14.4 S1/S2 doublet structure

`quick_check_chest.py` found evenly spaced envelope peaks at ~0.68 s in the
bare-capsule contact capture (interval CV 0.05, 100% of windows reliable), with
a visible smaller secondary bump after each. Its peak finder used
`distance=0.4 s`; systole is roughly 0.3 s, so **S2 was suppressed by
construction** and never counted.

`quick_check_doublet.py` lowers `distance` to 0.15 s, with prominence at 10% of
the median coarse prominence — data-derived but scale-free, so the same rule
works on the 50-count contact captures and the 4-count air control.

Pre-registered before the run: D1 peak count 1.6–2.4×; D2 lag-1 autocorrelation
of the interval series ≤ −0.30; D3 intervals bimodal with short fraction
0.25–0.45; D4 cycle rate within 5 BPM of the coarse rate. **The air capture was
the falsifier** — if it also alternated, the test measures nothing.

| | air | **bare on chest** | chamber on chest |
|---|---|---|---|
| D1 count ratio | 2.96× ✗ | **1.85× ✓** | 2.74× ✗ |
| D2 lag-1 autocorr | −0.082 ✗ | **−0.402 ✓** | −0.157 ✗ |
| D3 short fraction | 0.373 ✓ | **0.365 ✓** | 0.344 ✓ |
| D4 cycle vs coarse | Δ29.8 ✗ | **Δ0.8 ✓** | Δ12.3 ✗ |
| verdict | not established | **DOUBLET** | not established |

The falsifier behaved: lag-1 of −0.082 is essentially zero, and its interval
histogram is a single skewed lump.

The bare-capsule interval distribution is two cleanly separated clusters,
**short 0.252 s and long 0.438 s**, cycle 0.690 s, with the sequence forming a
visible sawtooth across 148 consecutive intervals.

**D2 is the load-bearing test.** It needs no threshold on interval length: a
strictly alternating short-long sequence has lag-1 correlation near −1, a
uniform sequence near 0.

**The systolic fraction is the strongest single number here.** 0.365 — systole
occupying just over a third of the cycle, which is the correct physiological
ratio. The 0.25–0.45 range was pre-registered before the run and the value
landed mid-range on the one capture that also passed the other three criteria.
Systole being shorter than diastole is close to diagnostic; little else that is
contact-dependent produces an alternating doublet at a plausible rate with the
right asymmetry.

### 14.5 The rate discrepancy — OPEN, unresolved

Three independent routes on the same capture agree closely:

| route | BPM |
|---|---|
| coarse peak intervals | 87.7 |
| envelope autocorrelation (adult parameters) | 88.2 |
| doublet cycle period (short + long) | 87.0 |

Against **two careful 60-second radial counts of 71 and 70 BPM**, taken later
the same day and agreeing with each other to 1 BPM.

**A 22% gap, and it is not closed.** The capture is internally consistent and
self-consistent across three methods. Its *identification with the subject's
heart rate* is not established. Two possibilities remain open:

- the periodicity is not cardiac, despite the doublet structure and correct
  systolic fraction; or
- the rate genuinely differed between the 13:36 capture and the 16:24 counts,
  which is possible but was not measured.

The counts were taken with the subject conscious of being measured, which
biases toward *lower* rates and therefore does not close a gap in this
direction.

**This is why the §14.4 result is stated as "surface-conducted heart sound
reaches the capsule" and not as a rate measurement.** The obvious resolution —
simultaneous AD8232 capture, using a QRS detector already validated at 99.69%
Se over 95,960 beats as the measuring instrument — was not performed, on time
grounds.

### 14.6 Chamber versus bare capsule — attempted and INVALIDATED

A controlled session was run on 2026-09-05: four captures, alternating
bare/cap/bare/cap to prevent drift aligning with condition, marked placement,
battery power, consistent grip.

**The session produced a hardware fault and no comparison was drawn.**

`ctrl_cap_B` reached **both ADC rails**: min 0, max 4095, p2p 4095. The trace
shows the signal collapsing toward zero for seconds at a time and recovering
instantly, at roughly 0.4, 6, 12, 16–22, 28 and 31 s. Long decays followed by
instant recovery are the signature of an **intermittent connection** — the
input floating, bleeding off, reconnecting — not acoustic content. `ctrl_cap_A`
shows one 4095 spike.

Independent evidence from the same session: DC bias moved 1857.7 → 1911.8 →
1987.0 → 1813.4 across twenty minutes, against 1907–1911 stable across every
prior session.

A `verify_bench` capture immediately afterwards cleared the assembly: bias
1911.1, no values near either rail, zero gaps, zero malformed lines. **The fault
was configuration-dependent, not a persistent wiring defect** — most likely
strain on the taped capsule leads from repeated pressing.

(That verification capture had std 35.3 against a 4-count baseline, with energy
piled below 30 Hz. Almost certainly handling. It clears the rails and the bias,
which is what it was for; it is not a noise-floor measurement and nothing is
indexed to it.)

The four control captures are retained as evidence of the fault. **No
cap-versus-bare conclusion is drawn from them.** The earlier chamber capture is
from 2026-09-01 and the bare from 2026-09-05 — different session, placement and
pressure — so those are two experiments, not a comparison either.

A mechanism was hypothesised and is **not tested**: pressing the metal capsule
can directly into skin is *contact* transduction, coupling vibration
mechanically into the capsule body and bypassing the tissue-to-air impedance
boundary the chamber exists to mitigate, whereas the chamber reintroduces an
air path through a taped and leaky cap.

### 14.7 Errors in this stage

Continuing the §11 numbering. The first is a fourth instance of the same
pattern.

**(e) A pre-registered criterion was not implemented.**
`quick_check_chest.py`'s docstring registered consistency as requiring **both**
a BPM match within ±10 of the counted pulse **and** ≥50% of windows clearing
threshold. The verdict code checked only the BPM. The chamber capture therefore
printed CONSISTENT on **21.2%** of windows — well under the 50% the criterion
demanded. Under the criterion as written, none of the three captures passed.
Fourth instance of §11's pattern: the criterion was right, the implementation
was blind to half of it.

**(f) The air control used the wrong transducer.** Specified as the chamber held
off the chest; executed as a bare capsule held off the chest. Weaker than
designed — it cannot rule out the chamber acting as a resonant cavity — though
it did answer contact-dependence.

**(g) The stage was estimated at "1–2 hours including building the cap".** It
consumed parts of two sessions across five days, produced a hardware fault, and
left its central comparison invalidated. The estimate was wrong by roughly an
order of magnitude and the error was in underestimating the physical work and
the number of controls a defensible claim needs.

---
