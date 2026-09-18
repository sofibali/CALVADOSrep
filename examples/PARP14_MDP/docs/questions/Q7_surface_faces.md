# Q7 — Which residues are exposed, and which face of each domain carries the active site?

**Status:** ANSWERED (2026-09-17). Three-face split + contact analysis; two
earlier numbers corrected in the process.

---

## The question

Two linked questions per domain:

1. Which residues are solvent-exposed vs buried?
2. For MD1 / MD2 / MD3 / ART, which residues sit on the **active-site face**
   versus the **back face** of the domain?

## Why it matters

Face assignment is what makes a binder campaign targetable — you want hotspots
on the face that competes with substrate, not on the opposite side ([[Q8]]). It
also gates the lysine analysis ([[Q6]]), which is restricted to exposed residues.

## How it was computed

`sim_analysis/figure_sasa_faces.py:4-20`:

- SASA from the AF2/AF3 structure, RSA normalised Chothia Gly-X-Gly (per Wu 2017)
- burial bins: buried < 30 Å², surface > 100 Å²
- **face assignment** by the sign of the dot product of
  `(residue_CA − domain_COM)` with `(active_site_COM − domain_COM)` —
  positive = active-site face, negative = back face

## What exists

`data/sasa_face_per_residue.csv` — **all 1801 residues**, columns:
`resid, resname, sasa_A2, max_asa_A2, rsa, burial, domain, face,
face_projection_A, is_catalytic_in, is_pocket_in`.

This is a complete, usable per-residue table. It is purely tabular — **no
summary text accompanies it.**

## It is already load-bearing

Even without a written conclusion, this table is a dependency:

- `bindcraft_md/README.md:73-75` — "Block hotspots = catalytic residues … +
  pocket residues that are **surface-exposed (RSA ≥ 0.15**, `data/sasa_face_per_residue.csv`)
  on the **active-site face**." So the entire binder hotspot definition ([[Q8]])
  rests on this file.
- `figure_lys_exposed_persistence.py` filters lysines at RSA ≥ 0.20 ([[Q6]]).

## Relationship to Q3

[[Q3]] answers a *different* accessibility question and should not be confused
with this one:

| | Q7 (this file) | [[Q3]] |
|---|---|---|
| quantity | per-residue RSA on a **static** structure | per-**pocket** solid-angle accessibility across the **MD ensemble** |
| asks | is this residue's surface exposed? | can a substrate-sized object **approach** this pocket? |
| source | AF2/AF3 model | 19 constructs × 215 replicates |

A residue can be solvent-exposed (high RSA) while its pocket is still
approach-blocked by a neighbouring domain — which is exactly what Q3 found for
MD1.

## The three-face split

A strict `projection > 0` test forced every residue onto one of two poles, which
over-claims: near the equator the **sign** of a near-zero projection carries no
geometric information. The classifier now emits three faces:

| face | meaning |
|---|---|
| `active` | the active-site pole |
| `rim` | the equatorial band — the lateral surface a domain presents to its neighbours while leaving both poles free |
| `back` | the opposite pole |

Split on the projection **normalised by the domain's own half-extent**, because
a fixed-Ångström band would be most of KH7a and a sliver of KH1-6.
`RIM_FRACTION = 0.35`, wide on purpose. Result: active ~27%, rim ~45-50%,
back ~25%.

### Two corrections this forced

**MD1L1's L1 linker was skewing the geometry.** Using `max|proj|` as the scale,
the protruding linker put **62%** of MD1L1 in the rim against ~35% elsewhere.
Fixed by taking the **95th percentile** instead of the max, and by computing the
face axis on the **structured core (794-978)** rather than the full MD1L1 range.
This is the same linker contamination that produced the retracted MD1-MD2
contact claim in [[Q5]].

That flipped a result: MD1L1's active face is **more** exposed than its back
(52.2 vs 50.2 Å²). The contaminated version reported the opposite
(56.9 vs 99.7) — it was measuring the linker, not the domain.

## Which face makes the inter-domain contacts?

`figure_face_contacts.py` counts per-residue inter-domain contacts (CA-CA
< 1.0 nm) across the MD ensemble and partitions them by face. Reported as
**enrichment**, since the rim holds ~half the residues and would dominate by
size alone:

```
enrichment(face) = (contacts on face / all contacts) ÷ (residues in face / all residues)
```

### Folded domains = the restrained residues

"Folded domain" here means **the residues the simulation actually restrains**
(`fl_optimized/input/domains.yaml`), not the full sequence unit. Those residues
are held rigid by harmonic restraints; everything else is free to flex, and a
flexible linker brushing past a domain is not a domain-domain interface.

| domain | full unit | FL_DOMAINS core | **restrained** |
|---|---|---|---|
| MD1L1 | 790–1004 | 794–978 | **800–968** (169) |
| MD2 | 1005–1193 | 1003–1190 | **1015–1183** (169) |
| MD3 | 1207–1388 | 1216–1384 | **1217–1378** (162) |
| ART | 1603–1801 | 1605–1801 | **1613–1791** (179) |
| WWE | 1534–1602 | 1537–1601 | **1549–1587** (39) |

WWE is the striking one — 39 restrained residues against 65 in the looser core,
reflecting its trim of 15 from the [[Q1]] sweep.

This tightened the face split to active ~48-51, rim ~71-83, back ~40-47 across
all four domains — far more uniform than either looser definition.

### The answer, across 6 constructs

| domain | pattern | enrichment range | consistency |
|---|---|---|---|
| **ART** | **back face** | 1.31 – 1.74 | 5/5 constructs |
| **MD3** | **back face** | 1.25 – 1.66 | 6/6 |
| MD1L1 | context-dependent | active 0.80 – 1.84 | 4/6 active |
| MD2 | weak / inconsistent | active 1.16 – 1.53 | 4/6 active |
| **rim** | **depleted everywhere** | 0.50 – 1.00 | 24/24 |

**ART is the clean result.** It consistently presents its *back* face, with its
catalytic face at or below parity (0.77 – 1.03). That is a structural mechanism
for the [[Q4]] finding that ART is the more available catalytic site: the
architecture keeps the writer's pocket face free.

**MD3 is equally consistent** in the same direction, and stronger than the
looser definition showed (1.17–1.30 → 1.25–1.66).

**The rim is depleted in all 24 domain×construct combinations** (0.50 – 1.00).
This is the opposite of the hypothesis that motivated naming it: the equatorial
band is *not* the preferred interface. Domains contact each other pole-to-pole,
not side-on.

### What restricting to restrained residues changed

Sharpening the domain definition strengthened the real signals and weakened a
marginal one — which is what a better definition should do:

| | looser (FL_DOMAINS core) | restrained only |
|---|---|---|
| ART back | 1.40 – 1.97 | 1.31 – 1.74 (all 5) |
| MD3 back | 1.17 – 1.30 (6/6) | **1.25 – 1.66 (6/6)** |
| rim depletion | 0.59 – 0.99 | **0.50 – 1.00, all 24** |
| MD2 active | 1.07 – 1.17 (5/6) | **flipped in 2 constructs (4/6)** |

**MD2's active-face preference does not survive.** It was partly carried by
flexible edge residues; on restrained residues only, `fl` and `mka_full` prefer
the back face. The earlier "MD2 packs through its pocket face" reading should be
treated as unsupported.

### MD1's context dependence — resolved, via a hypothesis that was refuted and then reinstated

Worth recording the arc, because the first answer was wrong for an instructive
reason.

**The observation.** MD1 prefers its *active* face in `fl` / `core` / `norrm` /
`kh1_art_full` and its *back* face in `mka_full` / `md_full` — a split that
tracked KH7a presence exactly (4/4 vs 2/2). With [[Q4]] showing that adding KH7a
drops MD1 accessibility, that looked like a mechanism.

**First test said no.** A per-partner breakdown on `core` gave KH7a only
**11.6%** of MD1's active-face contacts, behind MD3 (33%), ART (25%) and WWE
(17%). Recorded as refuted.

**That was a single-construct anecdote, and `core` is the worst construct to
test it on.** Running the breakdown across all six:

| construct | KH7a share of MD1 active-face contacts | preference |
|---|---|---|
| `fl` | **54%** | active |
| `kh1_art_full` | **36%** | active |
| `core` | 16% | active |
| `norrm` | 16% | active |
| `mka_full` | 0% (KH7a absent) | back |
| `md_full` | 0% (KH7a absent) | back |

`core` and `norrm` happen to be the two constructs where MD3 and ART crowd MD1
enough to mask KH7a. Averaged over all six: KH7a supplies **30.5%** of
active-face contacts in the active-preferring constructs and **0%** in the
back-preferring ones.

**But that comparison was confounded** — the back-preferring constructs do not
*contain* KH7a, so "0%" is trivially true, and they also have far fewer partners
able to reach MD1 at all (3.5 vs 7.5). Collective packing would fit equally well.

**The controlled test.** Two construct pairs differ by KH7a and nothing else:

| pair | construct | KH7a | preference | active-face contacts/frame |
|---|---|---|---|---|
| both **+ART** | `core` | ✓ | **active** | **105.1** |
| | `mka_full` | – | back | 4.7 |
| both **−ART** | `core_wwe_full_go` | ✓ | **active** | **16.2** |
| | `mka_wwe_full` | – | back | 4.2 |

Both pairs agree, with ART controlled in each. Removing KH7a collapses MD1's
active-face contact load **22×** and **4×** respectively, and flips the
preference. In `core_wwe_full_go` KH7a is the top partner at 46%.

**KH7a is the cause.** This closes the loop with [[Q4]], where adding KH7a
dropped MD1 accessibility 0.445 → 0.398: KH7a sits on the hydrolase's catalytic
face. That is a concrete structural mechanism for the eraser being the most
buried catalytic site in every construct ([[Q3]]).

**One alternative not excluded.** KH7a (738–789) is the immediate N-terminal
neighbour of MD1L1 (790–1004), so some contact is expected from chain
connectivity alone. If MD1's pocket happens to face N-terminally, adjacency
would produce active-face contact without any specific interaction. The
controlled pairs establish that KH7a *causes* the preference; they do not
establish that it is a specific interface rather than a proximity effect.
Distinguishing those needs either a construct where KH7a is present but
sequence-displaced, or a direct check of the pocket's orientation relative to
the domain's N-terminus.

## Figures

`figures/01_static_FL/<date>/face_contacts/`

| figure | shows |
|---|---|
| `face_composition` | residues per face per domain |
| `face_rsa` | mean RSA per face (static structure) |
| `face_contact_enrichment` | **the headline** — which face is preferred, per domain per construct |
| `active_face_contact_load` | absolute contacts per active-face residue |

Plus `data/face_contacts.csv` with every number behind them.

## Reproduce

```bash
cd sim_analysis
python figure_sasa_faces.py                 # three-face split -> sasa_face_per_residue.csv
python figure_face_contacts.py              # contacts + the four bar figures
```

## Still open

- Whether KH7a's effect is a **specific interface or a proximity effect** — it is
  the immediate N-terminal neighbour, so adjacency is a live alternative.
- The face assignment is still from a **single AF2 conformer** ([[Q2]] found the
  molecule is diffusive), so inter-domain burial in the RSA numbers is one
  arrangement's worth. The contact analysis above is ensemble-based and does not
  inherit that limitation.

## Outputs

- `data/sasa_face_per_residue.csv` (1801 rows)
- `figures/09_surface_gallery/2026-08-07/` — per-TICA-state renders:
  electrostatic / hydrophobic / face-classification surfaces (front + back) and
  per-pocket subdirectories (`pocket_MD1/MD2/MD3/ART`). No accompanying text;
  only a batch log at `tica_pipeline/surface_gallery_batch_log.txt`.
