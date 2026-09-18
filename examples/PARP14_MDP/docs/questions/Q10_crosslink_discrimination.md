# Q10 — Which constructs should we crosslink to map eraser vs reader vs writer exposure?

**Status:** ANSWERED (2026-09-16) · 19 constructs · the correlation holds and survives a size-confound check.

---

## The question

Different PARP14 constructs arrange their domains differently. If those
differences produce **different feasible crosslinks**, then XL-MS on a
well-chosen panel reads out the architecture — and specifically the exposure of
the eraser (MD1), the readers (MD2 / MD3 / WWE) and the writer (ART).

So: *can we predict, before doing the experiment, which constructs are worth
crosslinking?*

## Why it matters

This is the bridge between the computational and experimental halves of the
project. [[Q3]] and [[Q4]] give a *predicted* exposure hierarchy
(MD1 < MD2 < ART < MD3 < WWE) that is currently uncalibrated — nothing has
tested it at the bench. Crosslinking is the cheapest experiment that could,
because a crosslink is a distance restraint and the constructs differ in their
distances.

It also turns construct choice from intuition into a designed panel: rather than
crosslinking whatever is available, pick the minimal set whose predicted
signatures differ most.

## How it is tested

Three steps, in order, with the third being the one that decides whether the
whole idea works.

### 1. Which inter-domain crosslinks are feasible per construct?

From `analyze_lys_contacts.py --sasd` (see [[Q6]] for the method): DSS-reachable
lysine pairs, sequence-separation filtered at the derived threshold of 9, and
re-scored by solvent-accessible surface distance so buried pairs are dropped.

**Only inter-domain pairs count here.** An intra-domain crosslink reports on a
fold that does not change between constructs; only inter-domain ones report on
architecture.

### 2. Which of those discriminate between constructs?

A domain pair that is feasible in *every* construct containing both domains
carries no comparative information. The useful ones form in some constructs and
not others — and the ones that cross **role boundaries** (eraser ↔ reader ↔
writer) are the most valuable, because those directly report the geometry
between the three activities.

The script then computes a **greedy minimal panel**: the fewest constructs whose
combined crosslink sets cover every discriminating domain pair.

### 3. Does the crosslink signature actually track exposure?

**This is tested, not assumed.** For each catalytic site, its inter-domain
crosslink count is correlated against its solid-angle accessibility from
[[Q3]]/[[Q4]] across constructs.

| result | reading |
|---|---|
| \|r\| ≥ 0.6 | crosslink count is a usable proxy for exposure |
| 0.35 ≤ \|r\| < 0.6 | suggestive, not a proxy on its own |
| \|r\| < 0.35 | crosslinks report **contacts**, not pocket openness — the two must be measured separately |

A null result here is a real and useful answer. [[Q3]] found accessibility is
fold-determined and invariant across constructs, while [[Q5]] found inter-domain
contacts are rare and variable. Those two facts could easily mean crosslink
patterns vary *without* tracking exposure at all, in which case XL-MS would
constrain architecture but say nothing about which pocket is open.

## The flaw that had to be fixed first

The first attempt returned **zero** inter-domain crosslinks, which exposed a
real design error rather than an empty result:

> Measured on `fl`, inter-domain pairs average **0.021** persistence against
> **0.135** for intra-domain. They are systematically rarer, so ranking
> candidates globally by persistence buries every one of them — the first
> inter-domain pair appears only at rank 43, and the SASD budget never reached it.

The SASD pass now ranks **inter-domain pairs first**. On `fl` that changed the
candidate set from 0 to 400 inter-domain pairs (293 solvent-reachable) out of
34,047 available.

Two knobs worth knowing:

- `--sasd-top` is a **sample**, not an enumeration. 500 of ~34k inter-domain
  pairs is chosen for coverage of domain *pairs*, enough to say which interfaces
  are formable in which constructs — not enough to claim a complete crosslink
  inventory. Cost is ~0.75 s/pair.
- `--target-frames` (default 2000) picks the stride per trajectory. Without it
  the sweep is 9.46M frames ≈ 18 h; verified to reproduce the every-frame result
  exactly on `fl` (135/150 either way) at ~24× less cost.

## Reproduce

```bash
cd sim_analysis
bash run_crosslinks_all.sh                      # ~2-3 h, all 19 sets
python predict_crosslink_discrimination.py
```

## Outputs

- `data/{set}_lys_sasd.csv` — per-construct crosslinks, with `domain_i`,
  `domain_j`, `inter_domain` columns
- `data/crosslink_discrimination.csv` — per domain pair: which constructs form
  it, which do not, whether it discriminates
- `data/crosslink_shopping_list.csv` — the actual residue pairs to watch
- `figures/07_lysine_contacts/<date>/discrimination/` — construct × domain-pair
  heatmap, grey where a domain is absent from the construct

## Results

**Answer: yes — and the crosslink signature does track exposure.**

### 1. Constructs differ strongly in their inter-domain crosslinks

19 constructs, 500 inter-domain candidates each, SASD-filtered. Feasible
inter-domain crosslinks range from **0** (`md1_md2`, `md2_md3`) to **461**
(`norrm`), spanning 0-11 distinct domain pairs.

### 2. Nine domain pairs discriminate

Pairs that form in some constructs and not others, among constructs that contain
both domains:

| domain pair | role | forms in | absent in |
|---|---|---|---|
| **MD1-MD2** | eraser ↔ reader | 10 | 4 |
| **MD1-MD3** | eraser ↔ reader | 4 | 9 |
| **ART-MD1** | **writer ↔ eraser** | 2 | 5 |
| KHb-KH8-MD1 | eraser | 5 | 7 |
| MD2-MD3 | reader | 5 | 11 |
| KH7a-KHb-KH8 | scaffold | 6 | 4 |
| KH1-6-KHb-KH8 | scaffold | 4 | 3 |
| ART-KHb-KH8 | writer | 8 | 1 |
| KH1-6-RRM2 | scaffold | 2 | 1 |

The first three are the ones that matter: they cross role boundaries and so
report directly on the geometry between the three activities. **ART-MD1 is the
writer-eraser distance itself** — feasible only in `core` and `norrm`.

### 3. Crosslink count predicts exposure — and it is not a size artefact

Correlating each site's inter-domain crosslink count against its solid-angle
accessibility from [[Q3]]/[[Q4]]:

| site | n | r (XL vs SAA) | verdict |
|---|---|---|---|
| ART | 9 | **−0.911** | strong |
| MD2 | 17 | **−0.837** | strong |
| WWE | 16 | **−0.709** | strong |
| MD1 | 14 | **−0.706** | strong |
| MD3 | 18 | **−0.685** | strong |

Computed with SAA at the 3 nm probe, matching the crosslink ceiling. Every site
is stronger than at the earlier 5 nm probe (MD3 moved from −0.564 "weak" to
−0.685 "strong") — making the two measures use the same reachability criterion
improves their agreement, which is itself a small consistency check.

The correlation is **negative**, which is the physically sensible direction: a
site carrying many inter-domain crosslinks is packed against its neighbours, and
therefore occluded.

**The obvious confound is ruled out.** Bigger constructs could trivially have
both more crosslinks and more burial. They do not: `r(XL count, construct
length)` is between −0.07 and +0.09 for every site, and the partial correlation
controlling for length is *stronger* than the raw one in all five cases.

So crosslink count is a usable experimental proxy for pocket exposure. That is
the result that makes the whole approach worth doing at the bench.

### 4. Recommended panel — three constructs cover all nine discriminating pairs

| # | construct | covers | notable |
|---|---|---|---|
| 1 | `kh1_art_full` | 6 pairs | MD1-MD2, KHb-KH8-MD1, MD2-MD3, KH7a-KHb-KH8 |
| 2 | `core` | 5 pairs | **ART-MD1** (writer-eraser), MD1-MD3 |
| 3 | `fl` | 2 pairs | ART-KHb-KH8, KH1-6-RRM2 |

`core` earns its place specifically because it is one of only two constructs
where the direct writer-eraser crosslink is feasible.

### A prediction worth testing

The two-domain constructs `md1_md2` and `md2_md3` show **zero** inter-domain
crosslinks, while in larger constructs MD1-MD2 forms readily (111-147 in `core`,
`md_full`, `noart`, `norrm`). Read literally, the macrodomains only come into
contact when more of the chain is present — isolated pairs stay apart.

Caveat before trusting it: those two constructs are the smallest and have 5
replicates rather than 25, so a sampling difference cannot be fully excluded.
The effect size (0 vs 132) is large enough to be worth a direct test.

## Caveats

- `--sasd-top 500` samples the 500 most persistent inter-domain pairs out of
  ~25-34k per construct. Enough to say which *interfaces* are formable; not a
  complete crosslink inventory.
- Crosslink feasibility is predicted from a CG ensemble with a coarse surface
  ([[Q6]]). Treat SASD as a reachability filter, not a quantitative distance.
- The exposure correlation is across constructs, not a within-construct
  calibration. It says a construct with more MD1 crosslinks has a more buried
  MD1; it does not convert a crosslink count into an SAA value.
