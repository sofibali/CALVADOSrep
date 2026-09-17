# Q1 — Which domain boundaries and restraints reproduce the experimental structures?

**Status:** ANSWERED and adopted. One undocumented step (see *Open thread*).

---

## The question

CALVADOS holds folded domains rigid with harmonic restraints between residue
pairs inside a domain. Which residue ranges should be restrained? Restrain too
much and the linkers stiffen and the protein cannot flex; restrain too little
and folded domains inflate. The boundaries are a free parameter that has to be
calibrated against something real.

## Why it matters

Every downstream result inherits this choice. Inter-domain distances ([[Q5]]),
accessibility ([[Q3]]) and the conformational-states analysis ([[Q2]]) are all
sensitive to whether domains are the right size and rigidity.

## How it was tested

`archive_analysis/test_restraints.py` (note: in `archive_analysis/`, not the
repo root). A **trim sweep**: start from each domain's full extent, remove N
residues from each side for N ∈ {0, 5, 10, 15, 20, 25, 30}, each tested with and
without KH7a–KHb custom inter-domain restraints. 16 configs, each a 2 ns
full-length simulation.

Scoring (`test_restraints.py:508-520`):

```
score = 0.4*(1 − Rg_deviation) + 0.4*weighted_cmap_corr + 0.2*max(pLDDT_RMSF_corr, 0)
```

weighted per domain by `score_weight × mean_pLDDT/100`.

**Experimental references** (`compare_restraints.py:65-75`, with residue ranges):

| PDB | domain | FL residues |
|---|---|---|
| 3Q6Z | MD1L1 | 792–978 |
| 3VFQ | MD2 | 1005–1191 |
| 3GOY | WWE | 1534–1602 |
| 3GOY | ART | 1603–1720 |

Target Rg values (nm, `test_restraints.py:65-72`): MD1L1 1.509, MD2 1.533,
WWE 1.461, ART 1.636. KH domains have no PARP14 crystal structure, so homolog Rg
was used: single KH 1.072, tandem KH 1.485.

## The answer

Global ranking of the 16 configs (`restraint_tests/ranking.json`,
`summary_plots/simulation_summary.csv`):

| rank | config | score | Rg dev | CMap corr |
|---|---|---|---|---|
| 1 | `trim_00_kh7ab` | 0.7048 | 0.155 | 0.917 |
| 2 | `trim_00` | 0.7044 | 0.154 | 0.915 |
| 3 | `trim_05_kh7ab` | 0.7029 | 0.197 | 0.877 |
| 7 | `baseline_current` (the old 8-domain set) | 0.6123 | 0.524 | 0.743 |
| 16 | `no_restraints` | 0.1923 | 1.215 | 0.481 |

Three readings:

1. **Restraints matter enormously** — `no_restraints` scores 0.19 vs 0.70.
2. **Restrained trim choice barely matters.** The top configurations are
   separated by ~0.002 (0.7048 / 0.7044 / 0.7029) — far inside any meaningful
   resolution of this score. The sweep's real finding is a **null**: within the
   restrained regime, trim depth is not what determines quality.
3. **The old baseline was worse than every trim ≤ 10** (rank 7 of 16). Replacing
   it was justified.

### The decision that followed

Because trim depth made minimal difference, the choice was made on other
grounds:

> **Use the shortest boundaries that still match the experimental references,
> and add Go-model restraints for the KH domains.**

Shortest-possible boundaries leave the maximum amount of linker free to flex,
which matters for every downstream question about inter-domain geometry
([[Q5]], [[Q2]]) — a restraint that does not improve the score but does stiffen
a linker is a pure cost. The experimental Rg values are the constraint that
stops the trim going too far. The KH domains have no PARP14 crystal reference,
so they are held by Go-model restraints instead of being trimmed against a
reference they do not have.

### What actually shipped

`prepare_fl_optimized.py:1-31` — a **hand-picked per-domain mix**, not the
uniform winner:

| domain | trim |
|---|---|
| RRM1 / RRM2 / RRM3 | 10 |
| KH1 / KH2 / KH3 / KH5 | 5 |
| KH4 / KH6 | 10 |
| KH7a, KHb, KH8 | 0 (full extent — needed for the KH7a–KHb restraints) |
| MD1L1 / MD2 / MD3 | 10 |
| WWE | 15 |
| ART | 10 |

Plus KH7a–KHb custom inter-domain restraints: **141 CA–CA pairs within 0.9 nm**,
k = 350 kJ/mol/nm², generated from the AF2 structure at trim-0 boundaries.

This is why the shipped trims are non-zero even though `trim_00` topped the
ranking: the ranking difference was noise, so the tie was broken toward shorter
restrained cores and more free linker. WWE gets the deepest trim (15) as the
smallest domain, where a full-extent restraint would lock up proportionally the
most of its flanking sequence. KH7a / KHb / KH8 stay at trim 0 because their
full extents are what the Go-model KH7a–KHb pair list was generated against.

## Open thread — worth closing

**The per-domain rationale is not written down anywhere in the repo.** The
reasoning above is recorded here for the first time;
`restraint_tests/comparison_plots/best_per_domain.png` encodes part of it but
has **no accompanying text or CSV**. Worth a short note in
`prepare_fl_optimized.py` so the next person does not re-derive it.

**A second, parallel approach was never adopted.** `sim_analysis/analyze_kh_domains.py`
asked a related question — the FL `domains.yaml` omits the KH regions entirely
(315–789, 1389–1533), treating them as disordered; can AF3 PAE identify
structured KH sub-domains worth restraining? It produced six candidate boundary
sets (`data/domains_{conservative,moderate,moderate_merged,permissive,permissive_merged,aggressive}.yaml`,
15/13/11/11/9/9 blocks). **None is referenced by any script in the repo** — a
grep for those filenames returns zero hits. The KH restraints that shipped came
from the trim sweep, not from this PAE analysis. `ANALYSIS_README.md:125-127`
classifies the script as legacy.

## Outputs

- `restraint_tests/ranking.json`, `restraint_tests/summary_plots/simulation_summary.{csv,md}`
- `restraint_tests/comparison_plots/` (incl. `best_per_domain.png`)
- `input/custom_restraints.txt` — the 141 KH7a–KHb pairs
- `input/xtal_refs/` — 3VFQ, 3Q6Z, 3GOY, 1X4R
