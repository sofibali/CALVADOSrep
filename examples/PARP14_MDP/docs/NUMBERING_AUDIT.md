# PARP14 residue-numbering audit — 2026-09-18, resolved 2026-09-21

Checked every residue definition in the project against **UniProt Q460N5**
(PAR14_HUMAN, 1801 aa, SV=3), fetched live from the UniProt REST API.

Re-run any time with:

```bash
python sim_analysis/verify_numbering.py            # exits non-zero on failure
python sim_analysis/verify_numbering.py --offline  # skip the UniProt fetch
```

> ## STATUS: both defects fixed on 2026-09-21. All checks pass.
>
> | defect | resolution |
> |---|---|
> | catalytic residues didn't match their labels | corrected in `parp14/input/active_sites.yaml` from Đukić et al. 2023 (PMC10499325) + UniProt binding sites + canonical PARP motifs |
> | the numbers were duplicated in 8 modules | all now load `sim_registry.ACTIVE_SITES_FL`, which reads the YAML |
> | residue 1004 in two units | `md1l1` renamed to **`md1`** and narrowed to **790–1003**; MD2 still starts at 1004 |
>
> **Corrected catalytic residues** (UniProt Q460N5):
>
> | domain | residues | provenance |
> |---|---|---|
> | MD1 | N824, **G832**, D961 | G832E is Đukić's hydrolase-blocking mutant; N824/D961 are UniProt binding sites |
> | MD2 | S1034, **G1044**, G1135, H1176 | G1044 is Đukić's MD1-G832 equivalent (G1044E had no effect — MD2 is a reader); rest UniProt |
> | MD3 | S1247, V1258, G1334, F1371 | UniProt binding sites |
> | ART | **H1682**, **R1699**, **Y1714** | H/Y read off the canonical H-G-T and G-K-G-T-Y-F-A motifs; R1699A is Đukić's ADPr-deficient mutant |
> | WWE | unchanged | still unverified — RNF146 alignment, no PARP14 structure |
>
> PARP14 is a **mono**-ART, so there is no catalytic glutamate and no H-Y-E
> triad; the third triad position is hydrophobic and is deliberately not listed
> because the cited paper does not establish it.
>
> **`md1l1` still resolves** via `sim_registry.UNIT_ALIASES`, so the existing
> `metadata.json` files, fragment names and `representative_frames/` directory
> names keep working — nothing historical was rewritten.
>
> **Downstream:** Q3/Q4 never needed re-running (they use `pocket_residues`,
> which were already correct). The catalytic-specific analyses were re-run.
> BindCraft targets generated before this date used the old values; the shift is
> 1–4 residues and every MD1 hotspot was solvent-exposed either way, so those
> targets remain sane, but regenerate the settings before any new campaign.

---

## What follows is the original 2026-09-18 audit, kept as the record of how the
## defects were found.

## PASS — the simulated systems are correct

| check | result |
|---|---|
| local `parp14/PARP14.fasta` vs live UniProt Q460N5 | **byte-identical**, 1801 aa |
| 20/20 slab construct PDBs are exactly the FL slice they claim | **pass** |
| residue counts match `slab_meta.yaml` | **pass** (all 20) |
| `domains.yaml` = FL structured cores mapped to local numbering | **pass** (all 20) |
| FL cores vs UniProt Domain annotations | **exact match** |

The restraint domains are not approximations — they are UniProt's annotations:

| UniProt Domain | range | in `domains.yaml` |
|---|---|---|
| Macro 1 | 791–978 | yes |
| Macro 2 | 1003–1190 | yes |
| Macro 3 | 1216–1387 | yes |
| WWE | 1523–1601 | yes |
| PARP catalytic | 1605–1801 | yes |

(The RRM and KH unit boundaries are the project's own; UniProt does not annotate
them. Nothing contradicts them — they are simply unverifiable against UniProt.)

---

## FAIL 1 — catalytic residue numbers do not match their own labels

`parp14/input/active_sites.yaml` states *"Residue numbers are based on UniProt
Q460N5"*, but the residues it names are not the amino acids it says they are:

| labelled as | Q460N5 actually has |
|---|---|
| D831 | **G**831 |
| N923 | **S**923 |
| D962 | **V**962 |
| H1684 | **T**1684 |
| Y1705 | **N**1705 |
| E1706 | **A**1706 |
| I1722 | **S**1722 |

**No global offset reconciles them** (scanned −6…+6; the best shift matches only
1 of 7). So this is not a systematic frame-shift that could be fixed with a
constant — the numbers have mixed provenance.

Against UniProt's own annotated binding sites, most macrodomain values are
within ±1 — consistent with an off-by-one somewhere in whatever derived them —
but three are further out:

| site | yaml | nearest UniProt binding site | delta |
|---|---|---|---|
| MD1 | 923 | 922–926 | inside |
| MD1 | 962 | 961 | +1 |
| MD1 | 831 | 833 | **−2** |
| MD2 | 1035 | 1034 | +1 |
| MD2 | 1046 | 1046–1049 | inside |
| MD2 | 1134 | 1133–1137 | inside |
| MD2 | 1171 | 1175–1178 | **−4** |
| MD3 | 1248 | 1247 | +1 |
| MD3 | 1259 | 1258 | +1 |
| MD3 | 1330 | 1332–1336 | **−2** |
| MD3 | 1371 | 1371 | exact |

### The ART site is wrong in two separate ways

**Wrong numbers.** The canonical PARP catalytic motifs are unambiguous in the
sequence:

```
H-G-T      at 1682:  H1682 G1683 T1684
G-K-G-T-Y-F-A at 1710:  G1710 K1711 G1712 T1713 Y1714 F1715 A1716
```

So the catalytic histidine is **H1682** (not 1684) and the catalytic tyrosine is
**Y1714** (not 1705). The tyrosine is off by 9 residues.

**Wrong triad class.** CLAUDE.md and `active_sites.yaml` describe an
"H-Y-E catalytic triad". PARP14 is a **mono**-ADP-ribosyltransferase — UniProt's
own protein name is *"Protein mono-ADP-ribosyltransferase PARP14"*. H-Y-E is the
**poly**-ART signature (PARP1/2). Mono-ARTs like PARP14 carry H-Y-[I/L]; there is
no catalytic glutamate to point at. `E1706` is not merely mis-numbered, it is the
wrong kind of residue.

### What this does and does not affect

**Not affected — the accessibility and dominance results stand.**
`sim_analysis/analyze_accessibility.py` builds each site from the
**`pocket_residues`** list, not `catalytic_residues` (it takes the pocket COM and
casts rays from it). Those pocket lists are essentially correct:

| domain | UniProt binding sites covered by the pocket list |
|---|---|
| MD1 | 8 / 8 |
| MD2 | 12 / 16 (misses 1048, 1049, 1176, 1177) |
| MD3 | 10 / 10 |
| ART | contains both the true H1682 and Y1714 |

So **Q3 (active-site accessibility)** and **Q4 (writer-vs-eraser dominance,
`D = (SAA_ART − SAA_MD1)/(SAA_ART + SAA_MD1)`)** are computed over the right
pockets and do not need re-running on account of this.

**Affected — anything keyed on the catalytic list specifically:**

- `sim_analysis/analyze_active_sites.py` — `catalytic_rmsf` and inter-catalytic-site
  distances are computed over the wrong residues
- `sim_analysis/figure_sasa_faces.py` — the active-site-face vs back-face
  annotation uses `ACTIVE_SITES` (catalytic only)
- `sim_analysis/make_analysis_pse.py` — `cat_labels` render
  `D831 / N923 / D962` and `H1684 / Y1705 / E1706 / I1722` in the PyMOL session.
  These labels are simply false and should not go into a figure.
- `tica_pipeline/make_cluster_pse.py`, `tica_pipeline/cluster_states.py`

### Where the numbers live (6 hard-coded copies)

The same wrong values are duplicated rather than read from one source:

```
parp14/input/active_sites.yaml          (the nominal source of truth)
sim_registry.py:69,83
sim_analysis/analyze_all.py:217,229
sim_analysis/make_analysis_pse.py:76,99
sim_analysis/figure_sasa_faces.py:88,91
tica_pipeline/make_cluster_pse.py:57,60
CLAUDE.md                               (the domain/active-site table)
```

**Recommended fix:** decide the correct residues from the primary literature for
each site, put them in `active_sites.yaml` only, and have every consumer import
from `sim_registry`. Do not patch the six copies independently — that is how they
diverged. For the ART site the sequence itself settles it: H1682 and Y1714.

---

## FAIL 2 — residue 1004 belongs to two domain units

`sim_registry.DOMAIN_UNITS` (and the CLAUDE.md table it mirrors):

```python
'md1l1': (790, 1004),
'md2':   (1004, 1193),   # <- 1004 again
```

Residue 1004 is assigned to both units, so the 11 units sum to 1789 residues but
cover only 1788 distinct ones. Any per-residue → domain assignment double-counts
1004 or silently gives it to whichever unit is tested first.

Likely intent is `md1l1: (790, 1003)` — UniProt's Macro 2 domain starts at
**1003**, so 1003 already belongs to MD2's structured core, which makes the
`md1l1` upper bound of 1004 look like the error rather than `md2`'s lower bound.
Worth confirming against however the MD1L1 unit was originally drawn before
changing it.

The 13-residue gap at **1194–1206** is *not* a defect — it is the MD2/MD3 linker,
deliberately unassigned, and `sim_registry` documents it at length (it is exactly
what `CONTIGUOUS_FL_RANGE` exists to protect).

---

## Not checked here

- Whether the AF3 / AF2 input structures in `parp14/alphafold_outputs/` carry the
  same numbering (they are gitignored and were not on the critical path for the
  slab campaign). `verify_numbering.py` covers the slab tree; extending it to the
  fragment library and AF3 outputs would be straightforward.
- The literature provenance of each macrodomain catalytic residue. This audit
  establishes only that the numbers are inconsistent with their own labels and
  with UniProt's annotations — not which value is right.
