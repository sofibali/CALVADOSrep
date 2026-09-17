# Q4 — Is the writer (ART) or the eraser (MD1) favoured, and what controls it?

**Status:** ANSWERED as a *steric* ranking (2026-09-16). NOT calibrated to activity.

---

## The question

PARP14 carries both an ADP-ribosyltransferase (ART, the writer) and a
macrodomain hydrolase (MD1, the eraser). In a given construct, which of the two
has more opportunity to engage substrate — and which domains decide that?

## Why it matters

The net ADP-ribosylation output of a construct is writer minus eraser. If domain
content shifts that balance sterically, domain-deletion variants are not just
"less PARP14" — they are differently biased enzymes. This is the model-side
counterpart to the catalytic-dead controls in the puncta screen ([[Q9]]).

## How it was tested

From the SAA values in [[Q3]], a normalised replicate-wise contrast:

```
D = (SAA_ART − SAA_MD1) / (SAA_ART + SAA_MD1)

D > 0  writer site more sterically available
D < 0  eraser site more sterically available
```

Normalised rather than a raw difference so constructs whose sites are both open
and both closed stay comparable. Replicates are paired (same trajectory), so the
contrast is taken per replicate before averaging. 95 % CIs by percentile
bootstrap over replicates.

Two designed comparisons make it causal rather than correlative:

- **Nested series** — three constructs under one identical protocol, differing
  only in how much sequence sits N-terminal to MD1.
- **Matched ±ART pairs** — constructs differing by nothing but the presence of
  ART, isolating its steric footprint.

## The answer

### 1. Every construct is writer-leaning

D = **+0.127 to +0.232**, all bootstrap CIs excluding zero. Nothing is
eraser-leaning. `mka_full` +0.127 → `fl` +0.196 → `kh1_art_full` +0.232.

### 2. The KH region is what buries the eraser

| construct | adds | MD1 (eraser) | ART (writer) |
|---|---|---|---|
| `mka_full` | — (MD1–ART only) | 0.427 | 0.551 |
| `core_full_go` | + KH7a | 0.362 | 0.546 |
| `kh1_art_full` | + KH1–KH6 | **0.333** | 0.533 |

**MD1 falls 22 %, ART moves 3 %.** The occlusion is selective, not general
burial — the KH region packs against the hydrolase site specifically.

### 3. ART only reaches its neighbour

Matched ±ART pairs, change in each site's SAA when ART is present:

| site | pairs | mean Δ SAA |
|---|---|---|
| **WWE** | 6 | **−0.0928** |
| MD1 | 4 | −0.0133 |
| MD2 | 5 | −0.0099 |
| MD3 | 6 | −0.0060 |

ART occludes WWE **7–15× more** than any macrodomain — consistent with WWE
(1534–1602) being contiguous with ART (1603–1801) and folding against it
(cf. PDB 3GOY, which contains WWE+ART).

## Confidence and caveats

- **This is not an activity prediction.** It ranks steric opportunity. Turning
  it into a writer-vs-eraser activity claim requires an in-vitro ADPr
  transfer/hydrolysis readout on at least a few constructs to anchor the scale.
- The `fl` / `fl_optimized` ↔ `fl_wwe_full_go` ±ART pairs are **protocol-
  mismatched** (3.5 k-frame original/optimised restraints vs 100 k-frame Go
  model), so their contrast confounds ART removal with the restraint scheme.
  They are reported but excluded from the headline figure.
- `d/SD` in `art_steric_effect_saa.csv` divides by the replicate spread of
  identical starting structures, so it runs into the tens. It is a precision
  ratio, **not** a Cohen's d. Read `delta`.

## Reproduce

```bash
cd sim_analysis && python predict_enzyme_dominance.py --metric saa
```

## Outputs

- `data/enzyme_dominance_saa.csv`, `data/art_steric_effect_saa.csv`
- `figures/03_accessibility/latest/enzyme_dominance/writer_eraser_dominance.svg`
- `figures/03_accessibility/latest/enzyme_dominance/art_domain_steric_effect.svg`
