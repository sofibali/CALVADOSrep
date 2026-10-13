---
name: bindcraft-triage
description: "Diagnose where a BindCraft binder-design campaign is losing designs, and pre-screen a hotspot patch before spending GPU-days on it. Use when asked why a run produced few or no accepted binders, what the failure counts mean, which filter is the real bottleneck, whether a target site is worth running, or when reading failure_csv.csv / mpnn_design_stats.csv / trajectory_stats.csv."
---

# BindCraft triage

Reading a BindCraft campaign's output to say **where designs die** is easy to get
wrong: several of the counters do not mean what their names suggest. This skill
encodes the decode rules, the scripts that apply them, and a pre-screen that
predicts whether a site is designable at all.

## 1. Run the triage

```bash
python ~/.claude/skills/bindcraft-triage/scripts/triage.py <design_path> [<design_path> ...]
```

`<design_path>` is the folder holding `Trajectory/`, `Accepted/`,
`mpnn_design_stats.csv`, `failure_csv.csv`. Add `--filters <path>` if the campaign
used something other than `settings_filters/default_filters.json`.

It prints the true funnel, the per-sequence AF2 rejections, and a per-metric
failure table over the designs that reached Rosetta scoring, ending in a verdict.

## 2. Decode rules — the traps

Apply these whether or not you use the script. Each one caused a wrong conclusion
at least once.

**`failure_csv.csv` is a single flat row.** One header, one data row, cumulative
over the whole run. Read it with `pd.read_csv(...).iloc[0]`, not as a per-design table.

**Its counters are per-SEQUENCE events, not per-design.** Each relaxed trajectory
spawns `num_seqs` MPNN sequences (default 20). A value of 1000 against 50 relaxed
trajectories means 1000/1000 sequences — 100%, not "1000 designs failed". Always
normalise by `relaxed × num_seqs` before quoting a percentage.

**`Trajectory/LowConfidence/` conflates two different failures:** low final pLDDT
*and* too-few-contacts. Separate them with the `Trajectory_final_pLDDT` and
`Trajectory_Contacts` columns of `failure_csv.csv`. A site that fails on contacts
is a *reachability* problem; one that fails on pLDDT is a *foldability* problem.
They call for opposite fixes.

**Early aborts never produce a file.** `Trajectory_logits_pLDDT`,
`Trajectory_softmax_pLDDT` and `Trajectory_one-hot_pLDDT` abort during the 4-stage
hallucination before anything is written. They exist only in `failure_csv.csv`.
Counting PDBs in `Trajectory/*/` therefore **understates attempts by roughly 25%**
and overstates the trajectory yield. Always add them back.

**`InterfaceResidues` in the stats CSVs is indexed by BINDER residue, not target.**
From `functions/biopython_utils.py:171`, `hotspot_residues()` keys the dict on
`binder_residue.id[1]`. Do not use this column to ask which part of the *target*
a design bound — the numbers will look plausible and mean nothing. For target-side
contacts, run a KDTree over the complex directly (chain A atoms within 4.0 Å of
chain B), as `hotspot_residues()` does internally.

**Filters are evaluated per-model as well as on the average.**
`default_filters.json` carries `Average_*` **and** `1_*`–`5_*` keys (36 metrics × 6).
Checking only `Average_*` overcounts passing designs. With all 54 numeric gates the
"fails 0 gates" count matches `Accepted/` exactly; with averages alone it does not.
Residual gaps after that are the non-numeric checks (`InterfaceAAs` composition,
duplicate sequence).

**A metric failing 0% is not a bottleneck.** Report the gates that actually reject
designs and say explicitly which ones are not the problem, so nobody tunes a
threshold that was never binding.

## 3. Pre-screen a site before running it

```bash
python ~/.claude/skills/bindcraft-triage/scripts/pocket_chem.py <settings.json> [...]
# or: --pdb target.pdb --hotspots 112,115,119 [--chain A]
```

Reports apolar fraction, polar fraction, composition, net charge, H-bond
donor/acceptor burden per residue, and patch SASA.

**Apolar fraction of the hotspot patch is the strongest observed predictor of
designability**, and it behaves like a threshold, not a gradient:

| apolar % | patch | relaxed | scored | accepted |
|---|---|---|---|---|
| 43.8 | MD1 block | 12 | 51 | **9** |
| 42.9 | MD1–MD2 clamp state3 | 50 | 197 | **5** |
| 38.1 | MD2 block | 50 | 16 | 0 |
| 28.6 | MD1–MD2 clamp state1 | 28 | 49 | 0 |
| 28.6 | MD3 block | 50 | 0 | 0 |

Everything ≥42.9% produced binders; everything ≤38.1% produced none. The gap
between them is untested. A polar, charged, high-H-bond-burden patch is the
failure mode — it is not fixable by running longer.

If the script reports **zero hotspots matched**, that is a numbering mismatch
(UniProt vs construct), not an empty site. Check it before concluding anything.

## 4. Reading the verdict

The two walls seen so far, in order:

**i_pAE (threshold 0.35) is the upstream wall.** When `scored = 0` despite many
relaxed trajectories, the campaign never reached interface chemistry at all —
AF2 could not confidently re-predict the complex from sequence. The fix is
`predict_initial_guess` in the advanced settings, which re-predicts from the
trajectory's own pose rather than folding and docking from sequence alone. On
this project it took one target from 0 accepted (270 relaxed) to 9 accepted
(12 relaxed).

**Once the flag is on, i_pAE fails 0%** and the limiter moves downstream to the
interface-geometry gates — `ShapeComplementarity` (≥0.6) and
`n_InterfaceUnsatHbonds` (≤4) are the two that bind in practice. These are
properties of the *site*, so the honest answer is often "this patch is not
designable", not "run more trajectories".

Do not generalise a single run's limiter into a universal claim. Run the script
across every design_path in the campaign and compare before stating a pattern.

## 5. Comparing BindCraft 1 and BindCraft 2

BindCraft 2 drops PyRosetta entirely and ranks on `i_pDAE`. It has **no** `dG`,
`ShapeComplementarity`, `Surface_Hydrophobicity` or `n_InterfaceUnsatHbonds`
gates. Accepted counts are therefore **not comparable** between the two without
re-scoring BC2 output through the BC1 interface scorer. Say so rather than
putting the two numbers in the same column.

BC2 also treats hotspots as a *soft preference*. Without
`min_hotspot_contact_final`, it will happily accept designs with
`Hotspot_Contact_Fraction = 0.0` — binders to the wrong face that look like
successes. Check that column before reporting a BC2 yield.
