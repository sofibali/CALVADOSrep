# Campaign results

Committed outputs of the PARP14 binder campaigns. The full run directories live in
`../designs/` (BC1, 5.0 GB) and `/home/sbali/BindCraft2/results/` (BC2, 43 MB) and are
gitignored; what is here is everything needed to read, rank or re-score the work:
**every accepted structure, and every per-run statistics table**. The ~4.7 GB of
trajectory PDBs and the 35 MB JAX compile cache are excluded on purpose.

## BC1 (`bc1/`, BindCraft 1 + predict_initial_guess)

Per run: `trajectory_stats.csv`, `mpnn_design_stats.csv`, `final_design_stats.csv`,
`failure_csv.csv`, and `accepted/` where anything was accepted.

| run | relaxed | scored | accepted |
|---|---|---|---|
| `sweep_af3_guess` (MD1, the fix) | 12 | 51 | **9** |
| `guess_clamp_md1md2_state3` | 50 | 197 | **5** |
| `guess_clamp_md1md2_state1` | 28 | 49 | 0 |
| `guess_md2_block_af3` | 50 | 16 | 0 |
| `guess_md3_block_af3` | 50 | 0 | 0 |
| `guess_md3_block_sim` | 12 (partial) | 0 | 0 |
| `md1_block_af3` (baseline, no flag) | 270 | 10 | 0 |
| `md1_block_sim` | 91 | 0 | 0 |
| `clamp_md1md2_state3` (baseline) | 220 | 0 | 0 |

`accepted/` holds both the `Accepted/` copies and the rank-prefixed `Accepted/Ranked/`
copies, so `sweep_af3_guess` shows 18 files for 9 designs.

## BC2 (`bc2/`, BindCraft2)

`campaign_metadata.json` per campaign, plus `3_Ranked/` and `summary.csv` where a
campaign produced anything.

Only `parp14_md1_baseline` ran to completion: 60 trajectories, 2 accepted. **Those two
bind the wrong face** — `Hotspot_Contact_Fraction = 0.0`, zero overlap with the 16 MD1
hotspots — so they are kept as a record, not as candidates. The three selective
campaigns refused in ~23 s on the `forced_targeting` detarget bug (see `../bc2/README.md`).

## Caveat before comparing BC1 and BC2 numbers

BC2 applies **none** of the Rosetta biophysical gates — buried unsatisfied H-bonds,
shape complementarity, dG, dSASA, interface H-bond count, surface hydrophobicity — which
are precisely the gates that rejected every BC1 design. A higher accepted count under
BC2 does not by itself mean better binders. Re-score BC2 structures against the BC1
filters before treating the two as comparable.
