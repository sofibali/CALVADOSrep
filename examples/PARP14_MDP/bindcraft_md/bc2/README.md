# BindCraft2 arm (paused 2026-10-12)

BC2 (`PacesaLab/BindCraft2`) installed at `/home/sbali/BindCraft2` with its own `.venv`
(Python 3.12, jax 0.11.2, CUDA 13). It is a rewrite, not a version bump: **no PyRosetta**,
ranks on `i_pDAE`, and supports multi-target design with off-target penalties.

Everything here is the campaign definition. It lives in this repo rather than in the
upstream clone, which is read-only to us and would lose it on any re-clone.

## Two differences from BC1 that change what a result means

1. **Different target.** BC1's `md1_block_af3` is 404 residues (MD1 + MD2 scaffold);
   BC2's `PARP14_MD1` is 215 (MD1 alone), because selectivity needs the domains
   separable. This matters: of BC1's 270 relaxed MD1 trajectories, **170 contacted MD2**
   as well as MD1 — surface that does not exist in the BC2 target, and would not exist
   on an isolated macrodomain experimentally.
2. **Different acceptance bar.** BC2 drops every Rosetta biophysical gate — buried
   unsatisfied H-bonds, shape complementarity, dG, dSASA, interface H-bond count,
   surface hydrophobicity. Those are exactly the gates that rejected every design this
   project produced. **A higher accepted count under BC2 does not by itself mean better
   binders.** Re-score BC2 designs against the BC1 filters before trusting a comparison.

## Run state

| campaign | status |
|---|---|
| `md1_baseline` | RAN. 60 trajectories, 2 accepted — but both bind the WRONG face (`Hotspot_Contact_Fraction = 0.0`, zero overlap with our 16 hotspots). Not comparable to BC1's 9. Cause: no `min_hotspot_contact_final`, and BC2 treats hotspots as a soft preference where BC1 enforced them implicitly. |
| `md1_selective`, `md2_selective`, `md3_selective` | REFUSED after ~23 s. `forced_targeting` rebuilds non-epitope surface as lysine and applies to EVERY target, including the detargets, which deliberately carry no hotspots -> `epitope_targeting.py:55` raises. **Fix: drop `forced_targeting` from these three**; `min_hotspot_contact_final` is a separate filter and still enforces focus. |
| `md1_focused` | PREPARED, never launched. Arm 1 of the focused plan. |
| `md1_focused_strict` | PREPARED, never launched. Arm 2. |

## The focused MD1 plan (prepared 2026-10-07, not run)

Motivation: our 9 BC1 MD1 binders sit at i_pAE 0.21-0.34 against a 0.35 limit, median
0.26, with **3 of 9 within 0.02 of failing**. The goal is genuinely confident interfaces,
not marginal ones.

- **Arm 1 `md1_focused`** — `forced_targeting` (works here: single target WITH hotspots)
  plus `min_hotspot_contact_final 0.5`, default filters, 300 trajectories. Asks: does
  forced targeting put binders on the pocket, and at what yield?
- **Arm 2 `md1_focused_strict`** — same, plus `i_pAE <= 0.25` and `i_pTM >= 0.80`.
  Asks: can we beat our current median rather than merely pass?

Arm 1 runs first on purpose. Stacking focus and strict filters together risks the MD3
outcome — zero designs with no way to attribute the cause.

**`--initial-guess` is deliberately NOT used.** It re-predicts from the trajectory's own
pose, which inflates i_pAE rather than improving the binder. Since better i_pAE is the
goal, using it would defeat the measurement.

Weak signal worth testing: the two best BC1 designs (i_pAE 0.21) were both 80-mers; the
three worst (0.33-0.34) were 104 and 122. Nine points and not monotonic, so the 70-130
range was left open rather than narrowed on thin evidence.

## To run

    cd /home/sbali/BindCraft2 && source .venv/bin/activate
    cp <this repo>/bc2/campaigns/*.json parp14_campaigns/     # if the clone is fresh
    CUDA_VISIBLE_DEVICES=<n> bindcraft design parp14_campaigns/md1_focused.json

`bc2/scripts/wait_and_run_arm1.sh` waits for a GPU with zero compute processes from
anyone before starting, so it does not crowd shared cards.
