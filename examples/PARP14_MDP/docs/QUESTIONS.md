# PARP14 simulations — what we asked, and what we found

Question-oriented index. `README.md` describes the *infrastructure*;
`sim_analysis/ANALYSIS_README.md` describes the *scripts*. **This** describes the
*findings*.

One file per question in [`questions/`](questions/). Each states the question,
why it matters, how it was tested, the answer, and the caveats that limit it.

Last reviewed: **2026-09-16**.

---

## Status at a glance

| # | Question | Status | One-line answer |
|---|---|---|---|
| [1](questions/Q1_domain_restraints.md) | Which domain boundaries reproduce the experimental structures? | ✅ answered, adopted | Trim depth barely matters — so use the shortest boundaries matching experiment, plus Go restraints for the KHs. |
| [2](questions/Q2_conformational_states.md) | Does the macrodomain module occupy discrete states? | ✅ answered — **negative** | No. Timescales never plateau; it is a diffusive chain, not a state machine. |
| [3](questions/Q3_active_site_accessibility.md) | Are the catalytic sites sterically reachable? | ✅ answered | Invariant ordering, 0/19 violations: **MD1 < MD2 < ART < MD3 < WWE**. |
| [4](questions/Q4_writer_vs_eraser.md) | Is the writer or the eraser favoured? | ✅ answered (steric only) | Writer-leaning everywhere. The KH region selectively buries the eraser. |
| [5](questions/Q5_inter_domain_contacts.md) | Do the domains touch each other? | ✅ answered | Largely not — inter-domain contact is 20–40× rarer than intra-domain. |
| [6](questions/Q6_lysine_contacts.md) | Which lysine pairs can form a DSS crosslink? | ✅ answered (fl) | 135/150 reachable; SASD removed 15 buried pairs the distance test called confident. |
| [7](questions/Q7_surface_faces.md) | Which residues are exposed, and on which face? | ✅ answered | Three-face split on restrained residues. ART and MD3 present back faces; the rim is depleted in all 24 domain×construct combinations. |
| [8](questions/Q8_binder_design.md) | Can we block or clamp the macrodomains? | 🔄 ongoing | 11 targets + a 3-arm diagnostic sweep. First pass: 0 accepted, failing at interface confidence. |
| [9](questions/Q9_self_association_puncta.md) | Which constructs form puncta with RNA / ADPr? | 🔜 prepared | Slab sims built and validated; needs GPU + the sort-seq data. |
| [10](questions/Q10_crosslink_discrimination.md) | Which constructs should we crosslink to map eraser/reader/writer exposure? | ✅ answered | Yes. 9 discriminating domain pairs; XL count predicts exposure (r −0.56 to −0.89, size-independent). Panel: kh1_art_full + core + fl. |

---

## The through-line

Three answered questions tell one coherent story about what kind of molecule
PARP14 is:

**It is a flexible chain of domains that rarely touch** ([[Q5]]) — inter-domain
contact frequency ~0.003–0.006 against an intra-domain reference of ~0.12–0.13.

**So it has no metastable conformational states** ([[Q2]]). The implied
timescales rise linearly with lag and never plateau; going from 5 ns to 100 ns of
sampling did not produce a plateau, it just pushed the apparent timescale up in
lock-step. The right read is a diffusive multidomain chain, not a two- or
three-state switch. That is a genuine negative result, and it changed the
downstream plan: "representative states" means geometric diversity, not kinetic
basins.

**And yet the catalytic sites are ordered identically in every construct**
([[Q3]]): MD1 < MD2 < ART < MD3 < WWE, 0 violations in 19 constructs. Flexibility
at the linker level coexists with a rigid, fold-determined hierarchy of pocket
accessibility. The eraser is always the most buried site; the writer always has
more room ([[Q4]]), and it is specifically the **KH region** that buries the
eraser — adding KH7a then KH1–KH6 drops MD1 accessibility 22 % while ART moves 3 %.

The practical consequence runs through to design: because there is no stable
interface to stabilise, the binder campaign has to *enforce* a new one ([[Q8]]).

---

## Open threads, ranked by how much they block

1. **[Q8] Finish the diagnostic sweep.** Three arms are set up on a fixed target
   (`md1_block_af3` hotspots) to separate *target difficulty* from *filter
   strictness* — the first pass died at `i_pAE` / `i_pLDDT`. Runs are incomplete
   (`sweep_af3_diag.log` ends in a `KeyboardInterrupt`). Two READMEs still say
   BindCraft is not installed and nothing has been run; that is stale.
2. **[Q6] Run the remaining sets.** Only `fl` has been regenerated with the
   corrected method; the other sets still carry the superseded numbers.
3. **[Q7] Write the per-domain face summary.** The table supports it and two
   other analyses already depend on the file.
4. **[Q1] Note the trim rationale in the code.** The reasoning (trim depth is
   within noise → prefer shortest boundaries + Go restraints for the KHs) is now
   recorded in Q1 but not in `prepare_fl_optimized.py`.
5. **[Q9] Benchmark the GPU, then launch.** Cost estimates rest on an assumed
   50× speedup that has never been measured.

---

## Conventions that bite

Things that have already caused a wrong result once.

- **`MD1L1` is not MD1.** It is MD1 *plus* the linker to MD2. A contact/distance
  computed on that block is contaminated by the linker — this produced a since-
  retracted claim that MD1–MD2 were "permanently in contact at 0.38 nm"
  ([[Q5]]). Use trimmed cores for contact work.
- **Two sampling regimes, do not mix them.** The 25-replicate original sets
  (3–4 k frames) sample *structural* diversity across distinct AF3 seeds, spread
  ±0.04. The 5-replicate extensions (100 k frames) sample *conformational* time
  from similar starts, spread ±0.001. A tight CI in the second regime is
  precision, not accuracy.
- **N = 5 floor on significance.** With 5 replicates per group the smallest
  attainable Mann–Whitney p is 2/C(10,5) = 0.0079, so nothing at that depth
  survives correction at α = 0.001 (`README.md:200-222`).
- **Check the size confound before believing a cross-construct correlation.**
  The constructs differ ~5x in length, so anything that scales with length
  correlates with anything else that does. [[Q10]]'s XL-vs-exposure correlation
  only counts because `r(XL count, length)` is ~0 and the partial correlation
  survives.
- **Rank candidates by what the question needs, not by a global score.** [[Q10]]'s
  first pass returned zero inter-domain crosslinks because inter-domain pairs are
  ~6x rarer than intra-domain ones and a global persistence ranking buries them
  entirely. The budget has to be spent on the class that answers the question.
- **A filter threshold must be derived, not fitted to the data you are
  filtering.** [[Q6]]'s sequence-separation cutoff was first set to 20 by eye,
  to sit above an artefact band that had been observed rather than predicted —
  which silently discarded 245 informative pairs. The defensible value (9) falls
  out of chain geometry: `ceil(cutoff / bond_length) + 1`.
- **`latest` symlinks can point at empty directories** — seen under
  `figures/04_md_distances/` and `figures/07_lysine_contacts/`. Check the dated
  dir before concluding a figure is missing.
- **No GPU on this workstation.** `pollux.fraserlab.com` fails `nvidia-smi`;
  OpenMM sees only Reference/CPU, measured ~5.76e5 bead-steps/s. Multi-chain
  work and AF3 inference are GPU-only and target the SLURM GPU server that shares
  this filesystem.

---

## Primary sources

| Topic | File |
|---|---|
| TICA / MSM theory + the full Q2 derivation | [`../TICAresultDescreption.md`](../TICAresultDescreption.md) |
| Cross-set contact + distance prose | `../analysis_summary.txt` |
| Binder design rationale | `../bindcraft_md/README.md` |
| Simulation inventory + statistics conventions | [`../README.md`](../README.md) |
| Script reference | [`../sim_analysis/ANALYSIS_README.md`](../sim_analysis/ANALYSIS_README.md) |
