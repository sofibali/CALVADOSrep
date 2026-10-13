<!-- bindcraft-triage-skill.md -->
---
name: bindcraft-triage-skill
description: A bindcraft-triage skill exists that decodes BindCraft output counters, which are easy to misread
metadata:
  type: reference
---

`~/.claude/skills/bindcraft-triage/` (also version-controlled at
`CALVADOS/examples/PARP14_MDP/bindcraft_md/analysis/skill/`) decodes BindCraft
campaign output and pre-screens hotspot patches. Invoke it before hand-rolling
this analysis again.

The traps it encodes, each of which produced a wrong answer at least once:
`failure_csv.csv` is one flat row of **per-sequence** events (normalise by
relaxed x num_seqs); `Trajectory/LowConfidence/` conflates low pLDDT with
too-few-contacts; the three 4-stage early aborts never write a file so folder
counts understate attempts ~25%; `InterfaceResidues` is keyed by **binder**
residue, not target; and filters must be evaluated over the per-model `1_`-`5_`
keys as well as `Average_`, or passing designs are overcounted.

Apolar fraction of the hotspot patch predicts designability as a threshold:
everything >=42.9% yielded binders, everything <=38.1% yielded none.
See [[parp14-bindcraft-campaign-state]].

<!-- bindcraft-upstream-clone.md -->
---
name: bindcraft-upstream-clone
description: /home/sbali/BindCraft is an upstream clone, not a fork — nothing there is pushable, and params/ must stay untracked
metadata:
  type: project
---

`/home/sbali/BindCraft` has `origin = https://github.com/martinpacesa/BindCraft.git`
— **upstream, not the user's fork**. Nothing committed there can be pushed, so
project work belongs in the CALVADOS repo instead
(`CALVADOS/examples/PARP14_MDP/bindcraft_md/`), which the user owns.

`BindCraft/params/` is 5.3 GB of AlphaFold2 weights and must never be committed.

As of 2026-10-13 the local clone sits 12 commits behind upstream. The local
`functions/pyrosetta_utils.py` patch (PyRosetta 2026.x `set_interface` wanting a
`DockingPartners` object) is **superseded** — upstream now has a more robust
`getattr`-based version that handles both APIs, so a `git pull` is safe and the
local edit can be discarded.

<!-- calvados-repo-branch.md -->
---
name: calvados-repo-branch
description: In the CALVADOS repo the working branch is my-work, not main — push there
metadata:
  type: project
---

`/home/sbali/CALVADOS` has two remotes: `origin` = `sofibali/CALVADOSrep` (the
user's fork) and `upstream` = `KULL-Centre/CALVADOS`.

**All project work lives on the `my-work` branch**, which is also the fork's
default branch. The local `main` branch is a stale copy of upstream CALVADOS
(Jan 2026) with none of the user's work on it, and local `main` tracks
`upstream/main`.

Push with `git push origin my-work`. Do **not** assume `main`: on 2026-10-13 I ran
`git push origin main`, which reported success while pushing the stale branch and
leaving the real commit unpushed. Before claiming something is pushed, verify
against the actual branch, e.g.
`git merge-base --is-ancestor <sha> origin/my-work`.

See [[git-commit-scope]] and [[bindcraft-upstream-clone]].

<!-- git-commit-scope.md -->
---
name: git-commit-scope
description: Commit only my own changes; never rewrite the user's history or commit large binaries
metadata:
  type: feedback
---

Commit **only the changes I made**, not other people's work in progress, and
never rewrite the user's existing commits.

**Why:** these repos are shared with collaborators and the user wants to review
and attribute their own work. When a 4.3 GB tarball in one of their commits made
a push exceed GitHub's 2 GiB pack limit, the right move was to present options
and let them fix it themselves — which they did — not to rewrite history.

**How to apply:** stage explicit paths, never `git add -A`. Check `git status`
for files belonging to others before committing. Keep model weights, trajectory
tarballs and binaries out of git entirely — see [[bindcraft-upstream-clone]].

<!-- MEMORY.md -->
- [Sofia Bali — UCSF Fraser lab](user-sofia-fraser-lab.md) — who the user is and how she wants results delivered
- [Shared GPU etiquette](shared-gpu-etiquette.md) — one GPU, wait for others, no scheduler on the box
- [Git commit scope](git-commit-scope.md) — commit only my own changes, never rewrite her history
- [CALVADOS repo branch](calvados-repo-branch.md) — work lives on `my-work`, not `main`
- [BindCraft is an upstream clone](bindcraft-upstream-clone.md) — not pushable; params/ is 5.3 GB and stays untracked
- [PARP14 BindCraft campaign state](parp14-bindcraft-campaign-state.md) — finished and parked; what was left unfinished on purpose
- [bindcraft-triage skill](bindcraft-triage-skill.md) — decodes BindCraft counters that are easy to misread
- [PARP14 artifacts](parp14-artifacts.md) — report URL and how to rebuild it

<!-- parp14-artifacts.md -->
---
name: parp14-artifacts
description: URL of the BindCraft campaign report artifact and how to rebuild it
metadata:
  type: reference
---

BindCraft campaign report: https://claude.ai/artifact/DcdGTzfwaiib2iaTSDmGC5
("BindCraft Campaign - TMUX BC"). Private. Declares the `db` and `user`
capabilities so the user can add notes under any heading in the browser —
**omit `capabilities` when republishing** or those notes break.

Rebuild it from the repo alone (an earlier version was lost when a `/tmp`
scratchpad was cleared, so the pipeline is deliberately repo-relative):

```
cd CALVADOS/examples/PARP14_MDP/bindcraft_md
python analysis/make_figdata.py && python analysis/make_figures.py
python analysis/crop_figures.py && python analysis/build_report.py
```

Then publish `analysis/parp14_geometry_report.html` to the URL above.
See [[parp14-bindcraft-campaign-state]].

<!-- parp14-bindcraft-campaign-state.md -->
---
name: parp14-bindcraft-campaign-state
description: Where the PARP14 BindCraft campaign stands as of 2026-10-13 and what was deliberately left unfinished
metadata:
  type: project
---

As of **2026-10-13** the BindCraft 1 campaign is finished and parked; no jobs are
running and the GPU is free. Results, analysis scripts and the regeneration
pipeline are committed under
`CALVADOS/examples/PARP14_MDP/bindcraft_md/` (`results/`, `analysis/`).

Deliberately left unfinished, by the user's decision:

- **MD1-focused BC2 arm is prepped but NOT launched** — the user said "actually
  dont launch pause this arm". Configs are in `bc2/campaigns/`. Its goal is to
  improve i_pAE on the 9 BC1 MD1 binders (i_pAE 0.21–0.34, 3 of 9 within 0.02 of
  the 0.35 limit).
- **BC2 selective configs still contain a `forced_targeting` key that breaks
  them** (`epitope_targeting.py:55` raises, because the flag applies to every
  target including the detargets, which intentionally have no hotspots). This was
  documented rather than silently changed — do not "fix" it without asking.
- Remaining BC1 queue runs were dropped, left ready but not run.

Analysis caveat worth keeping: BC2 drops all the Rosetta biophysical gates, so
BC1 and BC2 accepted counts are **not comparable** without re-scoring.
See [[bindcraft-triage-skill]].

<!-- shared-gpu-etiquette.md -->
---
name: shared-gpu-etiquette
description: On the shared GPU box, use one GPU and wait for other users' jobs rather than crowding them
metadata:
  type: feedback
---

Use **one GPU at a time**, and when other users have jobs running, wait for them
to finish rather than starting alongside.

**Why:** the box is shared with the rest of the lab and has no scheduler (no
SLURM), so nothing enforces fair share. Crowding a GPU slows someone else's run
with no warning to them, and the user has to answer for it socially.

**How to apply:** check `nvidia-smi` for other users' processes before launching.
Queue work behind them instead of in parallel. Launch detached
(`setsid nohup ... &`, reparented to `ppid 1`) so jobs survive the session, and
cap trajectory counts rather than letting a run grow unbounded. Never edit a
running bash script — bash reads it incrementally and will execute the splice;
stop the runner, rewrite the whole file, relaunch. Watch out for `pgrep -f`
matching your own command line; use the `bindcraft[.]py` bracket trick or kill
by explicit PID.

<!-- user-sofia-fraser-lab.md -->
---
name: user-sofia-fraser-lab
description: Who the user is — UCSF Fraser lab, PARP14 de novo binder design, computational structural biology
metadata:
  type: user
---

Sofia Bali (sofia.bali@ucsf.edu), UCSF Fraser lab. Works on PARP14 macrodomains
(MD1/MD2/MD3): de novo binder design, MD simulation (CALVADOS), docking, and
wet-lab validation. Comfortable at the command line and reads raw stats CSVs
directly — give numbers and file paths, not reassurance.

Runs jobs manually on a shared GPU box; there is no SLURM, so job orchestration
is hand-rolled shell. See [[shared-gpu-etiquette]] and [[git-commit-scope]].
Likes results delivered as a single self-contained HTML artifact that can be
annotated in the browser — see [[parp14-artifacts]].

