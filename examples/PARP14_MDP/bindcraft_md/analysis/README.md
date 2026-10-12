# Report generator — regeneration state (2026-10-12)

`build_report.py` generates the campaign artifact
(claude.ai/artifact/DcdGTzfwaiib2iaTSDmGC5). The published Version 11 is intact.

**It cannot currently be re-run from this repo.** It was authored against a session
scratchpad under /tmp that has since been cleared, and two classes of input did not
survive:

- **figures** — every `crop_*.png` / `report_*_crop.png` it embeds (structure renders,
  pocket surfaces, binder hits, BC1-vs-BC2 targets). Regenerate with the `render_*.py`
  scripts beside this file; they read from `../designs/` and `../targets/`, which are intact.
- **figdata.json** — the committed copy holds only `_replay`, `_armA`, `_queue`.
  The generator also needs `_pockets`, `_clampcmp` and `_bccmp`, which were computed
  in-session from `../designs/*/mpnn_design_stats.csv`, `trajectory_stats.csv` and
  `failure_csv.csv` — all still present, so they are recomputable.

The analysis *inputs* are therefore all intact; only the derived figures and the
derived figdata entries need rebuilding. `build_report.py` itself is current and
contains every section of the published artifact.
