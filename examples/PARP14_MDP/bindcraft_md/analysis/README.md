# Report generator

`build_report.py` generates the campaign artifact
(claude.ai/artifact/DcdGTzfwaiib2iaTSDmGC5). **It rebuilds from this repository alone** —
everything it needs is committed, and no path points outside the repo.

    python analysis/make_figdata.py          # derived numbers  -> figdata.json
    pymol -cq analysis/make_figures.py       # renders          -> assets/
    python analysis/crop_figures.py          # trims/downscales -> assets/crop_*.png
    python analysis/build_report.py          # -> analysis/parp14_geometry_report.html

Run from `bindcraft_md/`. The first two need the `BindCraft` and `pymol-render` conda
envs respectively.

## Why these scripts exist

The figures and three `figdata.json` entries (`_pockets`, `_clampcmp`, `_bccmp`) were
originally computed in a session scratchpad under /tmp, which was later cleared. The
published artifact survived but could not be rebuilt. These scripts recompute all of it
from `../designs`, `../targets`, `../sweep` and `../bc2`, which are committed — so the
loss cannot recur. Regenerated values reproduce the originals exactly (e.g. pocket
apolar fractions 43.8 / 38.1 / 28.6 / 42.9 / 28.6).

## What is tracked

- `assets/crop_*.png`, `assets/report_*_crop.png` — the images the report embeds
- `figdata.json`, `hits_seqs.json`, `geometry_results.json`, `replay_md1_block_af3.csv`,
  `run_summary.csv` — derived data
- the four scripts above, plus `geometry_scan.py`, `replay_mpnn.py`, `render_*.py`

Not tracked: uncropped renders (`assets/.gitignore`) and the generated HTML — both
rebuild in minutes.
