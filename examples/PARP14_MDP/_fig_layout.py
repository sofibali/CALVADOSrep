"""Figure output layout: by-simulation (default) or the legacy dated layout.

Used by all analysis scripts. Each script picks its category (see below) and
gets back a Path object that's already created.

    figures/
      by_sim/<sim_name>/<category>/<subname?>/<filename>
          -- figures scoped to ONE simulation (the common case: any script
             invoked with a single --set / --sim-folder). Overwritten in
             place on re-run (no date folders) -- these are regenerable
             analysis products, not a lab notebook that needs history; the
             underlying data/*.npz cache already works the same way.
      comparisons/<category>/<sim1_sim2_.../<subname?>/<filename>
          -- figures spanning MULTIPLE simulations (built-in cross-set
             comparisons inside analyze_all.py, or compare_sims.py output).
      <category>/<YYYY-MM-DD>/<subname?>/<filename>   (legacy, sims=None)
          -- unmigrated scripts that aren't sim-scoped at all: they compare
             AF3 structure predictions, not CG-MD trajectories, so there's no
             single "simulation" to file them under (figure_sasa_faces.py,
             figure_domain_comparisons.py, figure_kh_domains.py,
             analyze_kh_domains.py, figure_schematic.py, make_analysis_pse.py).

Categories (unchanged meaning, now used as a subfolder under by_sim/<sim>/
or comparisons/ instead of a top-level dated root):
  01_static_FL, 02_main_analysis, 03_accessibility, 04_md_distances,
  05_clustering, 06_restraint_tests, 07_lysine_contacts, 08_kh_domains,
  09_surface_gallery, 99_misc

Call get_fig_dir(category, subname=None, sims=None):
  sims=['fl_go']              -> figures/by_sim/fl_go/<category>/<subname?>/
  sims=['fl', 'fl_optimized'] -> figures/comparisons/<category>/fl_fl_optimized/<subname?>/
  sims=None                   -> legacy dated layout (unchanged behavior)
"""
from datetime import date
from pathlib import Path

FIGURES_ROOT = Path(__file__).resolve().parent / 'figures'


def _sims_key(sims):
    """Stable, filesystem-safe folder name for a list of sim/set keys."""
    return '_'.join(sorted(str(s) for s in sims))


def get_fig_dir(category, subname=None, sims=None):
    """Return the output dir for a category (+ optional subname), creating it.

    sims: the set_key(s) this batch of figures belongs to.
      - None (default): legacy figures/<category>/<YYYY-MM-DD>/<subname?>/
        behavior, with a 'latest' symlink at the category root (unchanged).
      - one sim: figures/by_sim/<sim>/<category>/<subname?>/
      - multiple sims: figures/comparisons/<category>/<sims_joined>/<subname?>/
    """
    if sims:
        sims = list(sims)
        if len(sims) == 1:
            out_dir = FIGURES_ROOT / 'by_sim' / sims[0] / category
        else:
            out_dir = FIGURES_ROOT / 'comparisons' / category / _sims_key(sims)
        if subname:
            out_dir = out_dir / subname
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    # ---- legacy dated layout ----
    today = date.today().isoformat()
    cat_root = FIGURES_ROOT / category
    cat_root.mkdir(parents=True, exist_ok=True)
    out_dir = cat_root / today
    if subname:
        out_dir = out_dir / subname
    out_dir.mkdir(parents=True, exist_ok=True)

    latest = cat_root / 'latest'
    try:
        if latest.is_symlink() or latest.exists():
            latest.unlink()
    except OSError:
        pass
    try:
        latest.symlink_to(today)
    except OSError:
        pass

    return out_dir
