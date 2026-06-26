"""Tiny helper for the dated figures/<category>/<date>/ layout.

Used by all analysis scripts. Each script picks its category (see README.md
in figures/) and gets back a Path object that's already created. Also
updates a 'latest' symlink in the category root.

Categories (numbered for sort order):
  01_static_FL          (figure_sasa_faces.py)
  02_main_analysis      (analyze_all.py — per-set + cross-set)
  03_accessibility      (analyze_accessibility.py, figure_accessibility.py)
  04_md_distances       (figure_md_distances.py)
  05_clustering         (cluster_states.py)
  06_restraint_tests    (lives under restraint_tests/summary_plots/)
  07_lysine_contacts    (analyze_lys_contacts.py)
  08_kh_domains         (figure_kh_domains.py, analyze_kh_domains.py)
  99_misc               (figure_schematic.py, one-offs)
"""
from datetime import date
from pathlib import Path

FIGURES_ROOT = Path(__file__).resolve().parent / 'figures'


def get_fig_dir(category, subname=None):
    """Return figures/<category>/<YYYY-MM-DD>/<subname?>/, creating dirs.

    Also maintains figures/<category>/latest symlink.
    """
    today = date.today().isoformat()
    cat_root = FIGURES_ROOT / category
    cat_root.mkdir(parents=True, exist_ok=True)
    out_dir = cat_root / today
    if subname:
        out_dir = out_dir / subname
    out_dir.mkdir(parents=True, exist_ok=True)

    # 'latest' symlink at the category root
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
