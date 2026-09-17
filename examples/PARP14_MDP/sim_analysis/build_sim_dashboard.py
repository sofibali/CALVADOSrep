#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a per-simulation dashboard: one page with the sim's parameters at the
top and every figure that's been generated for it in a grid below.

Every analysis script in this project (analyze_all.py, analyze_lys_contacts.py,
analyze_active_sites.py, analyze_accessibility.py, figure_md_distances.py,
figure_lys_exposed_persistence.py, tica_pipeline/cluster_states.py) now writes
into figures/by_sim/<sim>/<category>/[<subname>/]<file>.png (+ matching .svg)
when invoked with a single --set/--sim-folder. This script doesn't generate
any new analysis -- it just scans that tree, reads the sim's own config to
build a parameter header, and renders:
  figures/by_sim/<sim>/dashboard.html     -- browsable page, PNGs inline,
                                              each with a link to its SVG
  figures/by_sim/<sim>/dashboard_montage.png -- single-file static export
                                              tiling every figure (for
                                              sharing/printing without a
                                              browser)

Usage:
    python build_sim_dashboard.py --sim fl_go
    python build_sim_dashboard.py --sim fl_optimized --sim-dir ../fl_optimized
"""
import argparse
import glob
import html
import os
import sys
from pathlib import Path

CWD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CWD))
import sim_registry as reg  # noqa: E402
from _fig_layout import FIGURES_ROOT  # noqa: E402

try:
    import yaml
except Exception:
    yaml = None


def find_sim_dir(sim, sim_dir_arg):
    if sim_dir_arg:
        return Path(sim_dir_arg).resolve()
    guess = CWD / sim
    if guess.is_dir():
        return guess
    if sim.startswith('frag_') and (CWD / 'fragments' / sim[5:]).is_dir():
        return CWD / 'fragments' / sim[5:]
    raise FileNotFoundError(
        f"can't locate the simulation folder for '{sim}' (tried {guess}); "
        f"pass --sim-dir explicitly")


def find_replicate_dirs(sim_dir):
    """Grid (seed-N_sample-M) or flat (any dir with a .dcd) replicate dirs."""
    grid = sorted(sim_dir.glob('seed-*_sample-*'))
    if grid:
        return grid
    return sorted(d for d in sim_dir.iterdir() if d.is_dir() and list(d.glob('*.dcd')))


def load_metadata(sim, sim_dir):
    """Collect whatever's available: metadata.json, units, replicate count,
    and the actual config.yaml/components.yaml of one representative
    replicate (sim parameters as actually run, not assumed defaults)."""
    meta = {}
    mfile = sim_dir / 'metadata.json'
    if mfile.is_file():
        import json
        meta.update(json.load(open(mfile)))

    reps = find_replicate_dirs(sim_dir)
    meta['n_replicates'] = len(reps)

    if reps and yaml is not None:
        cfg_file = reps[0] / 'config.yaml'
        comp_file = reps[0] / 'components.yaml'
        if cfg_file.is_file():
            cfg = yaml.safe_load(open(cfg_file)) or {}
            for k in ('box', 'steps', 'wfreq', 'temp', 'ionic', 'eps_lj',
                      'cutoff_lj', 'cutoff_yu', 'custom_restraints',
                      'custom_restraint_type', 'restart', 'sysname'):
                if k in cfg:
                    meta[k] = cfg[k]
        if comp_file.is_file():
            comp = yaml.safe_load(open(comp_file)) or {}
            defaults = comp.get('defaults', {})
            for k in ('restraint_type', 'k_harmonic', 'k_go', 'cutoff_restr',
                      'colabfold'):
                if k in defaults:
                    meta[k] = defaults[k]

    if 'units' not in meta:
        try:
            meta['units'] = reg._resolve_units(str(sim_dir))
        except Exception:
            pass

    return meta


def format_metadata_rows(meta):
    order = ['label', 'units', 'sysname', 'n_replicates', 'box', 'steps',
             'wfreq', 'temp', 'ionic', 'eps_lj', 'cutoff_lj', 'cutoff_yu',
             'restraint_type', 'k_harmonic', 'custom_restraints',
             'custom_restraint_type', 'k_go', 'cutoff_restr', 'restart']
    labels = {
        'label': 'Label', 'units': 'Domain units', 'sysname': 'System name',
        'n_replicates': 'Replicates found', 'box': 'Box (nm)',
        'steps': 'Steps', 'wfreq': 'Write frequency', 'temp': 'Temperature (K)',
        'ionic': 'Ionic strength (M)', 'eps_lj': 'AH epsilon (kJ/mol)',
        'cutoff_lj': 'AH cutoff (nm)', 'cutoff_yu': 'YU cutoff (nm)',
        'restraint_type': 'Restraint type', 'k_harmonic': 'k_harmonic (kJ/mol/nm^2)',
        'custom_restraints': 'Custom restraints', 'custom_restraint_type': 'Custom restraint type',
        'k_go': 'k_go (kJ/mol)', 'cutoff_restr': 'Restraint cutoff (nm)',
        'restart': 'Restart mode',
    }
    rows = []
    for k in order:
        if k in meta and meta[k] not in (None, '', []):
            v = meta[k]
            if isinstance(v, list):
                v = ', '.join(str(x) for x in v)
            rows.append((labels.get(k, k), v))
    return rows


def scan_figures(sim_dir_name):
    """{category: {subname_or_None: [png paths]}} under figures/by_sim/<sim>/."""
    root = FIGURES_ROOT / 'by_sim' / sim_dir_name
    tree = {}
    if not root.is_dir():
        return tree, root
    for png in sorted(root.rglob('*.png')):
        rel = png.relative_to(root)
        parts = rel.parts
        category = parts[0]
        subname = str(Path(*parts[1:-1])) if len(parts) > 2 else None
        tree.setdefault(category, {}).setdefault(subname, []).append(png)
    return tree, root


CATEGORY_TITLES = {
    '02_main_analysis': 'Main Analysis (Rg, distance/contact maps, energy, WCN, active sites, accessibility)',
    '03_accessibility': 'Steric Accessibility',
    '04_md_distances': 'Inter-Domain Distances',
    '05_clustering': 'TICA Clustering',
    '07_lysine_contacts': 'Lysine Contacts',
    '09_surface_gallery': 'Surface Gallery',
}


def build_html(sim, meta, tree, root, out_path):
    rows = format_metadata_rows(meta)
    meta_html = '\n'.join(
        f'<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>'
        for k, v in rows)

    sections = []
    for category in sorted(tree):
        title = CATEGORY_TITLES.get(category, category)
        sections.append(f'<h2>{html.escape(title)}</h2>')
        for subname in sorted(tree[category], key=lambda s: (s is None, s)):
            pngs = tree[category][subname]
            if subname:
                sections.append(f'<h3>{html.escape(subname)}</h3>')
            sections.append('<div class="grid">')
            for png in pngs:
                rel_png = png.relative_to(root)
                svg = png.with_suffix('.svg')
                card = [f'<figure><a href="{rel_png}" target="_blank">',
                        f'<img src="{rel_png}" loading="lazy" alt="{html.escape(png.stem)}">',
                        '</a><figcaption>', html.escape(png.stem)]
                if svg.is_file():
                    rel_svg = svg.relative_to(root)
                    card.append(f' &middot; <a href="{rel_svg}" download>SVG</a>')
                card.append(f' &middot; <a href="{rel_png}" download>PNG</a>')
                card.append('</figcaption></figure>')
                sections.append(''.join(card))
            sections.append('</div>')

    n_figs = sum(len(v) for d in tree.values() for v in d.values())
    html_doc = f"""<!doctype html>
<html><head><meta charset="utf-8">
<title>{html.escape(sim)} dashboard</title>
<style>
  body {{ font-family: -apple-system, Helvetica, Arial, sans-serif; margin: 24px;
         background: #fafafa; color: #1a1a1a; }}
  h1 {{ margin-bottom: 4px; }}
  .subtitle {{ color: #666; margin-top: 0; }}
  table {{ border-collapse: collapse; margin: 16px 0 32px; }}
  th, td {{ text-align: left; padding: 4px 12px; border-bottom: 1px solid #ddd; font-size: 14px; }}
  th {{ color: #555; font-weight: 600; white-space: nowrap; }}
  h2 {{ border-bottom: 2px solid #ccc; padding-bottom: 4px; margin-top: 40px; }}
  h3 {{ color: #444; margin-top: 24px; }}
  .grid {{ display: flex; flex-wrap: wrap; gap: 16px; margin-bottom: 16px; }}
  figure {{ margin: 0; background: white; border: 1px solid #ddd; border-radius: 6px;
           padding: 8px; width: 360px; }}
  figure img {{ width: 100%; height: auto; border-radius: 3px; }}
  figcaption {{ font-size: 12px; color: #555; margin-top: 6px; word-break: break-word; }}
  figcaption a {{ color: #06c; text-decoration: none; }}
  figcaption a:hover {{ text-decoration: underline; }}
  @media (prefers-color-scheme: dark) {{
    body {{ background: #16181c; color: #e6e6e6; }}
    figure {{ background: #22252b; border-color: #3a3d44; }}
    th, td {{ border-color: #3a3d44; }}
    th {{ color: #aaa; }} figcaption {{ color: #aaa; }}
    figcaption a {{ color: #6cf; }}
    h2 {{ border-color: #444; }}
  }}
</style></head>
<body>
<h1>{html.escape(sim)}</h1>
<p class="subtitle">{n_figs} figures &middot; static export: <a href="dashboard_montage.png" download>dashboard_montage.png</a></p>
<table>{meta_html}</table>
{"".join(sections) if sections else "<p><em>No figures found yet under figures/by_sim/" + html.escape(sim) + "/. Run an analysis script with --sim-folder first.</em></p>"}
</body></html>
"""
    out_path.write_text(html_doc)
    return n_figs


def build_montage(tree, root, out_path, thumb=(320, 240), cols=5):
    from PIL import Image, ImageDraw, ImageFont
    pngs = [p for d in tree.values() for lst in d.values() for p in lst]
    if not pngs:
        return None
    pad, label_h = 10, 20
    cell_w, cell_h = thumb[0] + pad, thumb[1] + label_h + pad
    n = len(pngs)
    rows = (n + cols - 1) // cols
    canvas = Image.new('RGB', (cols * cell_w + pad, rows * cell_h + pad), 'white')
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.load_default()
    except Exception:
        font = None
    for i, p in enumerate(pngs):
        r, c = divmod(i, cols)
        x, y = pad + c * cell_w, pad + r * cell_h
        try:
            im = Image.open(p).convert('RGB')
            im.thumbnail(thumb)
            ox = x + (thumb[0] - im.width) // 2
            oy = y + (thumb[1] - im.height) // 2
            canvas.paste(im, (ox, oy))
        except Exception:
            pass
        label = p.stem[:40]
        draw.text((x, y + thumb[1] + 2), label, fill='black', font=font)
    canvas.save(out_path)
    return out_path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sim', required=True, help='Sim/set key (folder basename under figures/by_sim/)')
    ap.add_argument('--sim-dir', default=None, help='Override: path to the simulation folder (default: <root>/<sim>)')
    args = ap.parse_args()

    try:
        sim_dir = find_sim_dir(args.sim, args.sim_dir)
        meta = load_metadata(args.sim, sim_dir)
    except FileNotFoundError as e:
        print(f"  WARN: {e}")
        meta = {}

    tree, root = scan_figures(args.sim)
    if not root.is_dir():
        root.mkdir(parents=True, exist_ok=True)

    out_html = root / 'dashboard.html'
    n_figs = build_html(args.sim, meta, tree, root, out_html)
    print(f"  {n_figs} figures found under {root}")
    print(f"  Saved: {out_html}")

    out_montage = root / 'dashboard_montage.png'
    result = build_montage(tree, root, out_montage)
    if result:
        print(f"  Saved: {out_montage}")
    else:
        print("  (no figures yet -- skipped montage)")


if __name__ == '__main__':
    main()
