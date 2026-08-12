#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Assemble the separate per-state, per-view PNGs produced by
make_surface_gallery.py / run_surface_gallery_batch.py into multi-panel
summary figures: one states x whole-molecule-views grid per run, plus one
states x pocket-views grid per active site present.

Reads whatever states/views actually exist on disk -- states or views a
prior batch run skipped (failed APBS prep, etc.) are just left blank with
a "missing" label, never silently dropped from the grid.

Usage:
    python make_surface_gallery_montage.py --run-dir figures/09_surface_gallery/<date>/<run_tag>
    python make_surface_gallery_montage.py --auto   # do every run_tag under the latest date
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

CWD = Path(__file__).resolve().parent
ROOT = CWD.parent
sys.path.insert(0, str(ROOT))
from _fig_layout import get_fig_dir  # noqa: E402

WHOLE_VIEWS = [
    ('domain_cartoon.png', 'Domain cartoon'),
    ('electrostatic_surface_front.png', 'Electrostatic potential'),
    ('hydrophobic_surface_front.png', 'Hydrophobicity'),
    ('face_classification_front.png', 'Face classification'),
]
POCKET_VIEWS = [
    ('electrostatic.png', 'Electrostatic potential'),
    ('hydrophobicity.png', 'Hydrophobicity'),
    ('residue_class.png', 'Residue chemistry'),
    ('hbond_donor_acceptor.png', 'H-bond donor/acceptor'),
]


def state_label_sort_key(p):
    try:
        return int(p.name.split('_')[1])
    except (IndexError, ValueError):
        return p.name


def find_states(run_dir):
    return sorted((d for d in run_dir.iterdir() if d.is_dir() and d.name.startswith('state_')),
                  key=state_label_sort_key)


def find_pocket_sites(state_dirs):
    sites = set()
    for sd in state_dirs:
        for p in sd.glob('pocket_*'):
            if p.is_dir():
                sites.add(p.name.replace('pocket_', ''))
    return sorted(sites)


def plot_grid(rows, row_labels, col_specs, title, out_path, panel_size=3.2):
    n_rows, n_cols = len(rows), len(col_specs)
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(panel_size * n_cols, panel_size * n_rows),
                              squeeze=False)
    for i, (row_dir, row_label) in enumerate(zip(rows, row_labels)):
        for j, (fname, col_label) in enumerate(col_specs):
            ax = axes[i][j]
            img_path = row_dir / fname
            if img_path.exists():
                ax.imshow(mpimg.imread(img_path))
            else:
                ax.text(0.5, 0.5, 'missing', ha='center', va='center',
                        fontsize=10, color='gray', transform=ax.transAxes)
                ax.set_facecolor('#f0f0f0')
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if i == 0:
                ax.set_title(col_label, fontsize=11)
            if j == 0:
                ax.set_ylabel(row_label, fontsize=11, rotation=90)
    fig.suptitle(title, fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    fig.savefig(str(out_path).replace('.png', '.svg'), bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out_path}")


def assemble_run(run_dir, out_dir):
    run_dir = Path(run_dir)
    state_dirs = find_states(run_dir)
    if not state_dirs:
        print(f"  no state_* dirs under {run_dir}, skipping")
        return
    row_labels = [d.name for d in state_dirs]
    run_tag = run_dir.name

    plot_grid(state_dirs, row_labels, WHOLE_VIEWS,
              f'{run_tag}: whole-molecule surface gallery',
              out_dir / f'{run_tag}_whole_molecule.png')

    for site in find_pocket_sites(state_dirs):
        pocket_dirs = [d / f'pocket_{site}' for d in state_dirs]
        plot_grid(pocket_dirs, row_labels, POCKET_VIEWS,
                  f'{run_tag}: {site} pocket',
                  out_dir / f'{run_tag}_pocket_{site}.png')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-dir', action='append', default=[],
                     help='figures/09_surface_gallery/<date>/<run_tag> (repeatable)')
    ap.add_argument('--auto', action='store_true',
                     help='assemble every run_tag under the latest 09_surface_gallery date')
    args = ap.parse_args()

    run_dirs = [Path(p) for p in args.run_dir]
    if args.auto:
        gallery_root = ROOT / 'figures' / '09_surface_gallery'
        latest = gallery_root / 'latest'
        if latest.is_dir():
            run_dirs += [d for d in latest.iterdir() if d.is_dir()]
    if not run_dirs:
        print("Nothing to do: pass --run-dir or --auto")
        return 1

    out_dir = get_fig_dir('09_surface_gallery', 'montages')
    for run_dir in run_dirs:
        print(f"Assembling {run_dir.name}...")
        assemble_run(run_dir, out_dir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
