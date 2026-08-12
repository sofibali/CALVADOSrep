#!/usr/bin/env python3
"""
Generate SVG figures for KH domain PAE/pLDDT analysis.

Figures produced (all in figures/):
  1. PAE heatmaps for each construct (kh1-kh6 isolated + norrm)
  2. Per-residue pLDDT line plots with domain block overlays per threshold
  3. Inter-domain PAE heatmap (norrm construct)
  4. Summary bar chart comparing threshold sets
  5. Per-residue pLDDT + local PAE dual-axis line plot for all KH regions

Requires: analyze_kh_domains.py output in data/kh_pae_analysis_*.npz
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

_HERE = Path('/home/sbali/CALVADOS/examples/PARP14_MDP')
DATA_DIR = _HERE / 'data'
# Dated category subdir: figures/08_kh_domains/<YYYY-MM-DD>/
import sys as _sys
_sys.path.insert(0, str(_HERE))
import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
from _fig_layout import get_fig_dir as _get_fig_dir
FIG_DIR = _get_fig_dir('08_kh_domains')

# Domain definitions
DOMAIN_DEFS = {
    'rrm1':    (1, 145),
    'rrm2':    (146, 224),
    'rrm3':    (225, 314),
    'kh1-kh6': (315, 737),
    'kh7a':    (738, 789),
    'md1l1':   (790, 1004),
    'md2':     (1004, 1193),
    'md3':     (1207, 1388),
    'khb-kh8': (1389, 1533),
    'wwe':     (1534, 1602),
    'art':     (1603, 1801),
}

# Domain colors (consistent with other PARP14 figures)
DOMAIN_COLORS = {
    'rrm1': '#1f77b4', 'rrm2': '#2ca02c', 'rrm3': '#ff7f0e',
    'kh1-kh6': '#d62728', 'kh7a': '#9467bd',
    'md1l1': '#8c564b', 'md2': '#e377c2', 'md3': '#7f7f7f',
    'khb-kh8': '#bcbd22', 'wwe': '#17becf', 'art': '#ff9896',
}

# Threshold sets (must match analyze_kh_domains.py)
THRESHOLD_SETS = {
    'conservative': {'pae': 4.0, 'plddt': 70.0, 'merge': 0, 'color': '#2166ac'},
    'moderate': {'pae': 5.0, 'plddt': 60.0, 'merge': 0, 'color': '#4393c3'},
    'moderate_merged': {'pae': 5.0, 'plddt': 60.0, 'merge': 5, 'color': '#92c5de'},
    'permissive': {'pae': 6.0, 'plddt': 50.0, 'merge': 0, 'color': '#f4a582'},
    'permissive_merged': {'pae': 6.0, 'plddt': 50.0, 'merge': 10, 'color': '#d6604d'},
    'aggressive': {'pae': 8.0, 'plddt': 50.0, 'merge': 15, 'color': '#b2182b'},
}

# Current (original) domains
CURRENT_DOMAINS = [
    (6, 88), (150, 223), (227, 301),
    (791, 978), (1003, 1190), (1216, 1387),
    (1523, 1601), (1605, 1801),
]


def load_npz(name):
    """Load a construct's analysis data."""
    f = DATA_DIR / f'kh_pae_analysis_{name}.npz'
    if not f.exists():
        return None
    d = np.load(f, allow_pickle=True)
    boundaries = json.loads(str(d['boundaries']))
    return {
        'mean_pae': d['mean_pae'],
        'mean_plddt': d['mean_plddt'],
        'boundaries': {k: tuple(v) for k, v in boundaries.items()},
    }


def load_domain_yaml(name):
    """Load a proposed domains.yaml file."""
    f = DATA_DIR / f'domains_{name}.yaml'
    if not f.exists():
        return []
    domains = []
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('- ['):
                parts = line.replace('- [', '').replace(']', '').split(',')
                domains.append((int(parts[0]), int(parts[1])))
    return domains


def pae_cmap():
    """PAE colormap: dark blue (0) -> white (15) -> red (30)."""
    return LinearSegmentedColormap.from_list(
        'pae', ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
                '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    )


# ================================================================
# Figure 1: PAE heatmaps
# ================================================================
def fig_pae_heatmaps():
    """PAE heatmaps for kh1-kh6 (isolated) and norrm constructs."""
    constructs = [
        ('kh1-kh6', 'KH1-KH6 (isolated)'),
        ('norrm', 'No-RRM construct'),
    ]

    for cname, title in constructs:
        data = load_npz(cname)
        if data is None:
            print(f"  Skip {cname}: no data")
            continue

        pae = data['mean_pae']
        bounds = data['boundaries']
        n = pae.shape[0]

        fig, ax = plt.subplots(figsize=(10, 8.5))
        im = ax.imshow(pae, cmap=pae_cmap(), vmin=0, vmax=30, origin='upper',
                        aspect='equal', interpolation='nearest')

        # Domain boundary lines and labels
        for dname, (s, e) in bounds.items():
            ax.axhline(s - 0.5, color='k', lw=0.5, alpha=0.5)
            ax.axvline(s - 0.5, color='k', lw=0.5, alpha=0.5)
            ax.axhline(e - 0.5, color='k', lw=0.5, alpha=0.5)
            ax.axvline(e - 0.5, color='k', lw=0.5, alpha=0.5)
            # Label at midpoint
            mid = (s + e) / 2
            fl_s, fl_e = DOMAIN_DEFS[dname]
            label = f"{dname}\n({fl_s}-{fl_e})"
            ax.text(mid, -n * 0.03, label, ha='center', va='bottom',
                    fontsize=7, rotation=0,
                    color=DOMAIN_COLORS.get(dname, 'k'))

        cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label('Predicted Aligned Error (Å)', fontsize=10)

        ax.set_title(f'AF3 PAE — {title}\n(mean over 25 models)', fontsize=12)
        ax.set_xlabel('Residue (construct index)', fontsize=10)
        ax.set_ylabel('Residue (construct index)', fontsize=10)

        fig.tight_layout()
        out = FIG_DIR / f'pae_heatmap_{cname}.svg'
        fig.savefig(out, format='svg', bbox_inches='tight')
        fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {out}")


# ================================================================
# Figure 2: Per-residue pLDDT + local PAE for KH regions
# ================================================================
def fig_plddt_line_plots():
    """Per-residue pLDDT and local PAE line plots for all KH regions,
    with domain boundary overlays for each threshold set."""
    data = load_npz('norrm')
    if data is None:
        print("  Skip: no norrm data")
        return

    plddt = data['mean_plddt']
    pae = data['mean_pae']
    bounds = data['boundaries']

    # Build FL residue numbers for each position
    fl_residues = []
    domain_labels = []
    for dname, (s, e) in bounds.items():
        fl_s, fl_e = DOMAIN_DEFS[dname]
        for i in range(s, e):
            offset = i - s
            fl_residues.append(fl_s + offset)
            domain_labels.append(dname)
    fl_residues = np.array(fl_residues)

    # Compute local PAE (window=10)
    n = len(plddt)
    local_pae = np.zeros(n)
    for i in range(n):
        lo = max(0, i - 10)
        hi = min(n, i + 10 + 1)
        local_pae[i] = np.mean(pae[i, lo:hi])

    # ---- Figure 2a: Full construct pLDDT + local PAE ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8), sharex=True,
                                     gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.08})

    # Domain background shading
    for dname, (s, e) in bounds.items():
        color = DOMAIN_COLORS.get(dname, '#cccccc')
        ax1.axvspan(fl_residues[s], fl_residues[min(e-1, n-1)],
                    alpha=0.12, color=color, zorder=0)
        ax2.axvspan(fl_residues[s], fl_residues[min(e-1, n-1)],
                    alpha=0.12, color=color, zorder=0)
        mid_idx = (s + e) // 2
        if mid_idx < n:
            ax1.text(fl_residues[mid_idx], 100, dname, ha='center', va='bottom',
                     fontsize=7, color=color, fontweight='bold')

    # pLDDT line
    ax1.plot(fl_residues[:n], plddt[:n], color='#2166ac', lw=0.8, alpha=0.9)
    ax1.axhline(70, color='#d6604d', ls='--', lw=0.8, alpha=0.6, label='pLDDT=70')
    ax1.axhline(60, color='#b2182b', ls=':', lw=0.8, alpha=0.6, label='pLDDT=60')
    ax1.axhline(50, color='#67001f', ls='-.', lw=0.8, alpha=0.6, label='pLDDT=50')
    ax1.set_ylabel('pLDDT', fontsize=11)
    ax1.set_ylim(20, 102)
    ax1.legend(loc='lower left', fontsize=8, framealpha=0.8)
    ax1.set_title('AF3 per-residue pLDDT and local PAE — No-RRM construct (mean of 25 models)',
                   fontsize=12)

    # Local PAE line
    ax2.plot(fl_residues[:n], local_pae[:n], color='#d6604d', lw=0.8, alpha=0.9)
    ax2.axhline(4.0, color='#2166ac', ls='--', lw=0.8, alpha=0.6, label='PAE=4')
    ax2.axhline(5.0, color='#4393c3', ls=':', lw=0.8, alpha=0.6, label='PAE=5')
    ax2.axhline(6.0, color='#f4a582', ls='-.', lw=0.8, alpha=0.6, label='PAE=6')
    ax2.axhline(8.0, color='#b2182b', ls='--', lw=0.8, alpha=0.4, label='PAE=8')
    ax2.set_ylabel('Local PAE (Å)', fontsize=11)
    ax2.set_xlabel('Full-length residue number', fontsize=11)
    ax2.set_ylim(0, 20)
    ax2.invert_yaxis()
    ax2.legend(loc='lower left', fontsize=8, framealpha=0.8)

    fig.tight_layout()
    out = FIG_DIR / 'plddt_localpae_norrm.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")

    # ---- Figure 2b: KH zoom with threshold overlays ----
    kh_regions = [
        ('KH1-KH6', 315, 737),
        ('KH7a', 738, 789),
        ('KHb-KH8', 1389, 1533),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(16, 12),
                              gridspec_kw={'height_ratios': [4, 1, 2], 'hspace': 0.25})

    for ax_idx, (region_name, fl_start, fl_end) in enumerate(kh_regions):
        ax = axes[ax_idx]
        # Find indices in our arrays for this FL range
        mask = (fl_residues >= fl_start) & (fl_residues <= fl_end)
        if not np.any(mask):
            continue
        x = fl_residues[mask]
        y_plddt = plddt[:n][mask]
        y_lpae = local_pae[:n][mask]

        # pLDDT on left axis
        ax.plot(x, y_plddt, color='#2166ac', lw=1.0, label='pLDDT', zorder=5)
        ax.set_ylabel('pLDDT', fontsize=10, color='#2166ac')
        ax.set_ylim(20, 100)
        ax.tick_params(axis='y', labelcolor='#2166ac')

        # Local PAE on right axis
        ax2r = ax.twinx()
        ax2r.plot(x, y_lpae, color='#d6604d', lw=1.0, alpha=0.7, label='Local PAE', zorder=4)
        ax2r.set_ylabel('Local PAE (Å)', fontsize=10, color='#d6604d')
        ax2r.set_ylim(0, 18)
        ax2r.tick_params(axis='y', labelcolor='#d6604d')

        # Overlay restrained blocks for each threshold
        y_offsets = np.linspace(15, 28, len(THRESHOLD_SETS))
        for i_set, (set_name, tparams) in enumerate(THRESHOLD_SETS.items()):
            domains = load_domain_yaml(set_name)
            y_pos = y_offsets[i_set]
            for ds, de in domains:
                if de >= fl_start and ds <= fl_end:
                    ds_clip = max(ds, fl_start)
                    de_clip = min(de, fl_end)
                    ax.plot([ds_clip, de_clip], [y_pos, y_pos],
                            color=tparams['color'], lw=4, solid_capstyle='butt',
                            alpha=0.8, zorder=6)
            # Label on left
            ax.text(fl_start - 2, y_pos, set_name, ha='right', va='center',
                    fontsize=6, color=tparams['color'], fontweight='bold')

        # Threshold lines
        ax.axhline(70, color='#d6604d', ls=':', lw=0.5, alpha=0.4)
        ax.axhline(60, color='#b2182b', ls=':', lw=0.5, alpha=0.4)
        ax.axhline(50, color='#67001f', ls=':', lw=0.5, alpha=0.4)

        ax.set_title(f'{region_name} (FL {fl_start}-{fl_end})', fontsize=11, fontweight='bold')
        ax.set_xlabel('Full-length residue number', fontsize=10)

        # Combined legend
        from matplotlib.lines import Line2D
        handles = [
            Line2D([0], [0], color='#2166ac', lw=1.5, label='pLDDT'),
            Line2D([0], [0], color='#d6604d', lw=1.5, label='Local PAE'),
        ]
        ax.legend(handles=handles, loc='lower right', fontsize=8, framealpha=0.8)

    fig.suptitle('KH Domain Analysis — pLDDT + Local PAE with Restraint Block Overlays',
                  fontsize=13, fontweight='bold', y=1.01)
    fig.tight_layout()
    out = FIG_DIR / 'kh_zoom_threshold_overlays.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ================================================================
# Figure 3: Inter-domain PAE heatmap (norrm construct)
# ================================================================
def fig_interdomain_pae():
    """Block-averaged inter-domain PAE matrix."""
    data = load_npz('norrm')
    if data is None:
        print("  Skip: no norrm data")
        return

    pae = data['mean_pae']
    bounds = data['boundaries']

    domain_names = list(bounds.keys())
    n_dom = len(domain_names)
    block_pae = np.zeros((n_dom, n_dom))

    for i, d1 in enumerate(domain_names):
        s1, e1 = bounds[d1]
        for j, d2 in enumerate(domain_names):
            s2, e2 = bounds[d2]
            block_pae[i, j] = np.mean(pae[s1:e1, s2:e2])

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(block_pae, cmap=pae_cmap(), vmin=0, vmax=32, aspect='equal')

    # Annotate cells
    for i in range(n_dom):
        for j in range(n_dom):
            val = block_pae[i, j]
            color = 'white' if val > 20 or val < 5 else 'black'
            ax.text(j, i, f'{val:.1f}', ha='center', va='center',
                    fontsize=7, color=color, fontweight='bold')

    ax.set_xticks(range(n_dom))
    ax.set_xticklabels(domain_names, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(n_dom))
    ax.set_yticklabels(domain_names, fontsize=9)

    # Color domain labels
    for i, dname in enumerate(domain_names):
        color = DOMAIN_COLORS.get(dname, 'k')
        ax.get_xticklabels()[i].set_color(color)
        ax.get_yticklabels()[i].set_color(color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Mean PAE (Å)', fontsize=10)

    ax.set_title('Inter-domain PAE — No-RRM construct\n(block-averaged, mean of 25 AF3 models)',
                  fontsize=12)

    # Highlight KH domains
    kh_indices = [i for i, d in enumerate(domain_names) if 'kh' in d.lower()]
    for idx in kh_indices:
        rect = mpatches.FancyBboxPatch(
            (idx - 0.5, idx - 0.5), 1, 1,
            boxstyle="round,pad=0", linewidth=2,
            edgecolor='red', facecolor='none', zorder=10
        )
        ax.add_patch(rect)

    fig.tight_layout()
    out = FIG_DIR / 'interdomain_pae_norrm.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ================================================================
# Figure 4: Summary comparison bar chart
# ================================================================
def fig_summary_comparison():
    """Bar chart comparing threshold sets: restrained residues, domains, coverage."""
    set_names = list(THRESHOLD_SETS.keys())
    n_kh_res = []
    n_total_res = []
    n_domains = []
    colors = []

    for sn in set_names:
        domains = load_domain_yaml(sn)
        if not domains:
            n_kh_res.append(0)
            n_total_res.append(0)
            n_domains.append(0)
            colors.append('#999999')
            continue
        # KH residues = total - current
        current_res = sum(e - s + 1 for s, e in CURRENT_DOMAINS)
        total = sum(e - s + 1 for s, e in domains)
        n_total_res.append(total)
        n_kh_res.append(total - current_res)
        n_domains.append(len(domains))
        colors.append(THRESHOLD_SETS[sn]['color'])

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: KH restrained residues
    ax = axes[0]
    bars = ax.barh(range(len(set_names)), n_kh_res, color=colors, edgecolor='k', lw=0.5)
    ax.set_yticks(range(len(set_names)))
    ax.set_yticklabels(set_names, fontsize=9)
    ax.set_xlabel('KH restrained residues', fontsize=10)
    ax.set_title('New KH Restraints', fontsize=11, fontweight='bold')
    for i, v in enumerate(n_kh_res):
        ax.text(v + 5, i, str(v), va='center', fontsize=9)
    ax.invert_yaxis()

    # Panel 2: Total restrained residues (% of 1801)
    ax = axes[1]
    pcts = [100 * v / 1801 for v in n_total_res]
    # Add the "current" bar
    current_pct = 100 * sum(e - s + 1 for s, e in CURRENT_DOMAINS) / 1801
    all_pcts = [current_pct] + pcts
    all_names = ['current\n(no KH)'] + set_names
    all_colors = ['#999999'] + colors

    bars = ax.barh(range(len(all_names)), all_pcts, color=all_colors,
                    edgecolor='k', lw=0.5)
    ax.set_yticks(range(len(all_names)))
    ax.set_yticklabels(all_names, fontsize=9)
    ax.set_xlabel('% of FL restrained', fontsize=10)
    ax.set_title('Total Coverage (% of 1801 res)', fontsize=11, fontweight='bold')
    ax.set_xlim(0, 100)
    for i, v in enumerate(all_pcts):
        ax.text(v + 0.5, i, f'{v:.1f}%', va='center', fontsize=9)
    ax.invert_yaxis()

    # Panel 3: Number of domains
    ax = axes[2]
    all_ndom = [len(CURRENT_DOMAINS)] + n_domains
    bars = ax.barh(range(len(all_names)), all_ndom, color=all_colors,
                    edgecolor='k', lw=0.5)
    ax.set_yticks(range(len(all_names)))
    ax.set_yticklabels(all_names, fontsize=9)
    ax.set_xlabel('Number of restrained domains', fontsize=10)
    ax.set_title('Domain Count', fontsize=11, fontweight='bold')
    for i, v in enumerate(all_ndom):
        ax.text(v + 0.1, i, str(v), va='center', fontsize=9)
    ax.invert_yaxis()

    fig.suptitle('KH Domain Restraint Boundary Comparison', fontsize=13, fontweight='bold')
    fig.tight_layout()
    out = FIG_DIR / 'kh_threshold_comparison.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ================================================================
# Figure 5: Domain architecture schematic with restraint blocks
# ================================================================
def fig_domain_schematic():
    """Linear domain schematic showing which regions are restrained under each threshold."""
    set_names = list(THRESHOLD_SETS.keys())

    fig, ax = plt.subplots(figsize=(16, 6))

    y_spacing = 1.2
    n_rows = len(set_names) + 2  # +1 for current, +1 for domain labels

    # Row 0: Domain architecture bar
    y = (n_rows - 1) * y_spacing
    for dname, (fl_s, fl_e) in DOMAIN_DEFS.items():
        color = DOMAIN_COLORS.get(dname, '#cccccc')
        ax.barh(y, fl_e - fl_s + 1, left=fl_s, height=0.8,
                color=color, edgecolor='k', lw=0.5, alpha=0.7)
        mid = (fl_s + fl_e) / 2
        ax.text(mid, y, dname, ha='center', va='center', fontsize=6,
                fontweight='bold', color='white' if dname not in ['kh7a', 'wwe'] else 'black')
    ax.text(-10, y, 'FL domains', ha='right', va='center', fontsize=9, fontweight='bold')

    # Row 1: Current restraints
    y = (n_rows - 2) * y_spacing
    ax.barh(y, 1801, left=1, height=0.8, color='#f0f0f0', edgecolor='#cccccc', lw=0.3)
    for s, e in CURRENT_DOMAINS:
        ax.barh(y, e - s + 1, left=s, height=0.8,
                color='#999999', edgecolor='k', lw=0.5)
    ax.text(-10, y, 'current\n(no KH)', ha='right', va='center', fontsize=8, color='#666666')

    # Rows 2+: Each threshold set
    for i_set, sname in enumerate(set_names):
        y = (n_rows - 3 - i_set) * y_spacing
        domains = load_domain_yaml(sname)
        color = THRESHOLD_SETS[sname]['color']

        # Gray background
        ax.barh(y, 1801, left=1, height=0.8, color='#f0f0f0', edgecolor='#cccccc', lw=0.3)

        for s, e in domains:
            # Check if this is a current (non-KH) domain
            is_current = (s, e) in CURRENT_DOMAINS
            fc = '#999999' if is_current else color
            ax.barh(y, e - s + 1, left=s, height=0.8,
                    color=fc, edgecolor='k', lw=0.5)

        ax.text(-10, y, sname, ha='right', va='center', fontsize=8,
                color=color, fontweight='bold')

    ax.set_xlim(-5, 1810)
    ax.set_ylim(-y_spacing, (n_rows) * y_spacing)
    ax.set_xlabel('Full-length residue number', fontsize=11)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # KH region highlight
    for fl_s, fl_e, label in [(315, 737, 'KH1-KH6'), (738, 789, 'KH7a'), (1389, 1533, 'KHb-KH8')]:
        ax.axvspan(fl_s, fl_e, alpha=0.06, color='red', zorder=0)
        ax.text((fl_s + fl_e) / 2, -0.5 * y_spacing, label,
                ha='center', va='top', fontsize=7, color='red', fontstyle='italic')

    ax.set_title('PARP14 Restraint Domain Boundaries — Threshold Comparison',
                  fontsize=13, fontweight='bold')

    # Legend
    handles = [
        mpatches.Patch(facecolor='#999999', edgecolor='k', label='Existing restraints'),
        mpatches.Patch(facecolor='#4393c3', edgecolor='k', label='New KH restraints'),
        mpatches.Patch(facecolor='#f0f0f0', edgecolor='#cccccc', label='Unrestrained'),
    ]
    ax.legend(handles=handles, loc='upper right', fontsize=8, framealpha=0.9)

    fig.tight_layout()
    out = FIG_DIR / 'kh_domain_schematic.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ================================================================
# Figure 6: PAE heatmap zoomed on KH1-KH6 with block boundaries
# ================================================================
def fig_pae_kh_zoom():
    """Zoomed PAE heatmap of KH1-KH6 region with block boundaries for each threshold."""
    data = load_npz('norrm')
    if data is None:
        data = load_npz('kh1-kh6')
    if data is None:
        print("  Skip: no data for KH zoom")
        return

    pae = data['mean_pae']
    bounds = data['boundaries']

    if 'kh1-kh6' not in bounds:
        print("  Skip: kh1-kh6 not in boundaries")
        return

    s, e = bounds['kh1-kh6']
    kh_pae = pae[s:e, s:e]
    n = kh_pae.shape[0]
    fl_offset = DOMAIN_DEFS['kh1-kh6'][0]

    fig, ax = plt.subplots(figsize=(10, 9))
    im = ax.imshow(kh_pae, cmap=pae_cmap(), vmin=0, vmax=20, origin='upper',
                    aspect='equal', interpolation='nearest')

    # Overlay block boundaries for moderate threshold
    moderate_domains = load_domain_yaml('moderate')
    for ds, de in moderate_domains:
        if ds >= fl_offset and de <= fl_offset + n:
            # Convert to local coordinates
            ls = ds - fl_offset
            le = de - fl_offset
            rect = mpatches.Rectangle(
                (ls - 0.5, ls - 0.5), le - ls + 1, le - ls + 1,
                linewidth=2, edgecolor='lime', facecolor='none', zorder=10,
                linestyle='-'
            )
            ax.add_patch(rect)

    # Also overlay conservative (dashed)
    conservative_domains = load_domain_yaml('conservative')
    for ds, de in conservative_domains:
        if ds >= fl_offset and de <= fl_offset + n:
            ls = ds - fl_offset
            le = de - fl_offset
            rect = mpatches.Rectangle(
                (ls - 0.5, ls - 0.5), le - ls + 1, le - ls + 1,
                linewidth=1.5, edgecolor='cyan', facecolor='none', zorder=9,
                linestyle='--'
            )
            ax.add_patch(rect)

    # Tick labels in FL numbering
    tick_positions = np.arange(0, n, 50)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([str(fl_offset + t) for t in tick_positions], fontsize=8)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels([str(fl_offset + t) for t in tick_positions], fontsize=8)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Predicted Aligned Error (Å)', fontsize=10)

    ax.set_xlabel('Full-length residue number', fontsize=10)
    ax.set_ylabel('Full-length residue number', fontsize=10)
    ax.set_title('PAE Heatmap — KH1-KH6 Region\n(No-RRM construct, mean of 25 AF3 models)',
                  fontsize=12)

    # Legend
    handles = [
        mpatches.Patch(edgecolor='lime', facecolor='none', lw=2, label='Moderate blocks'),
        mpatches.Patch(edgecolor='cyan', facecolor='none', lw=1.5, linestyle='--',
                       label='Conservative blocks'),
    ]
    ax.legend(handles=handles, loc='upper right', fontsize=9, framealpha=0.9)

    fig.tight_layout()
    out = FIG_DIR / 'pae_kh16_zoom.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ================================================================
def main():
    print("Generating KH domain analysis figures...")
    print()

    print("Figure 1: PAE heatmaps")
    fig_pae_heatmaps()
    print()

    print("Figure 2: pLDDT + local PAE line plots")
    fig_plddt_line_plots()
    print()

    print("Figure 3: Inter-domain PAE heatmap")
    fig_interdomain_pae()
    print()

    print("Figure 4: Summary comparison")
    fig_summary_comparison()
    print()

    print("Figure 5: Domain schematic")
    fig_domain_schematic()
    print()

    print("Figure 6: PAE KH1-KH6 zoom")
    fig_pae_kh_zoom()
    print()

    print("All figures saved to:", FIG_DIR)


if __name__ == '__main__':
    main()
