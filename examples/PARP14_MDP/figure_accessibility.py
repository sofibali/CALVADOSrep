#!/usr/bin/env python3
"""
Publication figures for PARP14 active site accessibility analysis.

Generates:
  Fig 1: Overview panel — SAA heatmap + cone angle + shell density (3-panel)
  Fig 2: Per-site SAA violin+swarm across constructs
  Fig 3: All metrics by construct — violin + points + SD
  Fig 4: Domain context — line plot (sites across constructs by chain length)
  Fig 5: Docking feasibility — cone angle vs SAA scatter
  Fig 6: Schematic — what the metrics measure
  Fig 7: Site ranking summary

Usage:
    conda run -n calvados python figure_accessibility.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle
from matplotlib.lines import Line2D
import os

CWD = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(CWD, 'data')
FIG_PATH = os.path.join(CWD, 'figures')
os.makedirs(FIG_PATH, exist_ok=True)

# ============================================================
# Load data
# ============================================================

data = np.load(os.path.join(DATA_PATH, 'accessibility_stats.npz'))

SETS = ['fl', 'norrm', 'noart', 'core', 'mka', 'md', 'md3art', 'fl_optimized']
SET_LABELS = {
    'fl':            'Full Length',
    'norrm':         'KH1-ART',
    'noart':         'KH1-WWE',
    'core':          'KH7-ART',
    'mka':           'MD1-ART',
    'md':            'MD1-MD3',
    'md3art':        'MD3-ART',
    'fl_optimized':  'FL (optimized)',
}
SET_LABELS_SHORT = {
    'fl':            'FL (1-1801)',
    'norrm':         'KH1-ART',
    'noart':         'KH1-WWE',
    'core':          'KH7-ART',
    'mka':           'MD1-ART',
    'md':            'MD1-MD3',
    'md3art':        'MD3-ART',
    'fl_optimized':  'FL optim.',
}
SET_COLORS = {
    'fl':            '#1f77b4',
    'norrm':         '#9467bd',
    'noart':         '#8c564b',
    'core':          '#2ca02c',
    'mka':           '#d62728',
    'md':            '#ff7f0e',
    'md3art':        '#17becf',
    'fl_optimized':  '#aec7e8',
}

SITES = ['MD1', 'MD2', 'MD3', 'ART']
SITE_COLORS = {'MD1': '#e6194b', 'MD2': '#3cb44b', 'MD3': '#4363d8', 'ART': '#f58231'}

CONSTRUCT_SITES = {
    'fl':            ['MD1', 'MD2', 'MD3', 'ART'],
    'norrm':         ['MD1', 'MD2', 'MD3', 'ART'],
    'noart':         ['MD1', 'MD2', 'MD3'],
    'core':          ['MD1', 'MD2', 'MD3', 'ART'],
    'mka':           ['MD1', 'MD2', 'MD3', 'ART'],
    'md':            ['MD1', 'MD2', 'MD3'],
    'md3art':        ['MD3', 'ART'],
    'fl_optimized':  ['MD1', 'MD2', 'MD3', 'ART'],
}


def get_vals(set_key, site, metric):
    """Retrieve data array, return empty if missing."""
    key = f'{set_key}_{site}_{metric}'
    if key in data:
        return data[key]
    return np.array([])


# ============================================================
# Shared style
# ============================================================

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.dpi': 200,
})


def add_panel_label(ax, label, x=-0.12, y=1.08, fontsize=16):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', va='top')


def swarmplot(ax, data_list, positions, colors, width=0.3, size=18, alpha=0.7):
    """Simple swarm-style jittered scatter."""
    rng = np.random.default_rng(42)
    for pos, vals, color in zip(positions, data_list, colors):
        if len(vals) == 0:
            continue
        jitter = rng.normal(0, width * 0.15, len(vals))
        ax.scatter(pos + jitter, vals, c=color, s=size, alpha=alpha,
                   edgecolors='black', linewidths=0.3, zorder=3)


def violin_with_points(ax, data_list, positions, colors, width=0.6, point_size=18):
    """Draw violin body + jitter points + mean line + SD error bar."""
    rng = np.random.default_rng(42)
    for pos, vals, color in zip(positions, data_list, colors):
        if len(vals) == 0:
            continue
        # Violin body
        parts = ax.violinplot([vals], positions=[pos], widths=width,
                              showmeans=False, showmedians=False, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(color)
            pc.set_alpha(0.35)
            pc.set_edgecolor(color)
            pc.set_linewidth(1)

        # Jitter points
        jitter = rng.normal(0, width * 0.08, len(vals))
        ax.scatter(pos + jitter, vals, c=color, s=point_size, alpha=0.7,
                   edgecolors='black', linewidths=0.3, zorder=3)

        # Mean line + SD bar
        m = np.mean(vals)
        sd = np.std(vals)
        ax.plot([pos - 0.15, pos + 0.15], [m, m], color='black', lw=2.5, zorder=5)
        ax.errorbar(pos, m, yerr=sd, color='black', lw=1.5, capsize=5,
                    capthick=1.5, zorder=4, fmt='none')


# ============================================================
# Figure 1: Overview summary (3-panel heatmap)
# ============================================================

def make_fig1():
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))

    metrics = [
        ('saa', 'Solid-Angle\nAccessibility', 'RdYlGn', 0, 0.8, '.2f', ''),
        ('cone', 'Max Approach\nCone Angle', 'RdYlGn', 30, 100, '.0f', '\u00b0'),
        ('shell', 'Shell Density\n(2\u20135 nm)', 'RdYlGn_r', 0, 0.04, '.3f', ''),
    ]

    for ax_idx, (metric, title, cmap, vmin, vmax, fmt, suffix) in enumerate(metrics):
        ax = axes[ax_idx]
        matrix = np.full((len(SETS), len(SITES)), np.nan)
        for si, sk in enumerate(SETS):
            for sj, sn in enumerate(SITES):
                v = get_vals(sk, sn, metric)
                if len(v) > 0:
                    matrix[si, sj] = np.mean(v)

        im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
        ax.set_yticks(range(len(SETS)))
        ax.set_yticklabels([SET_LABELS[s] for s in SETS])
        ax.set_xticks(range(len(SITES)))
        ax.set_xticklabels(SITES, fontweight='bold')

        for i in range(len(SETS)):
            for j in range(len(SITES)):
                v = matrix[i, j]
                if np.isnan(v):
                    ax.text(j, i, '\u2014', ha='center', va='center', fontsize=11, color='gray')
                else:
                    txt = f'{v:{fmt}}{suffix}'
                    norm_v = (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
                    if '_r' in cmap:
                        norm_v = 1 - norm_v
                    tc = 'white' if norm_v < 0.4 else 'black'
                    ax.text(j, i, txt, ha='center', va='center', fontsize=11,
                            fontweight='bold', color=tc)

        plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        ax.set_title(title, fontsize=12, fontweight='bold')
        add_panel_label(ax, chr(65 + ax_idx))

    fig.suptitle('Active Site Steric Accessibility Across PARP14 Constructs',
                 fontsize=14, fontweight='bold', y=1.06)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'fig1_accessibility_overview.png'))
    fig.savefig(os.path.join(FIG_PATH, 'fig1_accessibility_overview.svg'))
    plt.close()
    print("  Fig 1: fig1_accessibility_overview")


# ============================================================
# Figure 2: Per-site SAA violin+swarm across constructs
# ============================================================

def make_fig2():
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5), sharey=True)

    for si, sname in enumerate(SITES):
        ax = axes[si]
        set_data = []
        set_labels = []
        set_colors = []

        for sk in SETS:
            v = get_vals(sk, sname, 'saa')
            if len(v) > 0:
                set_data.append(v)
                set_labels.append(SET_LABELS_SHORT[sk])
                set_colors.append(SET_COLORS[sk])

        if not set_data:
            ax.set_title(f'{sname}\n(not in construct)', fontsize=11)
            add_panel_label(ax, chr(65 + si))
            continue

        positions = np.arange(len(set_data))
        violin_with_points(ax, set_data, positions, set_colors, width=0.6, point_size=22)

        ax.set_xticks(positions)
        ax.set_xticklabels(set_labels)
        ax.set_title(f'{sname} Active Site', fontsize=12, fontweight='bold',
                     color=SITE_COLORS[sname])
        if si == 0:
            ax.set_ylabel('Solid-Angle Accessibility')

        ax.axhline(0.5, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax.axhline(0.3, color='red', ls=':', lw=0.8, alpha=0.5)
        ax.set_ylim(0, 0.85)
        add_panel_label(ax, chr(65 + si))

    axes[-1].plot([], [], color='gray', ls='--', lw=0.8, label='50% open')
    axes[-1].plot([], [], color='red', ls=':', lw=0.8, label='30% open')
    axes[-1].legend(loc='upper right', fontsize=8)

    fig.suptitle('Solid-Angle Accessibility per Active Site\n'
                 '(fraction of approach directions unblocked by own chain)',
                 fontsize=13, fontweight='bold', y=1.08)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'fig2_saa_per_site.png'))
    fig.savefig(os.path.join(FIG_PATH, 'fig2_saa_per_site.svg'))
    plt.close()
    print("  Fig 2: fig2_saa_per_site")


# ============================================================
# Figure 3: All metrics by construct — violin + points + SD
# ============================================================

def make_fig3():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

    metric_info = [
        ('saa', 'Solid-Angle Accessibility', (0, 0.85)),
        ('cone', 'Max Cone Angle (\u00b0)', (0, 120)),
        ('shell', 'Shell Density (2\u20135 nm)', (0, 0.050)),
    ]

    for ax_idx, (metric, ylabel, ylim) in enumerate(metric_info):
        ax = axes[ax_idx]
        n_sets = len(SETS)
        n_sites = len(SITES)
        group_width = 0.8
        violin_w = group_width / n_sites * 0.85

        x_base = np.arange(n_sets)

        for j, sname in enumerate(SITES):
            offset = (j - n_sites / 2 + 0.5) * (group_width / n_sites)
            positions_j = []
            data_j = []
            colors_j = []

            for i, sk in enumerate(SETS):
                v = get_vals(sk, sname, metric)
                pos = x_base[i] + offset
                if len(v) > 0:
                    positions_j.append(pos)
                    data_j.append(v)
                    colors_j.append(SITE_COLORS[sname])

            if data_j:
                violin_with_points(ax, data_j, positions_j, colors_j,
                                   width=violin_w, point_size=12)

        ax.set_xticks(x_base)
        ax.set_xticklabels([SET_LABELS_SHORT[s] for s in SETS])
        ax.set_ylabel(ylabel)
        ax.set_ylim(ylim)

        # Site legend
        handles = [Line2D([0], [0], marker='o', color=SITE_COLORS[s], lw=0,
                          markersize=7, markeredgecolor='black', label=s)
                   for s in SITES]
        ax.legend(handles=handles, fontsize=9, ncol=2)

        if metric == 'saa':
            ax.axhline(0.5, color='gray', ls='--', lw=0.8, alpha=0.5)
            ax.axhline(0.3, color='red', ls=':', lw=0.8, alpha=0.4)
        elif metric == 'cone':
            ax.axhline(90, color='gray', ls='--', lw=0.8, alpha=0.5)

        add_panel_label(ax, chr(65 + ax_idx))

    fig.suptitle('Active Site Accessibility: All Metrics by Construct',
                 fontsize=14, fontweight='bold', y=1.03)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'fig3_metrics_by_construct.png'))
    fig.savefig(os.path.join(FIG_PATH, 'fig3_metrics_by_construct.svg'))
    plt.close()
    print("  Fig 3: fig3_metrics_by_construct")


# ============================================================
# Figure 4: Domain context effect — line plot with SD bands
# ============================================================

def make_fig4():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    metric_info = [
        ('saa', 'Solid-Angle Accessibility'),
        ('cone', 'Max Cone Angle (\u00b0)'),
        ('shell', 'Shell Density (2\u20135 nm)'),
    ]

    construct_order = ['md3art', 'md', 'mka', 'core', 'noart', 'norrm', 'fl', 'fl_optimized']
    construct_sizes = {'md3art': 595, 'md': 586, 'mka': 999, 'core': 1051,
                       'noart': 1275, 'norrm': 1474, 'fl': 1801, 'fl_optimized': 1801}
    x_labels = [f'{SET_LABELS_SHORT[s]}\n({construct_sizes[s]} res)' for s in construct_order]

    for ax_idx, (metric, ylabel) in enumerate(metric_info):
        ax = axes[ax_idx]

        for sname in SITES:
            means = []
            sds = []
            valid_x = []

            for i, sk in enumerate(construct_order):
                v = get_vals(sk, sname, metric)
                if len(v) > 0:
                    means.append(np.mean(v))
                    sds.append(np.std(v))
                    valid_x.append(i)

            if means:
                means = np.array(means)
                sds = np.array(sds)
                valid_x = np.array(valid_x)

                # Line with SD error bars
                ax.errorbar(valid_x, means, yerr=sds,
                           marker='o', markersize=8, capsize=5, lw=2,
                           color=SITE_COLORS[sname], label=sname,
                           markeredgecolor='black', markeredgewidth=0.5,
                           capthick=1.5)

                # SD shaded band
                ax.fill_between(valid_x, means - sds, means + sds,
                               color=SITE_COLORS[sname], alpha=0.1)

                # Individual replicate points (jittered)
                rng = np.random.default_rng(hash(sname) % 2**32)
                for xi, sk in zip(valid_x, [construct_order[ix] for ix in valid_x]):
                    v = get_vals(sk, sname, metric)
                    jitter = rng.normal(0, 0.06, len(v))
                    ax.scatter(xi + jitter, v, c=SITE_COLORS[sname], s=10,
                              alpha=0.3, edgecolors='none', zorder=2)

        ax.set_xticks(range(len(construct_order)))
        ax.set_xticklabels(x_labels, fontsize=9)
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=9)

        if metric == 'saa':
            ax.axhline(0.5, color='gray', ls='--', lw=0.8, alpha=0.4)
            ax.set_ylim(0.1, 0.8)
        elif metric == 'cone':
            ax.axhline(90, color='gray', ls='--', lw=0.8, alpha=0.4)
        if metric == 'shell':
            ax.annotate('more\ncrowded', xy=(0.95, 0.95), xycoords='axes fraction',
                       fontsize=8, ha='right', va='top', color='red', alpha=0.6)
            ax.annotate('more\nopen', xy=(0.95, 0.05), xycoords='axes fraction',
                       fontsize=8, ha='right', va='bottom', color='green', alpha=0.6)

        add_panel_label(ax, chr(65 + ax_idx))

    fig.suptitle('Domain Context Effect on Active Site Accessibility\n'
                 '(constructs ordered by increasing chain length)',
                 fontsize=13, fontweight='bold', y=1.08)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'fig4_domain_context.png'))
    fig.savefig(os.path.join(FIG_PATH, 'fig4_domain_context.svg'))
    plt.close()
    print("  Fig 4: fig4_domain_context")


# ============================================================
# Figure 5: Cone angle vs SAA scatter — docking feasibility
# ============================================================

def make_fig5():
    fig, ax = plt.subplots(figsize=(8, 6))

    # Quadrant annotations
    ax.axhline(90, color='gray', ls='--', lw=0.8, alpha=0.4)
    ax.axvline(0.5, color='gray', ls='--', lw=0.8, alpha=0.4)

    ax.fill_between([0.5, 1.0], 90, 180, alpha=0.06, color='green')
    ax.fill_between([0.0, 0.3], 0, 90, alpha=0.06, color='red')

    ax.text(0.75, 170, 'Readily\naccessible', ha='center', va='top',
            fontsize=10, color='green', fontweight='bold', alpha=0.6)
    ax.text(0.15, 20, 'Sterically\noccluded', ha='center', va='bottom',
            fontsize=10, color='red', fontweight='bold', alpha=0.6)

    markers = {'fl': 'o', 'norrm': 'p', 'noart': 'h', 'core': 'D', 'mka': '^',
               'md': 's', 'md3art': 'v', 'fl_optimized': '*'}
    rng = np.random.default_rng(42)

    for sk in SETS:
        for sname in CONSTRUCT_SITES[sk]:
            saa = get_vals(sk, sname, 'saa')
            cone = get_vals(sk, sname, 'cone')
            if len(saa) == 0:
                continue

            # Individual replicate points (faded)
            jx = rng.normal(0, 0.005, len(saa))
            jy = rng.normal(0, 0.5, len(cone))
            ax.scatter(saa + jx, cone + jy, marker=markers[sk], s=20,
                      color=SITE_COLORS[sname], alpha=0.25,
                      edgecolors='none', zorder=2)

            # Mean + SD crosshairs
            m_saa = np.mean(saa)
            m_cone = np.mean(cone)
            s_saa = np.std(saa)
            s_cone = np.std(cone)

            ax.errorbar(m_saa, m_cone, xerr=s_saa, yerr=s_cone,
                       marker=markers[sk], markersize=10, capsize=4,
                       color=SITE_COLORS[sname], markeredgecolor='black',
                       markeredgewidth=0.8, lw=1.5, capthick=1.2, zorder=5)

            # Label
            ax.annotate(f'{SET_LABELS_SHORT[sk]}',
                       (m_saa + 0.012, m_cone + 2.5),
                       fontsize=7, color=SITE_COLORS[sname], alpha=0.8)

    # Legend: sites by color
    site_handles = [Line2D([0], [0], marker='o', color=SITE_COLORS[s], lw=0,
                           markersize=8, markeredgecolor='black', label=s)
                    for s in SITES]
    set_handles = [Line2D([0], [0], marker=markers[s], color='gray', lw=0,
                          markersize=8, markeredgecolor='black',
                          label=SET_LABELS_SHORT[s])
                   for s in SETS]

    leg1 = ax.legend(handles=site_handles, title='Active Site',
                     loc='upper left', fontsize=9, title_fontsize=10)
    ax.add_artist(leg1)
    ax.legend(handles=set_handles, title='Construct',
              loc='lower right', fontsize=9, title_fontsize=10)

    ax.set_xlabel('Solid-Angle Accessibility (fraction of open directions)', fontsize=12)
    ax.set_ylabel('Max Approach Cone Half-Angle (\u00b0)', fontsize=12)
    ax.set_xlim(0.1, 0.8)
    ax.set_ylim(30, 115)
    ax.set_title('Docking Feasibility: Cone Angle vs Accessibility',
                 fontsize=14, fontweight='bold')

    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'fig5_docking_feasibility.png'))
    fig.savefig(os.path.join(FIG_PATH, 'fig5_docking_feasibility.svg'))
    plt.close()
    print("  Fig 5: fig5_docking_feasibility")


# ============================================================
# Figure 6: Schematic — what the metrics measure
# ============================================================

def make_fig6():
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for ax in axes:
        ax.set_xlim(-6, 6)
        ax.set_ylim(-6, 6)
        ax.set_aspect('equal')
        ax.axis('off')

    # --- Panel A: Solid-Angle Accessibility ---
    ax = axes[0]
    ax.set_title('Solid-Angle Accessibility (SAA)', fontsize=11, fontweight='bold', pad=12)

    blob_angles = np.linspace(0, 2*np.pi, 100)
    r_blob = 3.5 + 0.8*np.sin(3*blob_angles) + 0.5*np.cos(5*blob_angles)
    ax.fill(r_blob*np.cos(blob_angles), r_blob*np.sin(blob_angles),
            color='#ddd', edgecolor='#888', lw=1.5, zorder=1)
    ax.text(0, -1.5, 'protein\nchain', ha='center', va='center',
            fontsize=8, color='#666', style='italic')

    ax.plot(0, 0, 'o', color='red', markersize=12, zorder=5, markeredgecolor='black')
    ax.text(0, 0.6, 'pocket', ha='center', fontsize=8, fontweight='bold', color='red')

    n_show = 16
    for i in range(n_show):
        angle = 2 * np.pi * i / n_show
        dx, dy = np.cos(angle), np.sin(angle)
        r_at_angle = 3.5 + 0.8*np.sin(3*angle) + 0.5*np.cos(5*angle)

        if r_at_angle < 3.8:
            ax.annotate('', xy=(dx*r_at_angle*0.8, dy*r_at_angle*0.8), xytext=(0, 0),
                       arrowprops=dict(arrowstyle='->', color='red', lw=1, alpha=0.5))
        else:
            ax.annotate('', xy=(dx*5, dy*5), xytext=(0, 0),
                       arrowprops=dict(arrowstyle='->', color='green', lw=1.2, alpha=0.7))

    ax.text(0, -5.5, 'SAA = open rays / total rays', ha='center', fontsize=9,
            fontweight='bold', color='#333',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='gray'))
    add_panel_label(ax, 'A', x=-0.05, y=1.05)

    # --- Panel B: Max Cone Angle ---
    ax = axes[1]
    ax.set_title('Max Approach Cone', fontsize=11, fontweight='bold', pad=12)

    r_blob2 = 3.0 + 1.0*np.sin(2*blob_angles) + 0.3*np.cos(4*blob_angles)
    bx = r_blob2*np.cos(blob_angles) - 0.5
    by = r_blob2*np.sin(blob_angles) - 1.0
    ax.fill(bx, by, color='#ddd', edgecolor='#888', lw=1.5, zorder=1)

    site_x, site_y = 1.5, 1.5
    ax.plot(site_x, site_y, 'o', color='red', markersize=12, zorder=5, markeredgecolor='black')

    cone_dir = 50
    cone_half = 35
    wedge = Wedge((site_x, site_y), 5, cone_dir - cone_half, cone_dir + cone_half,
                  alpha=0.2, color='green', zorder=2)
    ax.add_patch(wedge)

    for sign in [-1, 1]:
        a = np.radians(cone_dir + sign * cone_half)
        ax.plot([site_x, site_x + 5*np.cos(a)],
                [site_y, site_y + 5*np.sin(a)],
                'g--', lw=1.5, alpha=0.7)

    arc_r = 2.0
    arc_angles = np.linspace(np.radians(cone_dir - cone_half),
                             np.radians(cone_dir + cone_half), 30)
    ax.plot(site_x + arc_r*np.cos(arc_angles),
            site_y + arc_r*np.sin(arc_angles), 'g-', lw=2)

    ax.annotate(f'{2*cone_half}\u00b0', xy=(site_x + 1.8, site_y + 2.2),
               fontsize=11, fontweight='bold', color='green')

    app_x = site_x + 4.2*np.cos(np.radians(cone_dir))
    app_y = site_y + 4.2*np.sin(np.radians(cone_dir))
    circle = Circle((app_x, app_y), 0.8, color='#f58231', alpha=0.6, zorder=3)
    ax.add_patch(circle)
    ax.text(app_x, app_y, 'ART', ha='center', va='center', fontsize=8, fontweight='bold')

    ax.text(0, -5.5, 'Widest unobstructed\napproach direction', ha='center', fontsize=9,
            fontweight='bold', color='#333',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='gray'))
    add_panel_label(ax, 'B', x=-0.05, y=1.05)

    # --- Panel C: Shell Density ---
    ax = axes[2]
    ax.set_title('Shell Density (2\u20135 nm)', fontsize=11, fontweight='bold', pad=12)

    ax.plot(0, 0, 'o', color='red', markersize=12, zorder=5, markeredgecolor='black')
    ax.text(0, 0.6, 'pocket', ha='center', fontsize=8, fontweight='bold', color='red')

    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(2*np.cos(theta), 2*np.sin(theta), 'b--', lw=1, alpha=0.5)
    ax.text(2.2, 0, '2 nm', fontsize=7, color='blue', alpha=0.6)
    ax.plot(5*np.cos(theta), 5*np.sin(theta), 'b--', lw=1, alpha=0.5)
    ax.text(5.2, 0, '5 nm', fontsize=7, color='blue', alpha=0.6)

    ax.fill_between(5*np.cos(theta), 5*np.sin(theta),
                    2*np.cos(theta), alpha=0.08, color='blue')

    rng = np.random.default_rng(123)
    n_beads = 40
    bead_r = rng.uniform(2.2, 4.8, n_beads)
    bead_a = rng.uniform(0, 2*np.pi, n_beads)
    ax.scatter(bead_r * np.cos(bead_a), bead_r * np.sin(bead_a),
              c='#888', s=30, alpha=0.6, edgecolors='#555', linewidths=0.5, zorder=3)

    ax.text(0, -5.5, 'Volume fraction occupied\nby own-chain beads', ha='center', fontsize=9,
            fontweight='bold', color='#333',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='gray'))
    add_panel_label(ax, 'C', x=-0.05, y=1.05)

    fig.suptitle('What Each Metric Measures', fontsize=14, fontweight='bold', y=1.04)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'fig6_schematic.png'))
    fig.savefig(os.path.join(FIG_PATH, 'fig6_schematic.svg'))
    plt.close()
    print("  Fig 6: fig6_schematic")


# ============================================================
# Figure 7: Site ranking summary
# ============================================================

def make_fig7():
    fig, ax = plt.subplots(figsize=(10, 5))

    construct_order = ['fl', 'norrm', 'noart', 'core', 'mka', 'md', 'md3art',
                       'fl_optimized']
    y_pos = 0
    y_ticks = []
    y_labels = []

    for sk in construct_order:
        site_means = {}
        for sname in CONSTRUCT_SITES[sk]:
            v = get_vals(sk, sname, 'saa')
            if len(v) > 0:
                site_means[sname] = (np.mean(v), np.std(v))

        sorted_sites = sorted(site_means.items(), key=lambda x: -x[1][0])

        y_ticks.append(y_pos)
        y_labels.append(SET_LABELS[sk])

        for rank, (sname, (mean, std)) in enumerate(sorted_sites):
            ax.barh(y_pos, mean, height=0.6, left=0,
                    color=SITE_COLORS[sname], alpha=0.7,
                    edgecolor='black', linewidth=0.5)
            ax.errorbar(mean, y_pos, xerr=std, capsize=3,
                       color='black', lw=1.2)
            label_x = mean - 0.02 if mean > 0.15 else mean + 0.02
            ha = 'right' if mean > 0.15 else 'left'
            ax.text(label_x, y_pos, f'{sname} ({mean:.2f}\u00b1{std:.2f})',
                   ha=ha, va='center', fontsize=9, fontweight='bold',
                   color='white' if mean > 0.15 else 'black')
            y_pos -= 0.8

        y_pos -= 0.5

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=11, fontweight='bold')
    ax.set_xlabel('Solid-Angle Accessibility', fontsize=12)
    ax.axvline(0.5, color='gray', ls='--', lw=0.8, alpha=0.5, label='50% open')
    ax.axvline(0.3, color='red', ls=':', lw=0.8, alpha=0.5, label='30% open')
    ax.legend(loc='lower right', fontsize=9)
    ax.set_xlim(0, 0.85)
    ax.set_title('Active Site Accessibility Ranking by Construct',
                 fontsize=14, fontweight='bold')
    ax.invert_yaxis()

    fig.tight_layout()
    fig.savefig(os.path.join(FIG_PATH, 'fig7_site_ranking.png'))
    fig.savefig(os.path.join(FIG_PATH, 'fig7_site_ranking.svg'))
    plt.close()
    print("  Fig 7: fig7_site_ranking")


# ============================================================
# Main
# ============================================================

print("=" * 60)
print("Generating accessibility figures")
print("=" * 60)

make_fig1()
make_fig2()
make_fig3()
make_fig4()
make_fig5()
make_fig6()
make_fig7()

print(f"\nAll figures saved to: {FIG_PATH}")
print("=" * 60)
