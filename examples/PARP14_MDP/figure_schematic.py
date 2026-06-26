#!/usr/bin/env python3
"""
Schematic figure: PARP14 domain architecture → CALVADOS CG-MD pipeline.
Simple, minimal text.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
from matplotlib.lines import Line2D
import os

CWD = os.path.dirname(os.path.abspath(__file__))
# Dated category subdir: figures/99_misc/<YYYY-MM-DD>/
import sys as _sys
_sys.path.insert(0, CWD)
from _fig_layout import get_fig_dir as _get_fig_dir
FIG_PATH = str(_get_fig_dir('99_misc'))

# ============================================================
# Domain architecture
# ============================================================

DOMAINS = [
    ('RRM1',    1,  145, '#A0A0A0'),
    ('RRM2',  146,  224, '#B0B0B0'),
    ('RRM3',  225,  314, '#C0C0C0'),
    ('KH1-6', 315,  737, '#4ECDC4'),
    ('KH7a',  738,  789, '#45B7AA'),
    ('MD1',   790, 1004, '#6B1F7A'),   # dark purple
    ('MD2',  1004, 1193, '#9450A8'),   # medium purple
    ('MD3',  1207, 1388, '#B87FCC'),   # light purple
    ('KHb-8',1389, 1533, '#45B7AA'),
    ('WWE',  1534, 1602, '#F7DC6F'),
    ('ART',  1603, 1801, '#E74C3C'),   # red
]

TOTAL_RES = 1801


def draw_domain_bar(ax, y, x0, x1, height, domains, total, label_size=7):
    """Draw a horizontal domain architecture bar."""
    bar_width = x1 - x0
    for name, start, end, color in domains:
        dx = (start - 1) / total * bar_width
        dw = (end - start + 1) / total * bar_width
        rect = FancyBboxPatch((x0 + dx, y - height/2), dw, height,
                               boxstyle='round,pad=0.003',
                               facecolor=color, edgecolor='black',
                               linewidth=0.5, zorder=3)
        ax.add_patch(rect)
        # Label
        cx = x0 + dx + dw/2
        if dw > 0.04:
            ax.text(cx, y, name, ha='center', va='center',
                    fontsize=label_size, fontweight='bold', color='white',
                    zorder=4)


def draw_bead_chain(ax, cx, cy, n_beads, radius, spread, colors, seed=42):
    """Draw a coarse-grained bead chain as connected spheres."""
    rng = np.random.default_rng(seed)
    # Random walk for chain path
    angles = np.cumsum(rng.normal(0, 0.4, n_beads))
    xs = np.cumsum(np.cos(angles) * spread) + cx
    ys = np.cumsum(np.sin(angles) * spread) + cy
    # Center
    xs -= (xs.mean() - cx)
    ys -= (ys.mean() - cy)

    # Bonds
    for i in range(n_beads - 1):
        ax.plot([xs[i], xs[i+1]], [ys[i], ys[i+1]],
                color='#555', lw=0.8, zorder=1)
    # Beads
    for i in range(n_beads):
        circle = Circle((xs[i], ys[i]), radius, facecolor=colors[i],
                        edgecolor='#333', linewidth=0.3, zorder=2)
        ax.add_patch(circle)
    return xs, ys


def draw_arrow(ax, x0, y0, x1, y1, color='black', lw=2):
    arrow = FancyArrowPatch((x0, y0), (x1, y1),
                            arrowstyle='->', mutation_scale=20,
                            color=color, lw=lw, zorder=5)
    ax.add_patch(arrow)


# ============================================================
# Figure
# ============================================================

fig, ax = plt.subplots(figsize=(14, 8))
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.set_aspect('equal')
ax.axis('off')

# --- Row 1: Full-length domain architecture ---
y_row1 = 0.92
ax.text(0.0, y_row1 + 0.05, 'PARP14 (1801 aa)', fontsize=13, fontweight='bold',
        va='bottom')
draw_domain_bar(ax, y_row1, 0.0, 1.0, 0.05, DOMAINS, TOTAL_RES, label_size=7)

# N/C terminus labels
ax.text(-0.02, y_row1, 'N', fontsize=9, ha='right', va='center', fontweight='bold')
ax.text(1.02, y_row1, 'C', fontsize=9, ha='left', va='center', fontweight='bold')

# --- Arrow down to constructs ---
draw_arrow(ax, 0.5, 0.87, 0.5, 0.82, color='#333', lw=2.5)
ax.text(0.54, 0.845, 'domain deletion\nconstructs', fontsize=8, va='center',
        ha='left', style='italic', color='#555')

# --- Row 2: Construct bars (6 constructs) ---
constructs = [
    ('FL (1801)',  list(range(11))),
    ('No-RRM (1474)', [3,4,5,6,7,8,9,10]),
    ('No-ART (1275)', [3,4,5,6,7,8,9]),
    ('Core (1051)',    [4,5,6,7,8,9,10]),
    ('MKA (999)',      [5,6,7,8,9,10]),
    ('MD only (586)',  [5,6,7]),
]

y_start = 0.78
dy = 0.045
for ci, (clabel, didx) in enumerate(constructs):
    y = y_start - ci * dy
    subset = [DOMAINS[i] for i in didx]
    # Compute width proportional to residue count
    n_res = sum(e - s + 1 for _, s, e, _ in subset)
    bar_w = n_res / TOTAL_RES * 1.0
    x0 = 0.15
    ax.text(x0 - 0.01, y, clabel, fontsize=7, ha='right', va='center',
            fontweight='bold')
    # Draw bar with internal scaling
    for name, start, end, color in subset:
        # Map to construct-relative position
        dx_total = 0
        for _, s2, e2, _ in subset:
            if s2 < start:
                dx_total += (e2 - s2 + 1)
        dw = (end - start + 1) / n_res * bar_w
        dx = dx_total / n_res * bar_w
        rect = FancyBboxPatch((x0 + dx, y - 0.015), dw, 0.03,
                               boxstyle='round,pad=0.002',
                               facecolor=color, edgecolor='black',
                               linewidth=0.4, zorder=3)
        ax.add_patch(rect)

# --- Arrow down to AF3 ---
y_af3 = 0.47
draw_arrow(ax, 0.5, y_start - len(constructs) * dy + 0.01, 0.5, y_af3 + 0.04,
           color='#333', lw=2.5)
ax.text(0.54, y_af3 + 0.065, 'AlphaFold3\n(25 models each)', fontsize=8,
        va='center', ha='left', style='italic', color='#555')

# --- Row 3: AF3 structure → CG beads ---
# Left: cartoon of folded structure
struct_cx, struct_cy = 0.22, y_af3 - 0.05
ax.text(struct_cx, struct_cy + 0.08, 'All-atom\nstructure', fontsize=8,
        ha='center', va='bottom', color='#555')
# Draw a simple folded protein cartoon (blobs)
rng = np.random.default_rng(77)
for i in range(25):
    angle = rng.uniform(0, 2*np.pi)
    r = rng.uniform(0.01, 0.06)
    bx = struct_cx + r * np.cos(angle)
    by = struct_cy + r * np.sin(angle) * 0.8
    c = Circle((bx, by), 0.008, facecolor='#888', edgecolor='#555',
              linewidth=0.3, alpha=0.6, zorder=2)
    ax.add_patch(c)

# Arrow to CG
draw_arrow(ax, struct_cx + 0.09, struct_cy, struct_cx + 0.18, struct_cy,
           color='#333', lw=2)
ax.text(struct_cx + 0.135, struct_cy + 0.02, 'CG', fontsize=9,
        ha='center', va='bottom', fontweight='bold', color='#333')

# Middle: CG bead model
cg_cx, cg_cy = 0.55, y_af3 - 0.05
ax.text(cg_cx, cg_cy + 0.08, '1 bead / residue', fontsize=8,
        ha='center', va='bottom', color='#555')

# Build color array for beads
n_show = 40
bead_colors = []
# Simplified: white with domain-colored regions
for i in range(n_show):
    frac = i / n_show
    if 0.4 < frac < 0.55:
        bead_colors.append('#6B1F7A')  # MD1
    elif 0.55 < frac < 0.7:
        bead_colors.append('#9450A8')  # MD2
    elif 0.7 < frac < 0.85:
        bead_colors.append('#B87FCC')  # MD3
    elif frac > 0.9:
        bead_colors.append('#E74C3C')  # ART
    else:
        bead_colors.append('#DDDDDD')  # white/linker

draw_bead_chain(ax, cg_cx, cg_cy, n_show, 0.008, 0.016, bead_colors, seed=42)

# Arrow to simulation
draw_arrow(ax, cg_cx + 0.14, cg_cy, cg_cx + 0.22, cg_cy,
           color='#333', lw=2)

# Right: simulation box
sim_cx, sim_cy = 0.88, y_af3 - 0.05
ax.text(sim_cx, sim_cy + 0.08, 'CALVADOS\nMD (20 ns)', fontsize=8,
        ha='center', va='bottom', color='#555')

# Draw box
box_w, box_h = 0.12, 0.12
rect = mpatches.Rectangle((sim_cx - box_w/2, sim_cy - box_h/2), box_w, box_h,
                            facecolor='#E8F4FF', edgecolor='#2980B9',
                            linewidth=1.5, linestyle='--', zorder=1)
ax.add_patch(rect)
# Chain inside box
box_beads = 20
box_colors = ['#DDDDDD'] * 5 + ['#6B1F7A'] * 4 + ['#9450A8'] * 4 + \
             ['#B87FCC'] * 3 + ['#E74C3C'] * 4
draw_bead_chain(ax, sim_cx, sim_cy, box_beads, 0.006, 0.012,
                box_colors[:box_beads], seed=99)
# Wavy lines for solvent
for _ in range(8):
    sx = sim_cx + rng.uniform(-0.04, 0.04)
    sy = sim_cy + rng.uniform(-0.04, 0.04)
    ax.text(sx, sy, '~', fontsize=5, color='#AED6F1', alpha=0.5,
            ha='center', va='center', zorder=0)

# --- Arrow down to analysis ---
y_analysis = 0.15
draw_arrow(ax, 0.5, y_af3 - 0.14, 0.5, y_analysis + 0.08, color='#333', lw=2.5)

# --- Row 4: Analysis outputs ---
ax.text(0.5, y_analysis + 0.06, 'Analysis', fontsize=10, fontweight='bold',
        ha='center', va='bottom')

# Analysis boxes
analysis_items = [
    (0.08, 'Rg / Ree'),
    (0.24, 'Energy\ndecomposition'),
    (0.40, 'Contact\nmaps'),
    (0.56, 'Active site\naccessibility'),
    (0.72, 'WCN\n(burial)'),
    (0.88, 'Domain\ndistances'),
]

for x, label in analysis_items:
    rect = FancyBboxPatch((x - 0.06, y_analysis - 0.05), 0.12, 0.06,
                           boxstyle='round,pad=0.01',
                           facecolor='#F0F0F0', edgecolor='#666',
                           linewidth=0.8, zorder=3)
    ax.add_patch(rect)
    ax.text(x, y_analysis - 0.02, label, fontsize=6.5, ha='center',
            va='center', fontweight='bold', color='#333', zorder=4)

# --- Domain legend ---
legend_y = 0.03
legend_items = [
    ('RRM1-3', '#A0A0A0'), ('KH', '#4ECDC4'), ('MD1', '#6B1F7A'),
    ('MD2', '#9450A8'), ('MD3', '#B87FCC'), ('WWE', '#F7DC6F'), ('ART', '#E74C3C'),
]
total_w = 0.7
item_w = total_w / len(legend_items)
x_start = 0.15
for i, (name, color) in enumerate(legend_items):
    x = x_start + i * item_w
    rect = mpatches.Rectangle((x, legend_y - 0.008), 0.025, 0.016,
                                facecolor=color, edgecolor='black',
                                linewidth=0.5, zorder=3)
    ax.add_patch(rect)
    ax.text(x + 0.03, legend_y, name, fontsize=7, va='center', fontweight='bold')

# --- Title ---
# (none — the figure speaks for itself)

fig.tight_layout(pad=0.5)
for ext in ['png', 'svg']:
    fig.savefig(os.path.join(FIG_PATH, f'schematic_pipeline.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: figures/schematic_pipeline.{{png,svg}}")
