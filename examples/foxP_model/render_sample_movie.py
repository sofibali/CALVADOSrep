#!/usr/bin/env python
"""
Render sample movie frames (first, middle, last) from all-atom trajectories.
Uses matplotlib 3D rendering - no PyMOL required.

Creates:
  - Per-trajectory animated GIF (3 frames: first, middle, last)
  - Per-trajectory panel PNG (3 frames side by side)

Usage:
    conda run -n CALVADOS python render_sample_movie.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

# Compatibility shim
import mdtraj
if not hasattr(mdtraj.Trajectory, 'bfactors'):
    mdtraj.Trajectory.bfactors = property(
        lambda self: getattr(self, '_bfactors', np.zeros((self.n_frames, self.n_atoms))),
        lambda self, val: setattr(self, '_bfactors', val),
    )
import mdtraj as md

BASE = os.path.dirname(os.path.abspath(__file__))
INDIVIDUAL = os.path.join(BASE, 'individual_chains')
OUTPUT = os.path.join(BASE, 'rendered_movies')
os.makedirs(OUTPUT, exist_ok=True)

# Domain definitions (1-indexed residue numbers within each chain)
DOMAINS = {
    'd1': (117, 198),   # 0-indexed: 118-199 in 1-indexed
    'd2': (298, 374),   # 299-375
    'd3': (450, 550),   # 451-551 (DBD)
}

# Colors
COLORS = {
    'FOXP4_idr': '#90EE90',   # light green
    'FOXP4_d1':  '#7CCD7C',   # medium green
    'FOXP4_d2':  '#66AA66',   # split pea
    'FOXP4_d3':  '#228B22',   # forest green
    'FOXP1_idr': '#FFDAB9',   # peach/light orange
    'FOXP1_d1':  '#FFA500',   # orange
    'FOXP1_d2':  '#FF8C00',   # dark orange
    'FOXP1_d3':  '#D2691E',   # chocolate
}

# Line widths
LW = {
    'idr': 0.8,
    'd1':  2.0,
    'd2':  2.0,
    'd3':  2.5,
}


def kabsch_align(mobile, target):
    """Align mobile to target using Kabsch algorithm. Returns aligned mobile."""
    mob_center = mobile.mean(axis=0)
    tgt_center = target.mean(axis=0)
    mob_c = mobile - mob_center
    tgt_c = target - tgt_center
    H = mob_c.T @ tgt_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ sign_matrix @ U.T
    return (mob_c @ R.T) + tgt_center


def get_ca_indices(top):
    """Get indices of CA atoms."""
    return top.select('name CA')


def classify_residue(res_idx, n_foxp4, is_two_chain):
    """Classify a residue index (0-based, within CA atoms) into chain+domain."""
    if is_two_chain:
        # 2-chain: res_idx 0..676 = FOXP4, 677..1356 = FOXP1
        if res_idx < n_foxp4:
            chain = 'FOXP4'
            local_idx = res_idx
        else:
            chain = 'FOXP1'
            local_idx = res_idx - n_foxp4
    else:
        # 1-chain merged: same logic
        if res_idx < n_foxp4:
            chain = 'FOXP4'
            local_idx = res_idx
        else:
            chain = 'FOXP1'
            local_idx = res_idx - n_foxp4

    for dname, (start, end) in DOMAINS.items():
        if start <= local_idx <= end:
            return f'{chain}_{dname}'
    return f'{chain}_idr'


def build_segments(n_residues, n_foxp4, is_two_chain):
    """Build colored line segments for the backbone trace."""
    segments = []  # list of (start_idx, end_idx, color, linewidth)
    current_class = classify_residue(0, n_foxp4, is_two_chain)
    seg_start = 0

    for i in range(1, n_residues):
        cls = classify_residue(i, n_foxp4, is_two_chain)
        if cls != current_class:
            # End current segment
            color = COLORS[current_class]
            domain_type = current_class.split('_')[1]
            lw = LW[domain_type]
            segments.append((seg_start, i, color, lw))
            seg_start = i
            current_class = cls

    # Final segment
    color = COLORS[current_class]
    domain_type = current_class.split('_')[1]
    lw = LW[domain_type]
    segments.append((seg_start, n_residues - 1, color, lw))

    return segments


def render_frame(ax, ca_coords, segments, title='', elev=20, azim=45):
    """Render one frame onto a 3D axes."""
    ax.clear()

    for seg_start, seg_end, color, lw in segments:
        idx = slice(seg_start, seg_end + 1)
        xs = ca_coords[idx, 0]
        ys = ca_coords[idx, 1]
        zs = ca_coords[idx, 2]
        ax.plot(xs, ys, zs, color=color, linewidth=lw, alpha=0.85)

    # Domain3 spheres (larger markers for DBD)
    n_foxp4 = 677
    for chain, offset in [('FOXP4', 0), ('FOXP1', n_foxp4)]:
        d3_start = DOMAINS['d3'][0] + offset
        d3_end = DOMAINS['d3'][1] + offset
        if d3_end < len(ca_coords):
            d3_coords = ca_coords[d3_start:d3_end+1]
            color = COLORS[f'{chain}_d3']
            # Plot sparse markers for domain3 to show volume
            ax.scatter(d3_coords[::5, 0], d3_coords[::5, 1], d3_coords[::5, 2],
                      c=color, s=15, alpha=0.3, edgecolors='none')

    ax.set_title(title, fontsize=11, fontweight='bold', pad=2)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_zlabel('')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # Remove tick marks
    ax.xaxis.set_tick_params(size=0)
    ax.yaxis.set_tick_params(size=0)
    ax.zaxis.set_tick_params(size=0)

    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect([1, 1, 1])


def set_consistent_limits(ax, all_coords_list, padding=0.5):
    """Set axis limits consistently across frames."""
    all_pts = np.concatenate(all_coords_list, axis=0)
    center = all_pts.mean(axis=0)
    max_range = (all_pts.max(axis=0) - all_pts.min(axis=0)).max() / 2 + padding
    ax.set_xlim(center[0] - max_range, center[0] + max_range)
    ax.set_ylim(center[1] - max_range, center[1] + max_range)
    ax.set_zlim(center[2] - max_range, center[2] + max_range)


def process_simulation(sim_name, sysname):
    """Process one simulation: render panels + GIF."""
    movie_dir = os.path.join(INDIVIDUAL, sim_name, 'viz', 'movie')
    top_file = os.path.join(movie_dir, 'allatom_top.pdb')
    dcd_file = os.path.join(movie_dir, 'allatom.dcd')

    if not os.path.exists(dcd_file):
        print(f"  SKIP: {dcd_file} not found")
        return

    print(f"\nProcessing: {sim_name}")
    traj = md.load(dcd_file, top=top_file)
    n_frames = traj.n_frames
    print(f"  Frames: {n_frames}")

    # Get CA atoms
    ca_idx = get_ca_indices(traj.top)
    ca_traj = traj.atom_slice(ca_idx)
    n_ca = ca_traj.n_atoms
    print(f"  CA atoms: {n_ca}")

    # Determine chain structure
    chains = list(ca_traj.top.chains)
    is_two_chain = len(chains) >= 2
    n_foxp4 = 677

    # Frame indices: first, middle, last
    frame_indices = [0, n_frames // 2, n_frames - 1]
    frame_labels = ['First frame', 'Middle frame', 'Last frame']

    # Extract CA coordinates (nm) for selected frames
    ca_coords_list = []
    for fi in frame_indices:
        coords = ca_traj.xyz[fi] * 10  # nm -> Angstrom
        ca_coords_list.append(coords)

    # Align all to first frame using FOXP4 d3
    d3_start = DOMAINS['d3'][0]
    d3_end = DOMAINS['d3'][1] + 1
    ref_d3 = ca_coords_list[0][d3_start:d3_end]
    for i in range(1, len(ca_coords_list)):
        mob_d3 = ca_coords_list[i][d3_start:d3_end]
        # Get rotation from d3, apply to all
        mob_center = mob_d3.mean(axis=0)
        ref_center = ref_d3.mean(axis=0)
        mob_c = mob_d3 - mob_center
        ref_c = ref_d3 - ref_center
        H = mob_c.T @ ref_c
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        sign_matrix = np.diag([1, 1, np.sign(d)])
        R = Vt.T @ sign_matrix @ U.T
        # Apply to full structure
        all_centered = ca_coords_list[i] - mob_center
        ca_coords_list[i] = (all_centered @ R.T) + ref_center

    # Build colored segments
    segments = build_segments(n_ca, n_foxp4, is_two_chain)

    # --- Panel figure (3 frames side by side) ---
    fig = plt.figure(figsize=(18, 6), dpi=150)
    fig.suptitle(sysname.replace('_', ' ').title(), fontsize=14, fontweight='bold', y=0.98)

    for i, (coords, label) in enumerate(zip(ca_coords_list, frame_labels)):
        ax = fig.add_subplot(1, 3, i + 1, projection='3d')
        render_frame(ax, coords, segments, title=label, elev=15, azim=45 + i * 20)
        set_consistent_limits(ax, ca_coords_list)
        ax.grid(False)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('lightgray')
        ax.yaxis.pane.set_edgecolor('lightgray')
        ax.zaxis.pane.set_edgecolor('lightgray')

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=COLORS['FOXP4_d3'], linewidth=3, label='FOXP4 DBD (d3)'),
        Line2D([0], [0], color=COLORS['FOXP4_d2'], linewidth=2, label='FOXP4 d2'),
        Line2D([0], [0], color=COLORS['FOXP4_d1'], linewidth=2, label='FOXP4 d1'),
        Line2D([0], [0], color=COLORS['FOXP4_idr'], linewidth=1, label='FOXP4 IDR'),
        Line2D([0], [0], color=COLORS['FOXP1_d3'], linewidth=3, label='FOXP1 DBD (d3)'),
        Line2D([0], [0], color=COLORS['FOXP1_d2'], linewidth=2, label='FOXP1 d2'),
        Line2D([0], [0], color=COLORS['FOXP1_d1'], linewidth=2, label='FOXP1 d1'),
        Line2D([0], [0], color=COLORS['FOXP1_idr'], linewidth=1, label='FOXP1 IDR'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=9,
              frameon=True, fancybox=True, shadow=False, borderpad=0.5)

    plt.tight_layout(rect=[0, 0.06, 1, 0.96])
    panel_path = os.path.join(OUTPUT, f'{sysname}_panels.png')
    fig.savefig(panel_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Panel PNG: {panel_path}")

    # --- Animated GIF (rotate + step through frames) ---
    fig_gif = plt.figure(figsize=(8, 8), dpi=120)
    ax_gif = fig_gif.add_subplot(111, projection='3d')
    ax_gif.grid(False)
    ax_gif.xaxis.pane.fill = False
    ax_gif.yaxis.pane.fill = False
    ax_gif.zaxis.pane.fill = False
    ax_gif.xaxis.pane.set_edgecolor('lightgray')
    ax_gif.yaxis.pane.set_edgecolor('lightgray')
    ax_gif.zaxis.pane.set_edgecolor('lightgray')

    gif_path = os.path.join(OUTPUT, f'{sysname}_movie.gif')
    writer = PillowWriter(fps=4)

    with writer.saving(fig_gif, gif_path, dpi=120):
        for fi, (coords, label) in enumerate(zip(ca_coords_list, frame_labels)):
            # Multiple rotation angles per frame for smooth rotation
            for azim in range(0, 360, 15):
                render_frame(ax_gif, coords, segments,
                           title=f'{sysname.replace("_", " ")}\n{label}',
                           elev=15, azim=azim)
                set_consistent_limits(ax_gif, ca_coords_list)
                ax_gif.grid(False)
                ax_gif.xaxis.pane.fill = False
                ax_gif.yaxis.pane.fill = False
                ax_gif.zaxis.pane.fill = False
                ax_gif.xaxis.pane.set_edgecolor('lightgray')
                ax_gif.yaxis.pane.set_edgecolor('lightgray')
                ax_gif.zaxis.pane.set_edgecolor('lightgray')
                writer.grab_frame()

    plt.close(fig_gif)
    gif_size = os.path.getsize(gif_path) / 1e6
    print(f"  Animated GIF: {gif_path} ({gif_size:.1f} MB)")


if __name__ == '__main__':
    sims = [
        ('option1_pull_interchain_d12', 'option1_pull_interchain_d12'),
        ('two_chains_coils_interchain', 'two_chains_coils_interchain'),
    ]

    for sim_name, sysname in sims:
        process_simulation(sim_name, sysname)

    print(f"\nAll outputs in: {OUTPUT}")
