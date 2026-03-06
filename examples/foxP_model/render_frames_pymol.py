#!/usr/bin/env python
"""
Render first/middle/last frames from FOXP trajectories using PyMOL headless.

Chain assignment: FOXP1 = chain A, FOXP4 = chain B (2-chain topology)
For single-chain topology: FOXP4 = resi 1-677, FOXP1 = resi 678-1357

After alignment, surface-showing domains (d2, d3) are extracted into
separate PyMOL objects to avoid surface abrasions at domain boundaries.

Usage:
    conda run -n pymol-render python render_frames_pymol.py
"""

import sys
import os
import numpy as np

# --- mdtraj compatibility shim ---
import mdtraj
if not hasattr(mdtraj.Trajectory, 'bfactors'):
    mdtraj.Trajectory.bfactors = property(
        lambda self: getattr(self, '_bfactors', np.zeros((self.n_frames, self.n_atoms))),
        lambda self, val: setattr(self, '_bfactors', val),
    )
import mdtraj as md

BASE = os.path.dirname(os.path.abspath(__file__))
INDIVIDUAL = os.path.join(BASE, 'individual_chains')
REF_PDB = os.path.join(INDIVIDUAL, '3g73_DBDA.pdb')
OUTPUT = os.path.join(BASE, 'rendered_movies')
os.makedirs(OUTPUT, exist_ok=True)

# Domain ranges (1-indexed residue numbers, within each chain)
D1 = (118, 199)
D2 = (299, 375)
D3 = (451, 551)

SIMS = {
    'option1_pull_interchain_d12': {
        'two_chain': True,
        'n_foxp4': 677,
    },
    'two_chains_coils_interchain': {
        'two_chain': False,
        'n_foxp4': 677,
    },
}


def extract_frames(sim_name):
    """Extract first, middle, last frames as individual PDB files via mdtraj."""
    movie_dir = os.path.join(INDIVIDUAL, sim_name, 'viz', 'movie')
    top_file = os.path.join(movie_dir, 'allatom_top.pdb')
    dcd_file = os.path.join(movie_dir, 'allatom.dcd')

    traj = md.load(dcd_file, top=top_file)
    n_frames = traj.n_frames
    indices = [0, n_frames // 2, n_frames - 1]
    labels = ['first', 'middle', 'last']

    tmp_dir = os.path.join(OUTPUT, f'{sim_name}_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    pdb_files = []
    for idx, label in zip(indices, labels):
        pdb_path = os.path.join(tmp_dir, f'frame_{label}.pdb')
        traj[idx].save_pdb(pdb_path)
        pdb_files.append((pdb_path, label))
        print(f"  Extracted frame {idx} -> {pdb_path}")

    return pdb_files


def setup_selections(sim_info):
    """Create chain/domain selections on the 'mol' object."""
    two_chain = sim_info['two_chain']
    n_foxp4 = sim_info['n_foxp4']
    OFF = n_foxp4

    if two_chain:
        # 2-chain: chain A = FOXP1, chain B = FOXP4
        cmd.select('FOXP1',     'mol and chain A')
        cmd.select('FOXP1_d1',  f'mol and chain A and resi {D1[0]}-{D1[1]}')
        cmd.select('FOXP1_d2',  f'mol and chain A and resi {D2[0]}-{D2[1]}')
        cmd.select('FOXP1_d3',  f'mol and chain A and resi {D3[0]}-{D3[1]}')
        cmd.select('FOXP1_idr', f'mol and chain A and not (resi {D1[0]}-{D1[1]} or resi {D2[0]}-{D2[1]} or resi {D3[0]}-{D3[1]})')

        cmd.select('FOXP4',     'mol and chain B')
        cmd.select('FOXP4_d1',  f'mol and chain B and resi {D1[0]}-{D1[1]}')
        cmd.select('FOXP4_d2',  f'mol and chain B and resi {D2[0]}-{D2[1]}')
        cmd.select('FOXP4_d3',  f'mol and chain B and resi {D3[0]}-{D3[1]}')
        cmd.select('FOXP4_idr', f'mol and chain B and not (resi {D1[0]}-{D1[1]} or resi {D2[0]}-{D2[1]} or resi {D3[0]}-{D3[1]})')
    else:
        # Single-chain: FOXP4 = resi 1-677, FOXP1 = resi 678-1357
        cmd.select('FOXP4',     f'mol and resi 1-{n_foxp4}')
        cmd.select('FOXP4_d1',  f'mol and resi {D1[0]}-{D1[1]}')
        cmd.select('FOXP4_d2',  f'mol and resi {D2[0]}-{D2[1]}')
        cmd.select('FOXP4_d3',  f'mol and resi {D3[0]}-{D3[1]}')
        cmd.select('FOXP4_idr', f'mol and resi 1-{n_foxp4} and not (resi {D1[0]}-{D1[1]} or resi {D2[0]}-{D2[1]} or resi {D3[0]}-{D3[1]})')

        cmd.select('FOXP1',     f'mol and resi {OFF+1}-{OFF+680}')
        cmd.select('FOXP1_d1',  f'mol and resi {D1[0]+OFF}-{D1[1]+OFF}')
        cmd.select('FOXP1_d2',  f'mol and resi {D2[0]+OFF}-{D2[1]+OFF}')
        cmd.select('FOXP1_d3',  f'mol and resi {D3[0]+OFF}-{D3[1]+OFF}')
        cmd.select('FOXP1_idr', f'mol and resi {OFF+1}-{OFF+680} and not (resi {D1[0]+OFF}-{D1[1]+OFF} or resi {D2[0]+OFF}-{D2[1]+OFF} or resi {D3[0]+OFF}-{D3[1]+OFF})')

    cmd.select('all_d1', 'FOXP4_d1 or FOXP1_d1')
    cmd.select('all_d2', 'FOXP4_d2 or FOXP1_d2')
    cmd.select('all_d3', 'FOXP4_d3 or FOXP1_d3')
    cmd.select('all_idr', 'FOXP4_idr or FOXP1_idr')
    cmd.select('ref_protein', 'ref and chain A and polymer.protein')
    cmd.select('ref_DNA', 'ref and (chain C or chain D) and polymer.nucleic')
    cmd.deselect()


def split_surface_domains():
    """Extract d2/d3 domains into separate objects for clean surface rendering.

    This avoids surface abrasions at domain boundaries. After extraction,
    'mol' retains IDR + d1 (cartoon only), while the extracted objects
    each get their own independent surface.
    """
    # Extract each surface domain into its own object
    # cmd.extract moves atoms OUT of 'mol' into the new object
    for protein, color in [('FOXP1', 'lightorange'), ('FOXP4', 'limegreen')]:
        for domain in ['d2', 'd3']:
            sel_name = f'{protein}_{domain}'
            obj_name = f'{protein}_{domain}_obj'
            cmd.extract(obj_name, sel_name)

            # Cartoon + surface on the extracted object
            cmd.show('cartoon', obj_name)
            cmd.show('surface', obj_name)
            cmd.color(color, obj_name)
            cmd.set('surface_color', color, obj_name)
            cmd.set('transparency', 0.65, obj_name)

    print("  Split d2/d3 into separate surface objects")


def setup_representations():
    """Set representations, coloring, and visual settings.

    Called AFTER split_surface_domains(), so 'mol' only contains IDR + d1.
    """
    cmd.hide('everything')

    # --- Reference DNA: gray cartoon ---
    cmd.show('cartoon', 'ref_DNA')
    cmd.color('gray60', 'ref_DNA')
    cmd.set('cartoon_ring_mode', 3)
    cmd.set('cartoon_ring_finder', 1)

    # Hide ref protein/ions (alignment only)
    cmd.hide('everything', 'ref and chain A')
    cmd.hide('everything', 'ref and chain B')

    # --- 'mol' now has IDR + d1 only (d2/d3 were extracted) ---
    cmd.show('cartoon', 'mol')
    # Color what remains in mol by chain
    cmd.color('lightorange', 'mol and FOXP1')
    cmd.color('limegreen', 'mol and FOXP4')
    # If the selections got invalidated by extract, recolor by residue range
    # (safety net for single-chain topology)
    if not cmd.count_atoms('mol and FOXP1'):
        # Selections on mol may be empty after extract; recolor mol broadly
        pass

    # d1: force coil secondary structure
    cmd.alter('mol', "ss='L'")  # mol only has IDR+d1 now, safe to set all to coil
    cmd.rebuild('mol')

    # IDR: thinner cartoon, slightly transparent
    cmd.set('cartoon_loop_radius', 0.10, 'mol')
    cmd.set('cartoon_transparency', 0.3, 'mol')

    # d1 is still in mol — give it normal cartoon thickness
    # (We can't easily distinguish d1 from IDR in mol after extract,
    #  but d1 is set to coil SS so it renders as loop anyway)

    # --- Show cartoon + surface on each extracted domain object ---
    for protein in ['FOXP1', 'FOXP4']:
        for domain in ['d2', 'd3']:
            obj = f'{protein}_{domain}_obj'
            cmd.show('cartoon', obj)
            cmd.show('surface', obj)

    # --- Global visual settings ---
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('antialias', 2)
    cmd.set('ray_trace_mode', 0)
    cmd.set('ray_shadows', 1)
    cmd.set('ray_shadow_decay_factor', 0.1)
    cmd.set('spec_reflect', 0.4)
    cmd.set('specular', 0.3)
    cmd.set('ambient', 0.35)
    cmd.set('direct', 0.6)
    cmd.set('cartoon_sampling', 7)
    cmd.set('cartoon_smooth_loops', 1)
    cmd.set('cartoon_oval_length', 1.0)
    cmd.set('surface_quality', 0)
    cmd.set('surface_color_smoothing', 1)


def compute_locked_view():
    """Compute a view with DNA at bottom, protein domains on top.

    Uses ONLY the reference structure (3g73) for orientation so the
    view is perfectly consistent across all frames and simulations.
    """
    import math

    cmd.center('ref_DNA')
    cmd.orient('ref_DNA')

    dna_center = cmd.centerofmass('ref_DNA')
    ref_prot_center = cmd.centerofmass('ref and chain A and polymer.protein')

    dy = ref_prot_center[1] - dna_center[1]
    dz = ref_prot_center[2] - dna_center[2]

    angle = -math.degrees(math.atan2(dz, dy))
    cmd.turn('x', angle)

    R = list(cmd.get_view())
    screen_y_dna  = R[1]*dna_center[0]      + R[4]*dna_center[1]      + R[7]*dna_center[2]
    screen_y_prot = R[1]*ref_prot_center[0]  + R[4]*ref_prot_center[1]  + R[7]*ref_prot_center[2]

    if screen_y_prot < screen_y_dna:
        cmd.turn('x', 180)

    cmd.move('y', 15)
    cmd.zoom('ref_DNA', buffer=60)

    return cmd.get_view()


def render_simulation(sim_name, sim_info):
    """Render first/middle/last frames for one simulation."""
    print(f"\n{'='*60}")
    print(f"Rendering: {sim_name}")
    print(f"{'='*60}")

    movie_dir = os.path.join(INDIVIDUAL, sim_name, 'viz', 'movie')
    if not os.path.exists(os.path.join(movie_dir, 'allatom.dcd')):
        print(f"  SKIP: no allatom.dcd")
        return

    pdb_files = extract_frames(sim_name)

    saved_view = None

    for pdb_path, label in pdb_files:
        # Fresh scene each frame
        cmd.delete('all')

        # Load reference
        cmd.load(REF_PDB, 'ref')

        # Load this frame
        cmd.load(pdb_path, 'mol')

        # Selections (on 'mol' before splitting)
        setup_selections(sim_info)

        # Align FOXP1 d3 → 3g73 chain A (reference FOXP1 DBD)
        cmd.align('FOXP1_d3 and name CA', 'ref_protein and name CA')

        # Split d2/d3 into separate objects (fixes surface abrasions)
        split_surface_domains()

        # Representations and coloring
        setup_representations()

        # Camera: compute once on first frame, reuse for others
        if saved_view is None:
            saved_view = compute_locked_view()
        else:
            cmd.set_view(saved_view)

        # Render
        fname = os.path.join(OUTPUT, f'{sim_name}_{label}.png')
        sys.stdout.write(f"  Ray-tracing {label}...")
        sys.stdout.flush()
        cmd.ray(1920, 1080)
        cmd.png(fname, dpi=300)
        sys.stdout.write(f" done -> {fname}\n")
        sys.stdout.flush()


if __name__ == '__main__':
    os.environ['PYMOL_FEEDBACK_LEVEL'] = '0'
    import pymol
    from pymol import cmd
    pymol.finish_launching(['pymol', '-cqQ'])
    cmd.feedback('disable', 'all', 'everything')

    for sim_name, sim_info in SIMS.items():
        render_simulation(sim_name, sim_info)

    cmd.quit()
    print(f"\nAll PNGs saved to: {OUTPUT}")
