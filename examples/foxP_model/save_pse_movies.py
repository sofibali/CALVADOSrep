#!/usr/bin/env python
"""
Save PyMOL session files (.pse) for FOXP trajectory movies.

Creates one PSE per simulation with:
  - All-atom trajectory loaded as multi-state object
  - FOXP1 d3 aligned to 3g73 reference across all states
  - d2/d3 split into separate objects (clean surfaces)
  - DNA at bottom, protein domains on top
  - View zoomed to fit all states (nothing goes off-screen)
  - Movie set up (mset) ready for playback/export

Usage:
    conda run -n pymol-render python save_pse_movies.py

To export movie locally in PyMOL (see bottom of script for full instructions).
"""

import sys
import os
import math

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
    'two_chains_coils_interchain': {
        'two_chain': False,
        'n_foxp4': 677,
    },
}


def setup_selections(sim_info):
    """Create chain/domain selections on the 'mol' object."""
    two_chain = sim_info['two_chain']
    n_foxp4 = sim_info['n_foxp4']
    OFF = n_foxp4

    if two_chain:
        # 2-chain topology: chain A = FOXP4 (677 res), chain B = FOXP1 (680 res)
        # (YAML system order: FOXP4 first, FOXP1 second)
        cmd.select('FOXP4',     'mol and chain A')
        cmd.select('FOXP4_d1',  f'mol and chain A and resi {D1[0]}-{D1[1]}')
        cmd.select('FOXP4_d2',  f'mol and chain A and resi {D2[0]}-{D2[1]}')
        cmd.select('FOXP4_d3',  f'mol and chain A and resi {D3[0]}-{D3[1]}')
        cmd.select('FOXP4_idr', f'mol and chain A and not (resi {D1[0]}-{D1[1]} or resi {D2[0]}-{D2[1]} or resi {D3[0]}-{D3[1]})')

        cmd.select('FOXP1',     'mol and chain B')
        cmd.select('FOXP1_d1',  f'mol and chain B and resi {D1[0]}-{D1[1]}')
        cmd.select('FOXP1_d2',  f'mol and chain B and resi {D2[0]}-{D2[1]}')
        cmd.select('FOXP1_d3',  f'mol and chain B and resi {D3[0]}-{D3[1]}')
        cmd.select('FOXP1_idr', f'mol and chain B and not (resi {D1[0]}-{D1[1]} or resi {D2[0]}-{D2[1]} or resi {D3[0]}-{D3[1]})')
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

    cmd.select('ref_protein', 'ref and chain A and polymer.protein')
    cmd.select('ref_DNA', 'ref and (chain C or chain D) and polymer.nucleic')
    cmd.deselect()


def split_surface_domains():
    """Extract d2/d3 domains into separate objects for clean surface rendering.

    cmd.extract moves atoms OUT of 'mol' into new objects (all states preserved).
    This avoids surface abrasions at domain boundaries.
    """
    for protein, color in [('FOXP1', 'lightorange'), ('FOXP4', 'limegreen')]:
        for domain in ['d2', 'd3']:
            sel_name = f'{protein}_{domain}'
            obj_name = f'{protein}_{domain}_obj'
            cmd.extract(obj_name, sel_name)
            cmd.show('cartoon', obj_name)
            cmd.show('surface', obj_name)
            cmd.color(color, obj_name)
            cmd.set('surface_color', color, obj_name)
            cmd.set('transparency', 0.65, obj_name)
    print("  Split d2/d3 into separate surface objects")


def setup_representations(sim_info):
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

    # --- mol: IDR + d1 only (d2/d3 were extracted) ---
    cmd.show('cartoon', 'mol')
    # Color by named selection first
    cmd.color('lightorange', 'mol and FOXP1')
    cmd.color('limegreen', 'mol and FOXP4')
    # Fallback: explicit residue-range coloring (robust after extract)
    if sim_info.get('two_chain'):
        cmd.color('limegreen', 'mol and chain A')   # FOXP4
        cmd.color('lightorange', 'mol and chain B')  # FOXP1
    else:
        n = sim_info.get('n_foxp4', 677)
        cmd.color('limegreen', f'mol and resi 1-{n}')        # FOXP4
        cmd.color('lightorange', f'mol and resi {n+1}-9999')  # FOXP1

    # d1: force coil secondary structure
    cmd.alter('mol', "ss='L'")
    cmd.rebuild('mol')

    # IDR: thinner cartoon, slightly transparent
    cmd.set('cartoon_loop_radius', 0.10, 'mol')
    cmd.set('cartoon_transparency', 0.3, 'mol')

    # --- Extracted domain objects: cartoon + surface ---
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


def compute_locked_view_all_states():
    """Compute a view with DNA at bottom, zoomed to fit ALL states.

    Uses the reference structure (3g73) for orientation so the view
    is consistent across simulations. Then zooms to fit all trajectory
    states so nothing goes off-screen during playback.
    """
    # Orient based on reference DNA
    cmd.center('ref_DNA')
    cmd.orient('ref_DNA')

    dna_center = cmd.centerofmass('ref_DNA')
    ref_prot_center = cmd.centerofmass('ref and chain A and polymer.protein')

    dy = ref_prot_center[1] - dna_center[1]
    dz = ref_prot_center[2] - dna_center[2]
    angle = -math.degrees(math.atan2(dz, dy))
    cmd.turn('x', angle)

    # Ensure protein is above DNA in screen-Y
    R = list(cmd.get_view())
    screen_y_dna  = R[1]*dna_center[0]  + R[4]*dna_center[1]  + R[7]*dna_center[2]
    screen_y_prot = R[1]*ref_prot_center[0] + R[4]*ref_prot_center[1] + R[7]*ref_prot_center[2]
    if screen_y_prot < screen_y_dna:
        cmd.turn('x', 180)

    # Save the rotation matrix (first 9 values of view tuple)
    oriented_view = list(cmd.get_view())
    rotation = oriented_view[:9]

    # Zoom to fit ALL states with generous buffer
    # state=0 means consider all states when computing bounding box
    cmd.zoom('all', state=0, buffer=15)

    # Shift DNA toward bottom of frame
    cmd.move('y', 12)

    # Get the zoomed view and replace rotation with our DNA-oriented rotation
    zoomed_view = list(cmd.get_view())
    final_view = list(rotation) + zoomed_view[9:]
    cmd.set_view(final_view)

    return tuple(final_view)


def build_pse(sim_name, sim_info):
    """Build a PSE session with trajectory movie."""
    print(f"\n{'='*60}")
    print(f"Building PSE: {sim_name}")
    print(f"{'='*60}")

    movie_dir = os.path.join(INDIVIDUAL, sim_name, 'viz', 'movie')
    top_file = os.path.join(movie_dir, 'allatom_top.pdb')
    dcd_file = os.path.join(movie_dir, 'allatom.dcd')

    if not os.path.exists(dcd_file):
        print(f"  SKIP: no {dcd_file}")
        return None

    cmd.delete('all')

    # Load reference structure
    cmd.load(REF_PDB, 'ref')
    print(f"  Loaded reference: {REF_PDB}")

    # Load topology then trajectory (state=1 overwrites PDB coords with DCD frame 1)
    cmd.load(top_file, 'mol')
    cmd.load_traj(dcd_file, 'mol', state=1)
    n_states = cmd.count_states('mol')
    print(f"  Loaded trajectory: {n_states} states, {cmd.count_atoms('mol')} atoms")

    # Define selections on mol
    setup_selections(sim_info)

    # Align state 1 of mol to the 3g73 reference via FOXP1 d3 CAs
    result = cmd.align('FOXP1_d3 and name CA', 'ref_protein and name CA', mobile_state=1)
    if result:
        print(f"  Aligned state 1 to reference (RMSD: {result[0]:.2f} A, {result[1]} atoms)")

    # Align all states to state 1 via FOXP1 d3 CAs
    # intra_fit uses the selection atoms for fitting but transforms ALL atoms in the object
    cmd.intra_fit('FOXP1_d3 and name CA')
    print(f"  Aligned all {n_states} states via intra_fit(FOXP1_d3)")

    # Split d2/d3 into separate objects (preserves all states)
    split_surface_domains()

    # Set up representations and coloring
    setup_representations(sim_info)

    # Compute view: DNA at bottom, zoomed to fit all states
    compute_locked_view_all_states()
    print("  Computed view (all states visible)")

    # Set up movie timeline
    cmd.mset(f'1 -{n_states}')
    cmd.set('movie_loop', 1)
    cmd.frame(1)
    print(f"  Movie set: {n_states} frames (looping)")

    # Save PSE
    pse_path = os.path.join(OUTPUT, f'{sim_name}.pse')
    cmd.save(pse_path)
    size_mb = os.path.getsize(pse_path) / 1e6
    print(f"  Saved: {pse_path} ({size_mb:.1f} MB)")

    return pse_path


if __name__ == '__main__':
    os.environ['PYMOL_FEEDBACK_LEVEL'] = '0'
    import pymol
    from pymol import cmd
    pymol.finish_launching(['pymol', '-cqQ'])
    cmd.feedback('disable', 'all', 'everything')

    pse_files = []
    for sim_name, sim_info in SIMS.items():
        pse = build_pse(sim_name, sim_info)
        if pse:
            pse_files.append(pse)

    cmd.quit()

    print(f"\n{'='*60}")
    print("PSE files saved:")
    for f in pse_files:
        print(f"  {f}")
    print(f"\n--- To export movie locally in PyMOL ---")
    print(f"")
    print(f"1. Open the PSE:")
    print(f"   pymol {pse_files[0] if pse_files else '<file>.pse'}")
    print(f"")
    print(f"2. Press the Play button to preview the movie")
    print(f"")
    print(f"3. Export ray-traced PNG frames (in PyMOL command line):")
    print(f"   set ray_trace_frames, 1")
    print(f"   mpng ~/Desktop/movie_frames/frame_, width=1920, height=1080")
    print(f"")
    print(f"   Or for quick non-ray-traced export:")
    print(f"   mpng ~/Desktop/movie_frames/frame_, width=1920, height=1080")
    print(f"")
    print(f"4. Convert PNG frames to MP4 (in terminal):")
    print(f"   ffmpeg -r 15 -i ~/Desktop/movie_frames/frame_%04d.png \\")
    print(f"     -c:v libx264 -pix_fmt yuv420p -crf 18 movie.mp4")
    print(f"{'='*60}")
