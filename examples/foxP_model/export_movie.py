#!/usr/bin/env python
"""
Export a movie from the FOXP two-chain coiled-coil simulation.

Loads the all-atom trajectory directly, sets up PyMOL scene
(same as save_pse_movies.py), and exports frames with 360-degree
rotation pauses at frames 1, 450, and 1000.

Movie structure:
  1. Frame 1: still + 360 rotation
  2. Frames 1-449: normal playback
  3. Frame 450: still + 360 rotation
  4. Frames 450-999: normal playback
  5. Frame 1000: still + 360 rotation

Usage:
    conda run -n pymol-render python export_movie.py
"""

import os
import sys
import time
import math
import subprocess
import tempfile
import glob

# --- Configuration ---
BASE = os.path.dirname(os.path.abspath(__file__))
INDIVIDUAL = os.path.join(BASE, 'individual_chains')
SIM_NAME = 'two_chains_coils_interchain'
MOVIE_DIR = os.path.join(INDIVIDUAL, SIM_NAME, 'viz', 'movie')
REF_PDB = os.path.join(INDIVIDUAL, '3g73_DBDA.pdb')
OUTPUT_DIR = os.path.join(BASE, 'rendered_movies')
os.makedirs(OUTPUT_DIR, exist_ok=True)

ALLATOM_TOP = os.path.join(MOVIE_DIR, 'allatom_top.pdb')
ALLATOM_DCD = os.path.join(MOVIE_DIR, 'allatom.dcd')

OUTPUT_MP4 = os.path.join(OUTPUT_DIR, f'{SIM_NAME}.mp4')
FPS = 60
WIDTH = 1920
HEIGHT = 1080

# Frames where we pause and do a 360-degree rotation (1-indexed PyMOL states)
ROTATION_FRAMES = [1, 450, 1000]
ROTATION_STEPS = 120  # frames per 360 rotation (= 2 seconds at 60fps)

# Domain ranges (1-indexed residue numbers, within each chain)
D1 = (118, 199)
D2 = (299, 375)
D3 = (451, 551)
N_FOXP4 = 677


def setup_scene(cmd):
    """Load trajectory, set up selections, representations, and view."""
    # --- Settings ---
    cmd.set('defer_builds_mode', 3)
    cmd.set('cache_frames', 0)
    cmd.set('async_builds', 0)
    cmd.set('cartoon_sampling', 7)
    cmd.set('cartoon_smooth_loops', 1)
    cmd.set('cartoon_oval_length', 1.0)
    cmd.set('surface_quality', 0)
    cmd.set('surface_color_smoothing', 1)
    cmd.set('ray_opaque_background', 1)
    cmd.set('antialias', 2)
    cmd.set('ray_trace_mode', 0)
    cmd.set('ray_shadows', 1)
    cmd.set('ray_shadow_decay_factor', 0.1)
    cmd.set('spec_reflect', 0.4)
    cmd.set('specular', 0.3)
    cmd.set('ambient', 0.35)
    cmd.set('direct', 0.6)
    cmd.bg_color('white')

    # --- Load reference ---
    if os.path.exists(REF_PDB):
        cmd.load(REF_PDB, 'ref')
        print(f"Loaded reference: {REF_PDB}")
    else:
        print(f"WARNING: Reference PDB not found: {REF_PDB}")

    # --- Load trajectory ---
    print(f"Loading trajectory...")
    cmd.load(ALLATOM_TOP, 'mol')
    cmd.load_traj(ALLATOM_DCD, 'mol', state=1)
    n_states = cmd.count_states('mol')
    print(f"Loaded {n_states} states, {cmd.count_atoms('mol')} atoms")

    # --- Selections (single-chain topology) ---
    OFF = N_FOXP4
    cmd.select('FOXP4',     f'mol and resi 1-{N_FOXP4}')
    cmd.select('FOXP4_d1',  f'mol and resi {D1[0]}-{D1[1]}')
    cmd.select('FOXP4_d2',  f'mol and resi {D2[0]}-{D2[1]}')
    cmd.select('FOXP4_d3',  f'mol and resi {D3[0]}-{D3[1]}')
    cmd.select('FOXP4_idr', f'mol and resi 1-{N_FOXP4} and not (resi {D1[0]}-{D1[1]} or resi {D2[0]}-{D2[1]} or resi {D3[0]}-{D3[1]})')

    cmd.select('FOXP1',     f'mol and resi {OFF+1}-{OFF+680}')
    cmd.select('FOXP1_d1',  f'mol and resi {D1[0]+OFF}-{D1[1]+OFF}')
    cmd.select('FOXP1_d2',  f'mol and resi {D2[0]+OFF}-{D2[1]+OFF}')
    cmd.select('FOXP1_d3',  f'mol and resi {D3[0]+OFF}-{D3[1]+OFF}')
    cmd.select('FOXP1_idr', f'mol and resi {OFF+1}-{OFF+680} and not (resi {D1[0]+OFF}-{D1[1]+OFF} or resi {D2[0]+OFF}-{D2[1]+OFF} or resi {D3[0]+OFF}-{D3[1]+OFF})')

    if os.path.exists(REF_PDB):
        cmd.select('ref_protein', 'ref and chain A and polymer.protein')
        cmd.select('ref_DNA', 'ref and (chain C or chain D) and polymer.nucleic')
    cmd.deselect()

    # --- Align to reference ---
    if os.path.exists(REF_PDB):
        result = cmd.align('FOXP1_d3 and name CA', 'ref_protein and name CA', mobile_state=1)
        if result:
            print(f"Aligned state 1 to reference (RMSD: {result[0]:.2f} A, {result[1]} atoms)")
        cmd.intra_fit('FOXP1_d3 and name CA')
        print(f"Aligned all {n_states} states via intra_fit(FOXP1_d3)")

    # --- Split d2/d3 into separate objects for clean surfaces ---
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

    # --- Representations ---
    cmd.hide('everything')

    # Reference DNA
    if os.path.exists(REF_PDB):
        cmd.show('cartoon', 'ref_DNA')
        cmd.color('gray60', 'ref_DNA')
        cmd.set('cartoon_ring_mode', 3)
        cmd.set('cartoon_ring_finder', 1)
        cmd.hide('everything', 'ref and chain A')
        cmd.hide('everything', 'ref and chain B')

    # mol: IDR + d1 only
    cmd.show('cartoon', 'mol')
    cmd.color('limegreen', f'mol and resi 1-{N_FOXP4}')
    cmd.color('lightorange', f'mol and resi {N_FOXP4+1}-9999')
    cmd.alter('mol', "ss='L'")
    cmd.rebuild('mol')
    cmd.set('cartoon_loop_radius', 0.10, 'mol')
    cmd.set('cartoon_transparency', 0.3, 'mol')

    # Extracted domain objects
    for protein in ['FOXP1', 'FOXP4']:
        for domain in ['d2', 'd3']:
            obj = f'{protein}_{domain}_obj'
            cmd.show('cartoon', obj)
            cmd.show('surface', obj)

    # --- Camera: DNA at bottom, zoomed to fit all states ---
    if os.path.exists(REF_PDB):
        cmd.center('ref_DNA')
        cmd.orient('ref_DNA')

        dna_center = cmd.centerofmass('ref_DNA')
        ref_prot_center = cmd.centerofmass('ref and chain A and polymer.protein')

        dy = ref_prot_center[1] - dna_center[1]
        dz = ref_prot_center[2] - dna_center[2]
        angle = -math.degrees(math.atan2(dz, dy))
        cmd.turn('x', angle)

        R = list(cmd.get_view())
        screen_y_dna  = R[1]*dna_center[0]  + R[4]*dna_center[1]  + R[7]*dna_center[2]
        screen_y_prot = R[1]*ref_prot_center[0] + R[4]*ref_prot_center[1] + R[7]*ref_prot_center[2]
        if screen_y_prot < screen_y_dna:
            cmd.turn('x', 180)

        oriented_view = list(cmd.get_view())
        rotation = oriented_view[:9]
        cmd.zoom('all', state=0, buffer=15)
        cmd.move('y', 12)
        zoomed_view = list(cmd.get_view())
        final_view = list(rotation) + zoomed_view[9:]
        cmd.set_view(final_view)
    else:
        cmd.zoom('all', state=0, buffer=15)

    cmd.viewport(WIDTH, HEIGHT)
    return n_states


def export_frames(cmd, n_states, tmpdir):
    """Export PNG frames with 360-degree rotation pauses at specified frames.

    Returns the total number of exported frames.
    """
    frame_num = 0
    base_view = list(cmd.get_view())

    # Sort rotation frames
    rot_frames = sorted(ROTATION_FRAMES)

    # Build segments: list of (start_state, end_state) with rotation pauses between
    # Frame 1 rotation, then 1->449, frame 450 rotation, then 450->999, frame 1000 rotation
    segments = []
    prev = 1
    for rf in rot_frames:
        if rf > prev:
            segments.append(('play', prev, rf - 1))
        segments.append(('rotate', rf))
        prev = rf + 1
    if prev <= n_states:
        segments.append(('play', prev, n_states))

    total_frames = 0
    for seg in segments:
        if seg[0] == 'play':
            total_frames += seg[2] - seg[1] + 1
        else:
            total_frames += ROTATION_STEPS
    print(f"Total frames to export: {total_frames}")

    for seg in segments:
        if seg[0] == 'rotate':
            state = seg[1]
            print(f"  Rotation at state {state} ({ROTATION_STEPS} frames)...")
            cmd.frame(state)
            cmd.refresh()
            # Save the current view, rotate around y-axis
            rot_base_view = list(cmd.get_view())
            for i in range(ROTATION_STEPS):
                angle = (360.0 / ROTATION_STEPS) * i
                cmd.set_view(rot_base_view)
                cmd.turn('y', angle)
                cmd.refresh()
                frame_num += 1
                fname = os.path.join(tmpdir, f"frame_{frame_num:06d}.png")
                cmd.png(fname, width=WIDTH, height=HEIGHT, ray=0, quiet=1)
                if frame_num % 30 == 0:
                    print(f"    Frame {frame_num}/{total_frames}")
            # Restore view after rotation
            cmd.set_view(rot_base_view)

        elif seg[0] == 'play':
            start_state, end_state = seg[1], seg[2]
            print(f"  Playing states {start_state}-{end_state}...")
            for state in range(start_state, end_state + 1):
                cmd.frame(state)
                cmd.refresh()
                frame_num += 1
                fname = os.path.join(tmpdir, f"frame_{frame_num:06d}.png")
                cmd.png(fname, width=WIDTH, height=HEIGHT, ray=0, quiet=1)
                if frame_num % 100 == 0:
                    print(f"    Frame {frame_num}/{total_frames}")

    print(f"Exported {frame_num} PNG frames.")
    return frame_num


def encode_mp4(tmpdir, total_frames):
    """Encode PNG frames to MP4 with ffmpeg."""
    pattern = os.path.join(tmpdir, "frame_%06d.png")

    ffmpeg_cmd = [
        "ffmpeg",
        "-y",
        "-framerate", str(FPS),
        "-start_number", "1",
        "-i", pattern,
        "-frames:v", str(total_frames),
        "-c:v", "libx264",
        "-preset", "medium",
        "-crf", "18",
        "-pix_fmt", "yuv420p",
        "-movflags", "+faststart",
        OUTPUT_MP4,
    ]
    print(f"Encoding MP4 at {FPS} fps...")
    print(f"Running: {' '.join(ffmpeg_cmd)}")

    t0 = time.time()
    result = subprocess.run(ffmpeg_cmd, capture_output=True, text=True)
    elapsed = time.time() - t0

    if result.returncode != 0:
        print("ffmpeg STDERR:")
        print(result.stderr)
        print("ERROR: ffmpeg failed.")
        sys.exit(1)

    print(f"Encoding complete in {elapsed:.1f} seconds.")
    sz = os.path.getsize(OUTPUT_MP4) / (1024 * 1024)
    print(f"Output: {OUTPUT_MP4} ({sz:.1f} MB)")


def main():
    # Verify input files exist
    for f in [ALLATOM_TOP, ALLATOM_DCD]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found")
            sys.exit(1)

    # Create temp directory for frames
    tmpdir = tempfile.mkdtemp(prefix="pymol_movie_frames_")
    print(f"Temporary frame directory: {tmpdir}")

    # Launch PyMOL headless
    print("Loading PyMOL...")
    import pymol
    from pymol import cmd
    pymol.finish_launching(["pymol", "-cqQ"])
    cmd.feedback('disable', 'all', 'everything')

    try:
        # Set up the full scene
        t0 = time.time()
        n_states = setup_scene(cmd)
        print(f"Scene setup in {time.time() - t0:.1f} seconds.")

        # Export frames with rotation pauses
        t0 = time.time()
        total_frames = export_frames(cmd, n_states, tmpdir)
        elapsed = time.time() - t0
        print(f"PNG export complete in {elapsed:.1f} seconds.")

    finally:
        cmd.quit()

    # Encode to MP4
    encode_mp4(tmpdir, total_frames)

    # Cleanup temp PNGs
    print("Cleaning up temporary PNG files...")
    for f in sorted(glob.glob(os.path.join(tmpdir, "*.png"))):
        os.remove(f)
    os.rmdir(tmpdir)
    print("Done.")


if __name__ == "__main__":
    main()
