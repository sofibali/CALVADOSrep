###############################################################################
#  pymol_movie_option1_pull.pml
#
#  Movie of the option1_pull_interchain_d12 simulation (full atom).
#  FOXP4 (chain A) + FOXP1 (chain B), interchain restraints on d1 & d2,
#  pulling force on domain3.
#
#  All states aligned to state 1 via FOXP4 domain3 (DBD).
#
#  Expects all-atom trajectory in:
#    individual_chains/option1_pull_interchain_d12/viz/movie/allatom_top.pdb
#    individual_chains/option1_pull_interchain_d12/viz/movie/allatom.dcd
#
#  Run with GUI:   pymol pymol_movie_option1_pull.pml
#  Render headless: pymol -cq pymol_movie_option1_pull.pml
###############################################################################

# --- 0. SETTINGS ---
set defer_builds_mode, 3
set cache_frames, 0
set async_builds, 0
set cartoon_sampling, 7
set cartoon_smooth_loops, 1
set cartoon_loop_radius, 0.15
set cartoon_oval_length, 1.0
set cartoon_rect_length, 1.2
set ray_opaque_background, 1
set antialias, 2
set spec_reflect, 0.3
set ambient, 0.25
set ray_trace_mode, 1
set ray_shadows, 0
bg_color white

python
import os, glob

base = os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else os.getcwd()
sim = "option1_pull_interchain_d12"
movie_dir = os.path.join(base, "individual_chains", sim, "viz", "movie")
frames_dir = os.path.join(base, "individual_chains", sim, "viz", "allatom_frames")
render_dir = os.path.join(base, "individual_chains", sim, "viz", "movie_frames")
os.makedirs(render_dir, exist_ok=True)

allatom_top = os.path.join(movie_dir, "allatom_top.pdb")
allatom_dcd = os.path.join(movie_dir, "allatom.dcd")

# Try loading DCD trajectory first, fall back to individual PDB frames
if os.path.exists(allatom_top) and os.path.exists(allatom_dcd):
    print(f"Loading all-atom trajectory from DCD...")
    cmd.load(allatom_top, "mol")
    cmd.load_traj(allatom_dcd, "mol")
    load_mode = "dcd"
elif os.path.exists(frames_dir):
    pdb_files = sorted(glob.glob(os.path.join(frames_dir, "frame_*_AA.pdb")))
    if pdb_files:
        print(f"Loading {len(pdb_files)} individual all-atom PDB frames...")
        for i, pf in enumerate(pdb_files):
            cmd.load(pf, "mol")
        load_mode = "pdbs"
    else:
        print("ERROR: No all-atom frames found!")
        print(f"  Run: conda run -n CALVADOS python prepare_movie_trajectories.py --sim {sim}")
        load_mode = "none"
else:
    print("ERROR: No all-atom trajectory found!")
    print(f"  Run: conda run -n CALVADOS python prepare_movie_trajectories.py --sim {sim}")
    load_mode = "none"

n_states = cmd.count_states("mol") if load_mode != "none" else 0
print(f"Loaded {n_states} states (mode: {load_mode})")
python end

# --- 1. SELECTIONS ---
# This simulation has 2 chains: A=FOXP4 (677 res), B=FOXP1 (680 res)
# Domain residue numbering is within each chain
select FOXP4,       mol and chain A
select FOXP4_d1,    mol and chain A and resi 118-199
select FOXP4_d2,    mol and chain A and resi 299-375
select FOXP4_d3,    mol and chain A and resi 451-551
select FOXP4_idr,   mol and chain A and not (resi 118-199 or resi 299-375 or resi 451-551)

select FOXP1,       mol and chain B
select FOXP1_d1,    mol and chain B and resi 118-199
select FOXP1_d2,    mol and chain B and resi 299-375
select FOXP1_d3,    mol and chain B and resi 451-551
select FOXP1_idr,   mol and chain B and not (resi 118-199 or resi 299-375 or resi 451-551)

select all_d3,      FOXP4_d3 or FOXP1_d3
select all_d2,      FOXP4_d2 or FOXP1_d2
select all_d1,      FOXP4_d1 or FOXP1_d1
deselect

# --- 2. ALIGN ALL STATES TO STATE 1 ---
python
n_states = cmd.count_states("mol")
if n_states > 1:
    # Align using FOXP4 domain3 (structured DBD)
    for state in range(2, n_states + 1):
        cmd.align(
            "FOXP4_d3 and name CA",
            "FOXP4_d3 and name CA",
            mobile_state=state,
            target_state=1,
            quiet=1
        )
    print(f"Aligned {n_states - 1} states to state 1 on FOXP4_d3 (DBD).")
python end

# --- 3. REPRESENTATIONS ---
hide everything

# Cartoon for everything
show cartoon, mol

# FOXP4 coloring (green palette)
color tv_green, FOXP4_idr
color chartreuse, FOXP4_d1
color splitpea, FOXP4_d2
color forest, FOXP4_d3

# FOXP1 coloring (orange palette)
color lightorange, FOXP1_idr
color tv_orange, FOXP1_d1
color orange, FOXP1_d2
color chocolate, FOXP1_d3

# Domain3 surfaces (semi-transparent)
show surface, FOXP4_d3
show surface, FOXP1_d3
set surface_quality, 1
set transparency, 0.65, FOXP4_d3
set transparency, 0.65, FOXP1_d3

# IDR representation: thinner cartoon
set cartoon_loop_radius, 0.10, FOXP4_idr
set cartoon_loop_radius, 0.10, FOXP1_idr

# --- 4. CAMERA ---
orient all_d3
zoom mol, 10

# --- 5. MOVIE ---
python
n_states = cmd.count_states("mol")
if n_states > 0:
    # Set up movie: one frame per state
    cmd.mset(f"1 -{n_states}")
    print(f"Movie set up: {n_states} frames")
    print("  GUI: press Play or use 'mplay' command")
    print("  Render: see section 6 below or run headless with -cq flag")
python end

# --- 6. RENDER (uncomment or run headless) ---
# To render PNG frames for ffmpeg, uncomment the block below
# or run: pymol -cq pymol_movie_option1_pull.pml
#
# Then combine with:
#   ffmpeg -framerate 15 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p movie.mp4

python
import sys

# Auto-render in headless mode
if not cmd.get_setting_boolean("auto_show_lines"):
    pass  # headless detection is tricky; use explicit flag instead

# Check if we should render (set RENDER=1 before running, or use -cq)
render_mode = os.environ.get("RENDER", "0") == "1"

# Simple check: if launched with -cq, the gui is off
try:
    gui_active = bool(cmd.get_setting_int("internal_gui"))
except:
    gui_active = True

if not gui_active or render_mode:
    n_states = cmd.count_states("mol")
    if n_states > 0:
        print(f"\nRendering {n_states} frames...")
        for i in range(1, n_states + 1):
            cmd.frame(i)
            cmd.refresh()
            fname = os.path.join(render_dir, f"frame_{i:04d}.png")
            cmd.png(fname, width=1920, height=1080, ray=1, dpi=150)
            print(f"  Frame {i}/{n_states} -> {fname}")
        print(f"\nRendered to: {render_dir}")
        print("Combine with: ffmpeg -framerate 15 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p option1_pull_movie.mp4")
python end
