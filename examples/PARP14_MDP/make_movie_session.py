#!/usr/bin/env python3
"""
Build a PyMOL session with 1000-frame FL CALVADOS movie (last 10 ns).

Colors:
  RRM1-3:   yellow
  KH1-6, KH7a, KHb-KH8: palegreen
  MD1:      purple (dark)
  MD2:      purple (medium)
  MD3:      purple (light)
  WWE:      gray
  ART:      red
  Linkers:  white

All states aligned on ART domain, ART centered in view.

Usage:
    conda run -n pymol-render python make_movie_session.py
"""

import os

CWD = os.path.dirname(os.path.abspath(__file__))
MOVIE_PDB = os.path.join(CWD, 'pymol_frames', 'fl_movie.pdb')
SESSION_OUT = os.path.join(CWD, 'fl_movie.pse')

# Domain residue ranges (FL 1-based)
DOMAINS = {
    'RRM':    [(1, 145), (146, 224), (225, 314)],
    'KH':     [(315, 737), (738, 789), (1389, 1533)],
    'MD1':    [(790, 1004)],
    'MD2':    [(1004, 1193)],
    'MD3':    [(1207, 1388)],
    'WWE':    [(1534, 1602)],
    'ART':    [(1603, 1801)],
}

# ART range for alignment
ART_START, ART_END = 1603, 1801


def main():
    import pymol
    from pymol import cmd

    pymol.finish_launching(['pymol', '-cq'])

    print("Loading multi-model PDB (1000 states)...")
    cmd.load(MOVIE_PDB, 'fl_movie')
    n_states = cmd.count_states('fl_movie')
    print(f"  Loaded: {n_states} states")

    # Show as surface
    cmd.hide('everything', 'fl_movie')
    cmd.show('cartoon', 'fl_movie')
    cmd.set('cartoon_trace_atoms', 1)
    cmd.set('cartoon_tube_radius', 0.3)
    cmd.set('ribbon_trace_atoms', 1)

    # Base color: white
    cmd.color('white', 'fl_movie')

    # Domain coloring
    print("Coloring domains...")

    # RRM1-3: yellow
    for s, e in DOMAINS['RRM']:
        cmd.color('yellow', f'fl_movie and resi {s}-{e}')

    # KH domains: palegreen
    for s, e in DOMAINS['KH']:
        cmd.color('palegreen', f'fl_movie and resi {s}-{e}')

    # MD1: dark purple
    cmd.set_color('md1_purple', [0.40, 0.10, 0.50])
    for s, e in DOMAINS['MD1']:
        cmd.color('md1_purple', f'fl_movie and resi {s}-{e}')

    # MD2: medium purple
    cmd.set_color('md2_purple', [0.58, 0.30, 0.68])
    for s, e in DOMAINS['MD2']:
        cmd.color('md2_purple', f'fl_movie and resi {s}-{e}')

    # MD3: light purple
    cmd.set_color('md3_purple', [0.72, 0.50, 0.80])
    for s, e in DOMAINS['MD3']:
        cmd.color('md3_purple', f'fl_movie and resi {s}-{e}')

    # WWE: gray
    for s, e in DOMAINS['WWE']:
        cmd.color('gray', f'fl_movie and resi {s}-{e}')

    # ART: red
    for s, e in DOMAINS['ART']:
        cmd.color('red', f'fl_movie and resi {s}-{e}')

    # Align all states to state 1 on the ART domain
    print("Aligning all states on ART domain...")
    art_sel = f'fl_movie and resi {ART_START}-{ART_END}'
    cmd.intra_fit(art_sel, state=0)
    print(f"  intra_fit done ({n_states} states)")

    # Center on ART
    cmd.zoom(art_sel, buffer=30, state=1)

    # Settings
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('ray_shadow', 0)
    cmd.set('state', 1)

    # Save session
    cmd.save(SESSION_OUT)
    print(f"\nSession saved: {SESSION_OUT}")
    print(f"  {n_states} states, aligned on ART (res {ART_START}-{ART_END})")
    print(f"  Open in PyMOL, press play to animate")

    cmd.quit()


if __name__ == '__main__':
    main()
