#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Render the vacuum-electrostatics surface (whole molecule + pocket close-ups)
in a freshly-launched, otherwise-untouched PyMOL process.

Why this is its own process: util.protein_vacuum_esp's derived "_e_pot"
object reproducibly fails PyMOL's selector registration -- silently, no
Python exception -- if ANY other cmd.set() call (any setting, not just
surface-related ones) has run in the same session first. Domain coloring,
ray_opaque_background, antialias... even one extra cmd.set() is enough.
Verified by direct A/B repro (see make_surface_gallery.py's ensure_esp()
docstring for the shared context). The only reliable fix found is to give
protein_vacuum_esp a session where NOTHING else has touched cmd.set()
before it runs.

Usage:
    conda run -n pymol-render python _render_esp_isolated.py \\
        --pdb state_1_allatom.pdb --out-dir OUTDIR \\
        --pockets '{"MD1": [822,823,...], "MD2": [...]}'
"""
import argparse
import json
from pathlib import Path

import pymol
from pymol import cmd, util
pymol.finish_launching(['pymol', '-cq'])

IMG_SIZE = (900, 900)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pdb', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--pockets', default='{}',
                     help='JSON {site_name: [pocket_resi, ...]}')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pockets = json.loads(args.pockets)

    # Nothing but these two calls before protein_vacuum_esp -- see module
    # docstring. Do not add cmd.set()/cmd.bg_color() calls above this line.
    cmd.set('surface_quality', 0)
    cmd.bg_color('white')

    cmd.load(args.pdb, 'mol')
    cmd.remove('solvent')
    cmd.hide('everything', 'mol')
    cmd.show('cartoon', 'mol')
    # Priming ray-trace: every A/B repro that succeeded had rendered
    # *something* with cmd.ray() before calling protein_vacuum_esp; every
    # attempt that skipped straight to protein_vacuum_esp with no prior
    # cmd.ray() failed. Not fully explained, but reproducible -- keep it.
    cmd.orient('mol')
    cmd.ray(300, 300)

    cmd.hide('everything')
    util.protein_vacuum_esp('mol', mode=2, quiet=1, _self=cmd)
    esp_obj = 'mol_e_pot'

    cmd.orient(esp_obj)
    cmd.zoom(esp_obj, buffer=3)
    cmd.ray(*IMG_SIZE)
    cmd.png(str(out_dir / 'electrostatic_surface_front.png'), dpi=200)
    cmd.turn('y', 180)
    cmd.ray(*IMG_SIZE)
    cmd.png(str(out_dir / 'electrostatic_surface_back.png'), dpi=200)
    print(f"saved electrostatic_surface_front/back.png -> {out_dir}", flush=True)

    for site, pocket in pockets.items():
        if not pocket:
            continue
        pocket_sel_str = "+".join(map(str, pocket))
        site_dir = out_dir / f'pocket_{site}'
        site_dir.mkdir(parents=True, exist_ok=True)
        sel = f'{esp_obj} and resi {pocket_sel_str}'
        cmd.orient(sel)
        cmd.zoom(sel, buffer=4)
        cmd.ray(*IMG_SIZE)
        cmd.png(str(site_dir / 'electrostatic.png'), dpi=200)
        print(f"saved pocket_{site}/electrostatic.png", flush=True)

    cmd.save(str(out_dir / '_esp.pse'))
    print("DONE", flush=True)
    cmd.quit()


if __name__ == '__main__':
    main()
