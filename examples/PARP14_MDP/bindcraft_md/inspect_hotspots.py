#!/usr/bin/env python3
"""
Render a target's hotspot residues on its cartoon so you can eyeball whether
they form one coherent, reachable patch (vs. buried or scattered).

Run headless with PyMOL (the pymol-render conda env already on this system):
    conda activate pymol-render
    pymol -cq inspect_hotspots.py                  # all 11 targets
    pymol -cq inspect_hotspots.py -- md1_block_af3  # just one target

Target cartoon in muted gray, hotspot residues as orange sticks + spheres,
labeled by residue number. PNGs land in figures/hotspot_renders/.
"""
import csv
import os
import sys

from pymol import cmd

HERE = "/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md"
SUMMARY_CSV = os.path.join(HERE, "targets_summary.csv")
OUT_DIR = os.path.join(HERE, "figures", "hotspot_renders")

COLOR_TARGET = "0x898781"   # muted ink
COLOR_HOTSPOT = "0xeb6834"  # slot 2 orange


def load_targets():
    with open(SUMMARY_CSV) as fh:
        rows = list(csv.DictReader(fh))
    return {r["target"]: r for r in rows}


def render_one(target, pdb_path, hotspots, out_path):
    cmd.reinitialize()
    cmd.load(pdb_path, target)
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.bg_color("white")
    cmd.color(COLOR_TARGET, "polymer")

    sel = " or ".join(f"resi {r}" for r in hotspots)
    cmd.select("hotspots", sel)
    cmd.show("sticks", "hotspots and not name C+N+O")
    cmd.show("spheres", "hotspots and name CA")
    cmd.set("sphere_scale", 0.4, "hotspots")
    cmd.color(COLOR_HOTSPOT, "hotspots")
    cmd.label("hotspots and name CA", '"%s" % resi')
    cmd.set("label_size", 14)
    cmd.set("label_color", "black")
    cmd.set("ray_opaque_background", 1)

    cmd.orient("hotspots")
    cmd.zoom("hotspots", buffer=8)
    cmd.ray(1200, 1000)
    cmd.png(out_path, dpi=200)
    print(f"Wrote {out_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_targets = load_targets()
    argv = sys.argv[1:]
    names = argv if argv else list(all_targets.keys())

    for name in names:
        row = all_targets.get(name)
        if row is None:
            print(f"{name} -- not in targets_summary.csv, skipping")
            continue
        hotspots = row["hotspots"].split()
        out_path = os.path.join(OUT_DIR, f"{name}_hotspots.png")
        render_one(name, row["pdb"], hotspots, out_path)


main()
