#!/usr/bin/env python3
"""
Render the top-ranked accepted binder for each PARP14 BindCraft target to a PNG.

Run headless with PyMOL (the pymol-render conda env already on this system):
    conda activate pymol-render
    pymol -cq render_top_binders.py                # all targets with an accepted design
    pymol -cq render_top_binders.py -- md1_block_af3  # just one target

Colors target chain A and the designed binder (all other chains) with the same
two categorical colors analyze_campaign.py uses for block/clamp, for visual
consistency across the campaign's figures.
"""
import csv
import glob
import os
import sys

from pymol import cmd

# __file__ does not resolve to this script under `pymol -cq` (it resolves
# into pymol's own package dir instead) -- hardcode the campaign dir.
HERE = "/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md"
DESIGNS_DIR = os.path.join(HERE, "designs")
SUMMARY_CSV = os.path.join(HERE, "targets_summary.csv")
OUT_DIR = os.path.join(HERE, "figures", "renders")

COLOR_TARGET = "0x898781"   # muted ink -- the fixed target
COLOR_BINDER = "0xeb6834"   # slot 2 orange -- the designed binder


def load_target_names():
    with open(SUMMARY_CSV) as fh:
        return [row["target"] for row in csv.DictReader(fh)]


def top_ranked_pdb(target):
    ranked_dir = os.path.join(DESIGNS_DIR, target, "Accepted", "Ranked")
    if not os.path.isdir(ranked_dir):
        return None
    hits = sorted(glob.glob(os.path.join(ranked_dir, "1_*.pdb")))
    if hits:
        return hits[0]
    # fall back to lowest-numbered rank present
    hits = sorted(glob.glob(os.path.join(ranked_dir, "*.pdb")))
    return hits[0] if hits else None


def render_one(target, pdb_path, out_path):
    cmd.reinitialize()
    cmd.load(pdb_path, target)
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.bg_color("white")
    cmd.color(COLOR_TARGET, "chain A")
    cmd.color(COLOR_BINDER, "not chain A")
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("ray_opaque_background", 1)
    cmd.orient()
    cmd.ray(1200, 1000)
    cmd.png(out_path, dpi=200)
    print(f"Wrote {out_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    # anything after `--` on the pymol -cq command line
    argv = sys.argv[1:]
    requested = argv if argv else None

    targets = load_target_names()
    if requested:
        targets = [t for t in targets if t in requested]

    for target in targets:
        pdb_path = top_ranked_pdb(target)
        if pdb_path is None:
            print(f"{target:24s} -- no accepted/ranked design yet, skipping")
            continue
        out_path = os.path.join(OUT_DIR, f"{target}_top1.png")
        render_one(target, pdb_path, out_path)


main()
