#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a permanent reference file of PyMOL-ready residue selections:

  1. Domain boundary ranges (energy-analysis FL_DOMAINS, remapped) for every
     finished CALVADOS simulation set that contains a macrodomain unit
     (MD1L1/MD2/MD3) or is full-length (fl/fl_optimized) -- i.e. all 8
     named sets in sim_registry.SETS, since every one of them qualifies.
  2. Active-site residues (catalytic + pocket, from ACTIVE_SITES_FL) for
     MD1/MD2/MD3/ART, remapped into each set's own construct numbering.
  3. Face-classification groups (active/back for MD1L1/MD2/MD3/ART,
     faceA/faceB k-means split for the rest) for fl_optimized and
     md_full_pose specifically, computed from their state_1 representative
     conformation (tica_pipeline/enhance_face_classification_pse.py bakes
     the same groups into each state's own .pse as named selections).

Run whenever sim_registry.py's domain/active-site tables change, or a new
named set finishes simulating.

Usage:
    conda run -n pymol-render python generate_residue_selections.py
"""
import sys
import time
from pathlib import Path

CWD = Path(__file__).resolve().parent
ROOT = CWD.parent
sys.path.insert(0, str(CWD))
sys.path.insert(0, str(ROOT))

import pymol
from pymol import cmd
pymol.finish_launching(['pymol', '-cq'])

import surface_gallery_lib as sgl                                              # noqa: E402
import sim_registry as sr                                                      # noqa: E402
from make_surface_gallery import classify_active_face, classify_generic_face, DOMAIN_TO_SITE  # noqa: E402

OUT_PATH = CWD / 'residue_selections.md'

# All 9 named sets qualify as "MDs or fl": each contains >=1 macrodomain
# unit (md1l1/md2/md3) and/or is the full-length numbering (fl/fl_optimized/fl_go).
ALL_SETS = ['fl', 'fl_optimized', 'fl_go', 'norrm', 'noart', 'core', 'mka', 'md', 'md3art']

REPRESENTATIVE = {
    'fl_optimized': ROOT / 'figures/09_surface_gallery/2026-08-07/fl_optimized_ca25_tica/state_1/_apbs_work/centered.pdb',
    'md': ROOT / 'figures/09_surface_gallery/2026-08-07/md_full_pose/state_1/_apbs_work/centered.pdb',
}


def plus_join(resids):
    return "+".join(str(r) for r in sorted(resids))


def domain_boundaries_section():
    lines = ["## 1. Domain boundaries per finished set\n",
             "Energy-analysis domain ranges (`sim_registry.FL_DOMAINS`), remapped into "
             "each set's own construct numbering. `fl`/`fl_optimized`/`norrm`/`noart`/"
             "`core`/`mka`/`md`/`md3art` are finished (25/25 replicates each); `fl_go` "
             "is a separate campaign of 11 Go-model-restraint extension runs (2us each, "
             "all complete) seeded from 5 fl_optimized TICA-clustered states -- same "
             "full-length (1801 res) numbering as `fl_optimized`, different restraint "
             "type (`custom_restraint_type: go`, see fl_go/input/custom_restraints_go.txt) "
             "not domain boundaries. Every set here contains a macrodomain unit and/or "
             "is full-length.\n"]
    for set_key in ALL_SETS:
        info = sr.SETS[set_key]
        lines.append(f"\n### `{set_key}` -- {info['label']} (sysname `{info['sysname']}`)\n")
        domains = sr.get_construct_domains_for_set(set_key)
        md_or_art = {d: r for d, r in domains.items() if d in ('MD1', 'MD2', 'MD3', 'ART')}
        if not md_or_art:
            lines.append("_(no MD1/MD2/MD3/ART unit in this construct)_\n")
            continue
        for dname, (s, e) in md_or_art.items():
            lines.append(f"- `{dname}`: resi {s}-{e}  ->  `select {set_key}_{dname}, resi {s}-{e}`\n")
    return lines


def active_site_section():
    lines = ["\n## 2. Active-site residues per finished set\n",
             "Catalytic + pocket residues (`sim_registry.ACTIVE_SITES_FL`), remapped "
             "into each set's own construct numbering via "
             "`sim_registry.get_active_sites_for_set()`.\n",
             "\n**Note:** a bug in `get_active_sites_for_set()`/`get_construct_domains_for_set()` "
             "was found and fixed while generating this file -- both used to compress out the "
             "13-residue MD2/MD3 linker gap (FL resi 1194-1206) for `fl_optimized`, which shifted "
             "every MD3/ART residue number by -13 for that set only (`fl` was already correct; "
             "`md`/`norrm`/`noart`/`core`/`mka`/`md3art` were unaffected since they are real "
             "sub-constructs where that gap is genuinely absent). Fixed by special-casing "
             "`fl_optimized` alongside `fl` (both are the full uncompressed 1801-residue sequence).\n"]
    for set_key in ALL_SETS:
        sites = sr.get_active_sites_for_set(set_key)
        if not sites:
            continue
        lines.append(f"\n### `{set_key}`\n")
        for site in sr.SITE_NAMES:
            if site not in sites:
                continue
            data = sites[site]
            lines.append(f"- **{site}** catalytic ({len(data['catalytic'])}): "
                          f"`select {set_key}_{site}_catalytic, resi {plus_join(data['catalytic'])}`")
            lines.append(f"  \n  pocket ({len(data['pocket'])}): "
                          f"`select {set_key}_{site}_pocket, resi {plus_join(data['pocket'])}`\n")
    return lines


def face_classification_section():
    lines = ["\n## 3. Face-classification residue groups\n",
             "Computed from each construct's `state_1` representative conformation "
             "(same method `make_surface_gallery.py` uses to color the face-classification "
             "surface, and that `enhance_face_classification_pse.py` bakes into each state's "
             ".pse as named selections). **Conformation-dependent** -- the k-means "
             "(`faceA`/`faceB`) splits especially will shift somewhat state to state; "
             "the active/back dot-product split for MD1L1/MD2/MD3/ART is more stable "
             "but still recomputed per state.\n"]
    for set_key, pdb_path in REPRESENTATIVE.items():
        if not pdb_path.exists():
            lines.append(f"\n### `{set_key}`: representative PDB not found ({pdb_path}), skipped\n")
            continue
        cmd.reinitialize()
        cmd.load(str(pdb_path), 'mol')
        cmd.remove('solvent')
        domains = sgl.get_domain_ranges_for_set(set_key)
        lines.append(f"\n### `{set_key}` (from state_1)\n")
        for i, (dname, (color_name, base_rgb, (cs, ce))) in enumerate(domains.items()):
            site = DOMAIN_TO_SITE.get(dname)
            faces = {}
            method = None
            if site:
                cat, _ = sgl.get_pocket_for_set(set_key, site)
                if cat:
                    faces = classify_active_face('mol', (cs, ce), cat)
                    method = 'active/back (dot-product toward pocket)'
            if not faces:
                faces = classify_generic_face('mol', (cs, ce), seed=42 + i)
                method = 'faceA/faceB (k-means)'
            if not faces:
                continue
            groups = {}
            for r, f in faces.items():
                groups.setdefault(f, []).append(r)
            lines.append(f"\n**{dname}** (resi {cs}-{ce}) -- {method}\n")
            # PyMOL selection/object names can't contain '-' (e.g. 'KH1-KH6'),
            # so sanitize for the actual `select` command; keep the
            # human-readable hyphenated form in the surrounding label text.
            safe_dname = dname.replace('-', '_')
            for fname, resids in sorted(groups.items()):
                lines.append(f"- `{dname}_{fname}` ({len(resids)} res): "
                              f"`select {set_key}_{safe_dname}_{fname}, resi {plus_join(resids)}`")
    return lines


def main():
    lines = [
        "# PARP14 residue selections (PyMOL-ready)\n",
        f"\nGenerated {time.strftime('%Y-%m-%d %H:%M')} by "
        "`tica_pipeline/generate_residue_selections.py`. Re-run after any change to "
        "`sim_registry.py`'s DOMAIN_UNITS/ACTIVE_SITES_FL/FL_DOMAINS tables, or after a "
        "new named set finishes simulating.\n",
    ]
    lines += domain_boundaries_section()
    lines += active_site_section()
    lines += face_classification_section()

    OUT_PATH.write_text("\n".join(lines) + "\n")
    print(f"Saved: {OUT_PATH}")
    cmd.quit()


if __name__ == '__main__':
    main()
