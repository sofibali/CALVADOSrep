"""
Generate a PyMOL PSE file for PARP14 with:
- Full-length AlphaFold2 structure (colored by domain, active sites labeled)
- Largest AF3 construct (25 models as states): rrm1_rrm2_rrm3_kh1-kh6_kh7a_md2_md3_khb-kh8_wwe_art
- MD1L1-to-end AF3 construct (25 models as states): md1l1_md2_md3_khb-kh8_wwe_art
"""
import os
import glob
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-cq'])

# =======================================================
# 1. Full-length AlphaFold2 structure
# =======================================================
pdb_path = "/home/sbali/CALVADOS/examples/PARP14_MDP/input/parp14.pdb"
cmd.load(pdb_path, "PARP14_AF2")

# Domain definitions (from PARP14_domains_atS.fasta, cut at serine residues)
domains = {
    "RRM1":   (1, 148),
    "RRM2-3": (149, 304),
    "KH_N":   (305, 516),
    "KH_Ca":  (517, 789),
    "MD1":    (790, 1006),
    "MD2":    (1007, 1188),
    "MD3":    (1189, 1398),
    "KH_b":   (1399, 1508),
    "WWE":    (1509, 1609),
    "ART":    (1610, 1801),
}

# Colors matching analysis figures (03_longest_construct_figure.py)
domain_colors_rgb = {
    "RRM1":   [0x1f/255, 0x77/255, 0xb4/255],  # #1f77b4
    "RRM2-3": [0x2c/255, 0xa0/255, 0x2c/255],  # #2ca02c
    "KH_N":   [0xff/255, 0x7f/255, 0x0e/255],  # #ff7f0e
    "KH_Ca":  [0xd6/255, 0x27/255, 0x28/255],  # #d62728
    "MD1":    [0xae/255, 0xc7/255, 0xe8/255],  # #aec7e8
    "MD2":    [0x94/255, 0x67/255, 0xbd/255],  # #9467bd
    "MD3":    [0x8c/255, 0x56/255, 0x4b/255],  # #8c564b
    "KH_b":   [0xe3/255, 0x77/255, 0xc2/255],  # #e377c2
    "WWE":    [0x7f/255, 0x7f/255, 0x7f/255],  # #7f7f7f
    "ART":    [0xbc/255, 0xbd/255, 0x22/255],  # #bcbd22
}

# Set custom colors in PyMOL
for name, rgb in domain_colors_rgb.items():
    color_name = f"dom_color_{name.replace('-', '_')}"
    cmd.set_color(color_name, rgb)

domain_colors = {name: f"dom_color_{name.replace('-', '_')}" for name in domain_colors_rgb}

# Active site definitions (from active_sites.yaml)
active_sites = {
    "MD1": {
        "residues": [831, 923, 962],
        "labels": ["D831", "N923", "D962"],
    },
    "MD2": {
        "residues": [1035, 1046, 1134, 1171],
        "labels": ["1035", "1046", "1134", "1171"],
    },
    "MD3": {
        "residues": [1248, 1259, 1330, 1371],
        "labels": ["1248", "1259", "1330", "1371"],
    },
    "ART": {
        "residues": [1684, 1705, 1706, 1722],
        "labels": ["H1684", "Y1705", "E1706", "1722"],
    },
}

# Style AF2 structure
cmd.hide("everything", "PARP14_AF2")
cmd.show("cartoon", "PARP14_AF2")
cmd.color("gray70", "PARP14_AF2")

# Color domains
for name, (start, end) in domains.items():
    sel_name = f"dom_{name.replace('-', '_')}"
    cmd.select(sel_name, f"PARP14_AF2 and resi {start}-{end}")
    cmd.color(domain_colors[name], sel_name)

# Active sites: spheres + labels
all_active = []
for domain, info in active_sites.items():
    for resi, label in zip(info["residues"], info["labels"]):
        sel = f"active_{domain}_{resi}"
        cmd.select(sel, f"PARP14_AF2 and resi {resi} and name CA")
        cmd.show("spheres", sel)
        cmd.set("sphere_scale", 1.2, sel)
        cmd.color("white", sel)
        cmd.label(sel, f'"{label}"')
        all_active.append(f"resi {resi}")

cmd.select("active_sites", f"PARP14_AF2 and ({' or '.join(all_active)}) and name CA")

# =======================================================
# 2. Load AF3 constructs (25 models each as states)
# =======================================================
af3_base = "/home/sbali/CALVADOS/parp14/alphafold_outputs"

af3_constructs = {
    "AF3_largest": "rrm1_rrm2_rrm3_kh1-kh6_kh7a_md2_md3_khb-kh8_wwe_art",
    "AF3_MD1L1_to_end": "md1l1_md2_md3_khb-kh8_wwe_art",
}

for obj_name, construct in af3_constructs.items():
    construct_dir = os.path.join(af3_base, construct)
    # Collect all seed/sample model files, sorted for reproducible ordering
    model_files = sorted(glob.glob(os.path.join(construct_dir, "seed-*_sample-*/model.cif")))
    print(f"Loading {len(model_files)} models for {obj_name} ({construct})")

    for i, cif_path in enumerate(model_files):
        tmp_name = f"_tmp_{i}"
        cmd.load(cif_path, tmp_name)
        if i == 0:
            cmd.create(obj_name, tmp_name, source_state=1, target_state=1)
        else:
            cmd.create(obj_name, tmp_name, source_state=1, target_state=i + 1)
        cmd.delete(tmp_name)

    # Style: cartoon, colored by domain where applicable
    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)
    cmd.color("gray70", obj_name)

    # Color domains that exist in this construct
    for dom_name, (start, end) in domains.items():
        sel = f"{obj_name} and resi {start}-{end}"
        if cmd.count_atoms(sel) > 0:
            cmd.color(domain_colors[dom_name], sel)

    # Mark active sites on AF3 constructs too
    for domain, info in active_sites.items():
        for resi, label in zip(info["residues"], info["labels"]):
            sel = f"{obj_name} and resi {resi} and name CA"
            if cmd.count_atoms(sel) > 0:
                cmd.show("spheres", sel)
                cmd.set("sphere_scale", 1.2, sel)
                cmd.color("white", sel)

    # Disable by default (user can toggle on)
    cmd.disable(obj_name)

# Enable AF2 reference by default
cmd.enable("PARP14_AF2")

# =======================================================
# 3. General display settings
# =======================================================
cmd.set("label_color", "black")
cmd.set("label_size", 16)
cmd.set("label_font_id", 7)
cmd.set("label_position", [0, 0, 3])
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)
cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")
cmd.set("antialias", 2)
cmd.set("ray_trace_mode", 1)

# Orient to AF2 reference
cmd.orient("PARP14_AF2")
cmd.zoom("PARP14_AF2", buffer=5)
cmd.deselect()

# =======================================================
# Save
# =======================================================
output_path = "/home/sbali/CALVADOS/parp14/PARP14_domains_active_sites.pse"
cmd.save(output_path)
print(f"\nSaved PSE to: {output_path}")
print("\nObjects in session:")
print("  PARP14_AF2         - Full-length AlphaFold2 (enabled, colored by domain)")
print("  AF3_largest        - Largest AF3 construct, 25 states (disabled)")
print("  AF3_MD1L1_to_end   - MD1L1-to-ART AF3 construct, 25 states (disabled)")
print("\nUse 'enable <object>' in PyMOL to toggle visibility.")
print("Use state slider or 'set state, N' to browse AF3 models.")

cmd.quit()
