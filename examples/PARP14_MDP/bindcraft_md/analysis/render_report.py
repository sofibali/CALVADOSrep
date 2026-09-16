import os
from pymol import cmd

OUT = "/tmp/claude-64170/-home-sbali-BindCraft/adecb958-709c-481a-8e03-db4ccc10a3b9/scratchpad"
TDIR = "/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md/targets"

# MD1 = 1-215, MD2 = 216-404
C_MD1 = "0x9aa4b2"   # cool grey  (MD1 body)
C_MD2 = "0xc9c2b6"   # warm grey  (MD2 body)
C_HS1 = "0x2a78d6"   # blue  = MD1 hotspots
C_HS2 = "0xeb6834"   # orange = MD2 hotspots

JOBS = [
    ("md1_block_sim", {"MD1": [33,37,39,41,42,43,44,46,47,134,135,136,138,172,173], "MD2": []}),
    ("md1_block_af3", {"MD1": [33,37,39,41,42,43,44,46,47,134,135,136,138,172,173], "MD2": []}),
    ("clamp_md1md2_state3", {"MD1": [71,102,103,106,107,109,110], "MD2": [351,354,357,395,396,399,400]}),
]

for name, hs in JOBS:
    cmd.reinitialize()
    cmd.load(os.path.join(TDIR, name + ".pdb"), name)
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_transparency", 0.25)

    cmd.color(C_MD1, "resi 1-215")
    cmd.color(C_MD2, "resi 216-404")

    for grp, col in (("MD1", C_HS1), ("MD2", C_HS2)):
        resis = hs[grp]
        if not resis:
            continue
        sel = grp + "_hs"
        cmd.select(sel, " or ".join(f"resi {r}" for r in resis))
        cmd.show("sticks", f"{sel} and not name C+N+O")
        cmd.show("spheres", f"{sel} and name CA")
        cmd.set("sphere_scale", 0.45, sel)
        cmd.color(col, sel)
        cmd.set("cartoon_transparency", 0.0, sel)

    cmd.orient(name)
    cmd.zoom(name, buffer=3)
    cmd.ray(1400, 1000)
    out = os.path.join(OUT, f"report_{name}.png")
    cmd.png(out, dpi=200)
    print("wrote", out)
