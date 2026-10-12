"""All PyMOL renders for the campaign artifact, written into analysis/assets/.

Run with:  pymol -cq analysis/make_figures.py      (from bindcraft_md/)

Repo-relative by design: the originals lived in a /tmp scratchpad that was cleared.
"""
import os, glob, re
from pymol import cmd

ROOT = "/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md"
OUT  = os.path.join(ROOT, "analysis", "assets")
os.makedirs(OUT, exist_ok=True)

PAL   = ["0x2a78d6","0xeb6834","0x1baf7a","0xeda100","0xe87ba4"]
GREY, GREY2, HS_C = "0xc2c8d0", "0xd9cfb8", "0x4a3aa7"
MD1HS = [33,35,37,39,41,42,43,44,46,47,134,135,136,138,172,173]
CLAMPHS = {"MD1":[71,102,103,106,107,109,110], "MD2":[351,354,357,395,396,399,400]}
CLS = {'ARG':'pos','LYS':'pos','HIS':'polar','ASP':'neg','GLU':'neg','ASN':'polar','GLN':'polar',
       'SER':'polar','THR':'polar','TYR':'arom','TRP':'arom','CYS':'polar','PHE':'arom',
       'ALA':'apol','VAL':'apol','LEU':'apol','ILE':'apol','MET':'apol','PRO':'apol','GLY':'gly'}
CCOL = {'apol':"0xd9cfb8",'arom':"0x9085e9",'polar':"0x1baf7a",'pos':"0x2a78d6",'neg':"0xe34948",'gly':"0xc2c8d0"}
sel_of = lambda rs: " or ".join(f"resi {r}" for r in rs)

def base():
    cmd.bg_color("white"); cmd.set("ray_opaque_background", 1)
    cmd.hide("everything"); cmd.show("cartoon")

# ---- 1. target overviews: MD1 af3 / MD1 sim / clamp state3 -------------------
for name, hs in [("md1_block_af3",{"MD1":MD1HS,"MD2":[]}),
                 ("md1_block_sim",{"MD1":MD1HS,"MD2":[]}),
                 ("clamp_md1md2_state3",CLAMPHS)]:
    cmd.reinitialize(); cmd.load(f"{ROOT}/targets/{name}.pdb", name); base()
    cmd.set("cartoon_transparency", 0.25)
    cmd.color("0x9aa4b2","resi 1-215"); cmd.color("0xc9c2b6","resi 216-404")
    for g, col in (("MD1","0x2a78d6"), ("MD2","0xeb6834")):
        if not hs[g]: continue
        cmd.select(g+"_hs", sel_of(hs[g]))
        cmd.show("sticks", f"{g}_hs and not name C+N+O")
        cmd.show("spheres", f"{g}_hs and name CA"); cmd.set("sphere_scale", 0.45, g+"_hs")
        cmd.color(col, g+"_hs"); cmd.set("cartoon_transparency", 0.0, g+"_hs")
    cmd.orient(name); cmd.zoom(name, 3)
    cmd.ray(1400,1000); cmd.png(f"{OUT}/report_{name}.png", dpi=200)
    print("target", name)

# ---- 2. MD1 binder hits: one panel per source backbone, shared camera --------
files = sorted(glob.glob(f"{ROOT}/designs/sweep_af3_guess/Accepted/Ranked/*.pdb"))
seen = {}
for f in files:
    seen.setdefault(re.search(r"(l\d+_s\d+)", f).group(1), f)
trajs = list(seen.items())
cmd.reinitialize(); cmd.load(trajs[0][1], "t0"); base()
cmd.remove("t0 and chain B"); cmd.select("hs", "t0 and (" + sel_of(MD1HS) + ")")
cmd.orient("hs"); cmd.zoom("t0", 2); cmd.turn("y",-25); cmd.turn("x",10)
VIEW = cmd.get_view()
for i,(traj,path) in enumerate(trajs):
    cmd.reinitialize(); base(); cmd.load(path, "cx")
    cmd.color(GREY, "cx and chain A")
    cmd.select("hs", "cx and chain A and (" + sel_of(MD1HS) + ")")
    cmd.show("sticks", "hs and not name C+N+O"); cmd.color(HS_C, "hs")
    cmd.color(PAL[i % len(PAL)], "cx and chain B")
    cmd.set_view(VIEW); cmd.ray(760,640); cmd.png(f"{OUT}/hit_{i}.png", dpi=180)
    print("hit", i, traj)

# ---- 3. clamp bridging both domains ----------------------------------------
cf = sorted(glob.glob(f"{ROOT}/designs/guess_clamp_md1md2_state3/Accepted/*.pdb"))[0]
cmd.reinitialize(); base(); cmd.load(cf, "cx")
cmd.color("0xa8b4c4","cx and chain A and resi 1-215")
cmd.color("0xd9cfb8","cx and chain A and resi 216-404")
cmd.color("0xeb6834","cx and chain B")
cmd.select("c1","byres (cx and chain A and resi 1-215 within 4 of (cx and chain B))")
cmd.select("c2","byres (cx and chain A and resi 216-404 within 4 of (cx and chain B))")
cmd.show("sticks","(c1 or c2) and not name C+N+O")
cmd.color("0x2a78d6","c1"); cmd.color("0x1baf7a","c2")
cmd.orient("cx and chain B"); cmd.turn("x",90); cmd.zoom("cx",2)
cmd.ray(1400,1000); cmd.png(f"{OUT}/hits_clamp.png", dpi=200); print("clamp")

# ---- 4. two MPNN variants on one backbone ----------------------------------
pair = [f for f in files if "l128_s514541" in f][:2]
cmd.reinitialize(); base(); cmd.load(pair[0],"v1"); cmd.load(pair[1],"v2")
cmd.color(GREY,"v1 and chain A"); cmd.remove("v2 and chain A")
cmd.color("0x2a78d6","v1 and chain B"); cmd.color("0xeb6834","v2")
cmd.orient("v1 and chain B"); cmd.zoom("v1 and chain B",3)
cmd.ray(1100,850); cmd.png(f"{OUT}/hits_variants.png", dpi=200); print("variants")

# ---- 5. pocket surfaces coloured by residue class ---------------------------
import json
for name, cfg in [("md1_block_af3", f"{ROOT}/settings/md1_block_af3.json"),
                  ("md2_block_af3", f"{ROOT}/sweep/queue/md2_block_af3.json"),
                  ("md3_block_af3", f"{ROOT}/sweep/queue/md3_block_af3.json")]:
    s = json.load(open(cfg)); hs = [int(x) for x in s["target_hotspot_residues"].split(",")]
    cmd.reinitialize(); cmd.bg_color("white"); cmd.set("ray_opaque_background",1)
    cmd.load(s["starting_pdb"], "t"); cmd.hide("everything")
    cmd.show("surface","t"); cmd.color("0xeeeeee","t"); cmd.set("transparency",0.10)
    cmd.select("hs", "t and (" + sel_of(hs) + ")")
    for r in hs:
        at = cmd.get_model(f"t and resi {r} and name CA").atom
        if at: cmd.color(CCOL.get(CLS.get(at[0].resn,'apol'),CCOL['apol']), f"t and resi {r}")
    cmd.orient("hs"); cmd.zoom("hs",5)
    cmd.ray(820,700); cmd.png(f"{OUT}/pocket_{name}.png", dpi=180)
    print("pocket", name)

# ---- 6. BC1 vs BC2 MD1 target, identical camera -----------------------------
BC1T, BC2T = f"{ROOT}/targets/md1_block_af3.pdb", f"{ROOT}/bc2/targets/PARP14_MD1.pdb"
cmd.reinitialize(); cmd.load(BC1T,"t"); base()
cmd.orient("t"); cmd.zoom("t",3); cmd.turn("y",15); VIEW2 = cmd.get_view()
for tag, path, md2 in (("bc1", BC1T, True), ("bc2", BC2T, False)):
    cmd.reinitialize(); cmd.load(path,"t"); base()
    cmd.color("0xa8b4c4","t and resi 1-215")
    if md2: cmd.color("0xeb6834","t and resi 216-404")
    cmd.select("hs", "t and resi 1-215 and (" + sel_of(MD1HS) + ")")
    cmd.show("sticks","hs and not name C+N+O"); cmd.color(HS_C,"hs")
    cmd.set_view(VIEW2); cmd.ray(900,760); cmd.png(f"{OUT}/bc_target_{tag}.png", dpi=180)
    print("bc target", tag)
print("ALL RENDERS DONE ->", OUT)
