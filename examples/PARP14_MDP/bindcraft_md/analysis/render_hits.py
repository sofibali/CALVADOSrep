import os, glob, re
from pymol import cmd
OUT="/tmp/claude-64170/-home-sbali-BindCraft/adecb958-709c-481a-8e03-db4ccc10a3b9/scratchpad"
B="/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md/designs"
PAL=["0x2a78d6","0xeb6834","0x1baf7a","0xeda100","0xe87ba4"]
GREY="0xc2c8d0"; GREY2="0xdbd5c8"; HS_C="0x4a3aa7"
HS=[33,37,39,41,42,43,44,46,47,134,135,136,138,172,173]

def base():
    cmd.bg_color("white"); cmd.set("ray_opaque_background",1)
    cmd.hide("everything"); cmd.show("cartoon")

files=sorted(glob.glob(f"{B}/sweep_af3_guess/Accepted/Ranked/*.pdb"))
seen={}
for f in files:
    seen.setdefault(re.search(r"(l\d+_s\d+)",f).group(1), f)
trajs=list(seen.items())

# --- fixed camera from the target, oriented on the hotspot patch ---
cmd.reinitialize(); cmd.load(trajs[0][1],"t0"); base()
cmd.remove("t0 and chain B")
cmd.select("hs","t0 and ("+" or ".join(f"resi {r}" for r in HS)+")")
cmd.orient("hs"); cmd.zoom("t0",2)
cmd.turn("y",-25); cmd.turn("x",10)
VIEW=cmd.get_view()

def panel(idx, traj, path, out):
    cmd.reinitialize(); base()
    cmd.load(path,"cx")
    cmd.color(GREY,"cx and chain A")
    cmd.select("hs","cx and chain A and ("+" or ".join(f"resi {r}" for r in HS)+")")
    cmd.show("sticks","hs and not name C+N+O"); cmd.color(HS_C,"hs")
    cmd.color(PAL[idx%len(PAL)],"cx and chain B")
    cmd.set_view(VIEW)
    cmd.ray(760,640); cmd.png(out,dpi=180)

for i,(traj,path) in enumerate(trajs):
    panel(i,traj,path,f"{OUT}/hit_{i}_{traj}.png")
    print("panel",i,traj)

# --- superposition, binders semi-transparent so the target stays readable ---
cmd.reinitialize(); base()
cmd.load(trajs[0][1],"scaf"); cmd.remove("scaf and chain B")
cmd.color(GREY,"scaf")
cmd.select("hs","scaf and ("+" or ".join(f"resi {r}" for r in HS)+")")
cmd.show("sticks","hs and not name C+N+O"); cmd.color(HS_C,"hs")
for i,(traj,path) in enumerate(trajs):
    o=f"b{i}"; cmd.load(path,o); cmd.remove(f"{o} and chain A")
    cmd.color(PAL[i%len(PAL)],o); cmd.set("cartoon_transparency",0.55,o)
cmd.set_view(VIEW); cmd.ray(1200,1000); cmd.png(f"{OUT}/hits_superpose.png",dpi=200)
print("wrote superposition")

# --- clamp bridge ---
cmd.reinitialize(); base()
cf=sorted(glob.glob(f"{B}/guess_clamp_md1md2_state3/Accepted/*.pdb"))
cmd.load(cf[0],"cx")
cmd.color(GREY,"cx and chain A and resi 1-215")
cmd.color(GREY2,"cx and chain A and resi 216-404")
cmd.color("0xeb6834","cx and chain B")
for g,res,col in (("MD1",[71,102,103,106,107,109,110],"0x2a78d6"),
                  ("MD2",[351,354,357,395,396,399,400],"0x1baf7a")):
    cmd.select(g,"cx and chain A and ("+" or ".join(f"resi {r}" for r in res)+")")
    cmd.show("sticks",f"{g} and not name C+N+O"); cmd.color(col,g)
cmd.orient("cx"); cmd.zoom("cx",3)
cmd.ray(1400,1050); cmd.png(f"{OUT}/hits_clamp.png",dpi=200)
print("wrote clamp")

# --- same-backbone MPNN variants ---
cmd.reinitialize(); base()
pair=[f for f in files if "l128_s514541" in f]
cmd.load(pair[0],"v1"); cmd.load(pair[1],"v2")
cmd.color(GREY,"v1 and chain A"); cmd.remove("v2 and chain A")
cmd.color("0x2a78d6","v1 and chain B"); cmd.color("0xeb6834","v2")
cmd.orient("v1 and chain B"); cmd.zoom("v1 and chain B",3)
cmd.ray(1100,850); cmd.png(f"{OUT}/hits_variants.png",dpi=200)
print("wrote variants")
