import glob
from pymol import cmd
OUT="/tmp/claude-64170/-home-sbali-BindCraft/adecb958-709c-481a-8e03-db4ccc10a3b9/scratchpad"
B="/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md/designs"
cf=sorted(glob.glob(f"{B}/guess_clamp_md1md2_state3/Accepted/*.pdb"))[0]
cmd.reinitialize(); cmd.bg_color("white"); cmd.set("ray_opaque_background",1)
cmd.load(cf,"cx"); cmd.hide("everything"); cmd.show("cartoon")
cmd.color("0xa8b4c4","cx and chain A and resi 1-215")     # MD1 body, cool grey
cmd.color("0xd9cfb8","cx and chain A and resi 216-404")   # MD2 body, warm grey
cmd.color("0xeb6834","cx and chain B")
# contact residues on each domain (within 4A of binder)
cmd.select("c1","byres (cx and chain A and resi 1-215 within 4 of (cx and chain B))")
cmd.select("c2","byres (cx and chain A and resi 216-404 within 4 of (cx and chain B))")
cmd.show("sticks","(c1 or c2) and not name C+N+O")
cmd.color("0x2a78d6","c1"); cmd.color("0x1baf7a","c2")
# view along the binder so both domains are visible either side
cmd.orient("cx and chain B")
cmd.turn("x",90)
cmd.zoom("cx",2)
cmd.ray(1400,1000); cmd.png(f"{OUT}/hits_clamp.png",dpi=200)
print("wrote clamp view")
print("MD1 contacts:",cmd.count_atoms("c1 and name CA"),"MD2 contacts:",cmd.count_atoms("c2 and name CA"))
