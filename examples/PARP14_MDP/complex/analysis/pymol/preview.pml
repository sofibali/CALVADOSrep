load p9_dtx3l_xl_h20_best.pdb, best
hide everything
show cartoon, best
color grey80, best
color 0x6A3FA0, best and chain A and resi 107-296
color 0x9169C7, best and chain A and resi 311-493
color 0x1CC497, best and chain A and resi 543-622
color 0x9C5C45, best and chain A and resi 644-823
color 0xB88C13, best and chain B and resi 11-87
color 0x117C60, best and chain B and resi 138-507
color 0xC8481A, best and chain B and resi 556-607
color 0x9C4090, best and chain B and resi 608-737
set cartoon_transparency, 0.1
bg_color white
orient
ray 1000, 800
png preview_best.png, dpi=120
