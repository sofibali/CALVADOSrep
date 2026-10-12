#!/usr/bin/env python3
"""Trim whitespace and downscale the renders into the names build_report.py embeds.

    python analysis/crop_figures.py       # run from bindcraft_md/
"""
import os
from PIL import Image, ImageChops
HERE = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(HERE, "assets")

def crop(src, dst, maxw, box=None):
    im = Image.open(os.path.join(A, src)).convert("RGB")
    bb = box or ImageChops.difference(im, Image.new("RGB", im.size, (255,255,255))).getbbox()
    if bb:
        p = 14
        bb = (max(0,bb[0]-p), max(0,bb[1]-p), min(im.size[0],bb[2]+p), min(im.size[1],bb[3]+p))
        im = im.crop(bb)
    im.thumbnail((maxw, maxw), Image.LANCZOS)
    im.save(os.path.join(A, dst), optimize=True)
    return bb, os.path.getsize(os.path.join(A, dst)) // 1024

for n in ("md1_block_af3","md1_block_sim","clamp_md1md2_state3"):
    _, kb = crop(f"report_{n}.png", f"report_{n}_crop.png", 950); print(f"  report_{n}_crop.png {kb}KB")
for i in range(5):
    _, kb = crop(f"hit_{i}.png", f"crop_hit_{i}.png", 620); print(f"  crop_hit_{i}.png {kb}KB")
_, kb = crop("hits_clamp.png", "crop_clamp.png", 900);    print(f"  crop_clamp.png {kb}KB")
_, kb = crop("hits_variants.png", "crop_variants.png", 760); print(f"  crop_variants.png {kb}KB")
for n in ("md1_block_af3","md2_block_af3","md3_block_af3"):
    _, kb = crop(f"pocket_{n}.png", f"crop_pocket_{n}.png", 560); print(f"  crop_pocket_{n}.png {kb}KB")
# BC1/BC2 share ONE crop box so the size difference stays honest
box, kb = crop("bc_target_bc1.png", "crop_bc_target_bc1.png", 560)
print(f"  crop_bc_target_bc1.png {kb}KB")
_, kb = crop("bc_target_bc2.png", "crop_bc_target_bc2.png", 560, box=box)
print(f"  crop_bc_target_bc2.png {kb}KB")
