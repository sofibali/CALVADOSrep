"""Place PARP9 and DTX3L separated in a 40 nm box, matching the CALVADOS
separated arm (chains ~17 nm apart, no predicted pose used).

Atom order MUST be parp9 then dtx3l, to match the PSF built in that order.
"""
import numpy as np

def read(fn):
    out = []
    for l in open(fn):
        if l.startswith('ATOM'):
            out.append([l, np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])])
    return out

def write(fn, chains, box_nm):
    L = box_nm * 10.0
    with open(fn, 'w') as g:
        g.write(f'CRYST1{L:9.3f}{L:9.3f}{L:9.3f}  90.00  90.00  90.00 P 1           1\n')
        n = 0
        for atoms in chains:
            for l, xyz in atoms:
                n += 1
                g.write(l[:6] + f'{n:5d}' + l[11:30] +
                        f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}' + l[54:])
        g.write('END\n')
    return n

BOX = 40.0
SEP = 17.0            # nm between chain centroids, as in the CALVADOS grid start
p9, dx = read('parp9_hyres.pdb'), read('dtx3l_hyres.pdb')
for atoms, shift in ((p9, np.array([-SEP/2, 0, 0])), (dx, np.array([+SEP/2, 0, 0]))):
    com = np.mean([a[1] for a in atoms], axis=0)
    for a in atoms:
        a[1] = a[1] - com + (np.array([BOX/2]*3) + shift) * 10.0
c1 = np.mean([a[1] for a in p9], axis=0); c2 = np.mean([a[1] for a in dx], axis=0)
n = write('complex_start.pdb', [p9, dx], BOX)
print(f'complex_start.pdb: {n} particles ({len(p9)} parp9 + {len(dx)} dtx3l)')
print(f'  box {BOX} nm, centroid separation {np.linalg.norm(c1-c2)/10:.1f} nm')
mn = np.min([a[1] for a in p9+dx], axis=0)/10; mx = np.max([a[1] for a in p9+dx], axis=0)/10
print(f'  extent {np.round(mn,1)} to {np.round(mx,1)} nm -- must sit inside 0..{BOX}')
