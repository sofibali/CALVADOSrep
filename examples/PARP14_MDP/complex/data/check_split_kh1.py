"""Is PARP9's split KH1 still one fold? 1/2KH1a (60-96) vs 1/2KH1b (516-539)."""
import sys, warnings
from pathlib import Path
import numpy as np, yaml, MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
warnings.filterwarnings('ignore')
HERE = Path('/home/sbali/CALVADOS/examples/PARP14_MDP/complex')
sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology
import analyze_convergence as ac

A = np.arange(60, 97) - 1          # 1/2KH1a
B = np.arange(516, 540) - 1        # 1/2KH1b

def stats(P):
    return float(np.linalg.norm(P[A].mean(0) - P[B].mean(0))), \
           float(distance_array(P[A], P[B]).min())

def calvados(s, root='binding', stride=40):
    out = []
    for rep in sorted((HERE/root/s).glob('rep-*'))[:5]:
        top = topology(rep); dcd = rep/f'{sysname(s)}.dcd'
        if top is None or not dcd.exists(): continue
        u = mda.Universe(str(top), str(dcd))
        for ts in u.trajectory[::stride]:
            out.append(stats(u.atoms.positions/10.0))
    return np.array(out)

def hyres(dirs, stride=20):
    out = []
    for d in dirs:
        run = HERE/d
        top = run/'complex_start.pdb'
        if not top.exists(): top = run/'start.pdb'
        u = mda.Universe(str(top), str(run/'system.dcd'))
        ca = u.select_atoms('name CA')
        for ts in u.trajectory[::stride]:
            out.append(stats(ca.positions/10.0))
    return np.array(out)

print(f"{'model':28s} {'COM-COM nm':>22s} {'min CA-CA nm':>20s}")
rows = [('CALVADOS unrestrained', calvados('p9_dtx3l')),
        ('CALVADOS +3 XL h20', calvados('p9_dtx3l_xl_h20', 'binding_go')),
        ('HyRes matched separated', hyres([f'hyres/runs/matched/separated_rep{i}' for i in (1,2,3)])),
        ('HyRes matched docked', hyres([f'hyres/runs/matched/docked_rep{i}' for i in (1,2,3)])),
        ('HyRes unrestrained', hyres(['hyres/runs/prod_rep1']))]
for name, a in rows:
    if not len(a): print(f'{name:28s}  no data'); continue
    print(f'{name:28s}  {a[:,0].mean():6.2f} (first {a[0,0]:.2f}, max {a[:,0].max():5.2f})'
          f'   {a[:,1].mean():6.2f} (max {a[:,1].max():5.2f})')
