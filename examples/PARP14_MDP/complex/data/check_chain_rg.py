"""Whole-chain Rg and max intra-chain domain separation: is HyRes collapsing the chains?"""
import sys, warnings
from pathlib import Path
import numpy as np, MDAnalysis as mda
warnings.filterwarnings('ignore')
HERE = Path('/home/sbali/CALVADOS/examples/PARP14_MDP/complex')
sys.path.insert(0, str(HERE))
from analyze_binding import sysname, topology
import analyze_convergence as ac
L9, LD = ac.LENGTHS['parp9'], ac.LENGTHS['dtx3l']

def rg(X):
    return float(np.sqrt(((X - X.mean(0))**2).sum(1).mean()))

def run(name, frames):
    a = np.array(frames)
    print(f'{name:30s} PARP9 Rg {a[:,0].mean():5.2f} nm   DTX3L Rg {a[:,1].mean():5.2f} nm'
          f'   (start {a[0,0]:.2f} / {a[0,1]:.2f})')

def calvados(s, root='binding', stride=50):
    out=[]
    for rep in sorted((HERE/root/s).glob('rep-*'))[:5]:
        top=topology(rep); dcd=rep/f'{sysname(s if root=="binding" else "p9_dtx3l")}.dcd'
        if top is None or not dcd.exists(): continue
        u=mda.Universe(str(top),str(dcd))
        for ts in u.trajectory[::stride]:
            P=u.atoms.positions/10.0
            out.append((rg(P[:L9]), rg(P[L9:L9+LD])))
    return out

def hyres(dirs, stride=25):
    out=[]
    for d in dirs:
        r=HERE/d; top=r/'complex_start.pdb'
        if not top.exists(): top=r/'start.pdb'
        u=mda.Universe(str(top),str(r/'system.dcd')); ca=u.select_atoms('name CA')
        for ts in u.trajectory[::stride]:
            P=ca.positions/10.0
            out.append((rg(P[:L9]), rg(P[L9:L9+LD])))
    return out

run('CALVADOS unrestrained', calvados('p9_dtx3l'))
run('CALVADOS +3 XL h20', calvados('p9_dtx3l_xl_h20','binding_go'))
run('HyRes matched separated', hyres([f'hyres/runs/matched/separated_rep{i}' for i in (1,2,3)]))
run('HyRes matched docked', hyres([f'hyres/runs/matched/docked_rep{i}' for i in (1,2,3)]))
run('HyRes unrestrained', hyres(['hyres/runs/prod_rep1']))
