#!/usr/bin/env python3
"""PyMOL time-course sessions: the box, the setup, and what happens after.

WHAT THIS SHOWS THAT THE STATE SESSIONS DO NOT
----------------------------------------------
make_crosslink_pse.py picks frames by how many crosslinks they satisfy, which
answers "what does agreement look like" but throws away time. This one walks the
trajectory instead: t = 0, 25, 50, 75 and 100% of each run, so you can see the
setup, whether the chains found each other, and whether they stayed together.

Each state is labelled with its time, the minimum inter-chain CA-CA distance and
how many published crosslinks it satisfies, and carries dashes for exactly the
crosslinks satisfied at that instant.

THE BOX IS DRAWN. MDAnalysis's PDBWriter emits no CRYST1 record, so PyMOL has no
unit cell to show and the 40 nm periodic box -- the thing that sets the 25.9 uM
concentration and makes encounter frequent -- is invisible. CRYST1 is written
explicitly here and `show cell` turns it on.

A caveat worth seeing rather than reading: chains are written whole, so a chain
that has diffused across a periodic boundary is drawn OUTSIDE the box rather
than wrapped. That is correct for measuring distances (minimum image is applied)
but it means a chain can appear to sit beyond the cell edge.

Usage:
    python make_timecourse_pse.py --root binding --set p9_dtx3l --rep 1
    python make_timecourse_pse.py --hyres hyres/runs/matched/separated_rep2
"""
import sys, warnings
from pathlib import Path
from argparse import ArgumentParser
import numpy as np, pandas as pd, yaml
import MDAnalysis as mda
from pbc_center import center_pair
from MDAnalysis.analysis.distances import distance_array
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis' / 'pymol'
OUT.mkdir(parents=True, exist_ok=True)
XL = HERE / 'data' / 'published_xlinks_v2.csv'
sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology
import analyze_convergence as ac

CUT_NM, BOUND_NM, BOX_NM = 3.0, 1.0, 40.0
FRACS = [0.0, 0.25, 0.50, 0.75, 1.0]
DOMS = {'parp9': ['1/2KH1a','MD1','MD2','1/2KH1b','KH2','ART'],
        'dtx3l': ['RRM','KH1','KH2','KH3','KH4','KH5','RING','DTC']}
COLOR = {'parp9': {'1/2KH1a':'0x16A07B','MD1':'0x6A3FA0','MD2':'0x9169C7',
                   '1/2KH1b':'0x16A07B','KH2':'0x1CC497','ART':'0x9C5C45'},
         'dtx3l': {'RRM':'0xB88C13','KH1':'0x117C60','KH2':'0x16A07B','KH3':'0x1CC497',
                   'KH4':'0x16A07B','KH5':'0x117C60','RING':'0xC8481A','DTC':'0x9C4090'}}


def load_xl():
    d = pd.read_csv(XL)
    return d[d.block == 'P9-Dtx3L'].drop_duplicates(subset=['res_a','res_b']).reset_index(drop=True)


def write_frame(u, sel, frame, out_pdb, n9):
    u.trajectory[frame]
    ag = sel
    if not hasattr(ag, 'chainIDs') or ag.chainIDs is None:
        pass
    try:
        ag.universe.add_TopologyAttr('chainID')
    except Exception:
        pass
    ag[:n9].chainIDs = 'A'; ag[n9:].chainIDs = 'B'
    # centre the pair and bring the partner into the nearest periodic image, so
    # the drawn separation equals the minimum-image separation that was measured
    ag.positions = center_pair(ag.positions, n9, BOX_NM)[0].astype('float32')
    tmp = out_pdb.with_suffix('.tmp.pdb')
    ag.write(str(tmp))
    L = BOX_NM * 10.0
    with open(out_pdb, 'w') as g:
        g.write(f'CRYST1{L:9.3f}{L:9.3f}{L:9.3f}  90.00  90.00  90.00 P 1           1\n')
        for line in open(tmp):
            if line.startswith(('ATOM', 'TER', 'END')):
                g.write(line)
    tmp.unlink()


def main():
    ap = ArgumentParser()
    ap.add_argument('--root', default='binding'); ap.add_argument('--set', default='p9_dtx3l')
    ap.add_argument('--rep', type=int, default=1); ap.add_argument('--hyres', default=None)
    a = ap.parse_args()
    xl = load_xl(); n9 = ac.LENGTHS['parp9']
    boxv = np.array([BOX_NM*10]*3 + [90.,90.,90.], np.float32)

    if a.hyres:
        run = Path(a.hyres); tag = run.name
        top = run/'complex_start.pdb'
        if not top.exists(): top = run/'start.pdb'
        u = mda.Universe(str(top), str(run/'system.dcd'))
        sel = u.select_atoms('name CA'); nsf = 0.05
    else:
        base = a.set if a.set in SETS else 'p9_dtx3l'
        rep = HERE/a.root/a.set/f'rep-{a.rep}'
        tag = f'{a.set}_rep{a.rep}'
        u = mda.Universe(str(topology(rep)), str(rep/f'{sysname(base)}.dcd'))
        sel = u.atoms; nsf = 0.05
    N = len(u.trajectory)
    ia = np.array([r.res_a-1 for _, r in xl.iterrows()])
    ib = np.array([n9+r.res_b-1 for _, r in xl.iterrows()])

    rows = []
    for f in FRACS:
        k = min(int(round(f*(N-1))), N-1)
        u.trajectory[k]
        P = sel.positions.astype(np.float32)
        mind = float(distance_array(P[:n9], P[n9:], box=boxv).min())/10.0
        d = np.array([distance_array(P[i:i+1], P[j:j+1], box=boxv)[0,0] for i, j in zip(ia, ib)])/10.0
        satmask = d <= CUT_NM
        name = f't{int(round(f*100)):03d}'
        pdb = OUT/f'{tag}_{name}.pdb'
        write_frame(u, sel, k, pdb, n9)
        rows.append(dict(state=name, pdb=pdb.name, frame=k, ns=round(k*nsf, 1),
                         min_nm=round(mind, 2), bound=bool(mind < BOUND_NM),
                         n_sat=int(satmask.sum()),
                         sat=';'.join(f'{r.res_a}_{r.res_b}' for (_, r), s
                                      in zip(xl.iterrows(), satmask) if s)))
        print(f'  {name}: {rows[-1]["ns"]:6.1f} ns  minD {mind:6.2f} nm  '
              f'{"BOUND " if rows[-1]["bound"] else "apart "} {satmask.sum():2d} XL')

    pd.DataFrame(rows).to_csv(OUT/f'{tag}_timecourse.csv', index=False)
    dom = yaml.safe_load(open(HERE/'binding'/'p9_dtx3l'/'input'/'domains.yaml'))
    L = ['from pymol import cmd', 'cmd.bg_color("white")',
         'cmd.set("cartoon_trace_atoms",1)', 'cmd.set("cartoon_tube_radius",1.0)',
         'cmd.set("dash_gap",0.4)', 'cmd.set("dash_width",2.5)',
         'cmd.set("sphere_scale",0.5)', 'cmd.set("label_size",18)',
         'cmd.set("cell_color","grey50")']
    for r in rows:
        ob = r['state']
        L += [f'cmd.load("{r["pdb"]}","{ob}")', f'cmd.hide("everything","{ob}")',
              f'cmd.show("cartoon","{ob}")', f'cmd.cartoon("tube","{ob}")',
              f'cmd.color("grey70","{ob}")', f'cmd.show("cell","{ob}")']
        for ch, cid in (('parp9','A'), ('dtx3l','B')):
            for nm, (lo, hi) in zip(DOMS[ch], dom[ch]):
                L.append(f'cmd.color("{COLOR[ch][nm]}","{ob} and chain {cid} and resi {lo}-{hi}")')
        for pair in (r['sat'].split(';') if r['sat'] else []):
            aa, bb = pair.split('_')
            nmobj = f'{ob}_xl_{aa}_{bb}'
            L += [f'cmd.distance("{nmobj}","{ob} and chain A and resi {aa} and name CA",'
                  f'"{ob} and chain B and resi {bb} and name CA")',
                  f'cmd.color("limon","{nmobj}")', f'cmd.set("dash_radius",0.35,"{nmobj}")',
                  f'cmd.hide("labels","{nmobj}")',
                  f'cmd.show("spheres","{ob} and chain A and resi {aa} and name CA")',
                  f'cmd.show("spheres","{ob} and chain B and resi {bb} and name CA")']
        L.append(f'cmd.group("{ob}_xlinks","{ob}_xl_*")')
        lab = f'{r["ns"]:.0f} ns | {r["min_nm"]:.2f} nm | {r["n_sat"]} XL'
        L.append(f'cmd.pseudoatom("{ob}_label",pos=[20,380,20],label="{lab}")')
    L += ['cmd.zoom("all",5)'] + \
         [f'cmd.disable("{r["state"]}")' for r in rows[1:]] + \
         [f'cmd.disable("{r["state"]}_label")' for r in rows[1:]] + \
         [f'cmd.save("{tag}_timecourse.pse")']
    (OUT/f'{tag}_timecourse.pml').write_text('\n'.join(L)+'\n')
    print(f'wrote {OUT}/{tag}_timecourse.pml')


if __name__ == '__main__':
    main()
