#!/usr/bin/env python3
"""PyMOL sessions of PARP9-DTX3L states, keyed to which published crosslinks they satisfy.

WHAT IT PICKS
-------------
The scoring showed CALVADOS visits every published crosslink geometry but never
weights any of them. That raises the obvious structural question: is there a
single conformation that satisfies MANY crosslinks at once, or does the ensemble
only ever satisfy them one at a time in mutually incompatible poses?

So the frames chosen are not arbitrary:
  best        the single frame satisfying the most crosslinks simultaneously
  patterns    representative frames for the distinct SATISFACTION PATTERNS --
              frames are described by the binary vector of which crosslinks they
              satisfy, those vectors are clustered, and one medoid per cluster is
              written. Different clusters = genuinely different binding modes.

A satisfied crosslink is CA-CA <= 3.0 nm (BS3 convention, the project's
LYS_LYS_CUTOFF). Crosslinks are drawn as dashes: green where satisfied in that
state, grey where not.

Usage:
    python make_crosslink_pse.py                                   # xl_h20 (restrained)
    python make_crosslink_pse.py --root binding --set p9_dtx3l     # unrestrained control
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

CUT_NM, EQ = 3.0, 500
DOMS = {'parp9': ['1/2KH1a','MD1','MD2','1/2KH1b','KH2','ART'],
        'dtx3l': ['RRM','KH1','KH2','KH3','KH4','KH5','RING','DTC']}
# role-based, from the project style guide extension
COLOR = {'parp9': {'1/2KH1a':'0x16A07B','MD1':'0x6A3FA0','MD2':'0x9169C7',
                   '1/2KH1b':'0x16A07B','KH2':'0x1CC497','ART':'0x9C5C45'},
         'dtx3l': {'RRM':'0xB88C13','KH1':'0x117C60','KH2':'0x16A07B','KH3':'0x1CC497',
                   'KH4':'0x16A07B','KH5':'0x117C60','RING':'0xC8481A','DTC':'0x9C4090'}}


def scan(root, set_name, stride):
    base = set_name if set_name in SETS else 'p9_dtx3l'
    chains = SETS[base]; off, o = {}, 0
    for c in chains: off[c] = o; o += ac.LENGTHS[c]
    d = pd.read_csv(XL)
    d = d[d.block == 'P9-Dtx3L'].drop_duplicates(subset=['res_a','res_b']).reset_index(drop=True)
    pairs = list(zip(d.res_a, d.res_b))
    ia = np.array([off['parp9']+a-1 for a,_ in pairs])
    ib = np.array([off['dtx3l']+b-1 for _,b in pairs])
    sat, src = [], []
    for rep in sorted((HERE/root/set_name).glob('rep-*'),
                      key=lambda p: int(p.name.split('-')[1])):
        dcd = rep/f'{sysname(base)}.dcd'; top = topology(rep)
        if top is None or not dcd.exists(): continue
        cfg = yaml.safe_load(open(rep/'config.yaml')); b = np.array(cfg['box'], float)
        boxv = np.array([b[0]*10,b[1]*10,b[2]*10,90.,90.,90.], np.float32)
        u = mda.Universe(str(top), str(dcd))
        if u.atoms.n_atoms != o: continue
        for k, ts in enumerate(u.trajectory[EQ::stride]):
            P = u.atoms.positions.astype(np.float32)
            dist = np.array([distance_array(P[i:i+1], P[j:j+1], box=boxv)[0,0]
                             for i, j in zip(ia, ib)])/10.0
            sat.append(dist <= CUT_NM); src.append((rep, EQ + k*stride))
    return np.array(sat), src, d, off


def medoids(S, k):
    """k medoids of the binary satisfaction vectors, by Hamming distance."""
    keep = np.where(S.sum(1) > 0)[0]
    if len(keep) < k: return list(keep)
    X = S[keep].astype(float)
    rng = np.random.default_rng(0)
    cen = X[rng.choice(len(X), k, replace=False)]
    for _ in range(50):
        lab = np.argmin(((X[:,None,:] != (cen[None]>0.5)).sum(2)), axis=1)
        new = np.array([X[lab==i].mean(0) if (lab==i).any() else cen[i] for i in range(k)])
        if np.allclose(new, cen): break
        cen = new
    out = []
    for i in range(k):
        m = np.where(lab == i)[0]
        if not len(m): continue
        dd = (X[m] != (cen[i]>0.5)).sum(1)
        out.append(keep[m[np.argmin(dd)]])
    return out


def write_frame(rep, frame, out_pdb, base):
    u = mda.Universe(str(topology(rep)), str(rep/f'{sysname(base)}.dcd'))
    u.trajectory[frame]
    n9 = ac.LENGTHS['parp9']
    u.add_TopologyAttr('chainID')
    u.atoms[:n9].chainIDs = 'A'; u.atoms[n9:].chainIDs = 'B'
    # centre and un-split across the periodic boundary before writing
    u.atoms.positions = center_pair(u.atoms.positions, n9, 40.0)[0].astype('float32')
    tmp = out_pdb.with_suffix('.tmp.pdb')
    u.atoms.write(str(tmp))
    with open(out_pdb, 'w') as g:
        g.write('CRYST1  400.000  400.000  400.000  90.00  90.00  90.00 P 1           1\n')
        for line in open(tmp):
            if line.startswith(('ATOM', 'TER', 'END')):
                g.write(line)
    tmp.unlink()


def main():
    ap = ArgumentParser()
    ap.add_argument('--root', default='binding_go')
    ap.add_argument('--set', default='p9_dtx3l_xl_h20')
    ap.add_argument('--stride', type=int, default=40)
    ap.add_argument('--k', type=int, default=4)
    a = ap.parse_args()
    base = a.set if a.set in SETS else 'p9_dtx3l'

    S, src, xl, off = scan(a.root, a.set, a.stride)
    n = S.sum(1)
    print(f'{len(S)} frames scanned, {len(xl)} unique crosslinks')
    print(f'  crosslinks satisfied per frame: max {n.max()}, mean {n.mean():.2f}, '
          f'frames with >=1: {100*(n>0).mean():.1f}%')

    picks = [('best', int(np.argmax(n)))]
    for i, m in enumerate(medoids(S, a.k), 1):
        if m != picks[0][1]: picks.append((f'mode{i}', int(m)))

    rows = []
    for name, idx in picks:
        rep, frame = src[idx]
        pdb = OUT / f'{a.set}_{name}.pdb'
        write_frame(rep, frame, pdb, base)
        which = [f'K{r.res_a}-K{r.res_b}' for (_, r), s in zip(xl.iterrows(), S[idx]) if s]
        rows.append(dict(state=name, pdb=pdb.name, rep=rep.name, frame=frame,
                         n_satisfied=int(S[idx].sum()), crosslinks=';'.join(which)))
        print(f'  {name}: {rep.name} frame {frame}, {S[idx].sum()} satisfied')

    pd.DataFrame(rows).to_csv(OUT / f'{a.set}_states.csv', index=False)

    # ---- the .pml ----
    dom = yaml.safe_load(open(HERE/'binding'/base/'input'/'domains.yaml'))
    L = [f'# PyMOL session: PARP9-DTX3L states by crosslink satisfaction',
         f'# set {a.set}; a crosslink is satisfied at CA-CA <= {CUT_NM} nm',
         'from pymol import cmd', 'cmd.bg_color("white")', 'cmd.set("dash_gap",0.4)',
         'cmd.set("dash_width",2.5)', 'cmd.set("label_size",14)',
         'cmd.set("ray_opaque_background",0)',
         'cmd.set("cartoon_trace_atoms",1)', 'cmd.set("cartoon_tube_radius",1.1)',
         'cmd.set("ribbon_trace_atoms",1)']
    for r in rows:
        ob = r['state']
        L += [f'cmd.load("{r["pdb"]}","{ob}")', f'cmd.hide("everything","{ob}")',
              f'cmd.show("cartoon","{ob}")', f'cmd.cartoon("tube","{ob}")',
              f'cmd.color("grey70","{ob}")']
        for ch, cid in (('parp9','A'), ('dtx3l','B')):
            for nm, (lo, hi) in zip(DOMS[ch], dom[ch]):
                L.append(f'cmd.color("{COLOR[ch][nm]}","{ob} and chain {cid} '
                         f'and resi {lo}-{hi}")')
        sat = set(r['crosslinks'].split(';')) if r['crosslinks'] else set()
        # Draw ONLY the crosslinks this state satisfies. Creating all 36 and
        # disabling the rest does not survive: cmd.group() re-enables its
        # members and the disable is silently undone, so every state ends up
        # showing every crosslink.
        for _, x in xl.iterrows():
            tag = f'K{x.res_a}-K{x.res_b}'
            if tag not in sat:
                continue
            nmobj = f'{ob}_xl_{x.res_a}_{x.res_b}'
            L.append(f'cmd.distance("{nmobj}","{ob} and chain A and resi {x.res_a} '
                     f'and name CA","{ob} and chain B and resi {x.res_b} and name CA")')
            L.append(f'cmd.color("limon","{nmobj}")')
            L.append(f'cmd.set("dash_radius",0.35,"{nmobj}")')
            L.append(f'cmd.hide("labels","{nmobj}")')
            L.append(f'cmd.show("spheres","{ob} and chain A and resi {x.res_a} and name CA")')
            L.append(f'cmd.show("spheres","{ob} and chain B and resi {x.res_b} and name CA")')
        L.append(f'cmd.group("{ob}_xlinks","{ob}_xl_*")')

    L += ['cmd.set("sphere_scale",0.5)', 'cmd.orient("best")',
          ] + [f'cmd.disable("{r["state"]}")' for r in rows if r['state'] != 'best'] + [
          f'cmd.save("{a.set}_crosslink_states.pse")']
    (OUT / f'{a.set}_crosslink_states.pml').write_text('\n'.join(L) + '\n')
    print(f'\nwrote {OUT}/{a.set}_crosslink_states.pml')


if __name__ == '__main__':
    main()
