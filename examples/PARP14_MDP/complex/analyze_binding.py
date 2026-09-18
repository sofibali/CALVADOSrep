#!/usr/bin/env python3
"""Stage-1 convergence analysis for the binding campaign (castor layout).

Adapter over analyze_convergence.py, which was written for Vincent's older tree
(ROOT=/home/vle/.../p14fl_test/complexes, system/variant/rep nesting, old system
names) and does NOT run unchanged here. Differences handled:

  * layout      binding/{set}/rep-N/            (flat, no variant level)
  * set names   p9_dtx3l / ternary / *_homo / *_docked  (not parp9_dtx3l etc.)
  * dcd name    {sysname}.dcd, where sysname drops the _docked suffix
  * topology    separated sets have no input/system.pdb -- CALVADOS writes
                top.pdb at startup; docked sets have both. Prefer system.pdb.

Everything else (observables, Rhat, tau/ESS, drift, contact maps) is reused from
analyze_convergence so the two campaigns stay comparable.

Usage:  python analyze_binding.py [set ...]      # default: every set with >=3 reps
"""
import sys, itertools, warnings
from pathlib import Path
import numpy as np, yaml
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).parent))
import analyze_convergence as ac

HERE = Path(__file__).parent
ROOT = HERE / 'binding'
OUT  = HERE / 'analysis'
(OUT / 'data').mkdir(parents=True, exist_ok=True)

# set -> ordered chain list (order = components.yaml `system` order = CALVADOS bead order)
SETS = {
    'p9_dtx3l':          ['parp9', 'dtx3l'],
    'p9_dtx3l_docked':   ['parp9', 'dtx3l'],
    'dtx3l_homo':        ['dtx3l', 'dtx3l'],
    'dtx3l_homo_docked': ['dtx3l', 'dtx3l'],
    'p9_homo':           ['parp9', 'parp9'],
    'p14_dtx3l':         ['parp14', 'dtx3l'],
    'p14_dtx3l_docked':  ['parp14', 'dtx3l'],
    'p14_p9':            ['parp14', 'parp9'],
    'p14_p9_docked':     ['parp14', 'parp9'],
    'ternary':           ['parp14', 'parp9', 'dtx3l'],
    'ternary_docked':    ['parp14', 'parp9', 'dtx3l'],
    'p14_homo':          ['parp14', 'parp14'],
}

MIN_REPS = 3


def sysname(s):
    return s[:-7] if s.endswith('_docked') else s


def topology(rep):
    """Docked sets ship input/system.pdb; separated sets only have runtime top.pdb."""
    for cand in (rep / 'input' / 'system.pdb', rep / 'top.pdb'):
        if cand.exists():
            return cand
    return None


def analyse_rep(s, rep, chains):
    """Per-replicate observables. Mirrors ac.analyse_rep with the two path fixes."""
    top, dcd = topology(rep), rep / f'{sysname(s)}.dcd'
    if top is None or not dcd.exists():
        return None
    cfg = yaml.safe_load(open(rep / 'config.yaml'))
    box = np.array(cfg['box'], float)
    u = mda.Universe(str(top), str(dcd))
    if u.atoms.n_atoms != sum(ac.LENGTHS[c] for c in chains):
        return None

    offs, o = [], 0
    for c in chains:
        offs.append((c, o, o + ac.LENGTHS[c])); o += ac.LENGTHS[c]
    pairs = list(itertools.combinations(range(len(chains)), 2))
    lab = {c: ac.labels_for(c) for c in set(chains)}
    rg   = {i: [] for i in range(len(chains))}
    mind = {p: [] for p in pairs}
    ncon = {p: [] for p in pairs}
    cmaps = {p: np.zeros((len(lab[chains[p[0]]][1]), len(lab[chains[p[1]]][1]))) for p in pairs}

    boxv = np.array([box[0]*10, box[1]*10, box[2]*10, 90., 90., 90.], np.float32)
    cut_A = ac.CUTOFF * 10.0
    nfr = 0
    for ts in u.trajectory[::ac.STRIDE]:
        P = u.atoms.positions.astype(np.float32)
        for i, (c, a, b) in enumerate(offs):
            X = P[a:b] / 10.0
            rg[i].append(float(np.sqrt(((X - X.mean(0))**2).sum(1).mean())))
        for p in pairs:
            i, j = p
            A = P[offs[i][1]:offs[i][2]]; B = P[offs[j][1]:offs[j][2]]
            D = distance_array(A, B, box=boxv)
            mind[p].append(float(D.min()) / 10.0)
            ii, jj = np.nonzero(D < cut_A)
            ncon[p].append(int(ii.size))
            if ii.size:
                np.add.at(cmaps[p], (lab[chains[i]][0][ii], lab[chains[j]][0][jj]), 1.0)
        nfr += 1
    for p in pairs:
        cmaps[p] /= max(nfr, 1)
    return dict(rg={i: np.array(v) for i, v in rg.items()},
                mind={p: np.array(v) for p, v in mind.items()},
                ncon={p: np.array(v) for p, v in ncon.items()},
                cmaps=cmaps, nframes=nfr, chains=chains,
                blocks={c: ac.labels_for(c)[1] for c in set(chains)})


def trim(res):
    """Truncate every replicate to the shortest one.

    Replicates are analysed as they finish, and a set can also contain a rep that
    died early. report2 stacks per-rep series into a 2-D array, so unequal lengths
    raise "inhomogeneous shape". Trimming to the common length keeps that stage
    working on partial data and costs nothing once all reps are complete.
    """
    n = min(a['nframes'] for _, a in res)
    if all(a['nframes'] == n for _, a in res):
        return res
    out = []
    for name, a in res:
        a = dict(a)
        a['rg'] = {k: v[:n] for k, v in a['rg'].items()}
        a['mind'] = {k: v[:n] for k, v in a['mind'].items()}
        a['ncon'] = {k: v[:n] for k, v in a['ncon'].items()}
        a['nframes'] = n
        out.append((name, a))
    return out


def main(only=None):
    store = {}
    for s, chains in SETS.items():
        if only and s not in only:
            continue
        sd = ROOT / s
        if not sd.is_dir():
            continue
        reps = [r for r in sorted(sd.glob('rep-*'), key=lambda p: int(p.name.split('-')[1]))
                if (r / f'{sysname(s)}.dcd').exists()]
        res = []
        for r in reps:
            try:
                a = analyse_rep(s, r, chains)
            except Exception as e:
                print(f'  ! {s}/{r.name}: {e}', flush=True); a = None
            if a:
                res.append((r.name, a))
        if len(res) < MIN_REPS:
            print(f'skip {s} ({len(res)} usable reps)', flush=True); continue
        res = trim(res)
        store[f'{s}|binding'] = res
        print(f'ok  {s}: {len(res)} reps x {res[0][1]["nframes"]} analysed frames', flush=True)

    # merge into any existing raw.npy so sets can be processed as they finish
    f = OUT / 'data' / 'raw.npy'
    if store:
        merged = {}
        if f.exists():
            try:
                merged = dict(np.load(f, allow_pickle=True)[0])
            except Exception as e:
                print(f'  ! could not read existing {f.name}: {e}', flush=True)
        merged.update(store)
        np.save(f, np.array([merged], dtype=object), allow_pickle=True)
        print(f'saved {f} ({len(merged)} sets)', flush=True)
    return store


if __name__ == '__main__':
    main(sys.argv[1:] or None)
