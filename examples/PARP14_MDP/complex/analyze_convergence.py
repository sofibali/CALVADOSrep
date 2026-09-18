#!/usr/bin/env python3
"""Convergence analysis + inter-chain domain contact maps for the PARP14/PARP9/DTX3L
CALVADOS complexes.

Stage 1: per-replicate observables (chain Rg, inter-chain min distance, inter-chain
contact count) -> convergence diagnostics (Gelman-Rubin Rhat across replicates,
integrated autocorrelation time / ESS, split-half contact-map correlation, drift).
Stage 2: domain x domain inter-chain contact-frequency maps, computed only for
(system, variant) combos that pass the convergence criteria.
Stage 3: difference maps between variants / systems.

Conventions follow analyze_complexes.py: 50 ps/frame, CA-CA contact cutoff 1.0 nm,
domains from prepare_complex.DOMAINS. Chains are written whole (no PBC wrapping of
bonds), but inter-chain distances use minimum image with the box from config.yaml.
"""
import os, sys, json, warnings, itertools
from pathlib import Path
import numpy as np, yaml
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
warnings.filterwarnings('ignore')

ROOT   = Path('/home/vle/Vincent/p14fl_test/complexes')
OUT    = Path('/tmp/claude-64170/-home-vle-Vincent/a5f0a8dd-a59a-412d-b58f-35aebe7f4fab/scratchpad/conv')
NS_PER_FRAME = 0.05          # 50 ps
CUTOFF  = 1.0                # nm, CA-CA inter-chain contact
STRIDE  = 10                 # 500 ps between analysed frames
EQ_NS   = 20.0               # discard first 20 ns as equilibration

LENGTHS = {'parp14':1801,'parp9':854,'dtx3l':740,'rrm1-kh1':190}
DOMAINS = {
 'parp14':{'RRM1':(7,88),'RRM2':(150,223),'RRM3':(229,304),'KH1':(322,397),'KH2':(398,461),
           'KH3':(467,529),'KH4':(532,595),'KH5':(596,668),'KH6':(670,736),'1/2KH7a':(739,776),
           'MD1':(792,978),'MD2':(1007,1190),'MD3':(1219,1386),'1/2KH7b':(1427,1453),
           'KH8':(1454,1533),'WWE':(1535,1599),'ART':(1618,1800)},
 'parp9':{'1/2KH1a':(60,96),'MD1':(107,296),'MD2':(311,493),'1/2KH1b':(516,539),
          'KH2':(543,622),'ART':(644,823)},
 'dtx3l':{'RRM':(11,87),'KH1':(138,201),'KH2':(233,302),'KH3':(305,362),'KH4':(378,447),
          'KH5':(449,507),'RING':(556,607),'DTC':(608,737)},
}
SYS_CHAINS = {
 'parp14_parp9':['parp14','parp9'], 'parp14_dtx3l':['parp14','dtx3l'],
 'parp9_dtx3l':['parp9','dtx3l'],   'parp14_parp9_dtx3l':['parp14','parp9','dtx3l'],
 'dtx3l_dtx3l':['dtx3l','dtx3l'],
}

def labels_for(chain):
    """Per-residue block id + ordered block names (domains, then 'linker')."""
    names=[n for n,_ in blocks_for(chain)]
    lab=np.full(LENGTHS[chain], -1, np.int64)
    for k,(n,idx) in enumerate(blocks_for(chain)): lab[idx]=k
    return lab, names


def blocks_for(chain):
    """Domain blocks + one 'linker' bucket, as 0-based index arrays within the chain."""
    n = LENGTHS[chain]; doms = DOMAINS[chain]
    covered = np.zeros(n, bool); out = []
    for name,(a,b) in doms.items():
        idx = np.arange(a-1, min(b, n)); covered[idx] = True
        out.append((name, idx))
    link = np.where(~covered)[0]
    if len(link): out.append(('linker', link))
    return out

# ---------- convergence statistics ----------
def autocorr_time(x):
    """Integrated autocorrelation time, initial-positive-sequence estimator."""
    x = np.asarray(x, float); x = x - x.mean()
    n = len(x)
    if n < 10 or x.std() == 0: return np.nan
    f = np.fft.rfft(x, 2*n); ac = np.fft.irfft(f*np.conj(f))[:n].real
    ac /= ac[0]
    tau = 1.0
    for k in range(1, n):
        if ac[k] <= 0: break
        tau += 2*ac[k]
        if k > 5*tau: break
    return tau

def gelman_rubin(chains):
    """Rhat over replicate chains (list of 1-D arrays), using each chain's 2nd half."""
    ch = [np.asarray(c,float)[len(c)//2:] for c in chains if len(c) >= 8]
    if len(ch) < 2: return np.nan
    n = min(len(c) for c in ch); ch = np.array([c[:n] for c in ch]); m = len(ch)
    means = ch.mean(1); W = ch.var(1, ddof=1).mean()
    if W == 0: return np.nan
    B = n*means.var(ddof=1)
    var = (n-1)/n*W + B/n
    return float(np.sqrt(var/W))

def drift_frac(x):
    """|slope*span| / std over the post-equilibration series (0 = flat)."""
    x = np.asarray(x,float); n=len(x)
    if n < 10 or x.std()==0: return np.nan
    t = np.arange(n); s = np.polyfit(t,x,1)[0]
    return float(abs(s*n)/x.std())

# ---------- per-replicate pass ----------
def analyse_rep(system, rep, chains):
    top = rep/'input'/'system.pdb'
    dcd = rep/f'{system}.dcd'
    if not top.exists() or not dcd.exists(): return None
    try: cfg = yaml.safe_load(open(rep/'config.yaml'))
    except Exception: return None
    box = np.array(cfg['box'], float)
    u = mda.Universe(str(top), str(dcd))
    if u.atoms.n_atoms != sum(LENGTHS[c] for c in chains): return None
    offs, o = [], 0
    for c in chains: offs.append((c,o,o+LENGTHS[c])); o += LENGTHS[c]
    pairs = list(itertools.combinations(range(len(chains)), 2))
    rg   = {i: [] for i in range(len(chains))}
    mind = {p: [] for p in pairs}
    ncon = {p: [] for p in pairs}
    cmaps= {p: None for p in pairs}
    nfr  = 0
    lab = {c: labels_for(c) for c in set(chains)}
    boxv = np.array([box[0]*10, box[1]*10, box[2]*10, 90., 90., 90.], np.float32)
    cut_A = CUTOFF*10.0
    for p in pairs:
        i,j = p
        cmaps[p] = np.zeros((len(lab[chains[i]][1]), len(lab[chains[j]][1])))
    for ts in u.trajectory[::STRIDE]:
        P = u.atoms.positions.astype(np.float32)          # Angstrom
        for i,(c,a,b) in enumerate(offs):
            X = P[a:b]/10.0
            rg[i].append(float(np.sqrt(((X-X.mean(0))**2).sum(1).mean())))
        for p in pairs:
            i,j = p
            A = P[offs[i][1]:offs[i][2]]; B = P[offs[j][1]:offs[j][2]]
            D = distance_array(A, B, box=boxv)
            mind[p].append(float(D.min())/10.0)
            ii, jj = np.nonzero(D < cut_A)
            ncon[p].append(int(ii.size))
            if ii.size:
                li = lab[chains[i]][0][ii]; lj = lab[chains[j]][0][jj]
                M = cmaps[p]
                np.add.at(M, (li, lj), 1.0)
        nfr += 1
    for p in pairs: cmaps[p] = cmaps[p]/max(nfr,1)
    return dict(rg={i:np.array(v) for i,v in rg.items()},
                mind={p:np.array(v) for p,v in mind.items()},
                ncon={p:np.array(v) for p,v in ncon.items()},
                cmaps=cmaps, nframes=nfr, chains=chains,
                blocks={c:labels_for(c)[1] for c in set(chains)})

def discover():
    jobs=[]
    for system, chains in SYS_CHAINS.items():
        sd = ROOT/system
        if not sd.is_dir(): continue
        for v in sorted(sd.iterdir()):
            if not v.is_dir(): continue
            if v.name.startswith('_quarantine'): continue
            if v.name in ('input','analysis'): continue
            reps = [r for r in sorted(v.iterdir()) if r.is_dir() and (r/f'{system}.dcd').exists()]
            if len(reps) >= 3: jobs.append((system, v.name, chains, reps))
        # flat layout (rep dirs directly under the system)
        flat=[r for r in sorted(sd.iterdir()) if r.is_dir() and (r/f'{system}.dcd').exists()]
        if len(flat) >= 3: jobs.append((system,'(flat)',chains,flat))
    return jobs

if __name__ == '__main__':
    jobs = discover()
    print(f'{len(jobs)} (system, variant) combos', flush=True)
    store={}
    for system, variant, chains, reps in jobs:
        res=[]
        for r in reps:
            try: a = analyse_rep(system, r, chains)
            except Exception as e: print(f'  ! {system}/{variant}/{r.name}: {e}', flush=True); a=None
            if a: res.append((r.name,a))
        if len(res) < 3: print(f'skip {system}/{variant} ({len(res)} reps)', flush=True); continue
        store[f'{system}|{variant}'] = res
        print(f'ok  {system}/{variant}: {len(res)} reps x {res[0][1]["nframes"]} frames', flush=True)
    np.save(OUT/'data'/'raw.npy', np.array([store], dtype=object), allow_pickle=True)
    print('saved raw.npy', flush=True)
