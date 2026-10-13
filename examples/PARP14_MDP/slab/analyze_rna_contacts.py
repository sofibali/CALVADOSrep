#!/usr/bin/env python
"""
What the slab run-2 panel says about RNA, domain compaction and condensation.

Three questions, one pass over each trajectory:

1. **Which domains touch the RNA?**  Protein-RNA contacts within `--cutoff`,
   resolved per PARP14 domain unit and normalised per residue, so a 215-residue
   macrodomain and a 52-residue KH domain are comparable.

2. **Which domains compact against each other?**  Protein-protein contacts
   BETWEEN DIFFERENT CHAINS, resolved per domain pair. Intra-chain contacts are
   excluded on purpose: those are largely fixed by the restraints and say
   nothing about assembly. This is the inter-molecular contact map that
   distinguishes a condensate's sticky surfaces from an inert one.

3. **Does it condense at all?**  Slab width and c_sat, read back from the
   per-run `*_ps_results.csv` / analysis logs that `analyze_slab.py` wrote.

BEAD LAYOUT (verified, not assumed)
-----------------------------------
Protein chains carry one bead per residue. RNA chains carry TWO beads per
nucleotide -- `rP` (phosphate) and `rN` (base) -- so a 40-mer is 80 beads, and
n_atoms exceeds n_residues by exactly the nucleotide count. Counting RNA
"residues" as beads would undercount its contact surface by half.

PERIODICITY
-----------
Boxes are 40 x 40 x ~520 nm and the slab sits in the middle of z, but chains do
wrap in x/y. Contacts are therefore computed with a periodic cKDTree
(`boxsize=`), after wrapping coordinates into [0, L). Ignoring that would drop
every contact across an x/y boundary.

USAGE
    python analyze_rna_contacts.py --all --workers 12
    python analyze_rna_contacts.py rnapanel/md_full__rna__seed23456

Writes into `rnapanel/analysis/`, one directory next to the runs it describes:

    rnapanel_condensation_summary.csv        slab width + c_sat, one row per run
    rnapanel_rna_contacts_by_domain.csv      +RNA runs ONLY (see below)
    rnapanel_interchain_domain_contacts.csv  one row per domain pair per run

Every table carries `contact_cutoff_nm` and `frames_used`, so a number can
always be traced back to how it was measured.
"""
import argparse
import glob
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import yaml

warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                      # examples/PARP14_MDP
PANEL = os.path.join(HERE, 'rnapanel')
OUT_DEFAULT = os.path.join(PANEL, 'analysis')

# FL domain units, UniProt Q460N5 numbering (CLAUDE.md).
DOMAIN_UNITS = [
    ('RRM1',    1,  145), ('RRM2',  146,  224), ('RRM3',  225,  314),
    ('KH1-KH6', 315, 737), ('KH7a',  738,  789), ('MD1',   790, 1003),
    ('MD2',    1004, 1193), ('MD3',  1207, 1388), ('KHb-KH8', 1389, 1533),
    ('WWE',    1534, 1602), ('ART',  1603, 1801),
]

# The construct spans, as prepare_slab_rna.py defines them.
SPANS = {
    'fl': (1, 1801), 'kh1_art_full': (315, 1801), 'core_full_go': (738, 1801),
    'md_full': (790, 1388), 'fl_khb': (1, 1533), 'kh1_khb': (315, 1533),
    'core_khb': (738, 1533),
}


def domain_map(construct, nres):
    """Local (1-based) residue ranges for each FL domain unit inside a construct."""
    lo, hi = SPANS[construct]
    out = []
    for name, a, b in DOMAIN_UNITS:
        s, e = max(a, lo), min(b, hi)
        if s > e:
            continue
        ls, le = s - lo + 1, e - lo + 1
        if le > nres:
            le = nres
        if ls <= le:
            out.append((name, ls, le))
    return out


def read_condensation(run_dir):
    """Slab width and c_sat, from the analyze_slab.py log for this run."""
    name = os.path.basename(run_dir.rstrip('/'))
    log = os.path.join(HERE, 'logs', f'analysis_{name}.log')
    if not os.path.isfile(log):
        return {}
    txt = open(log, errors='ignore').read().replace('\r', '\n')
    m = re.search(r'^rnapanel/\S+\s+(\d+)\s*->\s*(\d+)\s+(\S+)', txt, re.M)
    if not m:
        return {}
    csat = m.group(3)
    return dict(slab_eq_nm=int(m.group(1)), slab_end_nm=int(m.group(2)),
                c_sat_mM=(np.nan if not re.match(r'^[\d.]+$', csat) else float(csat)),
                phase=('dispersed' if csat.upper().startswith('NO') else 'condensed'))


def analyse(run_dir, cutoff=1.0, nframes=80, discard=0.5):
    import mdtraj as md
    from scipy.spatial import cKDTree

    run_dir = run_dir.rstrip('/')
    meta = yaml.safe_load(open(os.path.join(run_dir, 'slab_meta.yaml')))
    construct, arm, seed = meta['construct'], meta['arm'], meta['seed']
    nres, nchain = meta['nres'], meta['nchain']
    sysname = meta['sysname']

    dcd = os.path.join(run_dir, f'{sysname}.dcd')
    top = os.path.join(run_dir, 'checkpoint.pdb')
    if not (os.path.isfile(dcd) and os.path.isfile(top)):
        return None

    topo = md.load_frame(dcd, 0, top=top).topology
    n_prot_atoms = nchain * nres
    n_atoms = topo.n_atoms

    # RNA beads are 2 per nucleotide (rP + rN); everything after the protein.
    rna_idx = np.arange(n_prot_atoms, n_atoms)
    prot_idx = np.arange(n_prot_atoms)

    # atom -> chain, and atom -> domain id (-1 = linker / not in a unit)
    chain_of = prot_idx // nres
    local_res = prot_idx % nres + 1
    doms = domain_map(construct, nres)
    dom_of = np.full(n_prot_atoms, -1, dtype=np.int16)
    for k, (dname, ls, le) in enumerate(doms):
        dom_of[(local_res >= ls) & (local_res <= le)] = k
    dom_sizes = np.array([(le - ls + 1) for _, ls, le in doms], float)

    total = md.load(dcd, top=top).n_frames if False else None  # avoid full load
    with md.open(dcd) as fh:
        total = len(fh)
    start = int(total * discard)
    idx = np.unique(np.linspace(start, total - 1, min(nframes, total - start)).astype(int))

    nd = len(doms)
    rna_hits = np.zeros(nd)                      # protein-RNA contacts per domain
    pair_hits = np.zeros((nd, nd))               # inter-CHAIN protein-protein
    nf = 0
    for fi in idx:
        fr = md.load_frame(dcd, int(fi), top=top)
        box = fr.unitcell_lengths[0].astype(float)
        xyz = np.mod(fr.xyz[0].astype(float), box)      # wrap into [0,L)

        if len(rna_idx):
            tp = cKDTree(xyz[prot_idx], boxsize=box)
            tr = cKDTree(xyz[rna_idx], boxsize=box)
            for pi, lst in enumerate(tp.query_ball_tree(tr, cutoff)):
                if lst:
                    d = dom_of[pi]
                    if d >= 0:
                        rna_hits[d] += len(lst)

        tp = cKDTree(xyz[prot_idx], boxsize=box)
        pairs = tp.query_pairs(cutoff, output_type='ndarray')
        if len(pairs):
            a, b = pairs[:, 0], pairs[:, 1]
            inter = chain_of[a] != chain_of[b]          # different molecules only
            a, b = a[inter], b[inter]
            da, db = dom_of[a], dom_of[b]
            ok = (da >= 0) & (db >= 0)
            da, db = da[ok], db[ok]
            np.add.at(pair_hits, (da, db), 1)
            np.add.at(pair_hits, (db, da), 1)
        nf += 1

    rna_rows, pair_rows = [], []
    n_rna_nt = meta.get('n_rna', 0) * 40
    # Homotypic runs contain no RNA, so a "0 contacts" row there is a STRUCTURAL
    # zero, not a measurement. Emitting them invites averaging the two arms
    # together, which halves every per-residue number and inflates its SD to
    # ~100% of the mean. No RNA -> no rows.
    for k, (dname, ls, le) in enumerate(doms):
        if n_rna_nt == 0:
            break
        rna_rows.append(dict(
            run=os.path.basename(run_dir), construct=construct, arm=arm, seed=seed,
            domain=dname, dom_residues=int(dom_sizes[k]),
            rna_contacts_per_frame_allchains=rna_hits[k] / max(nf, 1),
            rna_contacts_per_residue_per_frame=rna_hits[k] / max(nf, 1) / (dom_sizes[k] * nchain),
            n_rna_nt=n_rna_nt, contact_cutoff_nm=cutoff, frames_used=nf))
    for i, (di, _, _) in enumerate(doms):
        for j, (dj, _, _) in enumerate(doms):
            if j < i:
                continue
            v = pair_hits[i, j] / max(nf, 1)
            if i == j:
                v /= 2.0           # self-pairs double counted by the symmetric add
            pair_rows.append(dict(
                run=os.path.basename(run_dir), construct=construct, arm=arm, seed=seed,
                domain_i=di, domain_j=dj,
                interchain_contacts_per_frame=v,
                interchain_contacts_per_residuepair_per_frame=(
                    v / (dom_sizes[i] * dom_sizes[j] * nchain)),
                contact_cutoff_nm=cutoff, frames_used=nf))

    cond = dict(run=os.path.basename(run_dir), construct=construct, arm=arm,
                seed=seed, nchain=nchain, nres=nres, n_rna=meta.get('n_rna', 0),
                frames_used=nf, **read_condensation(run_dir))
    return pd.DataFrame(rna_rows), pd.DataFrame(pair_rows), cond


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', nargs='?')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--cutoff', type=float, default=1.0)
    ap.add_argument('--frames', type=int, default=80)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--out', default=OUT_DEFAULT)
    a = ap.parse_args()

    if a.all:
        runs = sorted(d for d in glob.glob(os.path.join(PANEL, '*'))
                      if os.path.isfile(os.path.join(d, 'slab_meta.yaml'))
                      and not os.path.isfile(os.path.join(d, '.parked')))
    elif a.run_dir:
        runs = [os.path.abspath(a.run_dir)]
    else:
        ap.error('give a run dir or --all')

    os.makedirs(a.out, exist_ok=True)
    from concurrent.futures import ProcessPoolExecutor, as_completed
    rna, pairs, cond = [], [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(analyse, r, a.cutoff, a.frames): r for r in runs}
        for f in as_completed(futs):
            r = futs[f]
            try:
                res = f.result()
            except Exception as e:
                print(f'  FAILED {os.path.basename(r)}: {e}', flush=True)
                continue
            if res is None:
                print(f'  skipped {os.path.basename(r)} (no trajectory)', flush=True)
                continue
            rr, pp, cc = res
            rna.append(rr); pairs.append(pp); cond.append(cc)
            print(f'  done {os.path.basename(r)}', flush=True)

    if not cond:
        print('nothing analysed'); return
    rna = [d for d in rna if len(d)]
    if rna:
        pd.concat(rna).sort_values(['construct', 'seed', 'domain']).to_csv(
            os.path.join(a.out, 'rnapanel_rna_contacts_by_domain.csv'), index=False)
    pd.concat(pairs).sort_values(['construct', 'arm', 'seed']).to_csv(
        os.path.join(a.out, 'rnapanel_interchain_domain_contacts.csv'), index=False)
    pd.DataFrame(cond).sort_values(['construct', 'arm', 'seed']).to_csv(
        os.path.join(a.out, 'rnapanel_condensation_summary.csv'), index=False)
    print(f'\nwrote 3 tables to {a.out}')


if __name__ == '__main__':
    main()
