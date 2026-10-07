#!/usr/bin/env python3
"""Design an experimentally testable interface mutation panel from the contact data.

WHY MUTATIONS RATHER THAN A FORCE-FIELD NULL
--------------------------------------------
Zeroing the Ashbaugh-Hatch lambda asks whether CALVADOS's potential produces
the observed contacts. That validates the model and has no bench counterpart.
A point mutation is the opposite: the same construct can be made, purified and
measured, so the simulation makes a falsifiable prediction.

HOW THE PANEL IS BUILT
----------------------
1. Per-residue inter-chain contact frequency, pooled over replicates, using the
   same 1.0 nm CA-CA criterion as the rest of the analysis.
2. Keep only surface-exposed residues -- a buried residue cannot be at an
   interface, and mutating it tests folding, not binding.
3. Rank by contact frequency and propose a substitution whose direction is
   unambiguous in BOTH the model and the experiment:

     charged residue   -> charge reversal (R/K -> E, D/E -> K)
                          tests electrostatic steering; in CALVADOS this flips
                          q, so the predicted effect is large and signed.
     hydrophobic (high lambda, neutral) -> Ser
                          tests the hydrophobic contribution; drops lambda with
                          little size change.

4. Report dlambda and dq so the predicted effect size is explicit and the
   in-silico result can be checked against the measured one.

A mutation that the model says matters and the bench says does not is the most
informative outcome available here -- it bounds how far CALVADOS can be trusted
on folded-domain surfaces, which is the open question behind the HyRes
disagreement.

Usage:  python design_interface_mutations.py [sets ...] [--top 12] [--target-frames 500]
"""
import sys
from pathlib import Path
from argparse import ArgumentParser
from collections import defaultdict
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
sys.path.insert(0, str(HERE))
from figure_face_separation import analyse, SETS
import analyze_convergence as ac

RES = HERE / 'binding' / 'p9_dtx3l' / 'input' / 'residues_CALVADOS3.csv'
THREE = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN',
         'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
         'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
         'Y': 'TYR', 'V': 'VAL'}
DEFAULT_SETS = ['p14_homo', 'p14_dtx3l', 'p9_dtx3l']


def propose(aa, tbl):
    """(target aa, rationale) or None if this residue is not informative."""
    q = tbl.loc[aa, 'q']; lam = tbl.loc[aa, 'lambdas']
    if q > 0:
        return 'E', 'charge reversal (+ -> -): tests electrostatic steering'
    if q < 0:
        return 'K', 'charge reversal (- -> +): tests electrostatic steering'
    if lam >= 0.5:
        return 'S', 'hydrophobic -> polar: tests the hydrophobic contribution'
    return None


def main():
    ap = ArgumentParser()
    ap.add_argument('sets', nargs='*', default=DEFAULT_SETS)
    ap.add_argument('--top', type=int, default=12)
    ap.add_argument('--target-frames', type=int, default=500)
    a = ap.parse_args()

    tbl = pd.read_csv(RES).set_index('one')
    rows = []
    for s in a.sets:
        r = analyse(s, a.target_frames)
        if r is None:
            print(f'  {s}: no data'); continue
        # per-residue contacts, restricted to exposed residues in a restrained domain
        for cidx, c, lo, hi in r['offs']:
            fm = r['face_map'].get(cidx, {})
            cand = []
            for resid, v in fm.items():
                dom, face, exposed = v
                if not exposed:
                    continue
                cf = r['per_res'].get((cidx, resid), 0.0)
                if cf <= 0:
                    continue
                cand.append((cf, resid, dom, face))
            cand.sort(reverse=True)
            seq = ac  # placeholder; residue identity comes from the PDB below
            for cf, resid, dom, face in cand[:a.top]:
                rows.append(dict(set=s, chain=f'{c}#{cidx+1}', domain=dom,
                                 resid=resid, face=face,
                                 contacts_per_frame=round(cf, 4)))
        print(f'  {s}: ranked', flush=True)

    if not rows:
        print('no candidates'); return
    df = pd.DataFrame(rows)

    # residue identity from the input PDB (CALVADOS reads the sequence from there)
    import MDAnalysis as mda, warnings
    warnings.filterwarnings('ignore')
    from MDAnalysis.lib.util import convert_aa_code
    ident = {}
    for s in df.set.unique():
        for c in set(SETS[s]):
            p = HERE / 'binding' / s / 'input' / f'{c}.pdb'
            if not p.is_file():
                continue
            u = mda.Universe(str(p))
            ca = u.select_atoms('name CA')
            for rid, rn in zip(ca.resids, ca.resnames):
                try:
                    ident[(c, int(rid))] = convert_aa_code(rn)
                except Exception:
                    pass
    df['aa'] = [ident.get((ch.split('#')[0], int(rid))) for ch, rid in
                zip(df.chain, df.resid)]
    df = df[df.aa.notna()].copy()

    keep = []
    for _, r in df.iterrows():
        pr = propose(r.aa, tbl)
        if pr is None:
            continue
        to, why = pr
        keep.append(dict(**r.to_dict(), mutate_to=to,
                         mutation=f'{r.aa}{int(r.resid)}{to}',
                         pdb_rename=f'{THREE[r.aa]}->{THREE[to]}',
                         d_lambda=round(float(tbl.loc[to, 'lambdas'] - tbl.loc[r.aa, 'lambdas']), 3),
                         d_charge=round(float(tbl.loc[to, 'q'] - tbl.loc[r.aa, 'q']), 1),
                         rationale=why))
    out = pd.DataFrame(keep).sort_values(['set', 'contacts_per_frame'],
                                         ascending=[True, False])
    out.to_csv(OUT / 'interface_mutation_panel.csv', index=False)
    cols = ['set', 'chain', 'domain', 'face', 'mutation', 'contacts_per_frame',
            'd_lambda', 'd_charge', 'rationale']
    print('\n' + out[cols].to_string(index=False))
    print(f'\nsaved {OUT/"interface_mutation_panel.csv"}  ({len(out)} candidates)')
    print('\nTo build a mutant: rename that residue in binding/<set>/input/<chain>.pdb '
          '(CALVADOS reads the sequence from the PDB via seq_from_pdb).')


if __name__ == '__main__':
    main()
