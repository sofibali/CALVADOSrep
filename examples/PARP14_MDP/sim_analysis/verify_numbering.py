#!/usr/bin/env python
"""
Residue-numbering consistency checker for the PARP14 project.

Verifies, against the authoritative UniProt Q460N5 (PAR14_HUMAN, 1801 aa)
sequence, that everything downstream agrees on what residue N is:

  1. the local reference FASTA is byte-identical to UniProt Q460N5
  2. every slab construct's PDB is exactly the FL slice it claims to be,
     with the residue count `slab_meta.yaml` declares
  3. every construct's `domains.yaml` equals the FL structured cores,
     mapped into that construct's local numbering
  4. the FL structured cores match UniProt's own Domain annotations
  5. DOMAIN_UNITS tiles 1..1801 without overlaps
  6. the catalytic residues in `active_sites.yaml` have the amino-acid
     identity they are labelled with, none of the retired (mislabelled)
     values have been reinstated, and the literature/UniProt anchor residues
     still match

Checks 5 and 6 failed when this was written on 2026-09-18; both were fixed on
2026-09-21 (see docs/NUMBERING_AUDIT.md). All checks currently pass.

    python verify_numbering.py             # all checks, offline where possible
    python verify_numbering.py --offline   # skip the UniProt fetch

Exit status is non-zero if any check fails, so this is usable in CI.
"""
import argparse
import json
import os
import sys
import urllib.request

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                      # examples/PARP14_MDP
REPO = os.path.dirname(os.path.dirname(ROOT))     # CALVADOS
PARP14 = os.path.join(REPO, 'parp14')

UNIPROT = 'Q460N5'
FASTA_URL = f'https://rest.uniprot.org/uniprotkb/{UNIPROT}.fasta'
FEAT_URL = (f'https://rest.uniprot.org/uniprotkb/{UNIPROT}.json'
            '?fields=ft_act_site,ft_binding,ft_domain')

# The FL structured cores used for domain restraints (slab/*/fl/input/domains.yaml).
FL_CORES = [[6, 88], [150, 223], [227, 301], [791, 978],
            [1003, 1190], [1216, 1387], [1523, 1601], [1605, 1801]]

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}

# Residues whose identity is asserted somewhere (active_sites.yaml provenance
# comments, CLAUDE.md, the Dukic et al. 2023 mutants, or the canonical PARP
# motifs). If any of these stops matching Q460N5, a numbering error has crept
# back in.
CLAIMED_IDENTITY = {
    # MD1 -- conserved macrodomain Asn, the Dukic G832E site, the UniProt Asp
    824: 'N', 832: 'G', 961: 'D',
    # MD2 -- UniProt binding site, and the residue equivalent to MD1 G832
    1034: 'S', 1044: 'G',
    # MD3 -- UniProt binding sites
    1247: 'S', 1258: 'V', 1371: 'F',
    # ART -- canonical PARP H-G-T and G-K-G-T-Y-F-A motifs, plus Dukic R1699A
    1682: 'H', 1699: 'R', 1714: 'Y',
}

# Superseded values, kept so the audit fails loudly if anyone reinstates them.
RETIRED_IDENTITY = {831: 'D', 923: 'N', 962: 'D',
                    1684: 'H', 1705: 'Y', 1706: 'E', 1722: 'I'}

fails = []


def check(name, ok, detail=''):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ''))
    if not ok:
        fails.append(name)
    return ok


def read_fasta(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>'))


def pdb_ca_sequence(path):
    """One-letter sequence from the first chain's CA atoms, in file order."""
    seq, seen = [], set()
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            key = (line[21], line[22:27])
            if key in seen:
                continue
            seen.add(key)
            seq.append(AA3.get(line[17:20].strip(), 'X'))
    return ''.join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--offline', action='store_true',
                    help='skip the UniProt fetch; use the local FASTA as reference')
    args = ap.parse_args()

    sys.path.insert(0, ROOT)
    import sim_registry as reg

    local_fasta = os.path.join(PARP14, 'PARP14.fasta')
    fl = read_fasta(local_fasta)

    # ---- 1. local FASTA vs UniProt -------------------------------------
    print(f"\n1. Reference sequence ({UNIPROT})")
    check('local FASTA is 1801 aa', len(fl) == 1801, f'got {len(fl)}')
    uni_feats = None
    if not args.offline:
        try:
            with urllib.request.urlopen(FASTA_URL, timeout=30) as r:
                up = ''.join(l.strip() for l in r.read().decode().splitlines()
                             if not l.startswith('>'))
            check('local FASTA byte-identical to UniProt', up == fl)
            with urllib.request.urlopen(FEAT_URL, timeout=30) as r:
                uni_feats = json.load(r)
        except Exception as exc:                       # network optional
            print(f"  [SKIP] UniProt fetch unavailable ({exc})")

    # ---- 2/3. slab constructs ------------------------------------------
    print("\n2/3. Slab constructs: sequence, residue count, domain mapping")
    slab = os.path.join(ROOT, 'slab')
    cfl = reg.CONTIGUOUS_FL_RANGE
    n_checked = 0
    for arm in ('homotypic', 'rna'):
        armd = os.path.join(slab, arm)
        if not os.path.isdir(armd):
            continue
        for name in sorted(os.listdir(armd)):
            d = os.path.join(armd, name)
            fmeta = os.path.join(d, 'slab_meta.yaml')
            if not os.path.isfile(fmeta):
                continue
            meta = yaml.safe_load(open(fmeta))
            # full-length constructs are the identity map, not in CONTIGUOUS_FL_RANGE
            a, b = cfl.get(name, (1, 1801)) if name != 'fl' else (1, 1801)
            pdb = os.path.join(d, 'input', f"{meta['sysname']}.pdb")
            if not os.path.isfile(pdb):
                check(f'{arm}/{name}: PDB present', False, pdb)
                continue
            seq = pdb_ca_sequence(pdb)
            ok_seq = seq == fl[a - 1:b]
            ok_n = len(seq) == meta['nres'] == (b - a + 1)
            dom = [list(x) for x in
                   list(yaml.safe_load(open(os.path.join(d, 'input', 'domains.yaml'))).values())[0]]
            expect = [[max(s, a) - a + 1, min(e, b) - a + 1]
                      for s, e in FL_CORES if max(s, a) <= min(e, b)]
            ok_dom = dom == expect
            check(f'{arm}/{name} (FL {a}-{b}, {b - a + 1} res)',
                  ok_seq and ok_n and ok_dom,
                  '' if (ok_seq and ok_n and ok_dom) else
                  f"seq={ok_seq} count={ok_n} domains={ok_dom}")
            n_checked += 1
    print(f"  ({n_checked} construct directories checked)")

    # ---- 4. FL cores vs UniProt domains --------------------------------
    print("\n4. FL structured cores vs UniProt Domain annotations")
    if uni_feats:
        udom = {}
        for f in uni_feats.get('features', []):
            if f['type'] == 'Domain':
                udom[f['description']] = (f['location']['start']['value'],
                                          f['location']['end']['value'])
        for desc, (s, e) in sorted(udom.items(), key=lambda kv: kv[1]):
            hit = [c for c in FL_CORES if c == [s, e]]
            check(f'UniProt {desc} {s}-{e} present in domains.yaml', bool(hit))
    else:
        print("  [SKIP] needs UniProt features")

    # ---- 5. DOMAIN_UNITS tiling ----------------------------------------
    print("\n5. DOMAIN_UNITS tiling of 1..1801")
    units = sorted(reg.DOMAIN_UNITS.items(), key=lambda kv: kv[1][0])
    total = sum(e - s + 1 for _, (s, e) in units)
    cov = set()
    for _, (s, e) in units:
        cov |= set(range(s, e + 1))
    overlap = total - len(cov)
    check('no residue assigned to two units', overlap == 0,
          f'{overlap} residue(s) double-counted' if overlap else '')
    if overlap:
        prev = None
        for n, (s, e) in units:
            if prev and s <= prev[1]:
                print(f"         {prev[0]} ({prev[1]}) and {n} ({s}) share residue {s}")
            prev = (n, e)
    gap = sorted(set(range(1, 1802)) - cov)
    print(f"  [note] {len(gap)} residue(s) intentionally unassigned "
          f"(MD2/MD3 linker): {gap[0]}-{gap[-1]}" if gap else "")

    # ---- 6. catalytic residue identities -------------------------------
    print("\n6. Catalytic residue identities vs Q460N5")
    sites = yaml.safe_load(open(os.path.join(PARP14, 'input', 'active_sites.yaml')))
    bad = []
    for dom, info in sites.items():
        if not isinstance(info, dict) or 'catalytic_residues' not in info:
            continue
        for r in info['catalytic_residues']:
            aa = fl[r - 1]
            want = CLAIMED_IDENTITY.get(r)
            if want and aa != want:
                bad.append((dom, r, want, aa))
            if r in RETIRED_IDENTITY:
                bad.append((dom, r, 'RETIRED value reinstated', aa))
    check('every labelled catalytic residue has its labelled identity', not bad)

    # The anchors themselves must still be what the provenance says they are,
    # whether or not they appear in a catalytic_residues list.
    anchors = [(r, want, fl[r - 1]) for r, want in sorted(CLAIMED_IDENTITY.items())
               if fl[r - 1] != want]
    check('literature/UniProt anchor residues match Q460N5', not anchors)
    for r, want, aa in anchors:
        print(f'         anchor {want}{r} is actually {aa}')
    for dom, r, want, aa in bad:
        print(f"         {dom}: labelled {want}{r} but Q460N5 residue {r} is {aa}")
    if bad:
        print("\n  Canonical PARP motifs actually present in Q460N5:")
        for motif in ('HGT', 'GKGTYFA'):
            i = fl.find(motif)
            if i >= 0:
                print(f"         {motif} at {i + 1}: "
                      + ' '.join(f'{fl[i + k]}{i + 1 + k}' for k in range(len(motif))))
        print("         -> catalytic His is H1682 (H-G-T), catalytic Tyr is Y1714 (G-K-G-T-Y-F-A)")
        print("         -> PARP14 is a MONO-ART; its triad is H-Y-[I/L], not H-Y-E")

    print(f"\n{'=' * 64}")
    if fails:
        print(f"FAILED {len(fails)} check(s):")
        for f in fails:
            print(f"   - {f}")
    else:
        print("All checks passed.")
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
