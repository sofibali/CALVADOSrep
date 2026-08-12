#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Back-map CALVADOS CG (CA-only) state structures to full all-atom PDBs.

Strategy:
  Each restrained domain is internally rigid in CALVADOS (intra-domain
  harmonic restraints keep it close to the reference geometry). So we can
  rigid-body-align the AF2 all-atom domain onto each state's CG CA
  positions and inherit AF2 side chains.

For each input CG state PDB:
  1. Read CA positions per residue from CG state.
  2. For each domain in fdomains.yaml: Kabsch-align the AF2 domain (all-atom)
     onto the state's CA positions for the same residues → R, t.
  3. Apply R, t to ALL atoms in that domain → all-atom coords.
  4. Linker residues (outside any restraint domain): keep just CA (or
     reconstruct backbone via simple geometry if --reconstruct-backbone).
  5. Write all-atom PDB.

Usage:
    python backmap_states.py \\
        --states representative_frames/2026-05-19/fl_optimized_com_tica/state_*.pdb \\
        --reference input/parp14.pdb \\
        --fdomains fl_optimized/input/domains.yaml \\
        --out-dir backmapped_states/

    python backmap_states.py --auto --set fl_optimized --k 5    # auto-find states

Optionally chain through PULCHRA (if installed in PATH) for linker
reconstruction:
    python backmap_states.py ... --pulchra
"""
import os
import sys
import argparse
import subprocess
import yaml
import numpy as np
from pathlib import Path

CWD = Path(__file__).resolve().parent


# ─────────────────────────────────────────────────────────────────
# PDB I/O
# ─────────────────────────────────────────────────────────────────

def parse_pdb_atoms(pdb_path):
    """Parse a PDB or mmCIF file into a list of atom records.

    Returns list of dicts with keys:
      record, atom_serial, atom_name, alt_loc, res_name, chain_id,
      res_seq, i_code, x, y, z, occ, b, segment, element
    """
    pdb_path = str(pdb_path)
    if pdb_path.lower().endswith('.cif'):
        return _parse_cif_atoms(pdb_path)

    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            try:
                rec = {
                    'record': line[0:6].strip(),
                    'atom_serial': int(line[6:11]),
                    'atom_name': line[12:16].strip(),
                    'alt_loc': line[16:17],
                    'res_name': line[17:20].strip(),
                    'chain_id': line[21:22],
                    'res_seq': int(line[22:26]),
                    'i_code': line[26:27],
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'occ': float(line[54:60]) if line[54:60].strip() else 1.0,
                    'b': float(line[60:66]) if line[60:66].strip() else 0.0,
                    'segment': line[72:76].strip() if len(line) > 76 else '',
                    'element': line[76:78].strip() if len(line) > 76 else '',
                }
                atoms.append(rec)
            except (ValueError, IndexError):
                continue
    return atoms


def _parse_cif_atoms(cif_path):
    """Parse an mmCIF file (AlphaFold3 output) using BioPython."""
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('m', cif_path)
    atoms = []
    serial = 1
    for model in structure:
        for chain in model:
            chain_id = chain.id[:1] if chain.id else 'A'
            for residue in chain:
                if residue.id[0] != ' ':
                    continue
                for atom in residue:
                    coord = atom.get_vector().get_array()
                    atoms.append({
                        'record': 'ATOM',
                        'atom_serial': serial,
                        'atom_name': atom.get_name(),
                        'alt_loc': ' ',
                        'res_name': residue.resname.strip(),
                        'chain_id': chain_id,
                        'res_seq': residue.id[1],
                        'i_code': ' ',
                        'x': float(coord[0]),
                        'y': float(coord[1]),
                        'z': float(coord[2]),
                        'occ': 1.0,
                        'b': float(atom.bfactor),
                        'segment': '',
                        'element': atom.element if atom.element else '',
                    })
                    serial += 1
        break  # only use model 0
    return atoms


def write_pdb_atoms(atoms, out_path):
    """Write atom records to PDB."""
    with open(out_path, 'w') as f:
        for a in atoms:
            name = a['atom_name']
            # PDB column rules: atom name right-justified in cols 13-14 for
            # single-letter elements (C, N, O) with name in cols 14-15
            if len(name) < 4 and not a.get('element', '').strip() in ['', None]:
                name_field = f' {name:<3}'
            else:
                name_field = f'{name:<4}'

            line = (
                f"{a['record']:<6}"
                f"{a['atom_serial']:>5} "
                f"{name_field}"
                f"{a['alt_loc']:1}"
                f"{a['res_name']:<3} "
                f"{a['chain_id']:1}"
                f"{a['res_seq']:>4}"
                f"{a['i_code']:1}   "
                f"{a['x']:8.3f}"
                f"{a['y']:8.3f}"
                f"{a['z']:8.3f}"
                f"{a['occ']:6.2f}"
                f"{a['b']:6.2f}      "
                f"{a['segment']:<4}"
                f"{a['element']:>2}\n"
            )
            f.write(line)
        f.write('END\n')


# ─────────────────────────────────────────────────────────────────
# Kabsch alignment
# ─────────────────────────────────────────────────────────────────

def kabsch(P, Q):
    """Compute rotation R and translation t that align P onto Q.

    Inputs:
      P, Q : (N, 3) arrays
    Returns:
      R : (3, 3) rotation
      t : (3,) translation
      such that Q ≈ (P @ R) + t (when P, Q are centered with their COMs)
    """
    P_c = P.mean(axis=0)
    Q_c = Q.mean(axis=0)
    Pc = P - P_c
    Qc = Q - Q_c
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = Q_c - P_c @ R
    return R, t


# ─────────────────────────────────────────────────────────────────
# Back-mapping core
# ─────────────────────────────────────────────────────────────────

def index_atoms_by_residue(atoms):
    """Return dict {res_seq: [atom indices in atoms list]}."""
    by_res = {}
    for i, a in enumerate(atoms):
        by_res.setdefault(a['res_seq'], []).append(i)
    return by_res


def get_ca_coord(atoms, res_seq):
    """Return (x, y, z) of the CA atom of res_seq, or None."""
    for a in atoms:
        if a['res_seq'] == res_seq and a['atom_name'] == 'CA':
            return np.array([a['x'], a['y'], a['z']])
    return None


def backmap_state(state_pdb, ref_atoms, domains, out_pdb,
                   pulchra=False, verbose=True):
    """Back-map one CG state PDB to all-atom.

    domains : list of (start, end) FL residue ranges
    """
    state_atoms = parse_pdb_atoms(state_pdb)
    state_by_res = index_atoms_by_residue(state_atoms)
    ref_by_res = index_atoms_by_residue(ref_atoms)

    # Build a copy of the reference atom list that we'll transform per domain
    out_atoms = [dict(a) for a in ref_atoms]

    # Track which residues we've placed via domain alignment
    placed_residues = set()

    if verbose:
        print(f"  {state_pdb.name}: {len(state_atoms)} CG atoms, "
              f"{len(domains)} domains to align")

    for ds, de in domains:
        # Gather CA pairs (state, ref) for this domain
        sP = []
        qP = []
        atom_indices = []
        for r in range(ds, de + 1):
            ca_state = get_ca_coord(state_atoms, r)
            ca_ref = get_ca_coord(ref_atoms, r)
            if ca_state is None or ca_ref is None:
                continue
            sP.append(ca_state)
            qP.append(ca_ref)
            atom_indices.extend(ref_by_res.get(r, []))
        if len(sP) < 3:
            continue
        sP = np.array(sP)
        qP = np.array(qP)

        # Kabsch: align ref CAs (qP) onto state CAs (sP)
        R, t = kabsch(qP, sP)

        # Apply transformation to all atoms in this domain
        for ai in atom_indices:
            xyz = np.array([out_atoms[ai]['x'],
                            out_atoms[ai]['y'],
                            out_atoms[ai]['z']])
            new_xyz = xyz @ R + t
            out_atoms[ai]['x'] = float(new_xyz[0])
            out_atoms[ai]['y'] = float(new_xyz[1])
            out_atoms[ai]['z'] = float(new_xyz[2])
            placed_residues.add(out_atoms[ai]['res_seq'])

        rmsd_after = float(np.sqrt(np.mean(
            np.sum(((qP @ R + t) - sP) ** 2, axis=1))))
        if verbose:
            print(f"    domain {ds}-{de}: {len(sP)} CAs, "
                  f"post-alignment RMSD = {rmsd_after:.3f} Å")

    # For linker residues (not in any restraint domain): keep only the CA
    # at the state's CG position; ref side chains aren't representative
    # there since the chain has moved.
    linker_residues = []
    for r in sorted(state_by_res.keys()):
        if r in placed_residues:
            continue
        # Find the state CA for this residue
        ca = get_ca_coord(state_atoms, r)
        if ca is None:
            continue
        linker_residues.append(r)

        # Remove non-CA atoms from out_atoms for this residue, update CA pos
        for ai in ref_by_res.get(r, []):
            atom_name = out_atoms[ai]['atom_name']
            if atom_name == 'CA':
                out_atoms[ai]['x'] = float(ca[0])
                out_atoms[ai]['y'] = float(ca[1])
                out_atoms[ai]['z'] = float(ca[2])

    # Drop linker side-chain atoms (only keep CA for those residues)
    if linker_residues:
        keep = []
        for a in out_atoms:
            if a['res_seq'] in linker_residues and a['atom_name'] != 'CA':
                continue
            keep.append(a)
        out_atoms = keep

    # Renumber atom serials
    for i, a in enumerate(out_atoms, 1):
        a['atom_serial'] = i

    write_pdb_atoms(out_atoms, out_pdb)
    if verbose:
        print(f"    → {out_pdb}  ({len(out_atoms)} atoms, "
              f"{len(linker_residues)} linker residues kept as CA-only)")

    # Optional PULCHRA pass for backbone reconstruction in linker regions
    if pulchra:
        run_pulchra(out_pdb, verbose=verbose)

    return out_pdb


def run_pulchra(pdb_path, verbose=True):
    """Run PULCHRA to fill in missing backbone/side-chain atoms.

    Requires `pulchra` in PATH. Output overwrites the input PDB.
    """
    if subprocess.run(['which', 'pulchra'], capture_output=True).returncode != 0:
        if verbose:
            print(f"    PULCHRA not found in PATH — skipping")
        return False
    try:
        result = subprocess.run(['pulchra', str(pdb_path)],
                                 capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            # PULCHRA writes <input>.rebuilt.pdb — move to overwrite
            rebuilt = pdb_path.with_suffix('.rebuilt.pdb')
            if rebuilt.exists():
                rebuilt.replace(pdb_path)
            if verbose:
                print(f"    PULCHRA: rebuilt {pdb_path.name}")
            return True
    except subprocess.TimeoutExpired:
        if verbose:
            print(f"    PULCHRA: timed out")
    return False


# ─────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────

def auto_find_states(set_key, date_str=None, run_tag=None,
                      source='representative_frames'):
    """Auto-locate state PDBs.

    source='representative_frames': scan
        representative_frames/<date>/<run_tag>/state_*.pdb
    source='pymol_frames': scan
        pymol_frames/<set_key>_*.pdb
    """
    if source == 'pymol_frames':
        pf_dir = CWD / 'pymol_frames'
        if not pf_dir.is_dir():
            return []
        pdbs = sorted(pf_dir.glob(f'{set_key}_*.pdb'))
        # Skip the multi-state 'movie' files (they have multiple MODELs)
        return [p for p in pdbs if 'movie' not in p.name.lower()]

    rep_root = CWD / 'representative_frames'
    if not rep_root.is_dir():
        return []

    if date_str:
        date_dirs = [rep_root / date_str]
    else:
        date_dirs = sorted(rep_root.glob('20*'), reverse=True)

    for d in date_dirs:
        if not d.is_dir():
            continue
        if run_tag:
            tag_dirs = [d / run_tag]
        else:
            tag_dirs = [t for t in d.iterdir()
                        if t.is_dir() and set_key in t.name]
        for t in tag_dirs:
            if not t.is_dir():
                continue
            pdbs = sorted(t.glob('state_*.pdb'))
            if pdbs:
                return pdbs
    return []


def main():
    parser = argparse.ArgumentParser(
        description='Back-map CG state PDBs to all-atom via '
                    'per-domain rigid alignment with AF2.')
    parser.add_argument('--states', nargs='+', type=Path,
                        help='CG state PDB files to back-map')
    parser.add_argument('--auto', action='store_true',
                        help='Auto-find state PDBs (see --source)')
    parser.add_argument('--source',
                        choices=['representative_frames', 'pymol_frames'],
                        default='representative_frames',
                        help='Where to find state PDBs when --auto: '
                             "representative_frames (default; cluster_states.py "
                             "output) or pymol_frames (e.g. md_MD1_buried.pdb, "
                             "fl_compact.pdb).")
    parser.add_argument('--set', dest='set_key', default='fl_optimized',
                        help='Set prefix for --auto (default fl_optimized). '
                             "For pymol_frames source, matches the leading "
                             "set name in filenames (e.g. 'md', 'fl', 'core').")
    parser.add_argument('--run-tag', default=None,
                        help='Specific run_tag for --auto with '
                             'representative_frames source')
    parser.add_argument('--reference', '-r', type=Path,
                        default=None,
                        help='All-atom reference PDB or AF3 .cif. '
                             'Default: auto-pick based on --set. '
                             'fl/fl_optimized → input/parp14.pdb (AF2 FL). '
                             'md/core/mka/... → first AF3 CIF in '
                             'parp14/alphafold_outputs/<af3_dir>/seed-1_sample-0/model.cif')
    parser.add_argument('--fdomains', type=Path, default=None,
                        help='domains.yaml file. Default: '
                             '<set>/input/domains.yaml')
    parser.add_argument('--out-dir', type=Path,
                        default=CWD / 'backmapped_states',
                        help='Output directory')
    parser.add_argument('--pulchra', action='store_true',
                        help='Run PULCHRA after back-mapping to fill in '
                             'linker backbone atoms (requires pulchra in '
                             'PATH)')
    args = parser.parse_args()

    # Resolve state PDBs
    if args.auto:
        states = auto_find_states(args.set_key, run_tag=args.run_tag,
                                    source=args.source)
        if not states:
            print(f"ERROR: no state PDBs found for set={args.set_key}, "
                  f"source={args.source}")
            return 1
        print(f"Auto-found {len(states)} state PDBs ({args.source}):")
        for s in states:
            print(f"  {s.name}")
    elif args.states:
        states = [Path(s).resolve() for s in args.states]
    else:
        parser.error('Must specify --states or --auto')

    # Resolve reference auto-default
    if args.reference is None:
        if args.set_key in ('fl', 'fl_optimized'):
            args.reference = CWD / 'input' / 'parp14.pdb'
        else:
            # Try AF3 CIF for this set's underlying construct
            af3_dir_map = {
                'md':     'md1l1_md2_md3',
                'mka':    'md1l1_md2_md3_khb-kh8_wwe_art',
                'core':   'kh7a_md1l1_md2_md3_khb-kh8_wwe_art',
                'norrm':  'kh1-kh6_kh7a_md1l1_md2_md3_khb-kh8_wwe_art',
                'noart':  'kh1-kh6_kh7a_md1l1_md2_md3_khb-kh8_wwe',
                'md3art': 'md3_khb-kh8_wwe_art',
            }
            af3_dir = af3_dir_map.get(args.set_key, args.set_key)
            cif = (Path('/home/sbali/CALVADOS/parp14/alphafold_outputs') /
                   af3_dir / 'seed-1_sample-0' / 'model.cif')
            if cif.exists():
                args.reference = cif
            else:
                # Fallback: backbone PDB from prepare_all_fragments.py
                args.reference = (CWD / 'fragments' / af3_dir /
                                   'input' / 'parp14.pdb')
        print(f"  Auto-selected reference: {args.reference}")
    if not args.reference.exists():
        print(f"ERROR: reference not found: {args.reference}")
        return 1
    ref_atoms = parse_pdb_atoms(args.reference)
    if not ref_atoms:
        print(f"ERROR: no atoms parsed from {args.reference}")
        return 1
    print(f"\nReference: {args.reference} "
          f"({len(ref_atoms)} atoms, "
          f"{len(set(a['res_seq'] for a in ref_atoms))} residues)")

    # Find domains.yaml
    if args.fdomains is None:
        candidates = [
            CWD / args.set_key / 'input' / 'domains.yaml',
            CWD / 'fragments' / args.set_key / 'input' / 'domains.yaml',
            CWD / 'input' / 'domains.yaml',  # FL default
        ]
        for c in candidates:
            if c.exists():
                args.fdomains = c
                print(f"  Auto-selected fdomains: {c}")
                break
        else:
            print(f"ERROR: --fdomains not given and none of "
                  f"{candidates} found")
            return 1
    if not args.fdomains.exists():
        print(f"ERROR: domains file not found: {args.fdomains}")
        return 1

    with open(args.fdomains) as f:
        domains_yaml = yaml.safe_load(f)
    # The first key in the yaml maps to a list of [start, end] pairs
    first_key = next(iter(domains_yaml))
    domains = [tuple(d) for d in domains_yaml[first_key]]
    print(f"\nDomains ({len(domains)}):")
    for ds, de in domains:
        print(f"  {ds}-{de}  ({de - ds + 1} residues)")

    # Back-map each state
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nBack-mapping to {args.out_dir}/")
    for state_pdb in states:
        out_pdb = args.out_dir / f'{state_pdb.stem}_allatom.pdb'
        try:
            backmap_state(state_pdb, ref_atoms, domains, out_pdb,
                          pulchra=args.pulchra)
        except Exception as e:
            print(f"  FAILED: {state_pdb.name} — {e}")

    print(f"\nDone. {len(states)} states back-mapped to {args.out_dir}/")
    return 0


if __name__ == '__main__':
    sys.exit(main())
