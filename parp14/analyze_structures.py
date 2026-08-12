#!/usr/bin/env python
"""
Per-domain RMSD and inter-domain contact analysis across all PARP14 AF3 structures.

Inputs
------
- atS domain definitions: parp14/input/PARP14_domains_atS.fasta (11 domains)
- FL PARP14 sequence:     parp14/input/PARP14.fasta
- CSV combinatorial units: parp14/input/domain_boundaries.csv (used for AF3 input gen)
- AF3 outputs:            parp14/alphafold_outputs/<combo>/seed-{1-5}_sample-{0-4}/model.cif
- Crystal references:     examples/PARP14_MDP/input/xtal_refs/{1x4r,3goy,3q6z,3vfq}.pdb

Outputs (under parp14/analysis/structure_analysis/)
---------------------------------------------------
- reference/atS_map.json                            atS FL residue ranges, CSV unit residue ranges, FL ref CA coords (per atS domain)
- per_combo/<combo>__summary.parquet                model x domain summary: rmsd_to_FLref, n_res, mean_sasa, n_exposed
- per_combo/<combo>__dompairs.parquet               model x dom_i x dom_j contact counts (multi-cutoff, multi-type)
- per_combo/<combo>__kpairs.parquet                 every K-K and K-acidic contact pair (filtered post-hoc)
- per_combo/<combo>__residues.parquet               per residue: aa, sasa_abs, sasa_rel, exposed, ca_x/y/z, atS_domain
- crystal/                                           same layout, source='crystal'
- aggregated/                                        unioned tables for plotting

Each per_combo file's existence implies the combo is done (resumable).
"""

import argparse
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

ROOT = Path('/home/sbali/CALVADOS/parp14')
AF3_OUT = ROOT / 'alphafold_outputs'
ANALYSIS_DIR = ROOT / 'analysis' / 'structure_analysis'
REF_DIR = ANALYSIS_DIR / 'reference'
PER_COMBO_DIR = ANALYSIS_DIR / 'per_combo'
CRYSTAL_DIR = ANALYSIS_DIR / 'crystal'
AGG_DIR = ANALYSIS_DIR / 'aggregated'

FL_FASTA = ROOT / 'input' / 'PARP14.fasta'
ATS_FASTA = ROOT / 'input' / 'PARP14_domains_atS.fasta'
CSV_DOMAINS = ROOT / 'input' / 'domain_boundaries.csv'
XTAL_DIR = Path('/home/sbali/CALVADOS/examples/PARP14_MDP/input/xtal_refs')

ATS_DOMAIN_ORDER = ['RRM1', 'RRM2', 'RRM3', 'KH_N', 'KH_Ca',
                    'MD1', 'MD2', 'MD3', 'KH_b', 'WWE', 'ART']

# Crystal -> atS domain mapping
XTAL_DOMAIN = {'1x4r': 'WWE', '3goy': 'ART', '3q6z': 'MD1', '3vfq': 'MD2'}

# Contact cutoffs
SC_CUTOFFS_A = [6.0, 8.0, 10.0]
CA_CUTOFFS_A = [20.0, 25.0, 30.0]
EXPOSURE_THRESH = 0.20  # SASA relative > 20% = exposed
ACIDIC = {'D', 'E'}

# ============================================================================
# Reference setup
# ============================================================================

def read_fasta(path):
    out = {}
    cur = None
    for l in open(path):
        l = l.strip()
        if l.startswith('>'):
            cur = l[1:].split()[0]; out[cur] = ''
        elif cur:
            out[cur] += l
    return out


def load_atS_map():
    """Map each atS domain to its FL residue range."""
    fl_seqs = read_fasta(FL_FASTA)
    fl = list(fl_seqs.values())[0]
    ats = read_fasta(ATS_FASTA)
    ats_ranges = {}
    for d, s in ats.items():
        p = fl.find(s)
        if p < 0:
            raise RuntimeError(f'atS {d} not found in FL')
        ats_ranges[d] = (p + 1, p + len(s))  # 1-based inclusive
    return ats_ranges, fl


def load_csv_units():
    df = pd.read_csv(CSV_DOMAINS)
    # canonical name (matches AF3 dir token after lowercasing + md1l1->md1)
    units = {}
    for _, r in df.iterrows():
        token = r['Domain'].lower().replace('md1l1', 'md1')
        units[token] = (int(r['Start']), int(r['End']))
    # canonical order (file order)
    order = [r['Domain'].lower().replace('md1l1', 'md1') for _, r in df.iterrows()]
    return units, order


def parse_combo_name(name):
    """Lowercase combo dir name -> ordered list of CSV unit tokens."""
    units, order = load_csv_units()
    name_l = name.lower()
    # strip timestamp suffix
    parts = name_l.split('_')
    # rebuild tokens, accounting for multi-piece tokens like 'kh1-kh6' and 'khb-kh8'
    # since dir uses underscores between units, just match each unit token greedily
    present = []
    i = 0
    while i < len(parts):
        # try multi-token match (e.g., 'kh1-kh6' is single but already includes hyphen)
        token = parts[i]
        if token in units:
            present.append(token)
            i += 1
        else:
            # ignore (timestamps etc.)
            i += 1
    return present


def construct_residue_to_fl_from_seq(construct_seq, fl_seq, seed_k=30):
    """
    Build construct->FL residue mapping by sequence matching.
    The construct sequence is a concatenation of contiguous FL chunks.
    For each chunk, we find a unique seed k-mer in FL and extend the match.
    Returns array c2fl of length N (N = construct residues), c2fl[i] = FL residue (1-based).
    """
    n = len(construct_seq)
    c2fl = np.zeros(n, dtype=np.int32)
    pos = 0
    while pos < n:
        K = min(seed_k, n - pos)
        if K < 5:
            break
        # try shrinking seed if not unique
        for k_try in range(K, 4, -2):
            seed = construct_seq[pos:pos + k_try]
            occurrences = []
            start = 0
            while True:
                p = fl_seq.find(seed, start)
                if p < 0: break
                occurrences.append(p)
                start = p + 1
            if len(occurrences) == 1:
                fl_pos = occurrences[0]
                # extend as far as possible
                i = 0
                while (pos + i < n and fl_pos + i < len(fl_seq) and
                       construct_seq[pos + i] == fl_seq[fl_pos + i]):
                    c2fl[pos + i] = fl_pos + i + 1
                    i += 1
                pos += i
                break
        else:
            # No unique seed found; bail
            raise ValueError(f'Cannot uniquely match construct[{pos}:] starting with '
                             f'{construct_seq[pos:pos+seed_k]!r} in FL')
    if pos < n:
        # Trailing residues couldn't be matched; truncate
        c2fl = c2fl[:pos]
    return c2fl


def build_FL_reference():
    """
    Identify the FL all-domain combo, load its rank-0 model, extract per-atS-domain Cα coords.
    """
    units, order = load_csv_units()
    fl_combo = '_'.join(order)
    # Find the dir (try original, then any timestamped version)
    candidates = [AF3_OUT / fl_combo]
    candidates += sorted(AF3_OUT.glob(f'{fl_combo}_2*'))
    fl_dir = None
    for c in candidates:
        # rank-0 model is seed-1_sample-0/model.cif by AF3 convention
        if (c / 'seed-1_sample-0' / 'model.cif').exists():
            fl_dir = c; break
    if fl_dir is None:
        # fallback: any seed/sample present
        for c in candidates:
            for s in c.glob('seed-*_sample-*'):
                if (s / 'model.cif').exists():
                    fl_dir = c; break
            if fl_dir: break
    if fl_dir is None:
        raise RuntimeError(f'FL combo {fl_combo} not found in {AF3_OUT}')
    cif = next(p for p in [fl_dir / 'seed-1_sample-0' / 'model.cif'] +
                          sorted(fl_dir.glob('seed-*_sample-*/model.cif'))
               if p.exists())
    print(f'  FL reference: {cif}')
    return cif, fl_dir


# ============================================================================
# Structure parsing
# ============================================================================

import gemmi  # noqa: E402

# 3-letter to 1-letter amino acid lookup
AA3to1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M',
}

BACKBONE = {'N', 'CA', 'C', 'O'}


def parse_structure(path):
    """
    Returns dict:
        ca: (N,3) Cα coords (Å)
        sc: list of (M,3) side-chain heavy-atom coords per residue
        seq: str of length N (1-letter)
        resids: (N,) residue numbers as int (chain A)
    Uses chain A only, first model.
    """
    st = gemmi.read_structure(str(path))
    model = st[0]
    chain = model[0]  # AF3 outputs use chain A as first; crystals: chain A first too
    if chain.name not in ('A', 'A1'):
        # find chain 'A'
        for c in model:
            if c.name == 'A':
                chain = c; break
    ca, sc, seq, resids = [], [], [], []
    for res in chain:
        if res.name not in AA3to1:
            continue
        aa = AA3to1[res.name]
        ca_atom = None
        sc_xyz = []
        for atom in res:
            if atom.element.name == 'H':
                continue
            if atom.name == 'CA':
                ca_atom = atom.pos
            if atom.name not in BACKBONE:
                sc_xyz.append([atom.pos.x, atom.pos.y, atom.pos.z])
        if ca_atom is None:
            continue
        # If glycine (no SC), use CA as side-chain proxy
        if not sc_xyz:
            sc_xyz = [[ca_atom.x, ca_atom.y, ca_atom.z]]
        ca.append([ca_atom.x, ca_atom.y, ca_atom.z])
        sc.append(np.asarray(sc_xyz, dtype=np.float32))
        seq.append(aa)
        resids.append(int(res.seqid.num))
    return {
        'ca': np.asarray(ca, dtype=np.float32),
        'sc': sc,
        'seq': ''.join(seq),
        'resids': np.asarray(resids, dtype=np.int32),
    }


# ============================================================================
# SASA via freesasa
# ============================================================================

import freesasa  # noqa: E402

# Reference (Tien 2013 max ASA values, Å²)
MAX_ASA = {
    'A': 121, 'R': 265, 'N': 187, 'D': 187, 'C': 148, 'E': 214, 'Q': 214,
    'G':  97, 'H': 216, 'I': 195, 'L': 191, 'K': 230, 'M': 203, 'F': 228,
    'P': 154, 'S': 143, 'T': 163, 'W': 264, 'Y': 255, 'V': 165,
}


def compute_sasa(path):
    """Return per-residue absolute SASA (Å²) for chain A, in residue order."""
    # freesasa accepts PDB or struct; for CIF we need a temp PDB or use Structure
    # Use gemmi to convert to PDB string in memory
    st = gemmi.read_structure(str(path))
    pdb_str = st.make_pdb_string()
    # freesasa.Structure accepts file path; write temp
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.pdb', delete=False) as f:
        f.write(pdb_str)
        tmp_path = f.name
    try:
        struct = freesasa.Structure(tmp_path)
        result = freesasa.calc(struct)
        residues = result.residueAreas()
        # residues is dict[chain][resi_str] -> ResidueArea
        sasa_per_res = []
        if 'A' in residues:
            for resi_str, area in residues['A'].items():
                sasa_per_res.append((int(resi_str), area.total))
        sasa_per_res.sort()
        return dict(sasa_per_res)  # resi -> total SASA
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass


# ============================================================================
# Geometry: Kabsch alignment + RMSD
# ============================================================================

def kabsch_rmsd(P, Q):
    """RMSD after optimal Kabsch alignment of P onto Q. P, Q: (N,3)."""
    if P.shape != Q.shape or P.shape[0] < 3:
        return np.nan
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    Pa = Pc @ R.T
    return float(np.sqrt(np.mean(np.sum((Pa - Qc) ** 2, axis=1))))


# ============================================================================
# Per-domain processing for one structure
# ============================================================================

def find_atS_in_construct(c2fl, ats_effective):
    """
    Given construct->FL mapping array (1-based FL residues at each construct index)
    and atS effective residue sets (atS range ∩ FL combo coverage),
    return dict atS_name -> array of construct-residue indices (0-based) that
    fully cover the effective set, in FL order matching the FL reference.
    """
    fl_set = set(int(x) for x in c2fl)
    # Map FL residue -> construct index (0-based)
    fl2c = {int(c2fl[i]): i for i in range(len(c2fl))}
    out = {}
    for d, eff_set in ats_effective.items():
        if eff_set.issubset(fl_set):
            ordered = sorted(eff_set)
            out[d] = np.array([fl2c[r] for r in ordered], dtype=np.int32)
    return out


def build_atS_effective(ats_ranges, fl_c2fl):
    """For each atS domain, compute effective FL residue set = atS range ∩ FL combo coverage."""
    fl_set = set(int(x) for x in fl_c2fl)
    out = {}
    for d, (s, e) in ats_ranges.items():
        eff = {r for r in range(s, e + 1) if r in fl_set}
        out[d] = eff
    return out


def analyze_structure(cif_path, combo, model_label, c2fl, ats_in_construct,
                      fl_ref_atS_coords, source='af3', prebuilt_struct=None):
    """
    Process one model file. Returns:
        summary_rows: list of dicts
        dompair_rows: list of dicts
        kpair_rows:   list of dicts
        residue_rows: list of dicts
    """
    s = prebuilt_struct if prebuilt_struct is not None else parse_structure(cif_path)
    if len(s['seq']) != len(c2fl):
        # unexpected; pad/truncate as best-effort but flag
        n = min(len(s['seq']), len(c2fl))
        s = {k: (v[:n] if hasattr(v, '__getitem__') else v) for k, v in s.items()}
        s['ca'] = s['ca'][:n]
        s['sc'] = s['sc'][:n]
        s['resids'] = s['resids'][:n]
        s['seq'] = s['seq'][:n]
        c2fl = c2fl[:n]
    # SASA: compute on full structure file (use cif_path).
    # If prebuilt_struct masked away residues, the SASA is keyed on real residue numbers
    # which we filter by the kept resids below.
    sasa_dict = {}
    try:
        sasa_dict = compute_sasa(cif_path)
    except Exception:
        sasa_dict = {}
    sasa_abs = np.full(len(s['seq']), np.nan, dtype=np.float32)
    for i, r in enumerate(s['resids']):
        sasa_abs[i] = sasa_dict.get(int(r), np.nan)
    sasa_rel = np.full_like(sasa_abs, np.nan)
    for i, aa in enumerate(s['seq']):
        m = MAX_ASA.get(aa, np.nan)
        if not np.isnan(sasa_abs[i]) and m > 0:
            sasa_rel[i] = sasa_abs[i] / m
    exposed = sasa_rel > EXPOSURE_THRESH
    # Domain assignment per residue
    res_atS = np.array(['none'] * len(s['seq']), dtype=object)
    for d, idx in ats_in_construct.items():
        res_atS[idx] = d

    # ----- Summary (per atS domain) -----
    summary_rows = []
    for d, idx in ats_in_construct.items():
        if len(idx) < 3:
            continue
        ref = fl_ref_atS_coords.get(d)
        rmsd = np.nan
        if ref is not None and len(ref) == len(idx):
            rmsd = kabsch_rmsd(s['ca'][idx], ref)
        summary_rows.append({
            'source': source, 'combo': combo, 'model': model_label, 'atS_domain': d,
            'n_res': int(len(idx)),
            'rmsd_to_FLref': rmsd,
            'mean_sasa_abs': float(np.nanmean(sasa_abs[idx])),
            'mean_sasa_rel': float(np.nanmean(sasa_rel[idx])),
            'n_exposed': int(np.nansum(exposed[idx])),
            'frac_exposed': float(np.nanmean(exposed[idx])),
        })

    # ----- Inter-domain contacts -----
    # Per-residue side-chain min coord arrays for fast distance computation:
    # Use side-chain centroid? No — use per-atom distances. Build flat arrays per domain.
    dompair_rows = []
    kpair_rows = []
    domains_present = sorted(ats_in_construct.keys())
    # Pre-extract per-domain CA + side-chain atoms with residue indices
    dom_data = {}
    for d in domains_present:
        idx = ats_in_construct[d]
        # Flatten side-chain atoms with parent residue idx
        sc_atoms = []
        sc_owner = []
        for ridx in idx:
            atoms = s['sc'][ridx]
            sc_atoms.append(atoms)
            sc_owner.extend([ridx] * len(atoms))
        if sc_atoms:
            sc_flat = np.concatenate(sc_atoms, axis=0)
        else:
            sc_flat = np.zeros((0, 3), dtype=np.float32)
        dom_data[d] = {
            'idx': idx,
            'ca': s['ca'][idx],
            'sc_flat': sc_flat,
            'sc_owner': np.asarray(sc_owner, dtype=np.int32),
        }

    for i, di in enumerate(domains_present):
        for dj in domains_present[i + 1:]:
            A = dom_data[di]; B = dom_data[dj]
            # CA distance matrix (small, e.g. 200x200)
            ca_d = np.linalg.norm(A['ca'][:, None, :] - B['ca'][None, :, :], axis=-1)
            # SC distance: for each pair of residues, minimum atom-atom distance
            # Build by iterating residue pairs using sc_flat blocks
            # Faster: compute full atom-atom distance and reduce per residue pair
            # But that's O(NaA * NaB). For ~200x200 res, ~10 atoms each = 2k x 2k = 4M.
            # OK at small scale.
            sc_dmin = np.full(ca_d.shape, np.inf, dtype=np.float32)
            if A['sc_flat'].shape[0] > 0 and B['sc_flat'].shape[0] > 0:
                # full atom-atom matrix
                D2 = np.sum((A['sc_flat'][:, None, :] - B['sc_flat'][None, :, :]) ** 2, axis=-1)
                # reduce per residue pair via owner indices
                # remap owners to local indices
                a_owner_map = {r: k for k, r in enumerate(A['idx'])}
                b_owner_map = {r: k for k, r in enumerate(B['idx'])}
                a_local = np.array([a_owner_map[r] for r in A['sc_owner']])
                b_local = np.array([b_owner_map[r] for r in B['sc_owner']])
                # For each (a_atom, b_atom), update sc_dmin[a_local[a], b_local[b]]
                for ai in range(D2.shape[0]):
                    block = D2[ai]
                    al = a_local[ai]
                    # min with current row
                    np.minimum.at(sc_dmin[al], b_local, np.sqrt(block))
            # Counts at each cutoff
            row = {'source': source, 'combo': combo, 'model': model_label,
                   'dom_i': di, 'dom_j': dj}
            for cut in SC_CUTOFFS_A:
                mask = sc_dmin < cut
                row[f'sc_n_{int(cut)}A'] = int(mask.sum())
                # exposed-only: both residues exposed
                exA = exposed[A['idx']]; exB = exposed[B['idx']]
                exp_mask = mask & exA[:, None] & exB[None, :]
                row[f'sc_n_{int(cut)}A_exposed'] = int(exp_mask.sum())
            for cut in CA_CUTOFFS_A:
                row[f'ca_n_{int(cut)}A'] = int((ca_d < cut).sum())
            # K-acidic and K-K specific (sidechain 8Å, exposed both)
            seqA = np.array(list(s['seq']))[A['idx']]
            seqB = np.array(list(s['seq']))[B['idx']]
            kA = (seqA == 'K'); kB = (seqB == 'K')
            acA = np.isin(seqA, list(ACIDIC)); acB = np.isin(seqB, list(ACIDIC))
            for cut in (8.0,):
                m = sc_dmin < cut
                exp_m = m & exposed[A['idx']][:, None] & exposed[B['idx']][None, :]
                row['sc8_KE_or_KD'] = int((m & ((kA[:, None] & acB[None, :]) | (acA[:, None] & kB[None, :]))).sum())
                row['sc8_KE_or_KD_exposed'] = int((exp_m & ((kA[:, None] & acB[None, :]) | (acA[:, None] & kB[None, :]))).sum())
                row['sc8_KK'] = int((m & (kA[:, None] & kB[None, :])).sum())
                row['sc8_KK_exposed'] = int((exp_m & (kA[:, None] & kB[None, :])).sum())
            dompair_rows.append(row)

            # ----- Per-residue K pairs (filter K-K and K-acidic, sc 8Å OR ca 25Å) -----
            ax = np.where(kA[:, None] & (kB[None, :] | acB[None, :]))
            bx = np.where(acA[:, None] & kB[None, :])
            # Combine; but only emit pairs with sc<10 OR ca<30 (loose so caller can filter)
            cands = list(zip(ax[0].tolist(), ax[1].tolist())) + list(zip(bx[0].tolist(), bx[1].tolist()))
            cands = list({(a, b) for a, b in cands})
            for (a, b) in cands:
                sc_d = float(sc_dmin[a, b])
                cad = float(ca_d[a, b])
                if not (sc_d < 12.0 or cad < 30.0):
                    continue
                kpair_rows.append({
                    'source': source, 'combo': combo, 'model': model_label,
                    'dom_i': di, 'dom_j': dj,
                    'res_i_construct': int(A['idx'][a] + 1),
                    'res_j_construct': int(B['idx'][b] + 1),
                    'res_i_fl': int(c2fl[A['idx'][a]]),
                    'res_j_fl': int(c2fl[B['idx'][b]]),
                    'aa_i': str(seqA[a]), 'aa_j': str(seqB[b]),
                    'sc_dist_A': sc_d, 'ca_dist_A': cad,
                    'exposed_i': bool(exposed[A['idx'][a]]),
                    'exposed_j': bool(exposed[B['idx'][b]]),
                })

    # ----- Per-residue table -----
    residue_rows = []
    for i in range(len(s['seq'])):
        residue_rows.append({
            'source': source, 'combo': combo, 'model': model_label,
            'res_construct': int(i + 1),
            'res_fl': int(c2fl[i]),
            'aa': s['seq'][i],
            'atS_domain': str(res_atS[i]),
            'sasa_abs': float(sasa_abs[i]),
            'sasa_rel': float(sasa_rel[i]),
            'exposed': bool(exposed[i]),
            'ca_x': float(s['ca'][i, 0]),
            'ca_y': float(s['ca'][i, 1]),
            'ca_z': float(s['ca'][i, 2]),
        })
    return summary_rows, dompair_rows, kpair_rows, residue_rows


# ============================================================================
# Combo dispatcher
# ============================================================================

def find_combo_models(combo_dir):
    """Return list of (model_label, cif_path) for one combo dir."""
    out = []
    for sub in sorted(combo_dir.glob('seed-*_sample-*')):
        cif = sub / 'model.cif'
        if cif.exists():
            out.append((sub.name, cif))
    return out


def find_combo_dirs():
    """Return dict combo_name -> [list of dirs (preferred first)]. Picks best dir per combo."""
    # Walk top-level dirs that have at least one seed-*_sample-* subdir
    bins = {}
    for d in AF3_OUT.iterdir():
        if not d.is_dir():
            continue
        name = d.name
        # Strip timestamp suffix _2026... to recover canonical name
        base = name
        if '_2' in name:
            i = name.rfind('_2')
            tail = name[i + 1:]
            if tail.startswith('2') and len(tail) >= 8 and tail[:4].isdigit():
                base = name[:i]
        bins.setdefault(base, []).append(d)
    # Sort each list: prefer dir without timestamp first, then most recent timestamps
    out = {}
    for base, dirs in bins.items():
        dirs.sort(key=lambda p: (0 if p.name == base else 1, p.name), reverse=False)
        # Filter only those with at least one seed-*_sample-*
        valid = [d for d in dirs if any(d.glob('seed-*_sample-*/model.cif'))]
        if valid:
            out[base] = valid
    return out


def process_combo(combo_name, combo_dirs, ats_effective, fl_seq, fl_ref_atS_coords):
    """
    Process all models from one combo (concatenating across timestamped dirs, dedup by seed/sample).
    Skip if all per_combo files exist.
    """
    out_summary = PER_COMBO_DIR / f'{combo_name}__summary.parquet'
    out_dompair = PER_COMBO_DIR / f'{combo_name}__dompairs.parquet'
    out_kpair = PER_COMBO_DIR / f'{combo_name}__kpairs.parquet'
    out_res = PER_COMBO_DIR / f'{combo_name}__residues.parquet'
    if all(p.exists() for p in (out_summary, out_dompair, out_kpair, out_res)):
        return combo_name, 'skip', 0

    # Collect models from all timestamped dirs (use the first instance per seed/sample)
    seen = {}
    for d in combo_dirs:
        for label, cif in find_combo_models(d):
            if label not in seen:
                seen[label] = cif
    if not seen:
        return combo_name, 'skip-no-models', 0

    sums, dps, kps, ress = [], [], [], []
    c2fl = None  # built per model (fast); cache from first model and reuse across same construct
    ats_in_construct = None
    for label, cif in seen.items():
        try:
            if c2fl is None:
                # First model: parse and derive c2fl from sequence
                first_struct = parse_structure(cif)
                c2fl = construct_residue_to_fl_from_seq(first_struct['seq'], fl_seq)
                ats_in_construct = find_atS_in_construct(c2fl, ats_effective)
                if not ats_in_construct:
                    return combo_name, 'skip-no-atS', 0
            a, b, c, r = analyze_structure(cif, combo_name, label, c2fl,
                                           ats_in_construct, fl_ref_atS_coords)
            sums.extend(a); dps.extend(b); kps.extend(c); ress.extend(r)
        except Exception as e:
            print(f'  ERROR {combo_name}/{label}: {e}', file=sys.stderr)

    if not sums:
        return combo_name, 'no-output', 0

    pd.DataFrame(sums).to_parquet(out_summary, index=False)
    pd.DataFrame(dps).to_parquet(out_dompair, index=False)
    if kps:
        pd.DataFrame(kps).to_parquet(out_kpair, index=False)
    else:
        pd.DataFrame(columns=['source','combo','model','dom_i','dom_j',
                              'res_i_construct','res_j_construct','res_i_fl','res_j_fl',
                              'aa_i','aa_j','sc_dist_A','ca_dist_A',
                              'exposed_i','exposed_j']).to_parquet(out_kpair, index=False)
    pd.DataFrame(ress).to_parquet(out_res, index=False)
    return combo_name, 'ok', len(seen)


# ============================================================================
# Subcommands
# ============================================================================

def cmd_setup():
    PER_COMBO_DIR.mkdir(parents=True, exist_ok=True)
    REF_DIR.mkdir(parents=True, exist_ok=True)
    CRYSTAL_DIR.mkdir(parents=True, exist_ok=True)
    AGG_DIR.mkdir(parents=True, exist_ok=True)

    print('Loading atS domain ranges...')
    ats_ranges, fl_seq = load_atS_map()
    csv_units, csv_order = load_csv_units()
    print(f'  atS domains: {len(ats_ranges)}')
    print(f'  CSV units:   {len(csv_units)}')

    print('Building FL reference (atS Cα per domain)...')
    fl_cif, fl_dir = build_FL_reference()
    fl_struct_for_seq = parse_structure(fl_cif)
    fl_c2fl = construct_residue_to_fl_from_seq(fl_struct_for_seq['seq'], fl_seq)
    ats_effective = build_atS_effective(ats_ranges, fl_c2fl)
    for d in ATS_DOMAIN_ORDER:
        full = ats_ranges[d][1] - ats_ranges[d][0] + 1
        eff = len(ats_effective[d])
        if eff < full:
            print(f'  atS {d}: {eff}/{full} residues (linker/gap residues missing in FL CSV coverage)')
    fl_atS_in = find_atS_in_construct(fl_c2fl, ats_effective)
    fl_struct = fl_struct_for_seq
    fl_ref = {}
    fl_ref_residues = {}  # atS_domain -> list of FL residue numbers (the "effective" ordered set)
    for d, idx in fl_atS_in.items():
        fl_ref[d] = fl_struct['ca'][idx].astype(np.float32)
        fl_ref_residues[d] = sorted(ats_effective[d])
        print(f'  {d}: {len(idx)} residues')

    np.savez(REF_DIR / 'fl_ref_coords.npz', **fl_ref)
    json.dump({
        'atS_ranges': {k: list(v) for k, v in ats_ranges.items()},
        'atS_effective_residues': {k: sorted(list(v)) for k, v in ats_effective.items()},
        'csv_units': {k: list(v) for k, v in csv_units.items()},
        'csv_order': csv_order,
        'fl_combo_dir': str(fl_dir),
        'fl_ref_cif': str(fl_cif),
        'fl_ref_residues': fl_ref_residues,
    }, open(REF_DIR / 'atS_map.json', 'w'), indent=2)
    print('Saved reference data to', REF_DIR)


def _load_runtime_refs():
    """Load atS_effective sets, FL master sequence, and FL ref coords."""
    ats_map = json.load(open(REF_DIR / 'atS_map.json'))
    ats_effective = {k: set(v) for k, v in ats_map['atS_effective_residues'].items()}
    fl_seq = list(read_fasta(FL_FASTA).values())[0]
    fl_ref = dict(np.load(REF_DIR / 'fl_ref_coords.npz'))
    return ats_effective, fl_seq, fl_ref, ats_map


def cmd_test(combo=None):
    """Run on one combo to verify."""
    PER_COMBO_DIR.mkdir(parents=True, exist_ok=True)
    ats_effective, fl_seq, fl_ref, _ = _load_runtime_refs()
    if combo is None:
        combo = 'rrm1_md1_md2_md3_wwe_art'  # mid-size example
    combo_dirs = find_combo_dirs().get(combo)
    if not combo_dirs:
        print(f'No combo {combo} found.')
        return
    t0 = time.time()
    name, status, n = process_combo(combo, combo_dirs, ats_effective, fl_seq, fl_ref)
    print(f'{name}: {status} ({n} models) in {time.time()-t0:.1f}s')


def cmd_process(workers=16, only=None):
    PER_COMBO_DIR.mkdir(parents=True, exist_ok=True)
    ats_effective, fl_seq, fl_ref, _ = _load_runtime_refs()

    combos = find_combo_dirs()
    if only:
        combos = {k: v for k, v in combos.items() if k in set(only)}
    print(f'Combos to process: {len(combos)}')

    items = list(combos.items())
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(process_combo, name, dirs, ats_effective, fl_seq, fl_ref):
                name for name, dirs in items}
        n_done = 0; n_ok = 0; n_skip = 0; n_models = 0
        for fut in as_completed(futs):
            name = futs[fut]
            try:
                _, status, n = fut.result()
                n_done += 1
                if status == 'ok':
                    n_ok += 1; n_models += n
                else:
                    n_skip += 1
                if n_done % 25 == 0 or n_done == len(items):
                    el = time.time() - t0
                    rate = n_done / max(el, 1e-6)
                    eta = (len(items) - n_done) / rate
                    print(f'  [{n_done}/{len(items)}] ok={n_ok} skip={n_skip} models={n_models} '
                          f'elapsed={el/60:.1f}min eta={eta/60:.1f}min')
            except Exception as e:
                print(f'  FAIL {name}: {e}', file=sys.stderr)


def _xtal_build_c2fl(s, fl_seq, pdb_name_for_log='?'):
    """
    Build construct->FL residue mapping for a crystal structure.
    If residues use PARP14 numbering (resid > 200 generally), use directly.
    Else, align crystal sequence to FL via unique-seed search.
    Returns (c2fl array, mask of crystal residues to keep).
    """
    rids = s['resids']
    cseq = s['seq']
    # Heuristic: PARP14 numbering iff most residues > 200
    if (rids > 200).mean() > 0.8:
        # Verify by sequence: pick a unique-in-FL k-mer near middle and compute shift.
        mid = len(cseq) // 2
        for k in (30, 20, 15, 10):
            seed = cseq[mid:mid + k]
            occ = []
            start = 0
            while True:
                p = fl_seq.find(seed, start)
                if p < 0: break
                occ.append(p); start = p + 1
            if len(occ) == 1:
                # Shift = actual FL residue at mid - crystal numbering at mid
                shift = (occ[0] + 1) - int(rids[mid])
                if abs(shift) > 0:
                    print(f'    [warn] {pdb_name_for_log}: PDB numbering offset by {shift} '
                          f'(rids[{mid}]={rids[mid]}, FL match={occ[0]+1}); applying shift')
                fl_resids = rids.astype(np.int32) + shift
                return fl_resids, np.ones(len(rids), dtype=bool)
        # Fallback: use rids directly (with warning)
        print(f'    [warn] {pdb_name_for_log}: could not verify PDB numbering; using rids directly')
        return rids.astype(np.int32), np.ones(len(rids), dtype=bool)
    # Local numbering: search by sliding seeds
    fl_resids = np.zeros(len(cseq), dtype=np.int32)
    pos = 0
    while pos < len(cseq):
        # Find next unique k-mer match
        K = min(30, len(cseq) - pos)
        if K < 5: break
        for k_try in range(K, 4, -2):
            seed = cseq[pos:pos + k_try]
            occ = []
            start = 0
            while True:
                p = fl_seq.find(seed, start)
                if p < 0: break
                occ.append(p); start = p + 1
            if len(occ) == 1:
                fl_pos = occ[0]
                i = 0
                while (pos + i < len(cseq) and fl_pos + i < len(fl_seq) and
                       cseq[pos + i] == fl_seq[fl_pos + i]):
                    fl_resids[pos + i] = fl_pos + i + 1
                    i += 1
                pos += i
                break
        else:
            # Couldn't match here; mark as 0 (will be skipped) and advance
            pos += 1
    mask = fl_resids > 0
    return fl_resids, mask


def cmd_xtal():
    """Process crystal PDBs through the full analyze_structure pipeline (RMSD + SASA + contacts).
    Saves analogous parquet files under crystal/ for direct comparison.
    """
    CRYSTAL_DIR.mkdir(parents=True, exist_ok=True)
    ats_effective, fl_seq, fl_ref, _ = _load_runtime_refs()

    sums, dps, kps, ress = [], [], [], []
    for pdb in sorted(XTAL_DIR.glob('*.pdb')):
        name = pdb.stem.lower()
        s = parse_structure(pdb)
        c2fl, mask = _xtal_build_c2fl(s, fl_seq, pdb_name_for_log=name)
        if mask.sum() < 10:
            print(f'  {name}: too few residues ({mask.sum()}); skip')
            continue
        # Apply mask: keep only residues with valid FL position
        if not mask.all():
            keep = np.where(mask)[0]
            s = {
                'ca': s['ca'][keep],
                'sc': [s['sc'][i] for i in keep],
                'seq': ''.join(s['seq'][i] for i in keep),
                'resids': s['resids'][keep],
            }
            c2fl = c2fl[keep]
        ats_in_construct = find_atS_in_construct(c2fl, ats_effective)
        if not ats_in_construct:
            # Try partial: any atS domain with ≥50% of effective set present in crystal
            fl_set = set(int(x) for x in c2fl)
            for d, eff in ats_effective.items():
                inter = eff & fl_set
                if len(inter) >= max(20, 0.5 * len(eff)):
                    ordered = sorted(inter)
                    fl2c = {int(c2fl[i]): i for i in range(len(c2fl))}
                    ats_in_construct[d] = np.array([fl2c[r] for r in ordered], dtype=np.int32)
        if not ats_in_construct:
            print(f'  {name}: no atS domains map; skip')
            continue
        # For RMSD, the FL ref is keyed on the FULL effective set (e.g. 201 for ART).
        # If crystal has a SUBSET, align to the corresponding subset of fl_ref.
        # Build subset reference per domain
        fl_ref_subsets = {}
        for d, idx in ats_in_construct.items():
            full_eff = sorted(ats_effective[d])
            present_fl = [int(c2fl[i]) for i in idx]
            sub_idx = [full_eff.index(r) for r in present_fl]
            fl_ref_subsets[d] = fl_ref[d][sub_idx]
        # Run analyze_structure manually-ish: but we need cif path for SASA. Save tmp PDB from s.
        # Easier: write a small temp PDB representation. parse_structure was already done; for SASA, we need a path. Use the original PDB.
        try:
            a, b, c, r = analyze_structure(pdb, name, 'A', c2fl, ats_in_construct,
                                           fl_ref_subsets, source='crystal',
                                           prebuilt_struct=s)
            sums.extend(a); dps.extend(b); kps.extend(c); ress.extend(r)
            doms_str = ','.join(sorted(ats_in_construct.keys()))
            sub_n = {d: len(idx) for d, idx in ats_in_construct.items()}
            print(f'  {name}: {doms_str} [{sub_n}], len(ca)={len(s["seq"])}')
        except Exception as e:
            print(f'  {name}: ERROR {e}')

    if sums:
        pd.DataFrame(sums).to_parquet(CRYSTAL_DIR / 'crystal_summary.parquet', index=False)
        pd.DataFrame(dps).to_parquet(CRYSTAL_DIR / 'crystal_dompairs.parquet', index=False)
        pd.DataFrame(kps).to_parquet(CRYSTAL_DIR / 'crystal_kpairs.parquet', index=False)
        pd.DataFrame(ress).to_parquet(CRYSTAL_DIR / 'crystal_residues.parquet', index=False)
        # Print RMSD summary
        df = pd.DataFrame(sums)
        print('\nCrystal RMSDs to FL_AF3 reference:')
        for _, r in df.iterrows():
            print(f'  {r["combo"]:>6} -> {r["atS_domain"]:>6}: '
                  f'n_res={r["n_res"]:>4}  RMSD={r["rmsd_to_FLref"]:.2f} Å')


def cmd_aggregate():
    AGG_DIR.mkdir(parents=True, exist_ok=True)
    summary_files = sorted(PER_COMBO_DIR.glob('*__summary.parquet'))
    dompair_files = sorted(PER_COMBO_DIR.glob('*__dompairs.parquet'))
    kpair_files = sorted(PER_COMBO_DIR.glob('*__kpairs.parquet'))
    res_files = sorted(PER_COMBO_DIR.glob('*__residues.parquet'))
    print(f'Aggregating: {len(summary_files)} combos')
    pd.concat([pd.read_parquet(p) for p in summary_files], ignore_index=True) \
        .to_parquet(AGG_DIR / 'all_summary.parquet', index=False)
    pd.concat([pd.read_parquet(p) for p in dompair_files], ignore_index=True) \
        .to_parquet(AGG_DIR / 'all_dompairs.parquet', index=False)
    if kpair_files:
        pd.concat([pd.read_parquet(p) for p in kpair_files], ignore_index=True) \
            .to_parquet(AGG_DIR / 'all_kpairs.parquet', index=False)
    if res_files:
        # Residues table is large; write but consider per-domain split
        pd.concat([pd.read_parquet(p) for p in res_files], ignore_index=True) \
            .to_parquet(AGG_DIR / 'all_residues.parquet', index=False)
    # Add crystals
    xtal = CRYSTAL_DIR / 'crystal_summary.parquet'
    if xtal.exists():
        sa = pd.read_parquet(AGG_DIR / 'all_summary.parquet')
        cx = pd.read_parquet(xtal)
        pd.concat([sa, cx], ignore_index=True).to_parquet(AGG_DIR / 'all_summary.parquet', index=False)
    print('Aggregated tables under', AGG_DIR)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['setup', 'test', 'process', 'xtal', 'aggregate'])
    ap.add_argument('--combo', default=None)
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--only', nargs='*', default=None)
    args = ap.parse_args()

    if args.cmd == 'setup':
        cmd_setup()
    elif args.cmd == 'test':
        cmd_test(combo=args.combo)
    elif args.cmd == 'process':
        cmd_process(workers=args.workers, only=args.only)
    elif args.cmd == 'xtal':
        cmd_xtal()
    elif args.cmd == 'aggregate':
        cmd_aggregate()


if __name__ == '__main__':
    main()
