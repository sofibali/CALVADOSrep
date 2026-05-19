#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-residue SASA analysis and domain-face annotation for full-length PARP14.

For each residue:
  1. Compute SASA from the AF2 (or AF3) structure.
  2. Classify as buried (SASA <30 Å²) / partially exposed / surface (>100 Å²).
  3. For domains with active sites (MD1, MD2, MD3, ART), assign each residue
     to "active-site face" or "back face" based on the dot product between
     (residue_CA - domain_COM) and (active_site_COM - domain_COM).

Outputs:
    figures/sasa_per_residue_FL.png/svg            (1D track over FL sequence)
    figures/sasa_face_FL.png/svg                   (per-domain face annotation)
    data/sasa_face_per_residue.csv                 (per-residue annotations)
    figures/sasa_faces_FL.pse                      (PyMOL session if pymol available)

Usage:
    python figure_sasa_faces.py                    # use AF2 PDB
    python figure_sasa_faces.py --pdb input/parp14.pdb
    python figure_sasa_faces.py --no-pymol         # skip PyMOL session
"""

import os
import csv
import argparse
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import LinearSegmentedColormap

# ============================================================
# Configuration
# ============================================================

CWD = Path(__file__).resolve().parent
INPUT_DIR = CWD / 'input'
FIG_DIR = CWD / 'figures'
DATA_DIR = CWD / 'data'
FIG_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)

DEFAULT_PDB = INPUT_DIR / 'parp14.pdb'

# Domain boundaries (FL numbering)
DOMAINS = [
    ('RRM1',    1,    145, '#1f77b4'),
    ('RRM2',    146,  224, '#2ca02c'),
    ('RRM3',    225,  314, '#aec7e8'),
    ('KH1-KH6', 315,  737, '#ff7f0e'),
    ('KH7a',    738,  789, '#d62728'),
    ('MD1L1',   790,  1004, '#17becf'),
    ('MD2',     1005, 1193, '#9467bd'),
    ('MD3',     1207, 1388, '#8c564b'),
    ('KHb-KH8', 1389, 1533, '#e377c2'),
    ('WWE',     1534, 1602, '#7f7f7f'),
    ('ART',     1603, 1801, '#bcbd22'),
]

# Active sites (FL numbering, catalytic residues only)
ACTIVE_SITES = {
    'MD1L1': [831, 923, 962],
    'MD2':   [1035, 1046, 1134, 1171],
    'MD3':   [1248, 1259, 1330, 1371],
    'ART':   [1684, 1705, 1706, 1722],
}

# Pocket residues (more inclusive)
POCKET_RESIDUES = {
    'MD1L1': [822,823,824,825,826,827,828,829,830,831,832,833,834,835,836,
              919,920,921,922,923,924,925,926,927,961,962,966],
    'MD2':   [1021,1022,1023,1024,1034,1035,1036,1037,1038,1039,1040,1041,
              1042,1043,1044,1045,1046,1047,1130,1131,1132,1133,1134,1135,
              1136,1137,1138,1139,1140,1141,1170,1171,1175,1178],
    'MD3':   [1235,1236,1237,1247,1248,1249,1250,1251,1252,1253,1254,1255,
              1256,1257,1258,1259,1260,1261,1302,1303,1304,1324,1325,1326,
              1327,1328,1329,1330,1331,1332,1333,1334,1335,1336,1337,1369,
              1370,1371,1375],
    'ART':   [1681,1682,1683,1684,1685,1688,1701,1704,1705,1706,1707,1708,
              1709,1714,1715,1716,1721,1722,1726,1727,1781],
}

# SASA classification thresholds (Å²)
BURIED_MAX = 30
SURFACE_MIN = 80

# Chothia Gly-X-Gly extended-tripeptide max ASA (Å²), values from
# Wu et al. 2017 BioData Mining (https://doi.org/10.1186/s13040-016-0121-5)
# Used to compute RSA (relative solvent accessibility): RSA = SASA / max_ASA
# RSA ranges 0-1; standard threshold for "exposed" is RSA >= 0.20-0.25,
# "buried" is RSA < 0.05-0.10.
MAX_ASA_CHOTHIA = {
    'PHE': 210, 'ILE': 175, 'LEU': 170, 'VAL': 155, 'PRO': 145,
    'ALA': 115, 'GLY':  75, 'MET': 185, 'CYS': 135, 'TRP': 255,
    'TYR': 230, 'THR': 140, 'SER': 115, 'GLN': 180, 'ASN': 160,
    'GLU': 190, 'ASP': 150, 'HIS': 195, 'LYS': 200, 'ARG': 225,
}
RSA_BURIED_MAX = 0.05    # < 5% exposed = buried
RSA_PARTIAL_MAX = 0.20   # 5-20% = partial
# >20% = surface


# ============================================================
# Load structure & compute SASA
# ============================================================

def load_ca(pdb_path):
    """Return dict resid -> (resname, CA coord in Å)."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('p', str(pdb_path))
    ca = {}
    for r in structure[0].get_residues():
        if r.id[0] != ' ':
            continue
        cas = [a for a in r if a.name == 'CA']
        if cas:
            ca[r.id[1]] = (r.resname.strip(),
                           cas[0].get_vector().get_array())
    return ca


def compute_sasa(pdb_path):
    """Per-residue SASA in Å² using mdtraj."""
    import mdtraj as md
    traj = md.load(str(pdb_path))
    sasa = md.shrake_rupley(traj, mode='residue')[0]  # nm²
    return sasa * 100.0  # → Å²


# ============================================================
# Domain face assignment
# ============================================================

def classify_face(ca_coords_dict, domain_range, active_resids):
    """For residues in `domain_range`, classify as 'active' face or 'back'.

    Method: project residue offset (CA - domain_COM) onto unit vector from
    domain COM to active-site COM. Positive => active face.

    Returns: dict {resid: (dot_product, face)}
    """
    d_start, d_end = domain_range

    # Domain residues with coords
    d_coords = []
    d_resids = []
    for r in range(d_start, d_end + 1):
        if r in ca_coords_dict:
            d_resids.append(r)
            d_coords.append(ca_coords_dict[r][1])
    if not d_coords:
        return {}
    d_coords = np.array(d_coords)
    d_com = d_coords.mean(axis=0)

    # Active site COM
    act_coords = []
    for r in active_resids:
        if r in ca_coords_dict:
            act_coords.append(ca_coords_dict[r][1])
    if not act_coords:
        return {}
    act_com = np.array(act_coords).mean(axis=0)

    # Direction from domain COM to active site COM
    dir_vec = act_com - d_com
    dir_norm = np.linalg.norm(dir_vec)
    if dir_norm < 1e-6:
        return {}
    dir_unit = dir_vec / dir_norm

    # Per-residue projection
    result = {}
    for rid, coord in zip(d_resids, d_coords):
        proj = float(np.dot(coord - d_com, dir_unit))
        face = 'active' if proj > 0 else 'back'
        result[rid] = (proj, face)
    return result


def classify_sasa(sasa):
    """Categorize raw SASA (Å²) — kept for backward compat."""
    if sasa < BURIED_MAX:
        return 'buried'
    elif sasa < SURFACE_MIN:
        return 'partial'
    else:
        return 'surface'


def compute_rsa(sasa, resname):
    """Compute relative solvent accessibility (RSA) using Chothia Gly-X-Gly
    extended-tripeptide maximum SASA (Wu et al. 2017)."""
    max_asa = MAX_ASA_CHOTHIA.get(resname.upper())
    if not max_asa or max_asa == 0:
        return None
    return sasa / max_asa


def classify_rsa(rsa):
    """Categorize RSA into buried / partial / surface."""
    if rsa is None:
        return 'unknown'
    if rsa < RSA_BURIED_MAX:
        return 'buried'
    elif rsa < RSA_PARTIAL_MAX:
        return 'partial'
    else:
        return 'surface'


# ============================================================
# Main analysis
# ============================================================

def analyze(pdb_path):
    print(f"Loading structure: {pdb_path}")
    ca_dict = load_ca(pdb_path)
    print(f"  {len(ca_dict)} CA atoms loaded")

    print("Computing SASA...")
    sasa = compute_sasa(pdb_path)
    print(f"  Mean SASA: {sasa.mean():.1f} Å², "
          f"max: {sasa.max():.1f}, min: {sasa.min():.1f}")

    # Face assignment per active-site domain
    print("Assigning domain faces...")
    face_data = {}
    domain_lookup = {d[0]: (d[1], d[2]) for d in DOMAINS}

    for dname, sites in ACTIVE_SITES.items():
        d_range = domain_lookup[dname]
        faces = classify_face(ca_dict, d_range, sites)
        face_data[dname] = faces
        n_active = sum(1 for v in faces.values() if v[1] == 'active')
        n_back = sum(1 for v in faces.values() if v[1] == 'back')
        print(f"  {dname}: {n_active} active-face, {n_back} back-face residues")

    # Build per-residue annotation
    rows = []
    for rid in sorted(ca_dict.keys()):
        resname, coord = ca_dict[rid]
        s = sasa[rid - 1] if rid - 1 < len(sasa) else 0.0
        rsa = compute_rsa(s, resname)
        burial = classify_rsa(rsa)
        max_asa = MAX_ASA_CHOTHIA.get(resname.upper(), 0)

        # Domain
        domain_name = ''
        domain_color = '#cccccc'
        for dn, ds, de, dc in DOMAINS:
            if ds <= rid <= de:
                domain_name = dn
                domain_color = dc
                break

        # Active site?
        is_catalytic = ''
        for an, ar in ACTIVE_SITES.items():
            if rid in ar:
                is_catalytic = an
                break
        is_pocket = ''
        for an, ar in POCKET_RESIDUES.items():
            if rid in ar:
                is_pocket = an
                break

        # Face?
        face = ''
        face_proj = 0.0
        if domain_name in face_data and rid in face_data[domain_name]:
            face_proj, face = face_data[domain_name][rid]

        rows.append({
            'resid': rid,
            'resname': resname,
            'sasa_A2': round(s, 2),
            'max_asa_A2': max_asa,
            'rsa': round(rsa, 4) if rsa is not None else None,
            'burial': burial,
            'domain': domain_name,
            'face': face,
            'face_projection_A': round(face_proj, 2),
            'is_catalytic_in': is_catalytic,
            'is_pocket_in': is_pocket,
        })

    # Save CSV
    csv_path = DATA_DIR / 'sasa_face_per_residue.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved: {csv_path}")

    return ca_dict, sasa, face_data, rows


# ============================================================
# Plotting
# ============================================================

def plot_sasa_track(rows, outpath):
    """1D track of SASA along the sequence, colored by domain, with active
    sites and faces annotated."""
    fig, axes = plt.subplots(2, 1, figsize=(22, 8), sharex=True,
                              gridspec_kw={'height_ratios': [3, 1]})
    ax = axes[0]
    ax_face = axes[1]

    resids = np.array([r['resid'] for r in rows])
    sasas = np.array([r['sasa_A2'] for r in rows])

    # Domain background shading
    for dn, ds, de, dc in DOMAINS:
        ax.axvspan(ds, de, alpha=0.10, color=dc, zorder=0)
        ax_face.axvspan(ds, de, alpha=0.10, color=dc, zorder=0)
        ax.text((ds + de) / 2, ax.get_ylim()[1] * 0.95, dn,
                ha='center', va='top', fontsize=8, color=dc, fontweight='bold')

    # SASA bars colored by burial category
    colors = []
    for r in rows:
        if r['burial'] == 'buried':
            colors.append('#2c3e50')
        elif r['burial'] == 'partial':
            colors.append('#f39c12')
        else:
            colors.append('#e74c3c')

    ax.bar(resids, sasas, width=1.0, color=colors, edgecolor='none',
           alpha=0.7)

    # Mark catalytic residues
    for an, ar in ACTIVE_SITES.items():
        for rid in ar:
            ax.axvline(rid, color='red', linewidth=1.0, alpha=0.7, zorder=5)
            ax.plot([rid], [sasas[rid-1] if rid-1 < len(sasas) else 0],
                    marker='*', markersize=12, color='red', zorder=10,
                    markeredgecolor='black', markeredgewidth=0.5)

    ax.set_ylabel('SASA (Å²)', fontsize=12)
    ax.set_title('Per-Residue SASA along PARP14 Full-Length Sequence\n'
                 '(red stars = catalytic residues)', fontsize=13)
    ax.set_xlim(1, 1801)
    ax.set_ylim(0, max(sasas) * 1.05)

    # Bottom: face annotation track
    for r in rows:
        if r['face'] == 'active':
            color = '#27ae60'  # green for active face
        elif r['face'] == 'back':
            color = '#7f8c8d'  # gray for back
        else:
            color = '#f5f5f5'  # very light
        ax_face.bar(r['resid'], 1, width=1.0, color=color,
                    edgecolor='none', alpha=0.9)

    # Mark catalytic in face track
    for an, ar in ACTIVE_SITES.items():
        for rid in ar:
            ax_face.axvline(rid, color='red', linewidth=1.0, alpha=0.9,
                            zorder=5)

    ax_face.set_yticks([])
    ax_face.set_xlabel('Full-Length Residue Index', fontsize=12)
    ax_face.set_ylim(0, 1)
    ax_face.set_title('Domain face: green = active-site side, '
                       'gray = back side', fontsize=10)

    # Legend
    legend_elements = [
        Patch(facecolor='#2c3e50', label=f'Buried (SASA<{BURIED_MAX})'),
        Patch(facecolor='#f39c12', label=f'Partial ({BURIED_MAX}-{SURFACE_MIN})'),
        Patch(facecolor='#e74c3c', label=f'Surface (>{SURFACE_MIN})'),
        Patch(facecolor='#27ae60', label='Active face'),
        Patch(facecolor='#7f8c8d', label='Back face'),
        plt.Line2D([0], [0], color='red', linewidth=2, label='Catalytic'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=2)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.png")


def plot_face_per_domain(rows, outpath):
    """Per-domain bar chart: number of buried/exposed residues split by face."""
    # Active-site domains only
    domains_with_sites = list(ACTIVE_SITES.keys())

    categories = ['Active+Surface', 'Active+Partial', 'Active+Buried',
                  'Back+Surface', 'Back+Partial', 'Back+Buried']
    cat_colors = ['#27ae60', '#82c89f', '#1e6e3a',
                  '#7f8c8d', '#bdc3c7', '#34495e']

    data = {d: [0] * 6 for d in domains_with_sites}
    for r in rows:
        if not r['face'] or r['domain'] not in domains_with_sites:
            continue
        if r['face'] == 'active':
            if r['burial'] == 'surface':
                idx = 0
            elif r['burial'] == 'partial':
                idx = 1
            else:
                idx = 2
        else:
            if r['burial'] == 'surface':
                idx = 3
            elif r['burial'] == 'partial':
                idx = 4
            else:
                idx = 5
        data[r['domain']][idx] += 1

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(domains_with_sites))
    bottom = np.zeros(len(domains_with_sites))
    for i, (cat, col) in enumerate(zip(categories, cat_colors)):
        vals = [data[d][i] for d in domains_with_sites]
        ax.bar(x, vals, bottom=bottom, label=cat, color=col,
               edgecolor='black', linewidth=0.3)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(domains_with_sites, fontsize=11)
    ax.set_ylabel('Number of residues', fontsize=12)
    ax.set_title('Domain face composition: active-site face vs back face\n'
                 f'(RSA-based burial: buried<{RSA_BURIED_MAX}, '
                 f'partial<{RSA_PARTIAL_MAX}, surface>{RSA_PARTIAL_MAX})',
                 fontsize=12)
    ax.legend(fontsize=9, loc='upper right', ncol=2)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.png")


def plot_rsa_heatmap(rows, outpath):
    """RSA heatmap: full-length sequence as a 1D color strip showing the
    relative solvent accessibility (RSA = SASA / Chothia max).

    Top track:    RSA value (0-1, white→red gradient = buried→exposed)
    Middle track: domain colored band
    Bottom track: active-site face annotation
    """
    from matplotlib.colors import LinearSegmentedColormap

    fig, axes = plt.subplots(4, 1, figsize=(22, 7), sharex=True,
                              gridspec_kw={'height_ratios': [1, 1, 0.5, 0.5]})
    ax_rsa = axes[0]    # RSA heatmap strip
    ax_burial = axes[1] # burial category strip
    ax_dom = axes[2]    # domain colors
    ax_face = axes[3]   # active-site face

    n_fl = 1801
    rsa_arr = np.full(n_fl, np.nan)
    burial_arr = np.zeros(n_fl, dtype=int)  # 0=unknown, 1=buried, 2=partial, 3=surface
    for r in rows:
        rid = r['resid']
        if 1 <= rid <= n_fl and r['rsa'] is not None:
            rsa_arr[rid - 1] = r['rsa']
        b = r['burial']
        if 1 <= rid <= n_fl:
            burial_arr[rid - 1] = (
                1 if b == 'buried' else
                2 if b == 'partial' else
                3 if b == 'surface' else 0)

    # --- RSA strip (continuous colormap) ---
    rsa_cmap = LinearSegmentedColormap.from_list(
        'rsa', ['#08306b', '#deebf7', '#fff5f0', '#fb6a4a', '#67000d'])
    im = ax_rsa.imshow(rsa_arr[np.newaxis, :], aspect='auto', cmap=rsa_cmap,
                       vmin=0, vmax=1.0,
                       extent=[1, n_fl, 0, 1], interpolation='nearest')
    cbar = plt.colorbar(im, ax=ax_rsa, orientation='horizontal',
                        pad=0.05, fraction=0.08, aspect=60)
    cbar.set_label('RSA (SASA / Chothia Gly-X-Gly max ASA)', fontsize=10)
    ax_rsa.set_yticks([])
    ax_rsa.set_title('Per-Residue Relative Solvent Accessibility (RSA)\n'
                     'Buried (RSA<0.05) → Partial (0.05-0.20) → '
                     'Surface (>0.20), Chothia max ASA normalization',
                     fontsize=12)

    # Mark catalytic residues on RSA strip
    for an, ar in ACTIVE_SITES.items():
        for rid in ar:
            ax_rsa.plot([rid], [0.5], marker='v', markersize=10,
                        color='lime', markeredgecolor='black',
                        markeredgewidth=0.5, zorder=10)

    # --- Burial category strip (3 discrete colors) ---
    burial_cmap = LinearSegmentedColormap.from_list(
        'burial', ['#ffffff', '#2c3e50', '#f39c12', '#e74c3c'], N=4)
    ax_burial.imshow(burial_arr[np.newaxis, :], aspect='auto',
                     cmap=burial_cmap, vmin=0, vmax=3,
                     extent=[1, n_fl, 0, 1], interpolation='nearest')
    ax_burial.set_yticks([])
    ax_burial.set_title(
        f'Burial category: buried (RSA<{RSA_BURIED_MAX}), '
        f'partial ({RSA_BURIED_MAX}-{RSA_PARTIAL_MAX}), '
        f'surface (>{RSA_PARTIAL_MAX})', fontsize=10)

    # --- Domain colors strip ---
    domain_strip = np.full((1, n_fl, 3), 1.0)  # white default
    import matplotlib.colors as mcolors
    for dn, ds, de, dc in DOMAINS:
        rgb = mcolors.to_rgb(dc)
        for r in range(ds - 1, min(de, n_fl)):
            domain_strip[0, r] = rgb
    ax_dom.imshow(domain_strip, aspect='auto',
                  extent=[1, n_fl, 0, 1], interpolation='nearest')
    ax_dom.set_yticks([])
    ax_dom.set_title('Domain', fontsize=10)
    # Domain name labels
    for dn, ds, de, dc in DOMAINS:
        ax_dom.text((ds + de) / 2, 0.5, dn, ha='center', va='center',
                    fontsize=7, color='white', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.1',
                              facecolor=dc, edgecolor='none', alpha=0.8))

    # --- Active-site face strip ---
    face_strip = np.full((1, n_fl, 3), 0.96)  # near white
    for r in rows:
        rid = r['resid']
        if not (1 <= rid <= n_fl):
            continue
        if r['face'] == 'active':
            face_strip[0, rid - 1] = mcolors.to_rgb('#27ae60')
        elif r['face'] == 'back':
            face_strip[0, rid - 1] = mcolors.to_rgb('#7f8c8d')
    ax_face.imshow(face_strip, aspect='auto',
                   extent=[1, n_fl, 0, 1], interpolation='nearest')
    ax_face.set_yticks([])
    ax_face.set_xlabel('Full-Length Residue Index', fontsize=12)
    ax_face.set_title('Active-site face (green = active side, '
                       'gray = back side)', fontsize=10)

    # Mark catalytic on face strip too
    for an, ar in ACTIVE_SITES.items():
        for rid in ar:
            ax_face.axvline(rid, color='red', linewidth=1, alpha=0.9)

    for ax in axes:
        ax.set_xlim(1, n_fl)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.png")


def plot_rsa_per_domain_matrix(rows, outpath):
    """Per-domain RSA matrix heatmap: each row = a domain, each cell = a
    residue's RSA value, sorted within domain by position."""
    fig, ax = plt.subplots(figsize=(20, 6))

    # Group by domain and position
    domain_residues = {dn: [] for dn, _, _, _ in DOMAINS}
    for r in rows:
        if r['domain'] in domain_residues and r['rsa'] is not None:
            domain_residues[r['domain']].append((r['resid'], r['rsa']))

    # Build a matrix where each row is a domain, padded with NaN
    max_len = max(len(v) for v in domain_residues.values())
    matrix = np.full((len(DOMAINS), max_len), np.nan)
    labels = []
    for i, (dn, _, _, _) in enumerate(DOMAINS):
        sorted_pairs = sorted(domain_residues[dn])
        for j, (_, rsa) in enumerate(sorted_pairs):
            matrix[i, j] = rsa
        labels.append(f'{dn} (n={len(sorted_pairs)})')

    from matplotlib.colors import LinearSegmentedColormap
    rsa_cmap = LinearSegmentedColormap.from_list(
        'rsa', ['#08306b', '#deebf7', '#fff5f0', '#fb6a4a', '#67000d'])
    im = ax.imshow(matrix, aspect='auto', cmap=rsa_cmap, vmin=0, vmax=1.0,
                   interpolation='nearest')

    ax.set_yticks(range(len(DOMAINS)))
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('Residue position within domain (0-indexed)', fontsize=11)
    ax.set_title('Per-Domain RSA Heatmap\n'
                 '(rows = domains, columns = residue position from N to C)',
                 fontsize=12)
    cbar = plt.colorbar(im, ax=ax, label='RSA', shrink=0.8)

    # Mark catalytic residue columns
    for an, ar in ACTIVE_SITES.items():
        # Find domain for this active site
        d_idx = None
        d_start = None
        for i, (dn, ds, de, _) in enumerate(DOMAINS):
            if dn == an or (an == 'MD1' and dn == 'MD1L1'):
                d_idx = i
                d_start = ds
                break
            # fallback: catalytic residue falls within domain range
            if any(ds <= r <= de for r in ar) and d_idx is None:
                d_idx = i
                d_start = ds
        if d_idx is None or d_start is None:
            continue
        sorted_resids = sorted(rid for rid, _ in domain_residues[DOMAINS[d_idx][0]])
        for rid in ar:
            if rid in sorted_resids:
                col = sorted_resids.index(rid)
                ax.scatter([col], [d_idx], marker='*', s=80, c='black',
                           edgecolors='lime', linewidths=1, zorder=10)

    plt.tight_layout()
    fig.savefig(f'{outpath}.png', dpi=200, bbox_inches='tight')
    fig.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.png")


# ============================================================
# PyMOL session
# ============================================================

def write_pymol_script(rows, pdb_path, pml_path, pse_path):
    """Write a PyMOL .pml script that colors residues by SASA and creates
    objects per face."""
    with open(pml_path, 'w') as f:
        f.write(f"# PARP14 SASA + Face annotation\n")
        f.write(f"# Run: pymol {pml_path}\n\n")
        f.write(f"load {pdb_path}, parp14\n")
        f.write(f"hide everything\n")
        f.write(f"show cartoon\n")
        f.write(f"bg_color white\n")
        f.write(f"set cartoon_transparency, 0.3\n\n")

        # Color by SASA gradient (B-factor trick: write SASA as b-factor)
        # Create a duplicate object for B-factor coloring
        f.write("# --- Scene 1: SASA spectrum ---\n")
        f.write("copy parp14_sasa, parp14\n")
        # Assign SASA to alter b-factors
        f.write("# Set alter all to SASA values\n")
        for r in rows:
            if r['sasa_A2'] > 0:
                f.write(f"alter parp14_sasa and resi {r['resid']}, "
                        f"b={r['sasa_A2']:.2f}\n")
        f.write("spectrum b, blue_white_red, parp14_sasa, 0, 150\n\n")

        # Create face objects per active-site domain
        f.write("# --- Scene 2: Active-site faces per domain ---\n")
        domains_with_sites = list(ACTIVE_SITES.keys())
        for dname in domains_with_sites:
            active_resids = [str(r['resid']) for r in rows
                             if r['domain'] == dname and r['face'] == 'active']
            back_resids = [str(r['resid']) for r in rows
                           if r['domain'] == dname and r['face'] == 'back']
            if active_resids:
                f.write(f"select sel_{dname}_active, parp14 and resi "
                        f"{'+'.join(active_resids)}\n")
                f.write(f"create {dname}_active_face, sel_{dname}_active\n")
                f.write(f"color limegreen, {dname}_active_face\n")
            if back_resids:
                f.write(f"select sel_{dname}_back, parp14 and resi "
                        f"{'+'.join(back_resids)}\n")
                f.write(f"create {dname}_back_face, sel_{dname}_back\n")
                f.write(f"color gray60, {dname}_back_face\n")
            f.write("\n")

        # Catalytic residues
        f.write("# --- Catalytic residues (red spheres) ---\n")
        for an, ar in ACTIVE_SITES.items():
            f.write(f"select cat_{an}, parp14 and resi "
                    f"{'+'.join(map(str, ar))}\n")
            f.write(f"show spheres, cat_{an}\n")
            f.write(f"color red, cat_{an}\n")
            f.write(f"set sphere_scale, 0.6, cat_{an}\n")

        # Buried residues object
        f.write("\n# --- Buried residues ---\n")
        buried = [str(r['resid']) for r in rows if r['burial'] == 'buried']
        if buried:
            f.write(f"select buried_residues, parp14 and resi "
                    f"{'+'.join(buried)}\n")
            f.write(f"create buried, buried_residues\n")
            f.write(f"color slate, buried\n\n")

        # Surface (highly exposed) residues
        surface = [str(r['resid']) for r in rows if r['burial'] == 'surface']
        if surface:
            f.write(f"select surface_residues, parp14 and resi "
                    f"{'+'.join(surface)}\n")
            f.write(f"create surface, surface_residues\n")
            f.write(f"color firebrick, surface\n\n")

        # Scenes
        f.write("# --- Scenes ---\n")
        f.write("disable all\n")
        f.write("enable parp14_sasa\n")
        f.write("orient parp14_sasa\n")
        f.write("scene SASA_spectrum, store\n\n")

        for dname in domains_with_sites:
            f.write(f"disable all\n")
            f.write(f"enable {dname}_active_face\n")
            f.write(f"enable {dname}_back_face\n")
            f.write(f"enable cat_{dname}\n")
            f.write(f"orient {dname}_active_face\n")
            f.write(f"scene {dname}_faces, store\n\n")

        f.write("disable all\n")
        f.write("enable buried\n")
        f.write("enable surface\n")
        f.write("orient parp14\n")
        f.write("scene Buried_vs_Surface, store\n\n")

        # Save session
        f.write("scene SASA_spectrum, recall\n")
        f.write(f"save {pse_path}\n")


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', default=str(DEFAULT_PDB),
                        help='Input PDB (default: AF2 FL)')
    parser.add_argument('--no-pymol', action='store_true',
                        help='Skip PyMOL session generation')
    args = parser.parse_args()

    pdb_path = Path(args.pdb).resolve()
    if not pdb_path.exists():
        print(f"ERROR: PDB not found: {pdb_path}")
        return 1

    ca_dict, sasa, face_data, rows = analyze(pdb_path)

    # Plot
    print("\nGenerating plots...")
    plot_sasa_track(rows, FIG_DIR / 'sasa_per_residue_FL')
    plot_face_per_domain(rows, FIG_DIR / 'sasa_face_per_domain_FL')
    plot_rsa_heatmap(rows, FIG_DIR / 'rsa_heatmap_FL')
    plot_rsa_per_domain_matrix(rows, FIG_DIR / 'rsa_per_domain_matrix_FL')

    # PyMOL session
    if not args.no_pymol:
        print("\nWriting PyMOL session script...")
        pml_path = FIG_DIR / 'sasa_faces_FL.pml'
        pse_path = FIG_DIR / 'sasa_faces_FL.pse'
        write_pymol_script(rows, str(pdb_path),
                           str(pml_path), str(pse_path))
        print(f"  PyMOL script: {pml_path}")
        print(f"  To render: pymol -cq {pml_path}")

    # Summary stats
    print("\n=== Summary ===")
    n_buried = sum(1 for r in rows if r['burial'] == 'buried')
    n_partial = sum(1 for r in rows if r['burial'] == 'partial')
    n_surface = sum(1 for r in rows if r['burial'] == 'surface')
    print(f"  Buried  (<{BURIED_MAX} Å²): {n_buried}")
    print(f"  Partial:                    {n_partial}")
    print(f"  Surface (>{SURFACE_MIN} Å²): {n_surface}")

    for dname in ACTIVE_SITES:
        a_res = [r for r in rows if r['domain'] == dname and r['face'] == 'active']
        b_res = [r for r in rows if r['domain'] == dname and r['face'] == 'back']
        if not a_res or not b_res:
            continue
        a_mean = np.mean([r['sasa_A2'] for r in a_res])
        b_mean = np.mean([r['sasa_A2'] for r in b_res])
        print(f"  {dname}: active-face mean SASA={a_mean:.1f}, "
              f"back-face mean SASA={b_mean:.1f} "
              f"(active is {'more' if a_mean > b_mean else 'less'} exposed)")


if __name__ == '__main__':
    main()
