#!/usr/bin/env python3
"""
Analyze the longest AF3 construct vs AF2 full-length PARP14.

4-panel figure:
  A. PAE contact map (AF3 construct) with atS domain boundaries
  B. Per-residue pLDDT: AF3 construct vs AF2 full-length (full-length x-axis)
  C. Active site SASA: AF3 vs AF2 side-by-side
  D. Active site pLDDT comparison: AF3 vs AF2

Construct: rrm1_rrm2_rrm3_kh1-kh6_kh7a_md2_md3_khb-kh8_wwe_art
(1574 residues, 9 atS domains — longest available AF3 prediction)

Usage:
    python 03_longest_construct_figure.py
    python 03_longest_construct_figure.py --seed 1 --sample 0
    python 03_longest_construct_figure.py --all-models
"""

import json
import argparse
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# ============================================================
# Constants
# ============================================================

CONSTRUCT = 'rrm1_rrm2_rrm3_kh1-kh6_kh7a_md2_md3_khb-kh8_wwe_art'
AF3_BASE = Path(__file__).resolve().parent.parent / 'alphafold_outputs'
AF2_PDB = Path(__file__).resolve().parent.parent.parent / 'examples' / 'PARP14_MDP' / 'input' / 'parp14.pdb'
OUTPUT_DIR = Path(__file__).resolve().parent / 'figures'

# --- AF3 construct domain definitions (what the construct was actually built with) ---
CONSTRUCT_DOMAIN_DEFS = {
    'rrm1':    (1, 145),
    'rrm2':    (146, 224),
    'rrm3':    (225, 314),
    'kh1-kh6': (315, 737),
    'kh7a':    (738, 789),
    'md2':     (1004, 1193),
    'md3':     (1207, 1388),
    'khb-kh8': (1389, 1533),
    'wwe':     (1534, 1602),
    'art':     (1603, 1801),
}
CONSTRUCT_DOMAIN_ORDER = ['rrm1', 'rrm2', 'rrm3', 'kh1-kh6', 'kh7a',
                          'md2', 'md3', 'khb-kh8', 'wwe', 'art']

# --- Domain annotations from PARP14_domains_atS.fasta ---
# Full-length boundaries cut at serine residues (1-indexed, inclusive)
DOMAIN_DEFS_ATS = {
    'RRM1':   (1, 148),
    'RRM2-3': (149, 304),
    'KH_N':   (305, 516),
    'KH_Ca':  (517, 789),
    'MD1':    (790, 1006),
    'MD2':    (1007, 1188),
    'MD3':    (1189, 1398),
    'KH_b':   (1399, 1508),
    'WWE':    (1509, 1609),
    'ART':    (1610, 1810),
}

# All atS domains (for full-length annotation)
ALL_DOMAIN_ORDER = ['RRM1', 'RRM2-3', 'KH_N', 'KH_Ca', 'MD1',
                    'MD2', 'MD3', 'KH_b', 'WWE', 'ART']

# Domains present in the AF3 construct (MD1 is absent)
CONSTRUCT_ATS_ORDER = ['RRM1', 'RRM2-3', 'KH_N', 'KH_Ca',
                       'MD2', 'MD3', 'KH_b', 'WWE', 'ART']

# Active site residues (full-length numbering)
ACTIVE_SITES = {
    'MD1': [831, 923, 962],
    'MD2': [1035, 1046, 1134, 1171],
    'ART': [1684, 1705, 1706, 1722],
}

# Domain colors (all 10 atS domains)
DOMAIN_COLORS = {
    'RRM1':   '#1f77b4',
    'RRM2-3': '#2ca02c',
    'KH_N':   '#ff7f0e',
    'KH_Ca':  '#d62728',
    'MD1':    '#aec7e8',
    'MD2':    '#9467bd',
    'MD3':    '#8c564b',
    'KH_b':   '#e377c2',
    'WWE':    '#7f7f7f',
    'ART':    '#bcbd22',
}


# ============================================================
# Domain boundary computation
# ============================================================

def compute_construct_boundaries() -> Tuple[Dict[str, Tuple[int, int]], Dict[int, int], int]:
    """Compute 0-indexed atS domain boundaries on the AF3 construct.

    Returns: (boundaries, fl_to_construct, total_residues)
    """
    construct_to_fl = {}
    fl_to_construct = {}
    offset = 0
    prev_end_fl = None

    for dname in CONSTRUCT_DOMAIN_ORDER:
        fl_start, fl_end = CONSTRUCT_DOMAIN_DEFS[dname]
        n_residues = fl_end - fl_start + 1
        overlap = 0
        if prev_end_fl is not None and fl_start <= prev_end_fl:
            overlap = prev_end_fl - fl_start + 1
            n_residues -= overlap
        for i in range(n_residues):
            fl_res = fl_start + overlap + i
            construct_to_fl[offset + i] = fl_res
            fl_to_construct[fl_res] = offset + i
        offset += n_residues
        prev_end_fl = fl_end

    total_residues = offset

    boundaries = {}
    for dname in CONSTRUCT_ATS_ORDER:
        fl_start, fl_end = DOMAIN_DEFS_ATS[dname]
        positions = [c for c, fl in construct_to_fl.items()
                     if fl_start <= fl <= fl_end]
        if positions:
            boundaries[dname] = (min(positions), max(positions) + 1)

    return boundaries, fl_to_construct, total_residues


# ============================================================
# Data extraction
# ============================================================

def extract_plddt(conf: Dict) -> np.ndarray:
    """Extract per-residue pLDDT from AF3 confidences.json."""
    atom_plddts = np.array(conf['atom_plddts'])
    n_tokens = len(conf['token_chain_ids'])
    n_atoms = len(atom_plddts)
    atoms_per_token = n_atoms / n_tokens
    per_res = []
    for i in range(n_tokens):
        s = int(round(i * atoms_per_token))
        e = int(round((i + 1) * atoms_per_token))
        if s < n_atoms:
            per_res.append(float(np.mean(atom_plddts[s:e])))
    return np.array(per_res)


def compute_sasa_cif(cif_path: str) -> np.ndarray:
    """Compute per-residue SASA from CIF structure."""
    from Bio.PDB import MMCIFParser, ShrakeRupley
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure('model', cif_path)
    sr = ShrakeRupley()
    sr.compute(struct, level='R')
    return np.array([res.sasa for res in struct.get_residues()])


def compute_sasa_pdb(pdb_path: str) -> np.ndarray:
    """Compute per-residue SASA from PDB structure."""
    from Bio.PDB import PDBParser, ShrakeRupley
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('model', pdb_path)
    sr = ShrakeRupley()
    sr.compute(struct, level='R')
    return np.array([res.sasa for res in struct.get_residues()])


def load_af2_data() -> Dict:
    """Load AF2 full-length structure: pLDDT from B-factors, SASA from structure."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('af2', str(AF2_PDB))

    residues = list(struct.get_residues())
    plddt = []
    resids = []
    for res in residues:
        ca = res['CA'] if 'CA' in res else list(res.get_atoms())[0]
        plddt.append(ca.get_bfactor())
        resids.append(res.get_id()[1])

    print(f'  AF2: {len(residues)} residues, resids {resids[0]}-{resids[-1]}')
    sasa = compute_sasa_pdb(str(AF2_PDB))

    return {
        'plddt': np.array(plddt),
        'sasa': sasa,
        'resids': np.array(resids),
        'n_residues': len(residues),
    }


def load_model_data(seed: int, sample: int) -> Dict:
    """Load PAE, pLDDT, and SASA for one AF3 seed-sample model."""
    model_dir = AF3_BASE / CONSTRUCT / f'seed-{seed}_sample-{sample}'
    with open(model_dir / 'confidences.json') as f:
        conf = json.load(f)
    pae = np.array(conf['pae'])
    plddt = extract_plddt(conf)
    sasa = compute_sasa_cif(str(model_dir / 'model.cif'))
    with open(model_dir / 'summary_confidences.json') as f:
        summary = json.load(f)
    return {
        'pae': pae, 'plddt': plddt, 'sasa': sasa,
        'ptm': summary.get('ptm'),
        'ranking_score': summary.get('ranking_score'),
    }


def load_all_models() -> Dict:
    """Load and average across all 25 AF3 seed-sample models."""
    pae_list, plddt_list, sasa_list = [], [], []
    ptm_list, rs_list = [], []
    for seed in range(1, 6):
        for sample in range(5):
            model_dir = AF3_BASE / CONSTRUCT / f'seed-{seed}_sample-{sample}'
            if not model_dir.exists():
                continue
            print(f'  Loading seed-{seed}_sample-{sample}...')
            data = load_model_data(seed, sample)
            pae_list.append(data['pae'])
            plddt_list.append(data['plddt'])
            sasa_list.append(data['sasa'])
            if data['ptm'] is not None:
                ptm_list.append(data['ptm'])
            if data['ranking_score'] is not None:
                rs_list.append(data['ranking_score'])
    return {
        'pae': np.mean(pae_list, axis=0),
        'pae_sd': np.std(pae_list, axis=0),
        'plddt': np.mean(plddt_list, axis=0),
        'plddt_sd': np.std(plddt_list, axis=0),
        'sasa': np.mean(sasa_list, axis=0),
        'sasa_sd': np.std(sasa_list, axis=0),
        'ptm_mean': np.mean(ptm_list) if ptm_list else None,
        'ranking_score_mean': np.mean(rs_list) if rs_list else None,
        'n_models': len(pae_list),
    }


# ============================================================
# Plotting
# ============================================================

def make_figure(af3_data: Dict, af2_data: Dict, boundaries: Dict,
                fl_to_construct: Dict, n_residues: int,
                title_suffix: str = '', output_path: Path = None):
    """Generate the 4-panel figure comparing AF3 construct with AF2 full-length."""

    pae = af3_data['pae']
    af3_plddt = af3_data['plddt']
    af3_sasa = af3_data['sasa']
    af2_plddt = af2_data['plddt']
    af2_sasa = af2_data['sasa']
    af2_resids = af2_data['resids']

    fig = plt.figure(figsize=(18, 16))
    gs = GridSpec(3, 2, height_ratios=[1.0, 0.8, 0.8],
                  hspace=0.35, wspace=0.3)

    # ---- Panel A: PAE contact map (AF3 construct coordinates) ----
    ax_pae = fig.add_subplot(gs[0, 0])
    im = ax_pae.imshow(pae, cmap='Greens_r', vmin=0, vmax=30, aspect='equal',
                       origin='upper', interpolation='none')

    for dname, (s, e) in boundaries.items():
        for pos in [s, e]:
            ax_pae.axhline(pos, color='black', linewidth=0.5, alpha=0.6)
            ax_pae.axvline(pos, color='black', linewidth=0.5, alpha=0.6)

    for dname, (s, e) in boundaries.items():
        mid = (s + e) / 2
        ax_pae.text(mid, mid, dname, ha='center', va='center',
                    fontsize=6, fontweight='bold', color='white',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.5))

    for domain, residues in ACTIVE_SITES.items():
        for fl_res in residues:
            if fl_res in fl_to_construct:
                pos = fl_to_construct[fl_res]
                ax_pae.plot(pos, -15, 'v', color='red', markersize=4, clip_on=False)
                ax_pae.plot(-15, pos, '>', color='red', markersize=4, clip_on=False)

    cb = plt.colorbar(im, ax=ax_pae, shrink=0.8, pad=0.02)
    cb.set_label('PAE (\u00c5)', fontsize=10)
    ax_pae.set_xlabel('Residue (construct)', fontsize=10)
    ax_pae.set_ylabel('Residue (construct)', fontsize=10)
    ax_pae.set_title('A. AF3 Predicted Aligned Error (PAE)', fontsize=12, fontweight='bold')
    ax_pae.set_xlim(-0.5, n_residues - 0.5)
    ax_pae.set_ylim(n_residues - 0.5, -0.5)

    # ---- Panel B: Active site SASA comparison (AF3 vs AF2) ----
    ax_sasa = fig.add_subplot(gs[0, 1])

    active_labels = []
    af3_sasa_vals = []
    af3_sasa_sds = []
    af2_sasa_vals = []
    bar_domain_colors = []

    for domain, residues in ACTIVE_SITES.items():
        for fl_res in residues:
            active_labels.append(f'{domain}\n{fl_res}')
            # AF2 (always has full-length)
            af2_idx = fl_res - 1  # 0-indexed (AF2 resids start at 1)
            if af2_idx < len(af2_sasa):
                af2_sasa_vals.append(af2_sasa[af2_idx])
            else:
                af2_sasa_vals.append(0)
            # AF3 (only if in construct)
            if fl_res in fl_to_construct:
                pos = fl_to_construct[fl_res]
                af3_sasa_vals.append(af3_sasa[pos])
                if 'sasa_sd' in af3_data:
                    af3_sasa_sds.append(af3_data['sasa_sd'][pos])
                else:
                    af3_sasa_sds.append(0)
            else:
                af3_sasa_vals.append(np.nan)
                af3_sasa_sds.append(0)

            bar_domain_colors.append(DOMAIN_COLORS.get(domain, '#999999'))

    bar_x = np.arange(len(active_labels))
    width = 0.35

    # AF2 bars
    ax_sasa.bar(bar_x - width/2, af2_sasa_vals, width,
                color=[c for c in bar_domain_colors], alpha=0.4,
                edgecolor='black', linewidth=0.5, label='AF2 full-length')

    # AF3 bars
    af3_yerr = af3_sasa_sds if any(s > 0 for s in af3_sasa_sds) else None
    ax_sasa.bar(bar_x + width/2, af3_sasa_vals, width, yerr=af3_yerr,
                color=[c for c in bar_domain_colors], alpha=0.8,
                edgecolor='black', linewidth=0.5, capsize=2, label='AF3 construct')

    ax_sasa.set_xticks(bar_x)
    ax_sasa.set_xticklabels(active_labels, fontsize=7)
    ax_sasa.set_ylabel('SASA (\u00c5\u00b2)', fontsize=10)
    ax_sasa.set_title('B. Active Site SASA: AF2 vs AF3', fontsize=12, fontweight='bold')
    ax_sasa.axhline(50, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)
    ax_sasa.text(len(active_labels) - 0.5, 52, 'buried threshold', fontsize=7,
                 ha='right', color='gray')
    ax_sasa.legend(fontsize=8, loc='upper right')

    # ---- Panel C: Active site pLDDT comparison ----
    ax_plddt_bar = fig.add_subplot(gs[1, 1])

    af2_plddt_vals = []
    af3_plddt_vals = []
    af3_plddt_sds = []

    for domain, residues in ACTIVE_SITES.items():
        for fl_res in residues:
            af2_idx = fl_res - 1
            if af2_idx < len(af2_plddt):
                af2_plddt_vals.append(af2_plddt[af2_idx])
            else:
                af2_plddt_vals.append(0)
            if fl_res in fl_to_construct:
                pos = fl_to_construct[fl_res]
                af3_plddt_vals.append(af3_plddt[pos])
                if 'plddt_sd' in af3_data:
                    af3_plddt_sds.append(af3_data['plddt_sd'][pos])
                else:
                    af3_plddt_sds.append(0)
            else:
                af3_plddt_vals.append(np.nan)
                af3_plddt_sds.append(0)

    ax_plddt_bar.bar(bar_x - width/2, af2_plddt_vals, width,
                     color=[c for c in bar_domain_colors], alpha=0.4,
                     edgecolor='black', linewidth=0.5, label='AF2')
    af3_plddt_yerr = af3_plddt_sds if any(s > 0 for s in af3_plddt_sds) else None
    ax_plddt_bar.bar(bar_x + width/2, af3_plddt_vals, width, yerr=af3_plddt_yerr,
                     color=[c for c in bar_domain_colors], alpha=0.8,
                     edgecolor='black', linewidth=0.5, capsize=2, label='AF3')

    # Confidence thresholds
    for threshold, label, color in [(90, 'Very high', '#0053d6'),
                                     (70, 'Confident', '#65cbf3'),
                                     (50, 'Low', '#ffdb13')]:
        ax_plddt_bar.axhline(threshold, color=color, linewidth=0.5,
                              linestyle='--', alpha=0.5)

    ax_plddt_bar.set_xticks(bar_x)
    ax_plddt_bar.set_xticklabels(active_labels, fontsize=7)
    ax_plddt_bar.set_ylabel('pLDDT', fontsize=10)
    ax_plddt_bar.set_ylim(0, 105)
    ax_plddt_bar.set_title('C. Active Site pLDDT: AF2 vs AF3', fontsize=12, fontweight='bold')
    ax_plddt_bar.legend(fontsize=8, loc='lower right')

    # ---- Panel D: Full-length pLDDT trace (AF2 full + AF3 overlay) ----
    ax_plddt = fig.add_subplot(gs[1:, 0])

    fl_x = np.arange(1, len(af2_plddt) + 1)

    # AF2 trace (full-length, light gray background)
    ax_plddt.fill_between(fl_x, af2_plddt, alpha=0.15, color='gray')
    ax_plddt.plot(fl_x, af2_plddt, color='gray', linewidth=0.6, alpha=0.7, label='AF2 full-length')

    # AF3 construct trace, plotted at FL positions, colored by atS domain
    # Build AF3 pLDDT mapped to FL positions
    construct_to_fl_map = {}
    offset = 0
    prev_end_fl = None
    for dname in CONSTRUCT_DOMAIN_ORDER:
        fl_start, fl_end = CONSTRUCT_DOMAIN_DEFS[dname]
        n = fl_end - fl_start + 1
        overlap = 0
        if prev_end_fl is not None and fl_start <= prev_end_fl:
            overlap = prev_end_fl - fl_start + 1
            n -= overlap
        for i in range(n):
            construct_to_fl_map[offset + i] = fl_start + overlap + i
        offset += n
        prev_end_fl = fl_end

    # Plot AF3 colored by atS domain
    for dname in CONSTRUCT_ATS_ORDER:
        fl_s, fl_e = DOMAIN_DEFS_ATS[dname]
        xs = []
        ys = []
        for cpos, fl_pos in sorted(construct_to_fl_map.items()):
            if fl_s <= fl_pos <= fl_e:
                xs.append(fl_pos)
                ys.append(af3_plddt[cpos])
        if xs:
            color = DOMAIN_COLORS[dname]
            ax_plddt.fill_between(xs, ys, alpha=0.25, color=color)
            ax_plddt.plot(xs, ys, color=color, linewidth=1.0)

    # Domain region shading (all 10 atS domains on full-length axis)
    for dname in ALL_DOMAIN_ORDER:
        fl_s, fl_e = DOMAIN_DEFS_ATS[dname]
        color = DOMAIN_COLORS[dname]
        ax_plddt.axvspan(fl_s, fl_e, alpha=0.04, color=color)
        mid = (fl_s + fl_e) / 2
        ax_plddt.text(mid, 3, dname, ha='center', va='bottom', fontsize=6,
                      color=color, fontweight='bold', rotation=90 if (fl_e - fl_s) < 100 else 0)

    # Mark missing region (MD1)
    ax_plddt.axvspan(790, 1006, alpha=0.1, color='red', hatch='///')
    ax_plddt.text(898, 50, 'MD1\n(missing)', ha='center', va='center',
                  fontsize=8, color='red', fontweight='bold',
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # pLDDT thresholds
    for threshold, label, color in [(90, 'Very high', '#0053d6'),
                                     (70, 'Confident', '#65cbf3'),
                                     (50, 'Low', '#ffdb13')]:
        ax_plddt.axhline(threshold, color=color, linewidth=0.5, linestyle='--', alpha=0.4)

    # Mark active site residues
    for domain, residues in ACTIVE_SITES.items():
        for fl_res in residues:
            ax_plddt.axvline(fl_res, color='red', linewidth=0.5, alpha=0.3, linestyle=':')

    ax_plddt.set_xlim(1, max(len(af2_plddt), 1810))
    ax_plddt.set_ylim(0, 100)
    ax_plddt.set_xlabel('Full-length residue number', fontsize=10)
    ax_plddt.set_ylabel('pLDDT', fontsize=10)
    ax_plddt.set_title('D. Per-residue pLDDT: AF2 (gray) vs AF3 construct (colored)',
                        fontsize=12, fontweight='bold')

    # Legend
    handles = [mpatches.Patch(color='gray', alpha=0.4, label='AF2 full-length')]
    handles += [mpatches.Patch(color=DOMAIN_COLORS[d], label=d)
                for d in CONSTRUCT_ATS_ORDER]
    handles.append(mpatches.Patch(facecolor='red', alpha=0.15, label='MD1 (missing)'))
    ax_plddt.legend(handles=handles, loc='lower left', ncol=4, fontsize=6,
                    framealpha=0.8)

    # ---- Panel E: Domain-average pLDDT table ----
    ax_table = fig.add_subplot(gs[2, 1])
    ax_table.axis('off')

    # Compute domain-average pLDDT for AF2 and AF3
    table_data = []
    for dname in ALL_DOMAIN_ORDER:
        fl_s, fl_e = DOMAIN_DEFS_ATS[dname]
        # AF2
        af2_vals = af2_plddt[fl_s-1:fl_e]
        af2_mean = np.mean(af2_vals)
        # AF3
        af3_vals = []
        for cpos, fl_pos in construct_to_fl_map.items():
            if fl_s <= fl_pos <= fl_e:
                af3_vals.append(af3_plddt[cpos])
        if af3_vals:
            af3_mean = np.mean(af3_vals)
            delta = af3_mean - af2_mean
            table_data.append([dname, f'{af2_mean:.1f}', f'{af3_mean:.1f}',
                              f'{delta:+.1f}', f'{fl_e-fl_s+1}'])
        else:
            table_data.append([dname, f'{af2_mean:.1f}', 'N/A', 'N/A',
                              f'{fl_e-fl_s+1}'])

    col_labels = ['Domain', 'AF2 pLDDT', 'AF3 pLDDT', '\u0394pLDDT', 'Residues']
    table = ax_table.table(cellText=table_data, colLabels=col_labels,
                           loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.4)

    # Color the domain cells and highlight delta
    for i, row in enumerate(table_data):
        dname = row[0]
        color = DOMAIN_COLORS.get(dname, '#ffffff')
        table[i+1, 0].set_facecolor(color)
        table[i+1, 0].set_alpha(0.3)
        # Color delta cell
        if row[3] != 'N/A':
            delta = float(row[3])
            if delta > 2:
                table[i+1, 3].set_facecolor('#d4edda')
            elif delta < -2:
                table[i+1, 3].set_facecolor('#f8d7da')

    # Header styling
    for j in range(5):
        table[0, j].set_facecolor('#343a40')
        table[0, j].set_text_props(color='white', fontweight='bold')

    ax_table.set_title('E. Domain-average pLDDT comparison',
                        fontsize=12, fontweight='bold', pad=20)

    # Overall title
    ptm_val = af3_data.get('ptm_mean', af3_data.get('ptm'))
    rs_val = af3_data.get('ranking_score_mean', af3_data.get('ranking_score'))
    ptm_str = f'pTM={ptm_val:.2f}' if ptm_val else ''
    rs_str = f'Ranking={rs_val:.2f}' if rs_val else ''
    n_models = af3_data.get('n_models', 1)
    model_str = f'AF3 averaged over {n_models} models' if n_models > 1 else title_suffix

    fig.suptitle(
        f'PARP14: AF3 Construct vs AF2 Full-length\n'
        f'(Construct: {n_residues} res, {ptm_str} {rs_str}, {model_str}  |  '
        f'AF2: {af2_data["n_residues"]} res, mean pLDDT={np.mean(af2_plddt):.1f})',
        fontsize=13, fontweight='bold', y=0.99
    )

    if output_path is None:
        output_path = OUTPUT_DIR / 'longest_construct_vs_af2.png'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'Saved: {output_path}')
    plt.close(fig)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Analyze longest PARP14 AF3 construct vs AF2 full-length')
    parser.add_argument('--seed', type=int, default=1, help='Seed number (1-5)')
    parser.add_argument('--sample', type=int, default=0, help='Sample number (0-4)')
    parser.add_argument('--all-models', action='store_true',
                        help='Average over all 25 AF3 models')
    args = parser.parse_args()

    boundaries, fl_to_construct, n_residues = compute_construct_boundaries()
    print(f'Construct: {CONSTRUCT}')
    print(f'Residues: {n_residues}')
    print(f'Domains (atS): {", ".join(CONSTRUCT_ATS_ORDER)}')
    print()

    # Load AF2 data
    print('Loading AF2 full-length structure...')
    af2_data = load_af2_data()

    # Load AF3 data
    if args.all_models:
        print('Loading all AF3 models...')
        af3_data = load_all_models()
        title_suffix = f'averaged over {af3_data["n_models"]} models'
        out_name = 'longest_construct_vs_af2_all.png'
    else:
        print(f'Loading AF3 seed-{args.seed}_sample-{args.sample}...')
        af3_data = load_model_data(args.seed, args.sample)
        title_suffix = f'seed-{args.seed}_sample-{args.sample}'
        out_name = f'longest_construct_vs_af2_s{args.seed}m{args.sample}.png'

    output_path = OUTPUT_DIR / out_name
    print(f'\nGenerating figure...')
    make_figure(af3_data, af2_data, boundaries, fl_to_construct,
                n_residues, title_suffix=title_suffix, output_path=output_path)
    print('Done.')


if __name__ == '__main__':
    main()
