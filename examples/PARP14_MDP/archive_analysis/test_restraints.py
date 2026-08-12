#!/usr/bin/env python3
"""
Test framework for PARP14 domain restraint boundaries.

Tests systematic 5-residue boundary shifts (±15 residues from medium baseline)
for all domains at k=700, plus KH7a-KHb custom inter-domain restraints with
the same boundary sweep.

Each test is a 2 ns FL PARP14 simulation, analyzed for per-domain Rg, RMSF,
and contact map correlation vs crystal/AF2 reference.

Usage:
    python test_restraints.py --prepare           # generate all test configs
    python test_restraints.py --run --parallel 8  # run tests
    python test_restraints.py --analyze           # analyze + rank
    python test_restraints.py --all --parallel 8  # all three
    python test_restraints.py --list              # list configurations

Boundary shift scheme:
    shift_m15: medium start-15, end+15 (widest)
    shift_m10: medium start-10, end+10
    shift_m05: medium start-5,  end+5
    shift_000: medium baseline
    shift_p05: medium start+5,  end-5
    shift_p10: medium start+10, end-10
    shift_p15: medium start+15, end-15 (narrowest)

KH7a-KHb custom restraints:
    For each boundary shift, the residue ranges used to generate inter-domain
    pairs are shifted correspondingly. Pairs are all CA-CA within 0.9 nm in
    the AF2 structure between the (shifted) KH7a and KHb ranges.
"""

import os
import sys
import json
import yaml
import subprocess
import numpy as np
from pathlib import Path
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed

# ============================================================
# Configuration
# ============================================================

CWD = Path(__file__).resolve().parent
AF2_PDB = CWD / 'input' / 'parp14.pdb'
FL_RESIDUES = CWD / 'input' / 'residues_CALVADOS3.csv'
TEST_DIR = CWD / 'restraint_tests'

TEST_STEPS = 200_000    # 2 ns
TEST_WFREQ = 1_000      # 200 frames
EQUIL_FRAMES = 50       # 0.5 ns equilibration
TEMP = 293
IONIC = 0.19
PH = 7.0
PLATFORM = 'CPU'
K_HARMONIC = 700.0
K_CUSTOM = 350.0        # inter-domain (KH7a-KHb) restraint strength
CUTOFF_RESTR = 0.9      # nm, for custom restraint pair generation
N_FL = 1801

# Experimental Rg references (nm)
EXPERIMENTAL_RG = {
    'MD1L1': 1.509, 'MD2': 1.533, 'WWE': 1.461, 'ART': 1.636,
}
EXPERIMENTAL_WEIGHT = {
    'MD1L1': 2.0, 'MD2': 2.0, 'MD3': 1.0, 'WWE': 1.5, 'ART': 1.5,
}
KH_HOMOLOG_RG = {'single_kh': 1.072, 'tandem_kh': 1.485}

# Full domain extents (hard limits — boundaries cannot exceed these)
DOMAIN_EXTENTS = {
    'RRM1': (1, 145), 'RRM2': (146, 224), 'RRM3': (225, 314),
    'KH1': (315, 384), 'KH2': (385, 454), 'KH3': (455, 520),
    'KH4': (521, 593), 'KH5': (594, 665), 'KH6': (666, 737),
    'KH7a': (738, 789),
    'MD1L1': (790, 978),   # capped at 3VFQ boundary
    'MD2': (1005, 1193),   # capped at 3VFQ boundary
    'MD3': (1207, 1388),
    'KHb': (1389, 1461), 'KH8': (1462, 1533),
    'WWE': (1534, 1602), 'ART': (1603, 1801),
}

# Per-domain "structured extents" override DOMAIN_EXTENTS for the trim=0
# starting point. Use these for domains where the full extent includes a
# large unstructured segment (e.g., RRM1's disordered N-terminus).
# All other domains start trimming from their full DOMAIN_EXTENTS.
STRUCTURED_EXTENTS = {
    'RRM1': (6, 88),   # exclude disordered N-term (1-5) and tail (89-145)
}

# Domains EXCLUDED from restraints entirely
EXCLUDED_DOMAINS = set()

# Analysis domains (full extent for metric computation)
ANALYSIS_DOMAINS = dict(DOMAIN_EXTENTS)

DOMAIN_ORDER = list(DOMAIN_EXTENTS.keys())

# Boundary shifts to test (number of residues to TRIM from EACH side of every
# restrained domain). 0 = fully restrained at extent. Higher = more residues
# removed from both edges (i.e., narrower core).
SHIFTS = [0, 5, 10, 15, 20, 25, 30]
MIN_DOMAIN_SIZE = 10  # minimum residues for a domain to be included


# ============================================================
# AF2 reference
# ============================================================

def get_af2_plddt():
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('af2', str(AF2_PDB))
    return np.array([
        [a for a in r if a.name == 'CA'][0].bfactor
        for r in structure[0].get_residues() if r.id[0] == ' '
        and any(a.name == 'CA' for a in r)
    ])


def get_af2_ca_coords():
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('af2', str(AF2_PDB))
    return np.array([
        [a for a in r if a.name == 'CA'][0].get_vector().get_array()
        for r in structure[0].get_residues() if r.id[0] == ' '
        and any(a.name == 'CA' for a in r)
    ])


# ============================================================
# Boundary shifting
# ============================================================

def shift_boundaries(trim):
    """Apply a symmetric trim to each domain's starting extent.

    Starting extent = STRUCTURED_EXTENTS[name] if defined (for domains with
    disordered tails/heads), else DOMAIN_EXTENTS[name] (full domain).
    Each "shift" value = number of residues removed from EACH side
    of the starting extent (so total residues lost = 2 * trim).

    Domains shrunk below MIN_DOMAIN_SIZE are dropped.
    Domains in EXCLUDED_DOMAINS are skipped.
    """
    result = {}
    for dname, full_extent in DOMAIN_EXTENTS.items():
        if dname in EXCLUDED_DOMAINS:
            continue
        ext_s, ext_e = STRUCTURED_EXTENTS.get(dname, full_extent)
        new_s = ext_s + trim
        new_e = ext_e - trim
        if new_e - new_s + 1 >= MIN_DOMAIN_SIZE:
            result[dname] = (new_s, new_e)
    return result


def generate_custom_restraint_pairs(ca_coords, kh7a_range, khb_range,
                                     cutoff_nm=CUTOFF_RESTR):
    """Generate KH7a-KHb inter-domain restraint pairs from AF2 coords.

    ca_coords: (N, 3) array in Angstroms
    Returns list of (resid1, resid2, distance_nm)
    """
    pairs = []
    cutoff_ang = cutoff_nm * 10.0
    for r1 in range(kh7a_range[0], kh7a_range[1] + 1):
        if r1 - 1 >= len(ca_coords):
            continue
        c1 = ca_coords[r1 - 1]
        for r2 in range(khb_range[0], khb_range[1] + 1):
            if r2 - 1 >= len(ca_coords):
                continue
            c2 = ca_coords[r2 - 1]
            d = np.linalg.norm(c1 - c2)
            if d <= cutoff_ang:
                pairs.append((r1, r2, d / 10.0))
    return pairs


# ============================================================
# Generate test configurations
# ============================================================

def generate_test_configs():
    """Generate boundary-shift tests + custom restraint variants."""
    configs = []

    # Baseline: no restraints
    configs.append({
        'name': 'no_restraints',
        'shift': None,
        'custom_kh7a_khb': False,
        'description': 'No restraints',
        'domains': {},
    })

    # Baseline: current (no KH)
    configs.append({
        'name': 'baseline_current',
        'shift': None,
        'custom_kh7a_khb': False,
        'description': 'Current restraints (no KH domains)',
        'domains': {
            'RRM1': (6, 88), 'RRM2': (150, 223), 'RRM3': (227, 301),
            'MD1L1': (791, 978), 'MD2': (1003, 1190), 'MD3': (1216, 1387),
            'WWE': (1523, 1601), 'ART': (1605, 1801),
        },
    })

    # Trim tests: trim N residues from each side of every domain extent
    # (RRM1 always excluded). Tested with and without KH7a-KHb custom restraints.
    for trim in SHIFTS:
        label = f"{trim:02d}"
        domains = shift_boundaries(trim)

        configs.append({
            'name': f'trim_{label}',
            'shift': trim,
            'custom_kh7a_khb': False,
            'description': f'Trim {trim} residues each side, no custom restr',
            'domains': dict(domains),
        })

        configs.append({
            'name': f'trim_{label}_kh7ab',
            'shift': trim,
            'custom_kh7a_khb': True,
            'description': f'Trim {trim} residues each side + KH7a-KHb custom',
            'domains': dict(domains),
        })

    return configs


# ============================================================
# Prepare simulation directories
# ============================================================

def prepare_test(config, plddt, ca_coords):
    """Prepare one test simulation directory."""
    test_dir = TEST_DIR / config['name']
    input_dir = test_dir / 'input'
    input_dir.mkdir(parents=True, exist_ok=True)

    # Symlink shared files
    for fname in ['parp14.pdb', 'residues_CALVADOS3.csv']:
        dest = input_dir / fname
        src = CWD / 'input' / fname
        if not dest.exists() and src.exists():
            dest.symlink_to(src)

    domains = config['domains']
    # Always set restraint=True so CALVADOS reads sequence from PDB
    # (rather than from a placeholder FASTA file). With an empty domains.yaml,
    # no harmonic bonds get added — equivalent to "no restraints".
    use_restraints = True
    n_domains = len(domains)

    # Write domains.yaml
    domain_list = []
    domain_info = {}
    for dname in sorted(domains.keys(), key=lambda x: domains[x][0]):
        ds, de = domains[dname]
        domain_list.append([ds, de])
        domain_info[dname] = {
            'start': ds, 'end': de,
            'n_residues': de - ds + 1,
            'mean_plddt': float(plddt[ds-1:de].mean()),
        }

    with open(input_dir / 'domains.yaml', 'w') as f:
        yaml.dump({'parp14': domain_list}, f, default_flow_style=True)

    # Custom restraints (KH7a-KHb) — only if config requests AND has any domains
    use_custom = config['custom_kh7a_khb'] and n_domains > 0
    n_custom_pairs = 0
    if use_custom:
        # Get KH7a and KHb ranges from this config's domains
        kh7a_range = domains.get('KH7a')
        khb_range = domains.get('KHb')
        if kh7a_range and khb_range:
            pairs = generate_custom_restraint_pairs(
                ca_coords, kh7a_range, khb_range)
            n_custom_pairs = len(pairs)
            cres_path = input_dir / 'custom_restraints.txt'
            with open(cres_path, 'w') as f:
                for r1, r2, d_nm in sorted(pairs):
                    f.write(f'parp14 1 {r1} | parp14 1 {r2} | '
                            f'{d_nm:.3f} {K_CUSTOM:.1f}\n')
        else:
            use_custom = False

    # Write config.yaml from defaults
    default_config = (CWD.parent.parent / 'calvados' / 'data' /
                      'default_config.yaml')
    with open(default_config) as f:
        config_dict = yaml.safe_load(f)
    config_dict.update({
        'box': [300, 300, 300],
        'temp': TEMP, 'ionic': IONIC, 'pH': PH,
        'steps': TEST_STEPS, 'wfreq': TEST_WFREQ,
        'topol': 'center', 'platform': PLATFORM,
        'custom_restraints': use_custom,
    })
    if use_custom:
        config_dict['fcustom_restraints'] = 'input/custom_restraints.txt'
    with open(test_dir / 'config.yaml', 'w') as f:
        yaml.dump(config_dict, f, default_flow_style=False)

    # Write components.yaml
    components_dict = {
        'defaults': {
            'molecule_type': 'protein', 'nmol': 1,
            'charge_termini': 'both', 'alpha': 0,
            'ffasta': 'fastabib.fasta', 'kb': 8033.0,
            'ext_restraint': True,
            'restraint': use_restraints,
            'cutoff_restr': CUTOFF_RESTR,
            'pdb_folder': str(input_dir),
            'restraint_type': 'harmonic',
            'k_harmonic': K_HARMONIC,
            'fdomains': str(input_dir / 'domains.yaml'),
            'k_go': 15.0, 'use_com': True, 'periodic': False,
            'colabfold': 1,
            'bfac_shift': 0.8, 'bfac_width': 50.0,
            'pae_shift': 0.3, 'pae_width': 15.0,
            'rna_kb1': 8033.0, 'rna_kb2': 8033.0,
            'rna_ka': 7.24, 'rna_pa': 3.14,
            'rna_nb_sigma': 0.4, 'rna_nb_scale': 15,
            'rna_nb_cutoff': 0.6, 'n_ends': 1,
            'ptm_name': 'example_ptm', 'ptm_locations': [],
            'fresidues': str(FL_RESIDUES),
        },
        'system': {'parp14': {}},
    }
    with open(test_dir / 'components.yaml', 'w') as f:
        yaml.dump(components_dict, f, default_flow_style=False)

    # Write run.py
    with open(test_dir / 'run.py', 'w') as f:
        f.write("from calvados import sim\n"
                "from argparse import ArgumentParser\n\n"
                "if __name__ == '__main__':\n"
                "    parser = ArgumentParser()\n"
                "    parser.add_argument('--path', default='.')\n"
                "    parser.add_argument('--config', default='config.yaml')\n"
                "    parser.add_argument('--components', default='components.yaml')\n"
                "    args = parser.parse_args()\n"
                "    sim.run(path=args.path, fconfig=args.config, "
                "fcomponents=args.components)\n")

    # Save metadata
    metadata = {
        'name': config['name'],
        'shift': config['shift'],
        'custom_kh7a_khb': use_custom,
        'n_custom_pairs': n_custom_pairs,
        'description': config['description'],
        'domains': {k: list(v) for k, v in domains.items()},
        'k_values': {k: K_HARMONIC for k in domains},
        'k_harmonic_mean': K_HARMONIC,
        'n_domains': len(domain_list),
        'n_restrained_residues': sum(de-ds+1 for ds, de in domains.values()),
        'domain_info': domain_info,
    }
    with open(test_dir / 'test_config.json', 'w') as f:
        json.dump(metadata, f, indent=2)

    return config['name'], metadata


# ============================================================
# Run simulations
# ============================================================

def run_test(test_name):
    test_dir = TEST_DIR / test_name
    if not (test_dir / 'run.py').exists():
        return test_name, False, 'run.py not found'
    if list(test_dir.glob('*.dcd')):
        return test_name, True, 'already completed'

    python_exe = '/home/sbali/miniconda3/envs/CALVADOS/bin/python'
    if not os.path.exists(python_exe):
        python_exe = sys.executable

    try:
        result = subprocess.run(
            [python_exe, 'run.py'], cwd=str(test_dir),
            capture_output=True, text=True, timeout=14400)
        if result.returncode != 0:
            with open(test_dir / 'error.log', 'w') as f:
                f.write(result.stdout + '\n' + result.stderr)
            return test_name, False, result.stderr[-500:]
        return test_name, True, ''
    except subprocess.TimeoutExpired:
        return test_name, False, 'timeout (4h)'
    except Exception as e:
        return test_name, False, str(e)


# ============================================================
# Analysis
# ============================================================

def analyze_test(test_name, plddt, ca_coords_af2):
    import mdtraj as md

    test_dir = TEST_DIR / test_name
    dcd_files = list(test_dir.glob('*.dcd'))
    top_file = test_dir / 'top.pdb'
    if not dcd_files or not top_file.exists():
        return None

    with open(test_dir / 'test_config.json') as f:
        config = json.load(f)

    try:
        traj = md.load(str(dcd_files[0]), top=str(top_file))
    except Exception as e:
        return {'error': str(e)}

    traj = traj[EQUIL_FRAMES:]
    if len(traj) < 10:
        return {'error': 'too few frames'}

    ca_indices = traj.topology.select('name CA')
    if len(ca_indices) != N_FL:
        return {'error': f'got {len(ca_indices)} CA, expected {N_FL}'}

    ca_traj = traj.atom_slice(ca_indices)
    avg_pos = ca_traj.xyz.mean(axis=0)
    rmsf = np.sqrt(np.mean(np.sum((ca_traj.xyz - avg_pos)**2, axis=2), axis=0))

    domain_metrics = {}
    for dname, (ds, de) in ANALYSIS_DOMAINS.items():
        idx = list(range(ds-1, de))
        if len(idx) < 3:
            continue

        n_domain_res = len(idx)
        mean_plddt = plddt[ds-1:de].mean()

        # Rg
        domain_xyz = ca_traj.xyz[:, idx, :]
        com = domain_xyz.mean(axis=1, keepdims=True)
        rg_per_frame = np.sqrt(np.mean(
            np.sum((domain_xyz - com)**2, axis=2), axis=1))
        mean_rg = float(rg_per_frame.mean())
        std_rg = float(rg_per_frame.std())

        # Reference Rg
        af2_dc = ca_coords_af2[idx] / 10.0
        af2_com = af2_dc.mean(axis=0)
        af2_rg = float(np.sqrt(np.mean(np.sum((af2_dc - af2_com)**2, axis=1))))

        if dname in EXPERIMENTAL_RG:
            ref_rg = EXPERIMENTAL_RG[dname]
            ref_source = 'experimental'
        elif dname.startswith('KH') and n_domain_res < 100:
            ref_rg = KH_HOMOLOG_RG['single_kh'] * np.sqrt(n_domain_res / 65.0)
            ref_source = 'kh_homolog_scaled'
        else:
            ref_rg = af2_rg
            ref_source = 'af2'

        # Contact map correlation
        cmap_corr = 0.0
        if n_domain_res > 5:
            sim_dmap = np.zeros((n_domain_res, n_domain_res))
            for fxyz in domain_xyz:
                diff = fxyz[:, None, :] - fxyz[None, :, :]
                sim_dmap += (np.sqrt(np.sum(diff**2, axis=2)) < 1.2)
            sim_dmap /= len(domain_xyz)

            af2_diff = af2_dc[:, None, :] - af2_dc[None, :, :]
            af2_cmap = (np.sqrt(np.sum(af2_diff**2, axis=2)) < 1.2).astype(float)

            mask = np.triu(np.ones((n_domain_res, n_domain_res), dtype=bool), k=2)
            sf, af = sim_dmap[mask], af2_cmap[mask]
            if sf.std() > 0 and af.std() > 0:
                cmap_corr = float(np.corrcoef(sf, af)[0, 1])

        rg_ratio = mean_rg / ref_rg if ref_rg > 0 else 999.0
        weight = EXPERIMENTAL_WEIGHT.get(dname, 1.0)
        restrained = dname in config.get('domains', {})
        k_val = config.get('k_values', {}).get(dname, 0.0)

        domain_metrics[dname] = {
            'mean_plddt': float(mean_plddt),
            'mean_rg': mean_rg, 'std_rg': std_rg,
            'af2_rg': af2_rg, 'ref_rg': float(ref_rg),
            'ref_source': ref_source, 'rg_ratio': rg_ratio,
            'mean_rmsf': float(rmsf[idx].mean()),
            'cmap_corr': cmap_corr,
            'n_residues': n_domain_res,
            'restrained': restrained, 'k_value': k_val,
            'score_weight': weight,
        }

    # Global scores
    vals = list(domain_metrics.values())
    weights = np.array([v['score_weight'] * v['mean_plddt'] / 100 for v in vals])
    rg_dev = np.average([abs(v['rg_ratio'] - 1) for v in vals], weights=weights)
    w_cmap = np.average([v['cmap_corr'] for v in vals], weights=weights)

    from scipy.stats import spearmanr
    pl = np.array([v['mean_plddt'] for v in vals])
    rm = np.array([v['mean_rmsf'] for v in vals])
    plddt_rmsf_corr = float(spearmanr(pl, 1/(rm+0.001))[0]) if len(rm) > 3 else 0.0

    combined = (0.4 * (1 - min(float(rg_dev), 1.0)) +
                0.4 * max(float(w_cmap), 0.0) +
                0.2 * max(plddt_rmsf_corr, 0.0))

    result = {
        'test_name': test_name,
        'shift': config.get('shift'),
        'custom_kh7a_khb': config.get('custom_kh7a_khb', False),
        'n_domains': config['n_domains'],
        'n_restrained_residues': config['n_restrained_residues'],
        'n_custom_pairs': config.get('n_custom_pairs', 0),
        'domain_metrics': domain_metrics,
        'global_scores': {
            'rg_deviation': float(rg_dev),
            'weighted_cmap_corr': float(w_cmap),
            'plddt_rmsf_corr': plddt_rmsf_corr,
            'combined_score': float(combined),
        },
    }
    with open(test_dir / 'analysis_results.json', 'w') as f:
        json.dump(result, f, indent=2)
    return result


# ============================================================
# Reporting & plots
# ============================================================

def print_ranking(all_results):
    ranked = sorted(all_results,
                    key=lambda x: x['global_scores']['combined_score'],
                    reverse=True)

    print("\n" + "=" * 105)
    print("RESTRAINT TEST RANKING")
    print("=" * 105)
    print(f"{'Rk':>3} {'Test Name':<25} {'Shift':>5} {'KH7ab':>5} "
          f"{'nDom':>4} {'nRes':>5} {'nCust':>5} "
          f"{'RgDev':>6} {'CMap':>5} {'plRMSF':>6} {'Score':>6}")
    print("-" * 105)

    for rank, r in enumerate(ranked, 1):
        gs = r['global_scores']
        sh = r.get('shift', '-')
        sh = f"{sh:+d}" if isinstance(sh, int) else str(sh)
        cust = 'Y' if r.get('custom_kh7a_khb') else 'N'
        print(f"{rank:>3} {r['test_name']:<25} {sh:>5} {cust:>5} "
              f"{r['n_domains']:>4} {r['n_restrained_residues']:>5} "
              f"{r.get('n_custom_pairs', 0):>5} "
              f"{gs['rg_deviation']:>6.3f} {gs['weighted_cmap_corr']:>5.3f} "
              f"{gs['plddt_rmsf_corr']:>6.3f} {gs['combined_score']:>6.3f}")

    # Top 3 detail
    for rank, r in enumerate(ranked[:3], 1):
        print(f"\n--- #{rank}: {r['test_name']} "
              f"(score={r['global_scores']['combined_score']:.3f}) ---")
        print(f"  {'Domain':<8} {'pLDDT':>5} {'Rg':>5} {'Ref':>5} "
              f"{'Ratio':>5} {'RMSF':>5} {'CMap':>5} {'Src':<8} {'R':>1}")
        print("  " + "-" * 60)
        for dn in DOMAIN_ORDER:
            dm = r['domain_metrics'].get(dn)
            if not dm:
                continue
            print(f"  {dn:<8} {dm['mean_plddt']:>5.1f} {dm['mean_rg']:>5.2f} "
                  f"{dm['ref_rg']:>5.2f} {dm['rg_ratio']:>5.2f} "
                  f"{dm['mean_rmsf']:>5.3f} {dm['cmap_corr']:>5.3f} "
                  f"{dm['ref_source'][:8]:<8} "
                  f"{'Y' if dm['restrained'] else 'N':>1}")

    return ranked


def generate_summary_plots(all_results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot_dir = TEST_DIR / 'summary_plots'
    plot_dir.mkdir(exist_ok=True)

    ranked = sorted(all_results,
                    key=lambda x: x['global_scores']['combined_score'],
                    reverse=True)

    # Separate into no-custom and custom groups
    no_cust = [r for r in all_results
               if not r.get('custom_kh7a_khb') and r.get('shift') is not None]
    with_cust = [r for r in all_results if r.get('custom_kh7a_khb')]

    # Sort by shift
    no_cust.sort(key=lambda x: x.get('shift', 0))
    with_cust.sort(key=lambda x: x.get('shift', 0))

    # ── Plot 1: Combined score vs shift ──
    fig, ax = plt.subplots(figsize=(10, 6))
    if no_cust:
        shifts_nc = [r['shift'] for r in no_cust]
        scores_nc = [r['global_scores']['combined_score'] for r in no_cust]
        ax.plot(shifts_nc, scores_nc, 'o-', color='steelblue', linewidth=2,
                markersize=8, label='Domain restraints only')
    if with_cust:
        shifts_wc = [r['shift'] for r in with_cust]
        scores_wc = [r['global_scores']['combined_score'] for r in with_cust]
        ax.plot(shifts_wc, scores_wc, 's-', color='darkorange', linewidth=2,
                markersize=8, label='+ KH7a-KHb custom restraints')

    # Add baselines
    for r in all_results:
        if r['test_name'] == 'baseline_current':
            ax.axhline(r['global_scores']['combined_score'],
                       color='gray', linestyle='--', alpha=0.7,
                       label='Baseline (no KH)')
        elif r['test_name'] == 'no_restraints':
            ax.axhline(r['global_scores']['combined_score'],
                       color='red', linestyle=':', alpha=0.5,
                       label='No restraints')

    ax.set_xlabel('Boundary Shift (residues)', fontsize=12)
    ax.set_ylabel('Combined Score', fontsize=12)
    ax.set_title('Score vs Boundary Shift\n'
                 '(negative = wider, positive = narrower)', fontsize=14)
    ax.legend(fontsize=9)
    ax.set_xticks(SHIFTS)
    ax.set_xticklabels([f"{s:+d}" for s in SHIFTS])
    plt.tight_layout()
    fig.savefig(plot_dir / 'score_vs_shift.png', dpi=200, bbox_inches='tight')
    fig.savefig(plot_dir / 'score_vs_shift.pdf', bbox_inches='tight')
    plt.close()

    # ── Plot 2: Per-domain Rg ratio heatmap ──
    tests_to_show = no_cust + with_cust
    if not tests_to_show:
        tests_to_show = ranked

    fig, ax = plt.subplots(figsize=(18, max(4, len(tests_to_show) * 0.5)))
    rg_matrix = np.full((len(tests_to_show), len(DOMAIN_ORDER)), np.nan)
    for i, r in enumerate(tests_to_show):
        for j, dn in enumerate(DOMAIN_ORDER):
            dm = r['domain_metrics'].get(dn, {})
            rg_matrix[i, j] = dm.get('rg_ratio', np.nan)

    im = ax.imshow(rg_matrix, cmap='RdBu_r', vmin=0.5, vmax=2.0, aspect='auto')
    ax.set_xticks(range(len(DOMAIN_ORDER)))
    ax.set_xticklabels(DOMAIN_ORDER, rotation=45, ha='right', fontsize=9)
    labels = []
    for r in tests_to_show:
        sh = r.get('shift')
        cust = '+KH7ab' if r.get('custom_kh7a_khb') else ''
        labels.append(f"trim {sh}{cust}" if isinstance(sh, int)
                      else r['test_name'])
    ax.set_yticks(range(len(tests_to_show)))
    ax.set_yticklabels(labels, fontsize=8)
    plt.colorbar(im, ax=ax, label='Rg ratio (sim/ref)', shrink=0.8)

    for i in range(len(tests_to_show)):
        for j in range(len(DOMAIN_ORDER)):
            v = rg_matrix[i, j]
            if not np.isnan(v):
                c = 'white' if abs(v - 1.0) > 0.4 else 'black'
                ax.text(j, i, f'{v:.2f}', ha='center', va='center',
                        fontsize=5, color=c)

    ax.set_title('Per-Domain Rg Ratio (1.0 = ideal)', fontsize=14)
    plt.tight_layout()
    fig.savefig(plot_dir / 'rg_heatmap_shifts.png', dpi=200, bbox_inches='tight')
    fig.savefig(plot_dir / 'rg_heatmap_shifts.pdf', bbox_inches='tight')
    plt.close()

    # ── Plot 3: Per-domain metrics vs boundary shift (ALL domains) ──
    n_domains = len(DOMAIN_ORDER)
    n_cols = 3
    n_rows = (n_domains + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 4 * n_rows),
                             sharex=True)
    axes = axes.flatten()

    for idx, dn in enumerate(DOMAIN_ORDER):
        ax = axes[idx]
        ax2 = ax.twinx()

        # Without custom restraints
        if no_cust:
            shifts_nc = [r['shift'] for r in no_cust]
            rg_nc = [r['domain_metrics'].get(dn, {}).get('rg_ratio', np.nan)
                     for r in no_cust]
            cmap_nc = [r['domain_metrics'].get(dn, {}).get('cmap_corr', np.nan)
                       for r in no_cust]
            ax.plot(shifts_nc, rg_nc, 'o-', color='steelblue',
                    label='Rg', linewidth=1.5)
            ax2.plot(shifts_nc, cmap_nc, 's--', color='green', alpha=0.7,
                     label='CMap', linewidth=1)

        # With KH7a-KHb custom restraints
        if with_cust:
            shifts_wc = [r['shift'] for r in with_cust]
            rg_wc = [r['domain_metrics'].get(dn, {}).get('rg_ratio', np.nan)
                     for r in with_cust]
            cmap_wc = [r['domain_metrics'].get(dn, {}).get('cmap_corr', np.nan)
                       for r in with_cust]
            ax.plot(shifts_wc, rg_wc, '^-', color='darkorange',
                    label='Rg +KH7ab', linewidth=1.5)
            ax2.plot(shifts_wc, cmap_wc, 'v--', color='limegreen', alpha=0.7,
                     label='CMap +KH7ab', linewidth=1)

        ax2.set_ylim(0, 1.1)
        if idx % n_cols == n_cols - 1:
            ax2.set_ylabel('CMap corr', fontsize=8, color='green')
        else:
            ax2.set_yticklabels([])

        ax.axhline(1.0, color='red', linestyle='--', alpha=0.3)
        ds, de = ANALYSIS_DOMAINS.get(dn, (0, 0))
        ax.set_title(f'{dn} ({ds}-{de})', fontsize=10, fontweight='bold')
        ax.set_ylabel('Rg ratio', fontsize=8)
        ax.set_ylim(0.3, 2.5)
        if idx >= n_domains - n_cols:
            ax.set_xlabel('Boundary shift', fontsize=9)
        if idx == 0:
            # Combined legend from both axes
            lines1, labels1 = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax.legend(lines1 + lines2, labels1 + labels2,
                      fontsize=5, loc='upper left')

    # Hide unused subplots
    for idx in range(n_domains, len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle('All Domain Metrics vs Boundary Shift', fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(plot_dir / 'all_domain_metrics_vs_shift.png', dpi=200,
                bbox_inches='tight')
    fig.savefig(plot_dir / 'all_domain_metrics_vs_shift.pdf',
                bbox_inches='tight')
    plt.close()

    # ── Plot 4: Ranking bar chart ──
    fig, ax = plt.subplots(figsize=(12, max(4, len(ranked) * 0.4)))
    names = [r['test_name'] for r in ranked]
    scores = [r['global_scores']['combined_score'] for r in ranked]
    colors = ['darkorange' if r.get('custom_kh7a_khb') else 'steelblue'
              for r in ranked]
    ax.barh(range(len(ranked)), scores, color=colors, alpha=0.8)
    ax.set_yticks(range(len(ranked)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel('Combined Score', fontsize=12)
    ax.set_title('All Configurations Ranked', fontsize=14)
    ax.invert_yaxis()

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color='steelblue', label='Domain restraints only'),
        Patch(color='darkorange', label='+ KH7a-KHb custom'),
    ], fontsize=9)

    plt.tight_layout()
    fig.savefig(plot_dir / 'ranking.png', dpi=200, bbox_inches='tight')
    fig.savefig(plot_dir / 'ranking.pdf', bbox_inches='tight')
    plt.close()

    # ── Plot 5: Boundary schematic (sequence boxes) ──
    print("Generating boundary schematic...")
    render_boundary_schematic(plot_dir)

    # ── Summary table ──
    print("Writing simulation summary table...")
    write_simulation_summary_table(plot_dir)

    # ── Plot 6: Cartoon snapshot of every simulation ──
    print("Rendering cartoon snapshots of all simulations...")
    render_all_cartoons(all_results, plot_dir)

    print(f"\nPlots saved to {plot_dir}/")


def write_simulation_summary_table(plot_dir):
    """Write a summary table of every test configuration: boundaries,
    restraint strength, custom restraint pairs, and run/analysis status.

    Outputs:
      simulation_summary.csv  — flat row-per-config table
      simulation_summary.md   — markdown table for easy viewing
    """
    import csv

    configs = generate_test_configs()

    rows = []
    for c in configs:
        test_dir = TEST_DIR / c['name']
        meta_file = test_dir / 'test_config.json'
        results_file = test_dir / 'analysis_results.json'

        # Load metadata if prepared
        meta = {}
        if meta_file.exists():
            with open(meta_file) as f:
                meta = json.load(f)

        # Run status
        dcds = list(test_dir.glob('*.dcd')) if test_dir.exists() else []
        status = 'not_prepared'
        if test_dir.exists():
            if dcds:
                status = 'completed'
            elif (test_dir / 'error.log').exists():
                status = 'failed'
            elif meta_file.exists():
                status = 'prepared'

        # Analysis score
        score = None
        rg_dev = None
        cmap = None
        if results_file.exists():
            with open(results_file) as f:
                res = json.load(f)
            gs = res.get('global_scores', {})
            score = gs.get('combined_score')
            rg_dev = gs.get('rg_deviation')
            cmap = gs.get('weighted_cmap_corr')

        # Domain boundaries summary
        domains = meta.get('domains', c.get('domains', {}))
        domain_str_parts = []
        for dn in DOMAIN_ORDER:
            if dn in domains:
                ds, de = domains[dn]
                domain_str_parts.append(f"{dn}:{ds}-{de}")
        domain_summary = '; '.join(domain_str_parts)

        rows.append({
            'name': c['name'],
            'shift': c.get('shift') if c.get('shift') is not None else '',
            'custom_kh7a_khb': 'Y' if c.get('custom_kh7a_khb') else 'N',
            'k_harmonic': K_HARMONIC if domains else 0,
            'k_custom': K_CUSTOM if c.get('custom_kh7a_khb') else 0,
            'n_domains': len(domains),
            'n_restrained_residues': sum(de-ds+1 for ds, de in domains.values()),
            'n_custom_pairs': meta.get('n_custom_pairs', 0),
            'status': status,
            'combined_score': f"{score:.4f}" if score is not None else '',
            'rg_deviation': f"{rg_dev:.4f}" if rg_dev is not None else '',
            'cmap_corr': f"{cmap:.4f}" if cmap is not None else '',
            'domain_boundaries': domain_summary,
        })

    # Write CSV
    csv_path = plot_dir / 'simulation_summary.csv'
    fieldnames = list(rows[0].keys()) if rows else []
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    # Write markdown table (compact, for human reading)
    md_path = plot_dir / 'simulation_summary.md'
    with open(md_path, 'w') as f:
        f.write("# PARP14 Restraint Test Simulation Summary\n\n")
        f.write(f"Generated from {len(rows)} test configurations.\n\n")
        f.write("## Configuration Matrix\n\n")
        f.write("| Name | Shift | KH7ab | k_dom | k_cust | nDom | nRes | nPairs "
                "| Status | Score | RgDev | CMap |\n")
        f.write("|------|------:|:-----:|------:|-------:|-----:|-----:|-------:"
                "|--------|------:|------:|-----:|\n")
        for r in rows:
            sh = f"{r['shift']:+d}" if isinstance(r['shift'], int) else "-"
            f.write(f"| {r['name']} | {sh} | {r['custom_kh7a_khb']} | "
                    f"{r['k_harmonic']:.0f} | {r['k_custom']:.0f} | "
                    f"{r['n_domains']} | {r['n_restrained_residues']} | "
                    f"{r['n_custom_pairs']} | {r['status']} | "
                    f"{r['combined_score']} | {r['rg_deviation']} | "
                    f"{r['cmap_corr']} |\n")

        # Per-trim boundaries detail
        f.write("\n## Domain Boundaries per Trim Tier\n\n")
        f.write("Start point: every domain restrained at its **full extent**. "
                "RRM1 uses a custom **structured core (6-88)** that excludes "
                "the disordered N-terminus and tail. Trim N = remove N "
                "residues from EACH side (so total residues removed = 2N).\n\n")
        f.write("| Domain |")
        for s in SHIFTS:
            f.write(f" trim {s} |")
        f.write("\n|--------|")
        for _ in SHIFTS:
            f.write("-----|")
        f.write("\n")

        for dn in DOMAIN_ORDER:
            f.write(f"| **{dn}** |")
            for s in SHIFTS:
                bds = shift_boundaries(s)
                if dn in bds:
                    ds, de = bds[dn]
                    f.write(f" {ds}-{de} ({de-ds+1}) |")
                else:
                    f.write(" — |")
            f.write("\n")

        # Counts summary
        n_completed = sum(1 for r in rows if r['status'] == 'completed')
        n_failed = sum(1 for r in rows if r['status'] == 'failed')
        n_prepared = sum(1 for r in rows if r['status'] == 'prepared')
        f.write(f"\n## Status Summary\n\n")
        f.write(f"- **Total configs:** {len(rows)}\n")
        f.write(f"- **Completed (have trajectory):** {n_completed}\n")
        f.write(f"- **Failed:** {n_failed}\n")
        f.write(f"- **Prepared but not run:** {n_prepared}\n")
        f.write(f"- **Not prepared:** {len(rows) - n_completed - n_failed - n_prepared}\n")
        f.write(f"\n**Restraint strength:** k_harmonic = {K_HARMONIC} kJ/mol/nm² "
                f"(intra-domain), k_custom = {K_CUSTOM} kJ/mol/nm² (KH7a-KHb).\n")

    print(f"  Saved: {csv_path}")
    print(f"  Saved: {md_path}")


def render_boundary_schematic(plot_dir):
    """Draw a linear sequence schematic with domain extents as background
    boxes and each boundary shift tier as a colored bar. Serines are marked
    as ticks on the top row, with boundary serines (from PARP14_domains_atS.fasta)
    highlighted in a darker color."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, FancyBboxPatch
    from Bio.PDB import PDBParser
    from Bio.Data.IUPACData import protein_letters_3to1

    # Get all serine positions
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('af2', str(AF2_PDB))
    serines = []
    for r in structure[0].get_residues():
        if r.id[0] != ' ':
            continue
        rn = r.resname.strip()
        if protein_letters_3to1.get(rn.capitalize(), 'X') == 'S':
            serines.append(r.id[1])

    # Boundary serines from PARP14_domains_atS.fasta — these mark the domain
    # transitions used in the AF3 combinatorial library
    boundary_serines = {
        148:  'RRM1/RRM2',
        220:  'RRM2/RRM3',
        303:  'RRM3/KH_N',
        514:  'KH_N/KH_Ca',
        786:  'KH_Ca/MD1',
        1002: 'MD1/MD2',
        1183: 'MD2/MD3',
        1392: 'MD3/KH_b',
        1501: 'KH_b/WWE',
        1601: 'WWE/ART',
    }

    DOMAIN_COLORS = {
        'RRM1': '#1f77b4', 'RRM2': '#2ca02c', 'RRM3': '#aec7e8',
        'KH1': '#ff7f0e', 'KH2': '#ff9f40', 'KH3': '#ffbf70',
        'KH4': '#ffd700', 'KH5': '#e6c200', 'KH6': '#ccac00',
        'KH7a': '#d62728', 'MD1L1': '#17becf', 'MD2': '#9467bd',
        'MD3': '#8c564b', 'KHb': '#e377c2', 'KH8': '#f7b6d2',
        'WWE': '#7f7f7f', 'ART': '#bcbd22',
    }

    # Row labels: top row = sequence + serines, then one row per shift
    shift_labels = [f"trim {s}" for s in SHIFTS]
    row_labels = ['Sequence'] + shift_labels
    n_rows = len(row_labels)
    row_height = 0.8
    gap = 0.3

    fig_height = (n_rows + 1) * (row_height + gap) + 1.5
    fig, ax = plt.subplots(figsize=(24, fig_height))

    # Draw each row
    for row_idx, (label, shift) in enumerate(
            [('Sequence', None)] + list(zip(shift_labels, SHIFTS))):
        y = (n_rows - row_idx) * (row_height + gap)

        # Row label
        ax.text(-30, y + row_height / 2, label, ha='right', va='center',
                fontsize=9, fontweight='bold')

        if shift is None:
            # Top row: full sequence bar (light gray) with domain extent shading
            ax.add_patch(Rectangle((1, y), N_FL, row_height,
                                   facecolor='#f0f0f0', edgecolor='black',
                                   linewidth=0.5))

            # Domain extent boxes (light fill)
            for dname, (ds, de) in DOMAIN_EXTENTS.items():
                c = DOMAIN_COLORS.get(dname, '#cccccc')
                ax.add_patch(Rectangle(
                    (ds, y), de - ds + 1, row_height,
                    facecolor=c, edgecolor='none', alpha=0.25))
                # Domain name label
                mid = (ds + de) / 2
                ax.text(mid, y + row_height + 0.15, dname,
                        ha='center', va='bottom', fontsize=5.5,
                        rotation=45, color=c, fontweight='bold')

            # Serine tick marks (faint gray) for all serines
            for s_pos in serines:
                if s_pos in boundary_serines:
                    continue  # boundary serines drawn separately below
                ax.plot([s_pos, s_pos],
                        [y + row_height, y + row_height + 0.12],
                        color='gray', linewidth=0.4, alpha=0.4)

            # Boundary serines: darker, taller ticks with labels
            for s_pos, label in boundary_serines.items():
                ax.plot([s_pos, s_pos],
                        [y + row_height, y + row_height + 0.30],
                        color='red', linewidth=1.5, alpha=0.95)
                # Label above the tick
                ax.text(s_pos, y + row_height + 0.35, f'S{s_pos}',
                        ha='center', va='bottom', fontsize=4.5,
                        rotation=90, color='red', fontweight='bold')

            # Vertical guideline at boundary serines extending down all rows
            for s_pos in boundary_serines:
                ax.axvline(s_pos, color='red', linewidth=0.3, alpha=0.25,
                           ymin=0, ymax=1, zorder=0)

        else:
            # Shift row: gray background, colored boxes for restrained regions
            ax.add_patch(Rectangle((1, y), N_FL, row_height,
                                   facecolor='#f8f8f8', edgecolor='gray',
                                   linewidth=0.3))

            boundaries = shift_boundaries(shift)
            for dname, (ds, de) in boundaries.items():
                c = DOMAIN_COLORS.get(dname, '#999999')
                width = de - ds + 1
                ax.add_patch(FancyBboxPatch(
                    (ds, y + 0.05), width, row_height - 0.1,
                    boxstyle="round,pad=0.02",
                    facecolor=c, edgecolor='black', linewidth=0.3,
                    alpha=0.75))
                # Residue range inside box (use varying font / placement
                # depending on width)
                if width > 50:
                    ax.text((ds + de) / 2, y + row_height / 2,
                            f'{ds}-{de} ({width})',
                            ha='center', va='center', fontsize=4.5,
                            color='white', fontweight='bold')
                elif width > 20:
                    ax.text((ds + de) / 2, y + row_height / 2,
                            f'{ds}-{de}',
                            ha='center', va='center', fontsize=3.5,
                            color='white', fontweight='bold')

                # Always tag the start and end residue IDs above/below the box
                ax.text(ds, y - 0.05, str(ds),
                        ha='center', va='top', fontsize=3.5,
                        color=c, fontweight='bold', rotation=90)
                ax.text(de, y - 0.05, str(de),
                        ha='center', va='top', fontsize=3.5,
                        color=c, fontweight='bold', rotation=90)

    # Axis formatting
    ax.set_xlim(-80, N_FL + 30)
    ax.set_ylim(0, (n_rows + 1) * (row_height + gap) + 0.5)
    ax.set_xlabel('Full-Length Residue Position', fontsize=12)
    ax.set_title('Domain Restraint Boundaries per Trim Tier\n'
                 '(start = full extents [RRM1 = structured core 6-88]; '
                 'trim N = remove N residues from each side; '
                 'gray ticks = serines, red ticks = boundary serines)',
                 fontsize=14)

    # Tick marks every 100 residues
    ax.set_xticks(range(0, N_FL + 1, 100))
    ax.set_xticklabels([str(x) for x in range(0, N_FL + 1, 100)], fontsize=7)
    ax.set_yticks([])

    # Domain color legend at bottom
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=DOMAIN_COLORS[dn], alpha=0.75,
                             edgecolor='black', linewidth=0.3, label=dn)
                       for dn in DOMAIN_ORDER]
    legend_elements.append(plt.Line2D([0], [0], color='gray', linewidth=1,
                                       label='Serine'))
    legend_elements.append(plt.Line2D([0], [0], color='red', linewidth=2,
                                       label='Boundary serine'))
    ax.legend(handles=legend_elements, loc='lower center',
              ncol=10, fontsize=6, frameon=True,
              bbox_to_anchor=(0.5, -0.08))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.tight_layout()
    fig.savefig(plot_dir / 'boundary_schematic.png', dpi=250,
                bbox_inches='tight')
    fig.savefig(plot_dir / 'boundary_schematic.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {plot_dir / 'boundary_schematic.png'}")


def render_all_cartoons(all_results, plot_dir):
    """Render a cartoon-style 2D trace of each simulation's middle frame,
    colored by domain, in a single multi-panel figure."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    import mdtraj as md

    DOMAIN_COLORS = {
        'RRM1': '#1f77b4', 'RRM2': '#2ca02c', 'RRM3': '#aec7e8',
        'KH1': '#ff7f0e', 'KH2': '#ff9f40', 'KH3': '#ffbf70',
        'KH4': '#ffd700', 'KH5': '#e6c200', 'KH6': '#ccac00',
        'KH7a': '#d62728', 'MD1L1': '#17becf', 'MD2': '#9467bd',
        'MD3': '#8c564b', 'KHb': '#e377c2', 'KH8': '#f7b6d2',
        'WWE': '#7f7f7f', 'ART': '#bcbd22',
    }

    # Build residue->color map
    res_colors = ['#dddddd'] * N_FL
    for dname, (ds, de) in ANALYSIS_DOMAINS.items():
        c = DOMAIN_COLORS.get(dname, '#999999')
        for i in range(ds - 1, de):
            if i < N_FL:
                res_colors[i] = c

    # Collect test names that have trajectories
    tests_with_traj = []
    for r in sorted(all_results, key=lambda x: (
            x.get('shift') if x.get('shift') is not None else -999,
            0 if not x.get('custom_kh7a_khb') else 1)):
        td = TEST_DIR / r['test_name']
        dcds = list(td.glob('*.dcd'))
        top = td / 'top.pdb'
        if dcds and top.exists():
            tests_with_traj.append((r['test_name'], str(dcds[0]), str(top)))

    if not tests_with_traj:
        print("  No trajectories found for cartoon rendering")
        return

    # Also add AF2 reference
    af2_pdb = str(AF2_PDB)

    n_panels = len(tests_with_traj) + 1  # +1 for AF2
    n_cols = 4
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5 * n_cols, 5 * n_rows))
    axes = axes.flatten()

    def project_and_draw(ax, xyz, title):
        """Draw a 2D projection of the CA trace, colored by domain."""
        # PCA-like: use the two axes with most variance
        centered = xyz - xyz.mean(axis=0)
        cov = np.cov(centered.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        # Use the two largest eigenvectors
        proj = centered @ eigvecs[:, -2:]  # (N, 2)

        # Draw line segments colored by domain
        points = proj.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Color each segment by the first residue's domain
        seg_colors = [res_colors[i] for i in range(len(segments))]

        lc = LineCollection(segments, colors=seg_colors, linewidths=1.5,
                           alpha=0.8)
        ax.add_collection(lc)

        # Add domain-colored dots at every 50th residue for landmarks
        for i in range(0, len(proj), 50):
            ax.plot(proj[i, 0], proj[i, 1], 'o', color=res_colors[i],
                    markersize=3, zorder=5)

        ax.set_xlim(proj[:, 0].min() - 2, proj[:, 0].max() + 2)
        ax.set_ylim(proj[:, 1].min() - 2, proj[:, 1].max() + 2)
        ax.set_aspect('equal')
        ax.set_title(title, fontsize=8, fontweight='bold')
        ax.set_xticks([])
        ax.set_yticks([])

    # Panel 0: AF2 reference
    try:
        af2_traj = md.load(af2_pdb)
        af2_ca = af2_traj.atom_slice(af2_traj.topology.select('name CA'))
        project_and_draw(axes[0], af2_ca.xyz[0] * 10, 'AF2 Reference')  # nm->Å
    except Exception as e:
        axes[0].text(0.5, 0.5, f'AF2 error:\n{e}', transform=axes[0].transAxes,
                     ha='center', fontsize=7)

    # Remaining panels: each test simulation (middle frame)
    for panel_idx, (test_name, dcd, top) in enumerate(tests_with_traj, 1):
        ax = axes[panel_idx]
        try:
            traj = md.load(dcd, top=top)
            ca = traj.atom_slice(traj.topology.select('name CA'))
            mid = ca.n_frames // 2
            xyz = ca.xyz[mid] * 10  # nm -> Å
            project_and_draw(ax, xyz, test_name)
        except Exception as e:
            ax.text(0.5, 0.5, f'Error:\n{str(e)[:50]}',
                    transform=ax.transAxes, ha='center', fontsize=7)

    # Hide unused panels
    for idx in range(n_panels, len(axes)):
        axes[idx].set_visible(False)

    # Add a domain color legend at the bottom
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=DOMAIN_COLORS[dn], label=dn)
                       for dn in DOMAIN_ORDER if dn in DOMAIN_COLORS]
    fig.legend(handles=legend_elements, loc='lower center',
               ncol=9, fontsize=7, frameon=True,
               bbox_to_anchor=(0.5, -0.01))

    plt.suptitle('Cartoon Snapshots (middle frame, PCA projection, '
                 'colored by domain)', fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(plot_dir / 'all_cartoons.png', dpi=200, bbox_inches='tight')
    fig.savefig(plot_dir / 'all_cartoons.pdf', bbox_inches='tight')
    plt.close()

    print(f"  Saved: {plot_dir / 'all_cartoons.png'} "
          f"({len(tests_with_traj)} sims + AF2)")


# ============================================================
# Main
# ============================================================

def main():
    parser = ArgumentParser(
        description='Test PARP14 restraint boundary shifts')
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--analyze', action='store_true')
    parser.add_argument('--all', action='store_true')
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--parallel', type=int, default=1)
    parser.add_argument('--tests', nargs='+',
                        help='Only specific tests by name')
    args = parser.parse_args()

    if args.all:
        args.prepare = args.run = args.analyze = True
    if not (args.prepare or args.run or args.analyze or args.list):
        parser.print_help()
        return

    print("Loading AF2 reference...")
    plddt = get_af2_plddt()
    ca_coords = get_af2_ca_coords()
    print(f"  {len(plddt)} residues, mean pLDDT={plddt.mean():.1f}")

    configs = generate_test_configs()

    if args.list:
        print(f"\n{len(configs)} test configurations:")
        print(f"{'Name':<25} {'Shift':>5} {'KH7ab':>5} "
              f"{'nDom':>5} {'nRes':>5}")
        print("-" * 70)
        for c in configs:
            sh = c.get('shift')
            sh = f"{sh:+d}" if isinstance(sh, int) else '-'
            cust = 'Y' if c['custom_kh7a_khb'] else 'N'
            nd = len(c['domains'])
            nr = sum(e-s+1 for s, e in c['domains'].values()) if c['domains'] else 0
            print(f"{c['name']:<25} {sh:>5} {cust:>5} {nd:>5} {nr:>5}")
        return

    if args.tests:
        configs = [c for c in configs if c['name'] in args.tests]
        print(f"Filtered to {len(configs)} tests")

    if args.prepare:
        print(f"\nPreparing {len(configs)} tests...")
        TEST_DIR.mkdir(exist_ok=True)
        for c in configs:
            name, meta = prepare_test(c, plddt, ca_coords)
            cust_str = f", {meta.get('n_custom_pairs', 0)} custom pairs" \
                if meta.get('custom_kh7a_khb') else ""
            print(f"  {name}: {meta['n_domains']} domains, "
                  f"{meta['n_restrained_residues']} res{cust_str}")

    if args.run:
        test_names = [c['name'] for c in configs]
        to_run = [n for n in test_names
                  if not list((TEST_DIR / n).glob('*.dcd'))]
        skipped = len(test_names) - len(to_run)
        if skipped:
            print(f"  Skipping {skipped} completed")
        if not to_run:
            print("All done!")
        else:
            print(f"\nRunning {len(to_run)} tests ({args.parallel} parallel)...")
            if args.parallel > 1:
                with ProcessPoolExecutor(max_workers=args.parallel) as pool:
                    futs = {pool.submit(run_test, n): n for n in to_run}
                    for f in as_completed(futs):
                        nm, ok, msg = f.result()
                        print(f"  {nm}: {'OK' if ok else f'FAIL: {msg[:60]}'}")
            else:
                for n in to_run:
                    print(f"  {n}...", end=' ', flush=True)
                    _, ok, msg = run_test(n)
                    print('OK' if ok else f'FAIL: {msg[:60]}')

    if args.analyze:
        test_names = [c['name'] for c in configs]
        print(f"\nAnalyzing {len(test_names)} tests...")
        all_results = []
        for n in test_names:
            td = TEST_DIR / n
            if not list(td.glob('*.dcd')):
                continue
            rf = td / 'analysis_results.json'
            if rf.exists():
                with open(rf) as f:
                    r = json.load(f)
                print(f"  {n}: cached ({r['global_scores']['combined_score']:.3f})")
            else:
                print(f"  {n}...", end=' ', flush=True)
                r = analyze_test(n, plddt, ca_coords)
                if not r or 'error' in r:
                    print(f"FAIL: {(r or {}).get('error', '?')}")
                    continue
                print(f"{r['global_scores']['combined_score']:.3f}")
            all_results.append(r)

        if all_results:
            ranked = print_ranking(all_results)
            generate_summary_plots(all_results)
            with open(TEST_DIR / 'ranking.json', 'w') as f:
                json.dump([{
                    'rank': i+1, 'name': r['test_name'],
                    'shift': r.get('shift'),
                    'custom_kh7a_khb': r.get('custom_kh7a_khb', False),
                    'scores': r['global_scores'],
                } for i, r in enumerate(ranked)], f, indent=2)


if __name__ == '__main__':
    main()
