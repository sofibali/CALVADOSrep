#!/usr/bin/env python3
"""
Compute per-domain AF3 metrics for all PARP14 domain combinations.

For each of the ~782 domain combinations (25 AF3 models each = 5 seeds x 5 samples),
computes per-domain pLDDT and intra-domain PAE, plus global metrics (pTM, ranking_score).

Outputs a structured directory of per-combination JSON files for fast random access,
plus a summary CSV.

Parallelizes across combinations using multiprocessing.

Usage:
    python 01_compute_domain_metrics.py                     # 4 workers (default)
    python 01_compute_domain_metrics.py --workers 16        # 16 workers
    python 01_compute_domain_metrics.py --workers 16 --force  # recompute all
"""

import os
import sys
import json
import argparse
import numpy as np
from pathlib import Path
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Tuple, Optional

# ============================================================
# Domain definitions
# ============================================================

# Full-length PARP14 domain boundaries (1-indexed, inclusive)
DOMAIN_DEFS = {
    'rrm1':    (1, 145),
    'rrm2':    (146, 224),
    'rrm3':    (225, 314),
    'kh1-kh6': (315, 737),
    'kh7a':    (738, 789),
    'md1':     (790, 981),
    'md1l1':   (790, 1004),
    'md2':     (1004, 1193),
    'md3':     (1207, 1388),
    'khb-kh8': (1389, 1533),
    'wwe':     (1534, 1602),
    'art':     (1603, 1801),
}

# Canonical ordering for construct assembly (matches naming convention)
DOMAIN_ORDER = [
    'rrm1', 'rrm2', 'rrm3', 'kh1-kh6', 'kh7a',
    'md1', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art',
]

AF3_BASE = Path(__file__).resolve().parent.parent / 'alphafold_outputs'
OUTPUT_BASE = Path(__file__).resolve().parent / 'domain_metrics'


# ============================================================
# Parse construct composition from directory name
# ============================================================

def parse_construct_domains(name: str) -> List[str]:
    """Parse a construct directory name into ordered domain list.

    E.g. 'kh1-kh6_kh7a_md2_md3_art' -> ['kh1-kh6', 'kh7a', 'md2', 'md3', 'art']
    """
    # Domain names that contain underscores or hyphens need careful parsing.
    # Strategy: greedily match known domain names left-to-right.
    remaining = name
    domains = []
    while remaining:
        remaining = remaining.lstrip('_')
        if not remaining:
            break
        matched = False
        # Try longest domain names first to avoid partial matches
        for dname in sorted(DOMAIN_DEFS.keys(), key=len, reverse=True):
            if remaining.startswith(dname):
                rest = remaining[len(dname):]
                if rest == '' or rest.startswith('_'):
                    domains.append(dname)
                    remaining = rest
                    matched = True
                    break
        if not matched:
            # Unknown token — skip to next underscore
            parts = remaining.split('_', 1)
            remaining = parts[1] if len(parts) > 1 else ''
    return domains


def compute_domain_boundaries_in_construct(domains: List[str]) -> Dict[str, Tuple[int, int]]:
    """Compute 0-indexed (start, end_exclusive) residue ranges for each domain
    within a concatenated construct (no linkers between domains).

    Handles overlapping boundaries: if domain N ends at residue X and domain N+1
    starts at residue X, the shared residue is counted once.
    """
    boundaries = {}
    offset = 0
    prev_end_fl = None  # previous domain's full-length end residue

    for dname in domains:
        fl_start, fl_end = DOMAIN_DEFS[dname]
        n_residues = fl_end - fl_start + 1

        # If this domain's FL start overlaps with previous domain's FL end,
        # the shared residue was already counted
        if prev_end_fl is not None and fl_start <= prev_end_fl:
            overlap = prev_end_fl - fl_start + 1
            n_residues -= overlap
            # Don't adjust offset — the overlap residues are already in previous domain

        boundaries[dname] = (offset, offset + n_residues)
        offset += n_residues
        prev_end_fl = fl_end

    return boundaries


# ============================================================
# Extract metrics from AF3 outputs
# ============================================================

def extract_per_residue_plddt(conf: Dict) -> Optional[np.ndarray]:
    """Extract per-residue pLDDT from confidences.json (token-level averaging)."""
    if 'atom_plddts' not in conf:
        return None

    atom_plddts = np.array(conf['atom_plddts'])

    if 'token_chain_ids' in conf:
        n_tokens = len(conf['token_chain_ids'])
        if n_tokens == 0:
            return None
        n_atoms = len(atom_plddts)
        atoms_per_token = n_atoms / n_tokens
        per_res = []
        for i in range(n_tokens):
            s = int(round(i * atoms_per_token))
            e = int(round((i + 1) * atoms_per_token))
            if s < n_atoms:
                per_res.append(float(np.mean(atom_plddts[s:e])))
        return np.array(per_res)

    return atom_plddts


def analyze_single_model(seed_dir: Path, domain_boundaries: Dict[str, Tuple[int, int]]) -> Optional[Dict]:
    """Analyze one seed-sample model. Returns dict with per-domain and global metrics."""
    result = {}

    # Summary confidences
    sf = seed_dir / 'summary_confidences.json'
    if sf.exists():
        with open(sf) as f:
            sc = json.load(f)
        for k in ['ptm', 'ranking_score', 'fraction_disordered']:
            if k in sc and sc[k] is not None:
                result[k] = float(sc[k])

    # Full confidences
    ff = seed_dir / 'confidences.json'
    if not ff.exists():
        return result if result else None

    with open(ff) as f:
        fc = json.load(f)

    # Per-residue pLDDT
    plddt = extract_per_residue_plddt(fc)
    if plddt is None or len(plddt) < 2:
        return result if result else None

    result['global_mean_plddt'] = float(np.mean(plddt))
    result['global_median_plddt'] = float(np.median(plddt))
    result['n_residues'] = len(plddt)

    # Per-domain pLDDT
    for dname, (start, end) in domain_boundaries.items():
        if end <= len(plddt):
            domain_plddt = plddt[start:end]
            result[f'{dname}_mean_plddt'] = float(np.mean(domain_plddt))
            result[f'{dname}_median_plddt'] = float(np.median(domain_plddt))
            result[f'{dname}_min_plddt'] = float(np.min(domain_plddt))
            result[f'{dname}_n_residues'] = end - start

    # PAE matrix — compute intra-domain mean PAE
    if 'pae' in fc:
        pae = np.array(fc['pae'])
        if pae.ndim == 2:
            result['global_mean_pae'] = float(np.mean(pae))
            for dname, (start, end) in domain_boundaries.items():
                if end <= pae.shape[0]:
                    intra_pae = pae[start:end, start:end]
                    result[f'{dname}_intra_pae'] = float(np.mean(intra_pae))

            # Inter-domain PAE for all domain pairs
            domain_names = list(domain_boundaries.keys())
            for i, d1 in enumerate(domain_names):
                s1, e1 = domain_boundaries[d1]
                for d2 in domain_names[i+1:]:
                    s2, e2 = domain_boundaries[d2]
                    if e1 <= pae.shape[0] and e2 <= pae.shape[0]:
                        inter = pae[s1:e1, s2:e2]
                        result[f'{d1}_vs_{d2}_inter_pae'] = float(np.mean(inter))

    return result


def process_combination(args) -> Optional[Dict]:
    """Process one domain combination: analyze all 25 seed-sample models.

    Returns a dict with combination metadata and per-metric mean/sd across models.
    """
    combo_name, af3_base, output_base, force = args

    # Check cache
    out_file = output_base / f'{combo_name}.json'
    if not force and out_file.exists():
        try:
            with open(out_file) as f:
                return json.load(f)
        except (json.JSONDecodeError, KeyError):
            pass  # recompute if corrupt

    combo_dir = af3_base / combo_name

    # Parse domains
    domains = parse_construct_domains(combo_name)
    if not domains:
        return None

    domain_boundaries = compute_domain_boundaries_in_construct(domains)

    # Find seed-sample directories
    seed_dirs = []
    for d in sorted(combo_dir.iterdir()):
        if d.is_dir() and d.name.startswith('seed-') and '_sample-' in d.name:
            seed_dirs.append(d)

    if not seed_dirs:
        return None

    # Analyze each model
    model_results = []
    for sd in seed_dirs:
        r = analyze_single_model(sd, domain_boundaries)
        if r:
            model_results.append(r)

    if not model_results:
        return None

    # Aggregate: mean ± sd across models
    all_keys = set()
    for r in model_results:
        all_keys.update(r.keys())

    summary = {
        'combination': combo_name,
        'domains': domains,
        'domain_boundaries': {k: list(v) for k, v in domain_boundaries.items()},
        'n_models': len(model_results),
        'n_seed_dirs': len(seed_dirs),
        'metrics': {},
    }

    # Boolean flags for easy filtering
    summary['has_kh7a'] = 'kh7a' in domains
    summary['has_khb_kh8'] = 'khb-kh8' in domains
    summary['has_split_kh'] = 'kh7a' in domains and 'khb-kh8' in domains
    summary['has_md2'] = 'md2' in domains
    summary['has_md3'] = 'md3' in domains

    for key in sorted(all_keys):
        if key == 'n_residues':
            # Just take from first model (same for all)
            summary['n_residues'] = model_results[0].get('n_residues')
            continue
        if key.endswith('_n_residues'):
            summary['metrics'][key] = model_results[0].get(key)
            continue

        vals = [r[key] for r in model_results if key in r]
        if vals:
            summary['metrics'][key] = {
                'mean': float(np.mean(vals)),
                'sd': float(np.std(vals)),
                'min': float(np.min(vals)),
                'max': float(np.max(vals)),
                'n': len(vals),
            }

    # Save per-combination JSON
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, 'w') as f:
        json.dump(summary, f, indent=2)

    return summary


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Compute per-domain AF3 metrics for PARP14')
    parser.add_argument('--workers', type=int, default=4,
                        help='Number of parallel workers (default: 4)')
    parser.add_argument('--force', action='store_true',
                        help='Recompute even if cached results exist')
    parser.add_argument('--af3-dir', type=str, default=str(AF3_BASE),
                        help='AlphaFold3 outputs directory')
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_BASE),
                        help='Output directory for per-combination JSONs')
    args = parser.parse_args()

    af3_base = Path(args.af3_dir)
    output_base = Path(args.output_dir)
    output_base.mkdir(parents=True, exist_ok=True)

    # Discover all combinations
    combos = sorted([
        d.name for d in af3_base.iterdir()
        if d.is_dir() and (d / 'seed-1_sample-0').is_dir()
    ])
    print(f'Found {len(combos)} domain combinations with AF3 data')
    print(f'Output: {output_base}/')
    print(f'Workers: {args.workers}')

    # Count cached
    if not args.force:
        cached = sum(1 for c in combos if (output_base / f'{c}.json').exists())
        print(f'Cached: {cached}, to compute: {len(combos) - cached}')

    # Process in parallel
    work_args = [(c, af3_base, output_base, args.force) for c in combos]

    results = []
    with Pool(processes=args.workers) as pool:
        for i, result in enumerate(pool.imap_unordered(process_combination, work_args), 1):
            if result is not None:
                results.append(result)
            if i % 50 == 0 or i == len(combos):
                print(f'  [{i}/{len(combos)}] processed ({len(results)} successful)')

    print(f'\nCompleted: {len(results)}/{len(combos)} combinations')

    # Write summary CSV (flat, one row per combination, mean values only)
    write_summary_csv(results, output_base)

    # Write index file
    write_index(results, output_base)


def write_summary_csv(results: List[Dict], output_base: Path):
    """Write a flat summary CSV with one row per combination, mean values only."""
    import csv

    # Collect all metric keys
    all_metric_keys = set()
    for r in results:
        all_metric_keys.update(r.get('metrics', {}).keys())

    # Separate into _n_residues (int) vs normal metrics (mean/sd)
    count_keys = sorted(k for k in all_metric_keys if k.endswith('_n_residues'))
    metric_keys = sorted(k for k in all_metric_keys if not k.endswith('_n_residues'))

    header = [
        'combination', 'n_domains', 'n_models', 'n_residues',
        'has_kh7a', 'has_khb_kh8', 'has_split_kh', 'has_md2', 'has_md3',
    ]
    for k in metric_keys:
        header.extend([f'{k}_mean', f'{k}_sd'])
    for k in count_keys:
        header.append(k)

    csv_path = output_base / 'summary.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header, extrasaction='ignore')
        writer.writeheader()
        for r in sorted(results, key=lambda x: x['combination']):
            row = {
                'combination': r['combination'],
                'n_domains': len(r['domains']),
                'n_models': r['n_models'],
                'n_residues': r.get('n_residues'),
                'has_kh7a': r.get('has_kh7a', False),
                'has_khb_kh8': r.get('has_khb_kh8', False),
                'has_split_kh': r.get('has_split_kh', False),
                'has_md2': r.get('has_md2', False),
                'has_md3': r.get('has_md3', False),
            }
            metrics = r.get('metrics', {})
            for k in metric_keys:
                m = metrics.get(k, {})
                if isinstance(m, dict):
                    row[f'{k}_mean'] = m.get('mean')
                    row[f'{k}_sd'] = m.get('sd')
            for k in count_keys:
                row[k] = metrics.get(k)
            writer.writerow(row)

    print(f'Summary CSV: {csv_path} ({len(results)} rows)')


def write_index(results: List[Dict], output_base: Path):
    """Write an index JSON for quick lookup of what's available."""
    index = {}
    for r in results:
        name = r['combination']
        index[name] = {
            'domains': r['domains'],
            'n_models': r['n_models'],
            'n_residues': r.get('n_residues'),
            'has_split_kh': r.get('has_split_kh', False),
            'has_md2': r.get('has_md2', False),
            'has_md3': r.get('has_md3', False),
            'file': f'{name}.json',
        }

    index_path = output_base / 'index.json'
    with open(index_path, 'w') as f:
        json.dump(index, f, indent=2)
    print(f'Index: {index_path} ({len(index)} entries)')


if __name__ == '__main__':
    main()
