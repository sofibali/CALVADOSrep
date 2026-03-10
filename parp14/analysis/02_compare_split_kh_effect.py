#!/usr/bin/env python3
"""
Compare AF3 predictions: with vs without split KH (KH7a + KHb-KH8).

Reads per-combination JSON files produced by 01_compute_domain_metrics.py.
Finds matched pairs of domain combinations that differ ONLY by the presence
or absence of both KH7a and KHb-KH8 together, then computes deltas in
pLDDT, PAE, pTM, and ranking_score focused on MD2 and MD3.

Outputs:
  - matched_pairs.json: all matched pairs with their delta metrics
  - matched_pairs_summary.csv: flat table of key deltas
  - split_kh_effect_report.txt: human-readable summary

Usage:
    python 02_compare_split_kh_effect.py
    python 02_compare_split_kh_effect.py --metrics-dir path/to/domain_metrics
"""

import json
import csv
import argparse
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

METRICS_DIR = Path(__file__).resolve().parent / 'domain_metrics'
OUTPUT_DIR = Path(__file__).resolve().parent / 'split_kh_comparison'


# ============================================================
# Find matched pairs
# ============================================================

def load_index(metrics_dir: Path) -> Dict:
    """Load the index.json from 01_compute_domain_metrics.py."""
    index_path = metrics_dir / 'index.json'
    with open(index_path) as f:
        return json.load(f)


def load_combination(metrics_dir: Path, name: str) -> Dict:
    """Load per-combination JSON."""
    with open(metrics_dir / f'{name}.json') as f:
        return json.load(f)


def find_matched_pairs(index: Dict) -> List[Tuple[str, str]]:
    """Find pairs where one has both KH7a+KHb-KH8 and the other has neither.

    A valid pair: the domain lists are identical except for the presence/absence
    of BOTH kh7a AND khb-kh8. Both must be present or both absent.

    Returns list of (with_split_kh_name, without_split_kh_name) tuples.
    """
    split_kh_domains = {'kh7a', 'khb-kh8'}

    # Group combinations by their non-split-KH domain set
    groups = {}  # frozenset(other_domains) -> {'with': [names], 'without': [names]}
    for name, info in index.items():
        domains = set(info['domains'])
        has_kh7a = 'kh7a' in domains
        has_khb_kh8 = 'khb-kh8' in domains

        # We only want combinations where BOTH are present or BOTH are absent
        if has_kh7a != has_khb_kh8:
            continue  # skip combinations with only one of the two

        other_domains = frozenset(domains - split_kh_domains)
        if other_domains not in groups:
            groups[other_domains] = {'with': [], 'without': []}

        if has_kh7a and has_khb_kh8:
            groups[other_domains]['with'].append(name)
        else:
            groups[other_domains]['without'].append(name)

    # Build matched pairs
    pairs = []
    for other_domains, group in groups.items():
        for w in group['with']:
            for wo in group['without']:
                pairs.append((w, wo))

    return sorted(pairs)


# ============================================================
# Compute deltas
# ============================================================

def get_metric_mean(data: Dict, metric_key: str) -> Optional[float]:
    """Extract mean value for a metric from combination data."""
    m = data.get('metrics', {}).get(metric_key)
    if isinstance(m, dict):
        return m.get('mean')
    return None


def get_metric_sd(data: Dict, metric_key: str) -> Optional[float]:
    """Extract SD value for a metric from combination data."""
    m = data.get('metrics', {}).get(metric_key)
    if isinstance(m, dict):
        return m.get('sd')
    return None


def compute_pair_deltas(with_data: Dict, without_data: Dict) -> Dict:
    """Compute delta metrics between with-split-KH and without-split-KH.

    Positive delta = metric is HIGHER with split KH present.
    """
    result = {
        'with_combination': with_data['combination'],
        'without_combination': without_data['combination'],
        'with_domains': with_data['domains'],
        'without_domains': without_data['domains'],
        'with_n_models': with_data['n_models'],
        'without_n_models': without_data['n_models'],
        'with_n_residues': with_data.get('n_residues'),
        'without_n_residues': without_data.get('n_residues'),
        'shared_domains': sorted(
            set(without_data['domains']) & set(with_data['domains'])
        ),
    }

    # Global metrics
    for metric in ['global_mean_plddt', 'ptm', 'ranking_score',
                    'fraction_disordered', 'global_mean_pae']:
        w_val = get_metric_mean(with_data, metric)
        wo_val = get_metric_mean(without_data, metric)
        w_sd = get_metric_sd(with_data, metric)
        wo_sd = get_metric_sd(without_data, metric)

        result[f'with_{metric}'] = w_val
        result[f'without_{metric}'] = wo_val
        result[f'with_{metric}_sd'] = w_sd
        result[f'without_{metric}_sd'] = wo_sd

        if w_val is not None and wo_val is not None:
            result[f'delta_{metric}'] = w_val - wo_val

    # Per-domain metrics for shared domains (focus on MD2, MD3, but compute all)
    shared = set(without_data['domains']) & set(with_data['domains'])
    for domain in sorted(shared):
        for suffix in ['mean_plddt', 'median_plddt', 'min_plddt', 'intra_pae']:
            key = f'{domain}_{suffix}'
            w_val = get_metric_mean(with_data, key)
            wo_val = get_metric_mean(without_data, key)
            w_sd = get_metric_sd(with_data, key)
            wo_sd = get_metric_sd(without_data, key)

            result[f'with_{key}'] = w_val
            result[f'without_{key}'] = wo_val
            result[f'with_{key}_sd'] = w_sd
            result[f'without_{key}_sd'] = wo_sd

            if w_val is not None and wo_val is not None:
                result[f'delta_{key}'] = w_val - wo_val

    # Inter-domain PAE between shared domains
    shared_list = sorted(shared)
    for i, d1 in enumerate(shared_list):
        for d2 in shared_list[i+1:]:
            key = f'{d1}_vs_{d2}_inter_pae'
            w_val = get_metric_mean(with_data, key)
            wo_val = get_metric_mean(without_data, key)
            if w_val is not None and wo_val is not None:
                result[f'with_{key}'] = w_val
                result[f'without_{key}'] = wo_val
                result[f'delta_{key}'] = w_val - wo_val

    return result


# ============================================================
# Output
# ============================================================

def write_matched_pairs_json(pairs_data: List[Dict], output_dir: Path):
    """Write full matched pairs data as JSON."""
    out = output_dir / 'matched_pairs.json'
    with open(out, 'w') as f:
        json.dump(pairs_data, f, indent=2)
    print(f'  {out} ({len(pairs_data)} pairs)')


def write_summary_csv(pairs_data: List[Dict], output_dir: Path):
    """Write a flat CSV with key delta metrics per pair."""
    # Focus on the most important columns
    key_columns = [
        'with_combination', 'without_combination',
        'with_n_residues', 'without_n_residues',
        'with_n_models', 'without_n_models',
    ]

    # Global deltas
    for m in ['global_mean_plddt', 'ptm', 'ranking_score',
              'fraction_disordered', 'global_mean_pae']:
        key_columns.extend([f'delta_{m}', f'with_{m}', f'without_{m}',
                            f'with_{m}_sd', f'without_{m}_sd'])

    # MD2 and MD3 specific deltas
    for domain in ['md2', 'md3']:
        for suffix in ['mean_plddt', 'intra_pae']:
            key = f'{domain}_{suffix}'
            key_columns.extend([f'delta_{key}', f'with_{key}', f'without_{key}',
                                f'with_{key}_sd', f'without_{key}_sd'])

    # Filter columns to those that actually exist in the data
    available_cols = set()
    for pd in pairs_data:
        available_cols.update(pd.keys())
    columns = [c for c in key_columns if c in available_cols]

    csv_path = output_dir / 'matched_pairs_summary.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        for pd in sorted(pairs_data, key=lambda x: x['with_combination']):
            writer.writerow({k: pd.get(k) for k in columns})
    print(f'  {csv_path} ({len(pairs_data)} rows, {len(columns)} columns)')


def write_report(pairs_data: List[Dict], output_dir: Path):
    """Write a human-readable text report."""
    report_path = output_dir / 'split_kh_effect_report.txt'

    # Separate pairs by whether they contain MD2, MD3, or both
    md2_pairs = [p for p in pairs_data if 'md2' in p.get('shared_domains', [])]
    md3_pairs = [p for p in pairs_data if 'md3' in p.get('shared_domains', [])]
    md2_md3_pairs = [p for p in pairs_data
                     if 'md2' in p.get('shared_domains', [])
                     and 'md3' in p.get('shared_domains', [])]

    def _stats(values):
        v = [x for x in values if x is not None]
        if not v:
            return 'N/A'
        return f'mean={np.mean(v):+.2f}, median={np.median(v):+.2f}, n={len(v)}'

    with open(report_path, 'w') as f:
        f.write('=' * 70 + '\n')
        f.write('PARP14: Effect of Split KH (KH7a + KHb-KH8) on AF3 Predictions\n')
        f.write('=' * 70 + '\n\n')

        f.write(f'Total matched pairs: {len(pairs_data)}\n')
        f.write(f'  Pairs containing MD2: {len(md2_pairs)}\n')
        f.write(f'  Pairs containing MD3: {len(md3_pairs)}\n')
        f.write(f'  Pairs containing both MD2+MD3: {len(md2_md3_pairs)}\n\n')

        f.write('Positive delta = metric is HIGHER when split KH is PRESENT.\n')
        f.write('For pLDDT/pTM: positive = better with split KH.\n')
        f.write('For PAE: positive = worse with split KH (higher error).\n\n')

        f.write('-' * 70 + '\n')
        f.write('GLOBAL METRICS (all pairs)\n')
        f.write('-' * 70 + '\n')
        for metric in ['global_mean_plddt', 'ptm', 'ranking_score',
                        'fraction_disordered', 'global_mean_pae']:
            deltas = [p.get(f'delta_{metric}') for p in pairs_data]
            f.write(f'  delta_{metric}: {_stats(deltas)}\n')

        for label, subset in [
            ('MD2-CONTAINING PAIRS', md2_pairs),
            ('MD3-CONTAINING PAIRS', md3_pairs),
            ('MD2+MD3-CONTAINING PAIRS', md2_md3_pairs),
        ]:
            if not subset:
                continue
            f.write(f'\n{"-" * 70}\n')
            f.write(f'{label} (n={len(subset)})\n')
            f.write(f'{"-" * 70}\n')

            for metric in ['global_mean_plddt', 'ptm', 'ranking_score']:
                deltas = [p.get(f'delta_{metric}') for p in subset]
                f.write(f'  delta_{metric}: {_stats(deltas)}\n')

            # Per-domain metrics
            for domain in ['md2', 'md3']:
                if not any(domain in p.get('shared_domains', []) for p in subset):
                    continue
                f.write(f'\n  {domain.upper()} domain:\n')
                domain_subset = [p for p in subset if domain in p.get('shared_domains', [])]
                for suffix in ['mean_plddt', 'intra_pae']:
                    key = f'{domain}_{suffix}'
                    deltas = [p.get(f'delta_{key}') for p in domain_subset]
                    with_vals = [p.get(f'with_{key}') for p in domain_subset]
                    without_vals = [p.get(f'without_{key}') for p in domain_subset]
                    f.write(f'    delta_{key}: {_stats(deltas)}\n')
                    f.write(f'      with split KH:    {_stats(with_vals)}\n')
                    f.write(f'      without split KH: {_stats(without_vals)}\n')

                # SD comparison
                for suffix in ['mean_plddt']:
                    key = f'{domain}_{suffix}'
                    with_sds = [p.get(f'with_{key}_sd') for p in domain_subset]
                    without_sds = [p.get(f'without_{key}_sd') for p in domain_subset]
                    f.write(f'    SD of {key} (prediction variability):\n')
                    f.write(f'      with split KH:    {_stats(with_sds)}\n')
                    f.write(f'      without split KH: {_stats(without_sds)}\n')

        # Individual pair details (for pairs with MD2+MD3)
        f.write(f'\n{"=" * 70}\n')
        f.write(f'INDIVIDUAL PAIR DETAILS (MD2+MD3 pairs)\n')
        f.write(f'{"=" * 70}\n')

        for p in sorted(md2_md3_pairs, key=lambda x: x['with_combination']):
            f.write(f'\n  WITH:    {p["with_combination"]}\n')
            f.write(f'  WITHOUT: {p["without_combination"]}\n')
            f.write(f'  Shared domains: {", ".join(p.get("shared_domains", []))}\n')

            delta_plddt = p.get('delta_global_mean_plddt')
            delta_ptm = p.get('delta_ptm')
            delta_rs = p.get('delta_ranking_score')
            f.write(f'  Global: dpLDDT={_fmt(delta_plddt)} '
                    f'dpTM={_fmt(delta_ptm)} '
                    f'dRankScore={_fmt(delta_rs)}\n')

            for domain in ['md2', 'md3']:
                dp = p.get(f'delta_{domain}_mean_plddt')
                da = p.get(f'delta_{domain}_intra_pae')
                w_sd = p.get(f'with_{domain}_mean_plddt_sd')
                wo_sd = p.get(f'without_{domain}_mean_plddt_sd')
                f.write(f'  {domain.upper()}: dpLDDT={_fmt(dp)} '
                        f'dPAE={_fmt(da)} '
                        f'SD(with)={_fmt(w_sd)} SD(without)={_fmt(wo_sd)}\n')

        f.write(f'\n{"=" * 70}\n')

    print(f'  {report_path}')


def _fmt(val):
    if val is None:
        return 'N/A'
    return f'{val:+.2f}'


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Compare AF3 predictions with vs without split KH')
    parser.add_argument('--metrics-dir', type=str, default=str(METRICS_DIR),
                        help='Directory with per-combination JSONs from step 01')
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_DIR),
                        help='Output directory')
    args = parser.parse_args()

    metrics_dir = Path(args.metrics_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load index
    print('Loading index...')
    index = load_index(metrics_dir)
    print(f'  {len(index)} combinations available')

    # Find matched pairs
    print('\nFinding matched pairs (with vs without KH7a + KHb-KH8)...')
    pairs = find_matched_pairs(index)
    print(f'  Found {len(pairs)} matched pairs')

    if not pairs:
        print('No matched pairs found. Check that domain_metrics/ has been computed.')
        return

    # Compute deltas for each pair
    print('\nComputing deltas...')
    pairs_data = []
    for with_name, without_name in pairs:
        with_data = load_combination(metrics_dir, with_name)
        without_data = load_combination(metrics_dir, without_name)
        deltas = compute_pair_deltas(with_data, without_data)
        pairs_data.append(deltas)

    # Write outputs
    print(f'\nWriting outputs to {output_dir}/:')
    write_matched_pairs_json(pairs_data, output_dir)
    write_summary_csv(pairs_data, output_dir)
    write_report(pairs_data, output_dir)

    print('\nDone.')


if __name__ == '__main__':
    main()
