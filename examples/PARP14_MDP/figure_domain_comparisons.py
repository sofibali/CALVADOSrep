#!/usr/bin/env python3
"""
Compare AF3 PAE/pLDDT across PARP14 domain constructs.

Matched-pair comparisons:
  A) KH7a effect: with vs without KH7a+KHb-KH8
  B) Macrodomain removal: with vs without MD1L1, MD2, or MD3
  C) Combined: all macrodomains vs none

Reads pre-computed domain metrics from parp14/analysis/domain_metrics/*.json
(produced by 01_compute_domain_metrics.py).

Outputs 6 SVG figures to figures/domain_comparisons/.
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

# ============================================================
# Configuration
# ============================================================

METRICS_DIR = Path('/home/sbali/CALVADOS/parp14/analysis/domain_metrics')
# Dated category subdir: figures/02_main_analysis/<YYYY-MM-DD>/domain_comparisons/
_HERE = Path('/home/sbali/CALVADOS/examples/PARP14_MDP')
import sys as _sys
_sys.path.insert(0, str(_HERE))
from _fig_layout import get_fig_dir as _get_fig_dir
FIG_DIR = _get_fig_dir('02_main_analysis', subname='domain_comparisons')

DOMAIN_DEFS = {
    'rrm1': (1, 145), 'rrm2': (146, 224), 'rrm3': (225, 314),
    'kh1-kh6': (315, 737), 'kh7a': (738, 789),
    'md1l1': (790, 1004), 'md2': (1004, 1193), 'md3': (1207, 1388),
    'khb-kh8': (1389, 1533), 'wwe': (1534, 1602), 'art': (1603, 1801),
}

DOMAIN_ORDER = [
    'rrm1', 'rrm2', 'rrm3', 'kh1-kh6', 'kh7a',
    'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art',
]

DOMAIN_COLORS = {
    'rrm1': '#1f77b4', 'rrm2': '#2ca02c', 'rrm3': '#ff7f0e',
    'kh1-kh6': '#d62728', 'kh7a': '#9467bd',
    'md1l1': '#8c564b', 'md2': '#e377c2', 'md3': '#7f7f7f',
    'khb-kh8': '#bcbd22', 'wwe': '#17becf', 'art': '#ff9896',
}

COMPARISON_COLORS = {
    'kh7a': '#9467bd',
    'md1l1': '#8c564b',
    'md2': '#e377c2',
    'md3': '#7f7f7f',
    'all_md': '#2166ac',
}

# ============================================================
# Data loading (reused from 02_compare_split_kh_effect.py)
# ============================================================

def load_index(metrics_dir: Path) -> Dict:
    with open(metrics_dir / 'index.json') as f:
        return json.load(f)

def load_combination(metrics_dir: Path, name: str) -> Dict:
    with open(metrics_dir / f'{name}.json') as f:
        return json.load(f)

def get_metric_mean(data: Dict, key: str) -> Optional[float]:
    m = data.get('metrics', {}).get(key)
    if isinstance(m, dict):
        return m.get('mean')
    return None

def get_metric_sd(data: Dict, key: str) -> Optional[float]:
    m = data.get('metrics', {}).get(key)
    if isinstance(m, dict):
        return m.get('sd')
    return None


# ============================================================
# Matched-pair finding
# ============================================================

def find_pairs(index: Dict, target_domains: Set[str]) -> List[Tuple[str, str]]:
    """Find matched pairs differing only by presence/absence of ALL target_domains.

    For a valid pair, one construct has ALL target domains and the other has NONE.
    All other domains must be identical.

    Returns list of (with_targets, without_targets) name tuples.
    """
    groups = {}  # frozenset(other_domains) -> {'with': [], 'without': []}
    for name, info in index.items():
        domains = set(info['domains'])
        has_targets = target_domains & domains
        missing_targets = target_domains - domains

        # Must have ALL or NONE of the target domains
        if has_targets and missing_targets:
            continue

        other = frozenset(domains - target_domains)
        if other not in groups:
            groups[other] = {'with': [], 'without': []}

        if has_targets:
            groups[other]['with'].append(name)
        else:
            groups[other]['without'].append(name)

    pairs = []
    for other, group in groups.items():
        for w in group['with']:
            for wo in group['without']:
                pairs.append((w, wo))
    return sorted(pairs)


def compute_pair_deltas(with_data: Dict, without_data: Dict) -> Dict:
    """Compute delta metrics (with - without). Positive = higher with target present."""
    result = {
        'with_combination': with_data['combination'],
        'without_combination': without_data['combination'],
        'with_domains': with_data['domains'],
        'without_domains': without_data['domains'],
        'with_n_models': with_data['n_models'],
        'without_n_models': without_data['n_models'],
        'shared_domains': sorted(
            set(without_data['domains']) & set(with_data['domains'])
        ),
    }

    # Global metrics
    for metric in ['global_mean_plddt', 'ptm', 'ranking_score',
                    'fraction_disordered', 'global_mean_pae']:
        w_val = get_metric_mean(with_data, metric)
        wo_val = get_metric_mean(without_data, metric)
        result[f'with_{metric}'] = w_val
        result[f'without_{metric}'] = wo_val
        if w_val is not None and wo_val is not None:
            result[f'delta_{metric}'] = w_val - wo_val

    # Per-domain metrics for shared domains
    shared = set(without_data['domains']) & set(with_data['domains'])
    for domain in sorted(shared):
        for suffix in ['mean_plddt', 'median_plddt', 'min_plddt', 'intra_pae']:
            key = f'{domain}_{suffix}'
            w_val = get_metric_mean(with_data, key)
            wo_val = get_metric_mean(without_data, key)
            result[f'with_{key}'] = w_val
            result[f'without_{key}'] = wo_val
            if w_val is not None and wo_val is not None:
                result[f'delta_{key}'] = w_val - wo_val

    # Inter-domain PAE
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
# Aggregation
# ============================================================

def build_delta_matrix(pairs_data: List[Dict], domain_list: List[str]):
    """Build NxN matrix of mean inter-domain PAE deltas.

    Returns (matrix, count_matrix) where matrix[i,j] is the mean delta
    inter-PAE between domain_list[i] and domain_list[j], and count_matrix
    has the number of pairs contributing to each cell.
    """
    n = len(domain_list)
    sums = np.zeros((n, n))
    counts = np.zeros((n, n))

    for p in pairs_data:
        shared = set(p.get('shared_domains', []))
        for i, d1 in enumerate(domain_list):
            if d1 not in shared:
                continue
            # Diagonal: intra-domain PAE delta
            key = f'delta_{d1}_intra_pae'
            val = p.get(key)
            if val is not None:
                sums[i, i] += val
                counts[i, i] += 1
            # Off-diagonal: inter-domain PAE delta
            for j, d2 in enumerate(domain_list):
                if j <= i or d2 not in shared:
                    continue
                # Try both orderings of domain names
                key1 = f'delta_{d1}_vs_{d2}_inter_pae'
                key2 = f'delta_{d2}_vs_{d1}_inter_pae'
                val = p.get(key1) or p.get(key2)
                if val is not None:
                    sums[i, j] += val
                    sums[j, i] += val
                    counts[i, j] += 1
                    counts[j, i] += 1

    matrix = np.full((n, n), np.nan)
    mask = counts > 0
    matrix[mask] = sums[mask] / counts[mask]
    return matrix, counts.astype(int)


def collect_domain_deltas(pairs_data: List[Dict], metric_suffix: str) -> Dict[str, List[float]]:
    """Collect delta values per domain for a given metric suffix (e.g. 'mean_plddt', 'intra_pae')."""
    result = {}
    for p in pairs_data:
        for d in p.get('shared_domains', []):
            key = f'delta_{d}_{metric_suffix}'
            val = p.get(key)
            if val is not None:
                result.setdefault(d, []).append(val)
    return result


# ============================================================
# Statistical helpers
# ============================================================

def wilcoxon_p(values):
    """Return Wilcoxon signed-rank test p-value, or None if too few samples."""
    from scipy.stats import wilcoxon
    v = [x for x in values if x is not None and not np.isnan(x)]
    if len(v) < 10:
        return None
    try:
        _, p = wilcoxon(v)
        return p
    except Exception:
        return None


def sig_stars(p):
    if p is None:
        return ''
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'ns'


# ============================================================
# Figure 1: Delta PAE heatmap — KH7a effect
# ============================================================

def fig_delta_pae_heatmap_kh7a(kh7a_pairs):
    """Domain×domain delta inter-PAE heatmap for KH7a effect."""
    # Use only domains that appear in shared domains
    all_shared = set()
    for p in kh7a_pairs:
        all_shared.update(p.get('shared_domains', []))
    domain_list = [d for d in DOMAIN_ORDER if d in all_shared]

    if len(domain_list) < 2:
        print("  Skip: too few shared domains for KH7a heatmap")
        return

    matrix, counts = build_delta_matrix(kh7a_pairs, domain_list)
    n = len(domain_list)

    fig, ax = plt.subplots(figsize=(9, 8))
    vmax = np.nanmax(np.abs(matrix[np.isfinite(matrix)])) if np.any(np.isfinite(matrix)) else 5
    vmax = max(vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax.imshow(matrix, cmap='RdBu_r', norm=norm, aspect='equal')

    # Annotate cells
    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            cnt = counts[i, j]
            if np.isnan(val):
                ax.text(j, i, '—', ha='center', va='center', fontsize=7, color='#999999')
            else:
                color = 'white' if abs(val) > vmax * 0.6 else 'black'
                label = f'{val:+.1f}\n(n={cnt})'
                ax.text(j, i, label, ha='center', va='center', fontsize=6, color=color)

    ax.set_xticks(range(n))
    ax.set_xticklabels(domain_list, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(n))
    ax.set_yticklabels(domain_list, fontsize=9)
    for i, d in enumerate(domain_list):
        c = DOMAIN_COLORS.get(d, 'k')
        ax.get_xticklabels()[i].set_color(c)
        ax.get_yticklabels()[i].set_color(c)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Δ PAE (Å) — positive = worse with KH7a', fontsize=10)

    ax.set_title(f'Effect of KH7a+KHb-KH8 on Inter-Domain PAE\n'
                 f'({len(kh7a_pairs)} matched pairs, mean delta)',
                 fontsize=12)

    fig.tight_layout()
    out = FIG_DIR / 'delta_pae_heatmap_kh7a.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ============================================================
# Figure 2: Delta PAE heatmaps — macrodomain removal (3-panel)
# ============================================================

def fig_delta_pae_heatmap_macrodomains(md_pairs_dict):
    """3-panel figure: delta PAE for removing MD1L1, MD2, MD3."""
    targets = ['md1l1', 'md2', 'md3']
    labels = ['Remove MD1L1', 'Remove MD2', 'Remove MD3']

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Common domain list across all comparisons
    all_shared = set()
    for tgt in targets:
        for p in md_pairs_dict.get(tgt, []):
            all_shared.update(p.get('shared_domains', []))
    domain_list = [d for d in DOMAIN_ORDER if d in all_shared]

    if len(domain_list) < 2:
        print("  Skip: too few shared domains for macrodomain heatmaps")
        plt.close(fig)
        return

    # Find global vmax for consistent colorbar
    all_vals = []
    matrices = {}
    count_matrices = {}
    for tgt in targets:
        pairs = md_pairs_dict.get(tgt, [])
        if pairs:
            m, c = build_delta_matrix(pairs, domain_list)
            matrices[tgt] = m
            count_matrices[tgt] = c
            finite = m[np.isfinite(m)]
            if len(finite) > 0:
                all_vals.extend(finite.tolist())

    vmax = max(np.max(np.abs(all_vals)), 1.0) if all_vals else 5.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    for ax_idx, (tgt, label) in enumerate(zip(targets, labels)):
        ax = axes[ax_idx]
        matrix = matrices.get(tgt, np.full((len(domain_list), len(domain_list)), np.nan))
        counts = count_matrices.get(tgt, np.zeros((len(domain_list), len(domain_list)), dtype=int))
        n_pairs = len(md_pairs_dict.get(tgt, []))
        n = len(domain_list)

        im = ax.imshow(matrix, cmap='RdBu_r', norm=norm, aspect='equal')

        for i in range(n):
            for j in range(n):
                val = matrix[i, j]
                cnt = counts[i, j]
                if np.isnan(val):
                    ax.text(j, i, '—', ha='center', va='center', fontsize=5, color='#999999')
                else:
                    color = 'white' if abs(val) > vmax * 0.6 else 'black'
                    ax.text(j, i, f'{val:+.1f}', ha='center', va='center',
                            fontsize=5, color=color)

        ax.set_xticks(range(n))
        ax.set_xticklabels(domain_list, rotation=45, ha='right', fontsize=7)
        ax.set_yticks(range(n))
        ax.set_yticklabels(domain_list, fontsize=7)
        for i, d in enumerate(domain_list):
            c = DOMAIN_COLORS.get(d, 'k')
            ax.get_xticklabels()[i].set_color(c)
            ax.get_yticklabels()[i].set_color(c)

        # Highlight the removed domain row/col
        tgt_idx = domain_list.index(tgt) if tgt in domain_list else None
        if tgt_idx is not None:
            ax.axhline(tgt_idx - 0.5, color='red', lw=1.5, ls='--', alpha=0.5)
            ax.axhline(tgt_idx + 0.5, color='red', lw=1.5, ls='--', alpha=0.5)
            ax.axvline(tgt_idx - 0.5, color='red', lw=1.5, ls='--', alpha=0.5)
            ax.axvline(tgt_idx + 0.5, color='red', lw=1.5, ls='--', alpha=0.5)

        ax.set_title(f'{label}\n({n_pairs} pairs)', fontsize=10, fontweight='bold')

    # Shared colorbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Δ PAE (Å) — positive = worse when domain present', fontsize=9)

    fig.suptitle('Effect of Macrodomain Removal on Inter-Domain PAE', fontsize=13, fontweight='bold')

    out = FIG_DIR / 'delta_pae_heatmap_macrodomains.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ============================================================
# Figure 3: Delta pLDDT distributions (violin plots)
# ============================================================

def fig_plddt_distributions(all_comparisons):
    """Violin/box plots of per-domain delta pLDDT across matched pairs."""
    comp_names = list(all_comparisons.keys())
    n_comp = len(comp_names)

    fig, axes = plt.subplots(1, n_comp, figsize=(4 * n_comp, 6), sharey=True)
    if n_comp == 1:
        axes = [axes]

    for ax_idx, (comp_name, pairs) in enumerate(all_comparisons.items()):
        ax = axes[ax_idx]
        deltas = collect_domain_deltas(pairs, 'mean_plddt')

        # Order domains canonically, filter to those with data
        ordered = [d for d in DOMAIN_ORDER if d in deltas and len(deltas[d]) >= 3]
        if not ordered:
            ax.set_title(f'{comp_name}\n(no data)')
            continue

        data = [deltas[d] for d in ordered]
        positions = range(len(ordered))
        colors = [DOMAIN_COLORS.get(d, '#999999') for d in ordered]

        parts = ax.violinplot(data, positions=positions, showmeans=True,
                               showmedians=True, widths=0.7)
        for pc, color in zip(parts['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(0.5)
        for key in ['cmeans', 'cmedians', 'cbars', 'cmins', 'cmaxes']:
            if key in parts:
                parts[key].set_color('black')
                parts[key].set_linewidth(0.8)

        # Overlay individual points
        for i, (d, vals) in enumerate(zip(ordered, data)):
            jitter = np.random.default_rng(42).normal(0, 0.05, len(vals))
            ax.scatter(i + jitter, vals, c=colors[i], s=8, alpha=0.3, zorder=5,
                       edgecolors='none')

        ax.axhline(0, color='black', lw=0.8, ls='--', alpha=0.5)
        ax.set_xticks(positions)
        ax.set_xticklabels(ordered, rotation=45, ha='right', fontsize=8)
        for i, d in enumerate(ordered):
            ax.get_xticklabels()[i].set_color(DOMAIN_COLORS.get(d, 'k'))

        # Add significance stars
        for i, d in enumerate(ordered):
            p = wilcoxon_p(deltas[d])
            stars = sig_stars(p)
            if stars and stars != 'ns':
                ymax = max(deltas[d])
                ax.text(i, ymax + 0.5, stars, ha='center', va='bottom', fontsize=8,
                        fontweight='bold', color='red')

        comp_color = COMPARISON_COLORS.get(comp_name, '#333333')
        ax.set_title(f'+{comp_name}\n({len(pairs)} pairs)', fontsize=10,
                     fontweight='bold', color=comp_color)
        if ax_idx == 0:
            ax.set_ylabel('Δ pLDDT (with − without)', fontsize=10)

    fig.suptitle('Per-Domain pLDDT Change by Construct Comparison',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    out = FIG_DIR / 'delta_plddt_distributions.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ============================================================
# Figure 4: Grouped bar chart — mean delta per domain
# ============================================================

def fig_domain_bar_chart(all_comparisons):
    """Grouped bars: mean delta pLDDT and delta intra-PAE per domain per comparison."""
    comp_names = list(all_comparisons.keys())
    n_comp = len(comp_names)

    # Collect all domains that appear in any comparison
    all_domains = set()
    for pairs in all_comparisons.values():
        for p in pairs:
            all_domains.update(p.get('shared_domains', []))
    domain_list = [d for d in DOMAIN_ORDER if d in all_domains]
    n_dom = len(domain_list)

    if n_dom == 0:
        print("  Skip: no shared domains for bar chart")
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(14, n_dom * 1.5), 8), sharex=True)

    bar_width = 0.8 / n_comp
    x = np.arange(n_dom)

    for i_comp, (comp_name, pairs) in enumerate(all_comparisons.items()):
        plddt_deltas = collect_domain_deltas(pairs, 'mean_plddt')
        pae_deltas = collect_domain_deltas(pairs, 'intra_pae')
        color = COMPARISON_COLORS.get(comp_name, f'C{i_comp}')

        plddt_means = []
        plddt_sems = []
        pae_means = []
        pae_sems = []

        for d in domain_list:
            pv = plddt_deltas.get(d, [])
            av = pae_deltas.get(d, [])
            plddt_means.append(np.mean(pv) if pv else 0)
            plddt_sems.append(np.std(pv) / np.sqrt(len(pv)) if len(pv) > 1 else 0)
            pae_means.append(np.mean(av) if av else 0)
            pae_sems.append(np.std(av) / np.sqrt(len(av)) if len(av) > 1 else 0)

        offset = (i_comp - n_comp / 2 + 0.5) * bar_width
        ax1.bar(x + offset, plddt_means, bar_width, yerr=plddt_sems,
                color=color, alpha=0.7, edgecolor='k', lw=0.5,
                label=comp_name, capsize=2)
        ax2.bar(x + offset, pae_means, bar_width, yerr=pae_sems,
                color=color, alpha=0.7, edgecolor='k', lw=0.5,
                label=comp_name, capsize=2)

    ax1.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
    ax2.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)

    ax1.set_ylabel('Δ pLDDT', fontsize=11)
    ax2.set_ylabel('Δ intra-PAE (Å)', fontsize=11)
    ax2.set_xlabel('Domain', fontsize=11)

    ax2.set_xticks(x)
    ax2.set_xticklabels(domain_list, rotation=45, ha='right', fontsize=9)
    for i, d in enumerate(domain_list):
        ax2.get_xticklabels()[i].set_color(DOMAIN_COLORS.get(d, 'k'))

    ax1.legend(fontsize=8, loc='best', framealpha=0.9)
    ax1.set_title('Mean Δ pLDDT per Domain (with − without target domain)', fontsize=11)
    ax2.set_title('Mean Δ Intra-PAE per Domain', fontsize=11)

    fig.suptitle('Domain-Level Impact of Adding/Removing Domains',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    out = FIG_DIR / 'delta_domain_barchart.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ============================================================
# Figure 5: Scatter plots — pLDDT with vs without
# ============================================================

def fig_scatter_with_vs_without(all_comparisons):
    """2×2 scatter: pLDDT(with) vs pLDDT(without) per comparison."""
    comp_names = list(all_comparisons.keys())[:4]  # max 4 panels
    n_panels = len(comp_names)
    nrows = (n_panels + 1) // 2
    ncols = min(n_panels, 2)

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5.5 * nrows), squeeze=False)

    for idx, comp_name in enumerate(comp_names):
        ax = axes[idx // ncols][idx % ncols]
        pairs = all_comparisons[comp_name]

        with_vals = []
        without_vals = []
        n_domains_list = []

        for p in pairs:
            w = p.get('with_global_mean_plddt')
            wo = p.get('without_global_mean_plddt')
            if w is not None and wo is not None:
                with_vals.append(w)
                without_vals.append(wo)
                n_domains_list.append(len(p.get('with_domains', [])))

        if not with_vals:
            ax.set_title(f'{comp_name} (no data)')
            continue

        with_vals = np.array(with_vals)
        without_vals = np.array(without_vals)
        n_domains_arr = np.array(n_domains_list)

        sc = ax.scatter(without_vals, with_vals, c=n_domains_arr, cmap='viridis',
                        s=20, alpha=0.6, edgecolors='k', lw=0.3)

        # Diagonal
        lims = [min(without_vals.min(), with_vals.min()) - 2,
                max(without_vals.max(), with_vals.max()) + 2]
        ax.plot(lims, lims, 'k--', lw=0.8, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        ax.set_xlabel(f'pLDDT (without {comp_name})', fontsize=9)
        ax.set_ylabel(f'pLDDT (with {comp_name})', fontsize=9)
        ax.set_aspect('equal')

        # Count above/below diagonal
        n_above = np.sum(with_vals > without_vals)
        n_below = np.sum(with_vals < without_vals)
        comp_color = COMPARISON_COLORS.get(comp_name, '#333333')
        ax.set_title(f'+{comp_name} ({len(pairs)} pairs)\n'
                     f'↑{n_above} better, ↓{n_below} worse',
                     fontsize=10, fontweight='bold', color=comp_color)

        cbar = fig.colorbar(sc, ax=ax, shrink=0.7)
        cbar.set_label('# domains', fontsize=8)

    # Hide unused axes
    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle('Global pLDDT: With vs Without Target Domain(s)',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    out = FIG_DIR / 'scatter_plddt_with_vs_without.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ============================================================
# Figure 6: Combined MD effect (all MDs vs no MDs)
# ============================================================

def fig_combined_md_effect(all_md_pairs):
    """All-MDs vs no-MDs: bar chart + heatmap side by side."""
    if not all_md_pairs:
        print("  Skip: no all-MD pairs")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                     gridspec_kw={'width_ratios': [1, 1.2]})

    # Panel 1: bar chart of delta pLDDT per remaining domain
    plddt_deltas = collect_domain_deltas(all_md_pairs, 'mean_plddt')
    pae_deltas = collect_domain_deltas(all_md_pairs, 'intra_pae')
    domain_list = [d for d in DOMAIN_ORDER if d in plddt_deltas and len(plddt_deltas[d]) >= 3]

    if domain_list:
        x = np.arange(len(domain_list))
        means = [np.mean(plddt_deltas[d]) for d in domain_list]
        sems = [np.std(plddt_deltas[d]) / np.sqrt(len(plddt_deltas[d])) for d in domain_list]
        colors = [DOMAIN_COLORS.get(d, '#999999') for d in domain_list]

        ax1.bar(x, means, yerr=sems, color=colors, edgecolor='k', lw=0.5, capsize=3)
        ax1.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
        ax1.set_xticks(x)
        ax1.set_xticklabels(domain_list, rotation=45, ha='right', fontsize=9)
        for i, d in enumerate(domain_list):
            # Significance
            p = wilcoxon_p(plddt_deltas[d])
            stars = sig_stars(p)
            if stars and stars != 'ns':
                yval = means[i] + sems[i] if means[i] > 0 else means[i] - sems[i]
                va = 'bottom' if means[i] > 0 else 'top'
                ax1.text(i, yval, stars, ha='center', va=va, fontsize=9,
                         fontweight='bold', color='red')

        ax1.set_ylabel('Δ pLDDT (with all MDs − without)', fontsize=10)
        ax1.set_title(f'pLDDT Change (n={len(all_md_pairs)} pairs)', fontsize=11, fontweight='bold')

    # Panel 2: delta PAE heatmap
    all_shared = set()
    for p in all_md_pairs:
        all_shared.update(p.get('shared_domains', []))
    heatmap_domains = [d for d in DOMAIN_ORDER if d in all_shared]

    if len(heatmap_domains) >= 2:
        matrix, counts = build_delta_matrix(all_md_pairs, heatmap_domains)
        n = len(heatmap_domains)
        vmax = np.nanmax(np.abs(matrix[np.isfinite(matrix)])) if np.any(np.isfinite(matrix)) else 5
        vmax = max(vmax, 1.0)
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

        im = ax2.imshow(matrix, cmap='RdBu_r', norm=norm, aspect='equal')
        for i in range(n):
            for j in range(n):
                val = matrix[i, j]
                if not np.isnan(val):
                    color = 'white' if abs(val) > vmax * 0.6 else 'black'
                    ax2.text(j, i, f'{val:+.1f}', ha='center', va='center',
                             fontsize=7, color=color)

        ax2.set_xticks(range(n))
        ax2.set_xticklabels(heatmap_domains, rotation=45, ha='right', fontsize=8)
        ax2.set_yticks(range(n))
        ax2.set_yticklabels(heatmap_domains, fontsize=8)
        for i, d in enumerate(heatmap_domains):
            c = DOMAIN_COLORS.get(d, 'k')
            ax2.get_xticklabels()[i].set_color(c)
            ax2.get_yticklabels()[i].set_color(c)

        cbar = fig.colorbar(im, ax=ax2, shrink=0.8, pad=0.02)
        cbar.set_label('Δ PAE (Å)', fontsize=9)
        ax2.set_title('Δ Inter-Domain PAE', fontsize=11, fontweight='bold')

    fig.suptitle('Combined Effect: All Macrodomains (MD1L1+MD2+MD3) vs None',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    out = FIG_DIR / 'combined_md_effect.svg'
    fig.savefig(out, format='svg', bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ============================================================
# Main
# ============================================================

def main():
    if not (METRICS_DIR / 'index.json').exists():
        print(f"ERROR: {METRICS_DIR / 'index.json'} not found.")
        print("Run 01_compute_domain_metrics.py first:")
        print("  cd /home/sbali/CALVADOS/parp14/analysis")
        print("  python 01_compute_domain_metrics.py --workers 8")
        return

    print("Loading index...")
    index = load_index(METRICS_DIR)
    print(f"  {len(index)} combinations available")

    # ---- Find all matched pairs ----
    print("\nFinding matched pairs...")

    kh7a_pairs_raw = find_pairs(index, {'kh7a', 'khb-kh8'})
    md1l1_pairs_raw = find_pairs(index, {'md1l1'})
    md2_pairs_raw = find_pairs(index, {'md2'})
    md3_pairs_raw = find_pairs(index, {'md3'})
    all_md_pairs_raw = find_pairs(index, {'md1l1', 'md2', 'md3'})

    print(f"  KH7a+KHb-KH8: {len(kh7a_pairs_raw)} pairs")
    print(f"  MD1L1:         {len(md1l1_pairs_raw)} pairs")
    print(f"  MD2:           {len(md2_pairs_raw)} pairs")
    print(f"  MD3:           {len(md3_pairs_raw)} pairs")
    print(f"  All MDs:       {len(all_md_pairs_raw)} pairs")

    # ---- Compute deltas ----
    print("\nComputing deltas...")

    def compute_all_deltas(pairs_raw):
        results = []
        for w_name, wo_name in pairs_raw:
            w_data = load_combination(METRICS_DIR, w_name)
            wo_data = load_combination(METRICS_DIR, wo_name)
            results.append(compute_pair_deltas(w_data, wo_data))
        return results

    kh7a_pairs = compute_all_deltas(kh7a_pairs_raw)
    md1l1_pairs = compute_all_deltas(md1l1_pairs_raw)
    md2_pairs = compute_all_deltas(md2_pairs_raw)
    md3_pairs = compute_all_deltas(md3_pairs_raw)
    all_md_pairs = compute_all_deltas(all_md_pairs_raw)

    # ---- Print summary stats ----
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)

    for label, pairs in [
        ('KH7a+KHb-KH8', kh7a_pairs),
        ('MD1L1', md1l1_pairs),
        ('MD2', md2_pairs),
        ('MD3', md3_pairs),
        ('All MDs', all_md_pairs),
    ]:
        if not pairs:
            print(f"\n  {label}: no pairs")
            continue
        print(f"\n  {label} ({len(pairs)} pairs):")
        for metric in ['global_mean_plddt', 'ptm', 'ranking_score', 'global_mean_pae']:
            vals = [p.get(f'delta_{metric}') for p in pairs]
            vals = [v for v in vals if v is not None]
            if vals:
                p_val = wilcoxon_p(vals)
                p_str = f'p={p_val:.2e} {sig_stars(p_val)}' if p_val else 'n<10'
                print(f"    Δ {metric:25s}: mean={np.mean(vals):+.2f}, "
                      f"median={np.median(vals):+.2f}, sd={np.std(vals):.2f}, "
                      f"n={len(vals)}, {p_str}")

        # Per-domain pLDDT deltas for key neighbors
        plddt_deltas = collect_domain_deltas(pairs, 'mean_plddt')
        for d in DOMAIN_ORDER:
            if d in plddt_deltas and len(plddt_deltas[d]) >= 3:
                vals = plddt_deltas[d]
                p_val = wilcoxon_p(vals)
                p_str = f'p={p_val:.2e} {sig_stars(p_val)}' if p_val else 'n<10'
                print(f"      {d:12s} Δ pLDDT: mean={np.mean(vals):+.2f}, "
                      f"n={len(vals)}, {p_str}")

    # ---- Generate figures ----
    print("\n" + "=" * 70)
    print("GENERATING FIGURES")
    print("=" * 70)

    all_comparisons = {}
    if kh7a_pairs:
        all_comparisons['kh7a'] = kh7a_pairs
    if md1l1_pairs:
        all_comparisons['md1l1'] = md1l1_pairs
    if md2_pairs:
        all_comparisons['md2'] = md2_pairs
    if md3_pairs:
        all_comparisons['md3'] = md3_pairs

    md_pairs_dict = {
        'md1l1': md1l1_pairs,
        'md2': md2_pairs,
        'md3': md3_pairs,
    }

    print("\nFigure 1: Delta PAE heatmap — KH7a effect")
    if kh7a_pairs:
        fig_delta_pae_heatmap_kh7a(kh7a_pairs)
    else:
        print("  Skip: no KH7a pairs")

    print("\nFigure 2: Delta PAE heatmaps — macrodomain removal")
    fig_delta_pae_heatmap_macrodomains(md_pairs_dict)

    print("\nFigure 3: Delta pLDDT distributions")
    if all_comparisons:
        fig_plddt_distributions(all_comparisons)

    print("\nFigure 4: Domain bar chart")
    if all_comparisons:
        fig_domain_bar_chart(all_comparisons)

    print("\nFigure 5: Scatter with vs without")
    if all_comparisons:
        fig_scatter_with_vs_without(all_comparisons)

    print("\nFigure 6: Combined MD effect")
    fig_combined_md_effect(all_md_pairs)

    print(f"\nAll figures saved to: {FIG_DIR}/")
    print("Done.")


if __name__ == '__main__':
    main()
