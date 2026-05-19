#!/usr/bin/env python3
"""
Analyze AF3 PAE matrices and pLDDT to determine stable domain boundaries
for KH domains (KH1-KH6, KH7a, KHb-KH8) in PARP14.

The current FL domains.yaml omits KH domains entirely (residues 315-789
and 1389-1533), treating them as disordered. This script uses AF3
predicted aligned error (PAE) to identify structured sub-domains within
the KH regions that should be restrained in CALVADOS simulations.

Strategy:
1. Load PAE matrices from multiple AF3 constructs containing KH domains
2. Compute per-residue pLDDT and intra-block PAE
3. Identify tightly-coupled residue blocks (low PAE) within KH regions
4. Compare with existing restrained domains for validation
5. Output proposed domain boundaries

Constructs analyzed:
- kh1-kh6 (isolated, 423 res)
- kh1-kh6_kh7a (475 res)
- khb-kh8 (isolated, 145 res)
- kh1-kh6_kh7a_md1l1_md2_md3_khb-kh8_wwe_art (norrm-equivalent, 1474 res)
"""

import json
import numpy as np
import os
from pathlib import Path

AF3_BASE = Path('/home/sbali/CALVADOS/parp14/alphafold_outputs')
OUTPUT_DIR = Path('/home/sbali/CALVADOS/examples/PARP14_MDP/data')

# Full-length PARP14 domain boundaries (1-indexed, inclusive)
DOMAIN_DEFS = {
    'rrm1':    (1, 145),
    'rrm2':    (146, 224),
    'rrm3':    (225, 314),
    'kh1-kh6': (315, 737),
    'kh7a':    (738, 789),
    'md1l1':   (790, 1004),
    'md2':     (1004, 1193),
    'md3':     (1207, 1388),
    'khb-kh8': (1389, 1533),
    'wwe':     (1534, 1602),
    'art':     (1603, 1801),
}

# Current FL restraint domains (from input/domains.yaml) - 0-indexed
CURRENT_DOMAINS = [
    (6, 88),       # RRM1 core
    (150, 223),    # RRM2 core
    (227, 301),    # RRM3 core
    (791, 978),    # MD1 core
    (1003, 1190),  # MD2
    (1216, 1387),  # MD3
    (1523, 1601),  # WWE
    (1605, 1801),  # ART
]

# Constructs to analyze
CONSTRUCTS = {
    'kh1-kh6': {
        'dir': 'kh1-kh6',
        'domains': ['kh1-kh6'],
        'description': 'Isolated KH1-KH6 (423 residues)',
    },
    'kh1-kh6_kh7a': {
        'dir': 'kh1-kh6_kh7a',
        'domains': ['kh1-kh6', 'kh7a'],
        'description': 'KH1-KH6 + KH7a (475 residues)',
    },
    'khb-kh8': {
        'dir': 'khb-kh8',
        'domains': ['khb-kh8'],
        'description': 'Isolated KHb-KH8 (145 residues)',
    },
    'norrm': {
        'dir': 'kh1-kh6_kh7a_md1l1_md2_md3_khb-kh8_wwe_art',
        'domains': ['kh1-kh6', 'kh7a', 'md1l1', 'md2', 'md3', 'khb-kh8', 'wwe', 'art'],
        'description': 'No-RRM construct (1474 residues)',
    },
}


def compute_construct_boundaries(domains):
    """Compute 0-indexed (start, end_exclusive) for each domain in a concatenated construct."""
    boundaries = {}
    offset = 0
    prev_end_fl = None
    for dname in domains:
        fl_start, fl_end = DOMAIN_DEFS[dname]
        n_residues = fl_end - fl_start + 1
        if prev_end_fl is not None and fl_start <= prev_end_fl:
            overlap = prev_end_fl - fl_start + 1
            n_residues -= overlap
        boundaries[dname] = (offset, offset + n_residues)
        offset += n_residues
        prev_end_fl = fl_end
    return boundaries, offset


def extract_plddt(conf):
    """Extract per-residue pLDDT from confidences.json."""
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


def load_all_models(construct_dir):
    """Load PAE and pLDDT from all seed-sample models."""
    models = []
    for seed in range(1, 6):
        for sample in range(5):
            sd = construct_dir / f'seed-{seed}_sample-{sample}'
            cf = sd / 'confidences.json'
            if not cf.exists():
                continue
            with open(cf) as f:
                conf = json.load(f)
            plddt = extract_plddt(conf)
            pae = np.array(conf['pae']) if 'pae' in conf else None
            if plddt is not None and pae is not None:
                models.append({'plddt': plddt, 'pae': pae, 'seed': seed, 'sample': sample})
    return models


def find_structured_blocks(pae_matrix, plddt, window=10, pae_threshold=5.0, plddt_threshold=60.0):
    """Identify contiguous blocks of structured residues using PAE and pLDDT.

    A residue is 'structured' if:
    1. Its pLDDT > plddt_threshold
    2. Mean intra-block PAE within a sliding window is < pae_threshold

    Returns list of (start, end) 0-indexed inclusive ranges.
    """
    n = len(plddt)

    # Per-residue mean local PAE (average PAE to nearby residues within window)
    local_pae = np.zeros(n)
    for i in range(n):
        lo = max(0, i - window)
        hi = min(n, i + window + 1)
        local_pae[i] = np.mean(pae_matrix[i, lo:hi])

    # Structured mask
    structured = (plddt > plddt_threshold) & (local_pae < pae_threshold)

    # Find contiguous blocks
    blocks = []
    in_block = False
    start = 0
    for i in range(n):
        if structured[i] and not in_block:
            start = i
            in_block = True
        elif not structured[i] and in_block:
            if i - start >= 10:  # minimum block size
                blocks.append((start, i - 1))
            in_block = False
    if in_block and n - start >= 10:
        blocks.append((start, n - 1))

    return blocks, structured, local_pae


def analyze_pae_blocks(pae_matrix, block_ranges):
    """Compute mean PAE within and between identified blocks."""
    n_blocks = len(block_ranges)
    intra = []
    inter = np.zeros((n_blocks, n_blocks))

    for i, (s1, e1) in enumerate(block_ranges):
        intra_pae = pae_matrix[s1:e1+1, s1:e1+1]
        intra.append(np.mean(intra_pae))
        for j, (s2, e2) in enumerate(block_ranges):
            inter[i, j] = np.mean(pae_matrix[s1:e1+1, s2:e2+1])

    return intra, inter


def map_construct_to_fl(block_start, block_end, domain_name, domain_offset_in_construct):
    """Map construct-local block indices to full-length residue numbers."""
    fl_start, fl_end = DOMAIN_DEFS[domain_name]
    # block positions are relative to construct, domain starts at domain_offset_in_construct
    fl_block_start = fl_start + (block_start - domain_offset_in_construct)
    fl_block_end = fl_start + (block_end - domain_offset_in_construct)
    return fl_block_start, fl_block_end


def print_separator(char='=', width=80):
    print(char * width)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = {}

    for name, info in CONSTRUCTS.items():
        construct_dir = AF3_BASE / info['dir']
        if not construct_dir.exists():
            print(f"SKIP {name}: directory not found")
            continue

        print_separator()
        print(f"CONSTRUCT: {name}")
        print(f"  {info['description']}")
        print(f"  Domains: {info['domains']}")

        boundaries, total_res = compute_construct_boundaries(info['domains'])
        print(f"  Total residues: {total_res}")
        for dname, (s, e) in boundaries.items():
            fl_s, fl_e = DOMAIN_DEFS[dname]
            print(f"    {dname}: construct [{s}:{e}) = FL residues {fl_s}-{fl_e}")

        # Load all models
        models = load_all_models(construct_dir)
        print(f"  Models loaded: {len(models)}")

        if not models:
            print("  NO DATA - skipping")
            continue

        # Average PAE and pLDDT across models
        pae_matrices = [m['pae'] for m in models]
        plddt_arrays = [m['plddt'] for m in models]

        # Check consistent sizes
        sizes = set(len(p) for p in plddt_arrays)
        pae_sizes = set(p.shape for p in pae_matrices)
        print(f"  pLDDT sizes: {sizes}, PAE shapes: {pae_sizes}")

        # Use most common size
        target_size = max(sizes, key=lambda s: sum(1 for p in plddt_arrays if len(p) == s))
        valid = [(m['pae'], m['plddt']) for m in models if len(m['plddt']) == target_size]

        mean_pae = np.mean([v[0] for v in valid], axis=0)
        mean_plddt = np.mean([v[1] for v in valid], axis=0)
        std_plddt = np.std([v[1] for v in valid], axis=0)

        print(f"\n  Global metrics (averaged over {len(valid)} models):")
        print(f"    Mean pLDDT: {np.mean(mean_plddt):.1f}")
        print(f"    Mean PAE: {np.mean(mean_pae):.1f}")

        # Per-domain metrics
        print(f"\n  Per-domain metrics:")
        for dname, (s, e) in boundaries.items():
            if e <= len(mean_plddt):
                d_plddt = mean_plddt[s:e]
                d_pae = mean_pae[s:e, s:e]
                print(f"    {dname:12s}: pLDDT={np.mean(d_plddt):5.1f} ± {np.std(d_plddt):4.1f}, "
                      f"intra-PAE={np.mean(d_pae):5.1f}, "
                      f"min_pLDDT={np.min(d_plddt):5.1f}, "
                      f"FL range: {DOMAIN_DEFS[dname][0]}-{DOMAIN_DEFS[dname][1]}")

        # Find structured blocks within KH regions
        kh_domains = [d for d in info['domains'] if 'kh' in d.lower()]

        print(f"\n  === Structured block analysis for KH domains ===")

        construct_results = {}

        for dname in kh_domains:
            s, e = boundaries[dname]
            if e > len(mean_plddt):
                continue

            d_plddt = mean_plddt[s:e]
            d_pae = mean_pae[s:e, s:e]

            print(f"\n  --- {dname} (construct [{s}:{e}), FL {DOMAIN_DEFS[dname][0]}-{DOMAIN_DEFS[dname][1]}) ---")

            # Find structured blocks at different thresholds
            for pae_thresh in [4.0, 5.0, 6.0, 8.0]:
                blocks, structured, local_pae = find_structured_blocks(
                    d_pae, d_plddt, window=10, pae_threshold=pae_thresh, plddt_threshold=60.0
                )

                if blocks:
                    print(f"    PAE < {pae_thresh}, pLDDT > 60:")
                    for bs, be in blocks:
                        fl_bs = DOMAIN_DEFS[dname][0] + bs
                        fl_be = DOMAIN_DEFS[dname][0] + be
                        block_plddt = d_plddt[bs:be+1]
                        block_pae = d_pae[bs:be+1, bs:be+1]
                        print(f"      Block [{bs}-{be}] = FL [{fl_bs}-{fl_be}] "
                              f"({be-bs+1} res): pLDDT={np.mean(block_plddt):.1f}, "
                              f"intra-PAE={np.mean(block_pae):.1f}")

            # Also try with pLDDT > 50 for more permissive
            blocks_permissive, _, _ = find_structured_blocks(
                d_pae, d_plddt, window=10, pae_threshold=6.0, plddt_threshold=50.0
            )
            if blocks_permissive:
                print(f"    PAE < 6.0, pLDDT > 50 (permissive):")
                for bs, be in blocks_permissive:
                    fl_bs = DOMAIN_DEFS[dname][0] + bs
                    fl_be = DOMAIN_DEFS[dname][0] + be
                    block_plddt = d_plddt[bs:be+1]
                    block_pae = d_pae[bs:be+1, bs:be+1]
                    print(f"      Block [{bs}-{be}] = FL [{fl_bs}-{fl_be}] "
                          f"({be-bs+1} res): pLDDT={np.mean(block_plddt):.1f}, "
                          f"intra-PAE={np.mean(block_pae):.1f}")

            # Store the moderate threshold results
            blocks_moderate, _, local_pae_arr = find_structured_blocks(
                d_pae, d_plddt, window=10, pae_threshold=5.0, plddt_threshold=60.0
            )
            construct_results[dname] = {
                'mean_plddt': d_plddt,
                'mean_pae': d_pae,
                'local_pae': local_pae_arr,
                'blocks': blocks_moderate,
            }

        # Inter-domain PAE between KH and non-KH domains
        print(f"\n  === Inter-domain PAE matrix ===")
        domain_names = list(boundaries.keys())
        header = f"{'':>12s}"
        for d in domain_names:
            header += f" {d:>10s}"
        print(f"    {header}")
        for d1 in domain_names:
            s1, e1 = boundaries[d1]
            row = f"    {d1:>12s}"
            for d2 in domain_names:
                s2, e2 = boundaries[d2]
                if e1 <= mean_pae.shape[0] and e2 <= mean_pae.shape[1]:
                    val = np.mean(mean_pae[s1:e1, s2:e2])
                    row += f" {val:10.1f}"
            print(row)

        all_results[name] = {
            'boundaries': boundaries,
            'mean_pae': mean_pae,
            'mean_plddt': mean_plddt,
            'kh_results': construct_results,
        }

        # Save per-residue data
        np.savez(OUTPUT_DIR / f'kh_pae_analysis_{name}.npz',
                 mean_pae=mean_pae,
                 mean_plddt=mean_plddt,
                 boundaries=json.dumps({k: list(v) for k, v in boundaries.items()}))

    # ============================================================
    # Cross-construct comparison for KH domains
    # ============================================================
    print_separator()
    print("CROSS-CONSTRUCT COMPARISON")
    print_separator()

    # Compare KH1-KH6 blocks across isolated vs norrm context
    for dname in ['kh1-kh6', 'kh7a', 'khb-kh8']:
        print(f"\n  {dname} structured blocks across constructs:")
        for cname, cdata in all_results.items():
            if dname in cdata.get('kh_results', {}):
                blocks = cdata['kh_results'][dname]['blocks']
                if blocks:
                    for bs, be in blocks:
                        fl_bs = DOMAIN_DEFS[dname][0] + bs
                        fl_be = DOMAIN_DEFS[dname][0] + be
                        d_plddt = cdata['kh_results'][dname]['mean_plddt']
                        d_pae = cdata['kh_results'][dname]['mean_pae']
                        print(f"    {cname:20s}: FL [{fl_bs:4d}-{fl_be:4d}] ({be-bs+1:3d} res) "
                              f"pLDDT={np.mean(d_plddt[bs:be+1]):.1f} "
                              f"PAE={np.mean(d_pae[bs:be+1, bs:be+1]):.1f}")
                else:
                    print(f"    {cname:20s}: no structured blocks found")

    # ============================================================
    # Propose MULTIPLE boundary sets at different thresholds
    # ============================================================
    print_separator()
    print("PROPOSED DOMAIN BOUNDARY SETS (multiple thresholds to test)")
    print_separator()

    print("\nCurrent FL domains.yaml (0-indexed):")
    for s, e in CURRENT_DOMAINS:
        print(f"  [{s}, {e}]")

    print("\nGap analysis (unrestrained regions):")
    gaps = []
    sorted_domains = sorted(CURRENT_DOMAINS)
    for i in range(len(sorted_domains) - 1):
        gap_start = sorted_domains[i][1] + 1
        gap_end = sorted_domains[i+1][0] - 1
        if gap_end > gap_start:
            region = "unknown"
            if gap_start >= 301 and gap_end <= 791:
                region = "KH1-KH6 + KH7a"
            elif gap_start >= 1190 and gap_end <= 1216:
                region = "MD2-MD3 linker"
            elif gap_start >= 1387 and gap_end <= 1523:
                region = "KHb-KH8"
            print(f"  Gap: FL [{gap_start}-{gap_end}] ({gap_end-gap_start+1} res) = {region}")
            gaps.append((gap_start, gap_end, region))

    best_construct = 'norrm' if 'norrm' in all_results else 'kh1-kh6'
    if best_construct not in all_results:
        best_construct = list(all_results.keys())[0] if all_results else None

    if not best_construct:
        print("ERROR: no AF3 data available")
        return

    print(f"\nUsing construct: {best_construct}")

    # ---- Define threshold levels ----
    THRESHOLD_SETS = {
        'conservative': {
            'label': 'Conservative (tight cores only)',
            'pae_threshold': 4.0,
            'plddt_threshold': 70.0,
            'min_block': 15,
            'merge_gap': 0,  # no merging
        },
        'moderate': {
            'label': 'Moderate (well-structured blocks)',
            'pae_threshold': 5.0,
            'plddt_threshold': 60.0,
            'min_block': 10,
            'merge_gap': 0,
        },
        'moderate_merged': {
            'label': 'Moderate + merge nearby blocks (gap <= 5 res)',
            'pae_threshold': 5.0,
            'plddt_threshold': 60.0,
            'min_block': 10,
            'merge_gap': 5,
        },
        'permissive': {
            'label': 'Permissive (include marginal regions)',
            'pae_threshold': 6.0,
            'plddt_threshold': 50.0,
            'min_block': 10,
            'merge_gap': 0,
        },
        'permissive_merged': {
            'label': 'Permissive + merge nearby blocks (gap <= 10 res)',
            'pae_threshold': 6.0,
            'plddt_threshold': 50.0,
            'min_block': 10,
            'merge_gap': 10,
        },
        'aggressive': {
            'label': 'Aggressive (maximize restrained residues)',
            'pae_threshold': 8.0,
            'plddt_threshold': 50.0,
            'min_block': 10,
            'merge_gap': 15,
        },
    }

    def merge_blocks(blocks, max_gap):
        """Merge blocks separated by <= max_gap residues."""
        if not blocks or max_gap <= 0:
            return blocks
        merged = [list(blocks[0])]
        for bs, be in blocks[1:]:
            if bs - merged[-1][1] <= max_gap + 1:
                merged[-1][1] = be
            else:
                merged.append([bs, be])
        return [tuple(b) for b in merged]

    def get_kh_blocks_at_threshold(construct_data, pae_thresh, plddt_thresh, min_block, merge_gap):
        """Extract KH domain blocks at given thresholds from construct data."""
        kh_blocks_fl = []
        for dname in ['kh1-kh6', 'kh7a', 'khb-kh8']:
            if dname not in construct_data.get('kh_results', {}):
                continue
            kr = construct_data['kh_results'][dname]
            d_plddt = kr['mean_plddt']
            d_pae = kr['mean_pae']

            blocks, _, _ = find_structured_blocks(
                d_pae, d_plddt, window=10,
                pae_threshold=pae_thresh,
                plddt_threshold=plddt_thresh,
            )
            # Filter by min block size
            blocks = [(s, e) for s, e in blocks if e - s + 1 >= min_block]
            # Merge nearby blocks
            blocks = merge_blocks(blocks, merge_gap)
            # Map to FL numbering
            fl_offset = DOMAIN_DEFS[dname][0]
            for bs, be in blocks:
                fl_bs = fl_offset + bs
                fl_be = fl_offset + be
                # Compute quality metrics for this block
                block_plddt = d_plddt[bs:be+1]
                block_pae = d_pae[bs:be+1, bs:be+1]
                kh_blocks_fl.append({
                    'fl_start': fl_bs,
                    'fl_end': fl_be,
                    'n_res': be - bs + 1,
                    'source_domain': dname,
                    'mean_plddt': float(np.mean(block_plddt)),
                    'mean_pae': float(np.mean(block_pae)),
                })
        return kh_blocks_fl

    # ---- Generate all boundary sets ----
    boundary_sets = {}
    cdata = all_results[best_construct]

    for set_name, params in THRESHOLD_SETS.items():
        kh_blocks = get_kh_blocks_at_threshold(
            cdata,
            params['pae_threshold'],
            params['plddt_threshold'],
            params['min_block'],
            params['merge_gap'],
        )
        # Combine with existing non-KH domains
        all_domains = sorted(
            CURRENT_DOMAINS + [(b['fl_start'], b['fl_end']) for b in kh_blocks]
        )
        # Check for overlaps and merge if needed
        cleaned = []
        for s, e in all_domains:
            if cleaned and s <= cleaned[-1][1] + 1:
                cleaned[-1] = (cleaned[-1][0], max(cleaned[-1][1], e))
            else:
                cleaned.append((s, e))

        boundary_sets[set_name] = {
            'params': params,
            'kh_blocks': kh_blocks,
            'all_domains': cleaned,
        }

    # ---- Print comparison table ----
    print_separator('─')
    for set_name, bset in boundary_sets.items():
        params = bset['params']
        kh_blocks = bset['kh_blocks']
        all_domains = bset['all_domains']

        n_kh_res = sum(b['n_res'] for b in kh_blocks)
        n_total_res = sum(e - s + 1 for s, e in all_domains)
        n_domains = len(all_domains)

        print(f"\n  === {set_name.upper()}: {params['label']} ===")
        print(f"  Thresholds: PAE < {params['pae_threshold']}, pLDDT > {params['plddt_threshold']}, "
              f"min_block={params['min_block']}, merge_gap={params['merge_gap']}")
        print(f"  KH blocks: {len(kh_blocks)}, KH restrained residues: {n_kh_res}")
        print(f"  Total domains: {n_domains}, total restrained residues: {n_total_res}/1801 "
              f"({100*n_total_res/1801:.1f}%)")

        print(f"\n  New KH restraint blocks:")
        if kh_blocks:
            for b in kh_blocks:
                print(f"    [{b['fl_start']:4d}, {b['fl_end']:4d}]  "
                      f"({b['n_res']:3d} res, {b['source_domain']}) "
                      f"pLDDT={b['mean_plddt']:.1f}, PAE={b['mean_pae']:.1f}")
        else:
            print(f"    (none)")

        print(f"\n  Complete domains.yaml:")
        print(f"  parp14:")
        for s, e in all_domains:
            is_current = (s, e) in CURRENT_DOMAINS
            marker = "" if is_current else "  # NEW"
            print(f"   - [{s},{e}]{marker}")

        # Save each set
        yaml_lines = ["parp14:\n"]
        for s, e in all_domains:
            yaml_lines.append(f" - [{s},{e}]\n")
        proposed_path = OUTPUT_DIR / f'domains_{set_name}.yaml'
        with open(proposed_path, 'w') as f:
            f.writelines(yaml_lines)
        print(f"  Saved: {proposed_path}")

    # ---- Summary comparison table ----
    print_separator()
    print("SUMMARY COMPARISON TABLE")
    print_separator()
    print(f"\n  {'Set':<22s} {'PAE':<6s} {'pLDDT':<7s} {'merge':<7s} "
          f"{'#KH_blk':<8s} {'KH_res':<8s} {'#dom':<6s} {'tot_res':<9s} {'%restr':<7s}")
    print(f"  {'─'*22} {'─'*6} {'─'*7} {'─'*7} {'─'*8} {'─'*8} {'─'*6} {'─'*9} {'─'*7}")
    for set_name, bset in boundary_sets.items():
        p = bset['params']
        kb = bset['kh_blocks']
        ad = bset['all_domains']
        n_kh_res = sum(b['n_res'] for b in kb)
        n_tot = sum(e - s + 1 for s, e in ad)
        print(f"  {set_name:<22s} {p['pae_threshold']:<6.1f} {p['plddt_threshold']:<7.1f} "
              f"{p['merge_gap']:<7d} {len(kb):<8d} {n_kh_res:<8d} {len(ad):<6d} "
              f"{n_tot:<9d} {100*n_tot/1801:<7.1f}")

    # ---- Residue coverage comparison ----
    print(f"\n  Unrestrained gap sizes by set:")
    for set_name, bset in boundary_sets.items():
        ad = sorted(bset['all_domains'])
        gap_sizes = []
        for i in range(len(ad) - 1):
            gap = ad[i+1][0] - ad[i][1] - 1
            if gap > 0:
                gap_sizes.append((ad[i][1]+1, ad[i+1][0]-1, gap))
        total_gap = sum(g[2] for g in gap_sizes)
        big_gaps = [(s, e, sz) for s, e, sz in gap_sizes if sz > 10]
        print(f"    {set_name:<22s}: {total_gap:4d} unrestrained res, "
              f"{len(big_gaps)} gaps > 10 res: ", end="")
        if big_gaps:
            print(", ".join(f"[{s}-{e}]({sz})" for s, e, sz in big_gaps))
        else:
            print("(none)")

    # ============================================================
    # Detailed per-residue pLDDT profile for KH region
    # ============================================================
    print_separator()
    print("PER-RESIDUE pLDDT PROFILE (KH region, 20-residue windows)")
    print_separator()

    if best_construct and 'kh1-kh6' in all_results[best_construct].get('kh_results', {}):
        d_plddt = all_results[best_construct]['kh_results']['kh1-kh6']['mean_plddt']
        d_local_pae = all_results[best_construct]['kh_results']['kh1-kh6']['local_pae']
        fl_offset = DOMAIN_DEFS['kh1-kh6'][0]

        print(f"  {'FL_resid':>10s} {'pLDDT':>8s} {'local_PAE':>10s} {'status':>10s}")
        for i in range(0, len(d_plddt), 20):
            end = min(i + 20, len(d_plddt))
            avg_plddt = np.mean(d_plddt[i:end])
            avg_lpae = np.mean(d_local_pae[i:end])
            status = "STABLE" if avg_plddt > 60 and avg_lpae < 5.0 else "FLEX" if avg_plddt > 50 else "DISORDERED"
            print(f"  {fl_offset+i:4d}-{fl_offset+end-1:4d}  {avg_plddt:8.1f} {avg_lpae:10.1f} {status:>10s}")

    # Same for KHb-KH8
    if best_construct and 'khb-kh8' in all_results[best_construct].get('kh_results', {}):
        print()
        d_plddt = all_results[best_construct]['kh_results']['khb-kh8']['mean_plddt']
        d_local_pae = all_results[best_construct]['kh_results']['khb-kh8']['local_pae']
        fl_offset = DOMAIN_DEFS['khb-kh8'][0]

        print(f"  {'FL_resid':>10s} {'pLDDT':>8s} {'local_PAE':>10s} {'status':>10s}")
        for i in range(0, len(d_plddt), 20):
            end = min(i + 20, len(d_plddt))
            avg_plddt = np.mean(d_plddt[i:end])
            avg_lpae = np.mean(d_local_pae[i:end])
            status = "STABLE" if avg_plddt > 60 and avg_lpae < 5.0 else "FLEX" if avg_plddt > 50 else "DISORDERED"
            print(f"  {fl_offset+i:4d}-{fl_offset+end-1:4d}  {avg_plddt:8.1f} {avg_lpae:10.1f} {status:>10s}")

    # Same for KH7a
    if best_construct and 'kh7a' in all_results[best_construct].get('kh_results', {}):
        print()
        d_plddt = all_results[best_construct]['kh_results']['kh7a']['mean_plddt']
        d_local_pae = all_results[best_construct]['kh_results']['kh7a']['local_pae']
        fl_offset = DOMAIN_DEFS['kh7a'][0]

        print(f"  {'FL_resid':>10s} {'pLDDT':>8s} {'local_PAE':>10s} {'status':>10s}")
        for i in range(0, len(d_plddt), 10):
            end = min(i + 10, len(d_plddt))
            avg_plddt = np.mean(d_plddt[i:end])
            avg_lpae = np.mean(d_local_pae[i:end])
            status = "STABLE" if avg_plddt > 60 and avg_lpae < 5.0 else "FLEX" if avg_plddt > 50 else "DISORDERED"
            print(f"  {fl_offset+i:4d}-{fl_offset+end-1:4d}  {avg_plddt:8.1f} {avg_lpae:10.1f} {status:>10s}")

    print("\nDone.")


if __name__ == '__main__':
    main()
