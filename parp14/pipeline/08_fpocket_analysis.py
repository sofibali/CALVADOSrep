#!/usr/bin/env python
"""
08_fpocket_analysis.py - Run fpocket pocket detection on AlphaFold3 structures

This script:
1. Converts mmCIF to PDB format
2. Runs fpocket for pocket detection
3. Extracts pocket metrics (volume, druggability, etc.)
4. Maps pockets to active sites (MD1, MD2, MD3, ART)
5. Outputs comprehensive pocket analysis CSV

Usage:
    python 08_fpocket_analysis.py --input ../alphafold_outputs --output ../results/pocket_analysis
"""

import os
import sys
import argparse
import subprocess
import json
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile
import shutil

# Active site definitions
ACTIVE_SITES = {
    'MD1': {
        'domain_range': (790, 981),
        'catalytic': [831, 923, 962],
        'pocket': list(range(822, 837)) + list(range(919, 928)) + [961, 966]
    },
    'MD2': {
        'domain_range': (1005, 1193),
        'catalytic': [1035, 1046, 1134, 1171],
        'pocket': list(range(1021, 1025)) + list(range(1034, 1048)) +
                  list(range(1130, 1142)) + [1170, 1171, 1175, 1178]
    },
    'MD3': {
        'domain_range': (1207, 1388),
        'catalytic': [1248, 1259, 1330, 1371],
        'pocket': list(range(1235, 1238)) + list(range(1247, 1262)) +
                  list(range(1302, 1305)) + list(range(1324, 1338)) +
                  list(range(1369, 1372)) + [1375]
    },
    'ART': {
        'domain_range': (1603, 1801),
        'catalytic': [1684, 1705, 1706, 1722],
        'pocket': list(range(1681, 1686)) + [1688, 1701] + list(range(1704, 1710)) +
                  list(range(1714, 1717)) + [1721, 1722, 1726, 1727, 1781]
    }
}


def cif_to_pdb(cif_path, pdb_path):
    """Convert mmCIF to PDB format using BioPython."""
    from Bio.PDB import MMCIFParser, PDBIO

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_path)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)

    return pdb_path


def run_fpocket(pdb_path, output_dir):
    """Run fpocket on a PDB file."""
    try:
        # Run fpocket
        cmd = ['fpocket', '-f', pdb_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            return None

        # fpocket creates output in same directory as input
        base_name = Path(pdb_path).stem
        fpocket_dir = Path(pdb_path).parent / f"{base_name}_out"

        if not fpocket_dir.exists():
            return None

        return fpocket_dir

    except subprocess.TimeoutExpired:
        return None
    except FileNotFoundError:
        print("ERROR: fpocket not found. Install with: conda install -c conda-forge fpocket")
        return None


def parse_fpocket_info(info_file):
    """Parse fpocket info file to extract pocket descriptors."""
    pockets = []

    if not os.path.exists(info_file):
        return pockets

    with open(info_file, 'r') as f:
        content = f.read()

    # Split by pocket
    pocket_blocks = content.split('Pocket')[1:]  # Skip header

    for i, block in enumerate(pocket_blocks):
        pocket = {'pocket_id': i + 1}

        lines = block.strip().split('\n')
        for line in lines:
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip().lower().replace(' ', '_').replace('-', '_')
                value = value.strip()

                # Try to convert to float
                try:
                    value = float(value)
                except ValueError:
                    pass

                pocket[key] = value

        pockets.append(pocket)

    return pockets


def parse_pocket_pdb(pocket_pdb):
    """Extract residue numbers from pocket PDB file."""
    residues = set()

    if not os.path.exists(pocket_pdb):
        return residues

    with open(pocket_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    resnum = int(line[22:26].strip())
                    residues.add(resnum)
                except ValueError:
                    continue

    return residues


def map_pocket_to_active_site(pocket_residues, min_overlap=3):
    """Map a detected pocket to known active sites."""
    mappings = []

    for site_name, site_def in ACTIVE_SITES.items():
        site_residues = set(site_def['pocket'])
        overlap = pocket_residues & site_residues

        if len(overlap) >= min_overlap:
            catalytic_overlap = pocket_residues & set(site_def['catalytic'])
            mappings.append({
                'active_site': site_name,
                'overlap_count': len(overlap),
                'overlap_fraction': len(overlap) / len(site_residues),
                'catalytic_overlap': len(catalytic_overlap),
                'catalytic_covered': len(catalytic_overlap) / len(site_def['catalytic'])
            })

    return mappings


def check_domain_in_structure(structure_name):
    """Check which domains are present in structure based on name."""
    name_lower = structure_name.lower()
    present = {}

    domain_patterns = {
        'MD1': ['md1', 'md1l1'],
        'MD2': ['md2'],
        'MD3': ['md3'],
        'ART': ['art']
    }

    for domain, patterns in domain_patterns.items():
        present[domain] = any(p in name_lower for p in patterns)

    return present


def analyze_structure(args):
    """Analyze a single structure with fpocket."""
    structure_dir, seed_dir, output_base = args

    structure_name = os.path.basename(structure_dir)
    seed_name = os.path.basename(seed_dir)

    cif_path = os.path.join(seed_dir, 'model.cif')
    if not os.path.exists(cif_path):
        return None

    results = []

    try:
        # Create temporary directory for processing
        with tempfile.TemporaryDirectory() as tmpdir:
            # Convert CIF to PDB
            pdb_path = os.path.join(tmpdir, 'structure.pdb')
            cif_to_pdb(cif_path, pdb_path)

            # Run fpocket
            fpocket_dir = run_fpocket(pdb_path, tmpdir)

            if fpocket_dir is None:
                return None

            # Parse results
            info_file = fpocket_dir / f"structure_info.txt"
            pockets = parse_fpocket_info(info_file)

            # Check which domains are in this structure
            domains_present = check_domain_in_structure(structure_name)

            # Analyze each pocket
            for pocket in pockets:
                pocket_id = pocket.get('pocket_id', 0)
                pocket_pdb = fpocket_dir / f"pockets" / f"pocket{pocket_id}_atm.pdb"

                pocket_residues = parse_pocket_pdb(pocket_pdb)
                active_site_mappings = map_pocket_to_active_site(pocket_residues)

                # Base result
                result = {
                    'structure': structure_name,
                    'seed': seed_name,
                    'pocket_id': pocket_id,
                    'volume': pocket.get('volume', None),
                    'druggability_score': pocket.get('druggability_score', None),
                    'score': pocket.get('score', None),
                    'hydrophobicity_score': pocket.get('hydrophobicity_score', None),
                    'polarity_score': pocket.get('polarity_score', None),
                    'apolar_sasa': pocket.get('apolar_alpha_sphere_sasa', None),
                    'polar_sasa': pocket.get('polar_alpha_sphere_sasa', None),
                    'n_residues': len(pocket_residues),
                    **{f'has_{d}': v for d, v in domains_present.items()}
                }

                # Add active site mappings
                for mapping in active_site_mappings:
                    result_with_mapping = result.copy()
                    result_with_mapping.update({
                        'mapped_active_site': mapping['active_site'],
                        'overlap_count': mapping['overlap_count'],
                        'overlap_fraction': mapping['overlap_fraction'],
                        'catalytic_overlap': mapping['catalytic_overlap'],
                        'catalytic_covered': mapping['catalytic_covered']
                    })
                    results.append(result_with_mapping)

                # If no mapping, still record the pocket
                if not active_site_mappings:
                    result['mapped_active_site'] = None
                    results.append(result)

            # Also save full fpocket output
            output_struct_dir = os.path.join(output_base, structure_name, seed_name)
            os.makedirs(output_struct_dir, exist_ok=True)

            if fpocket_dir.exists():
                shutil.copytree(fpocket_dir, os.path.join(output_struct_dir, 'fpocket_out'),
                               dirs_exist_ok=True)

    except Exception as e:
        print(f"Error processing {structure_name}/{seed_name}: {e}")
        return None

    return results


def main():
    parser = argparse.ArgumentParser(description='Run fpocket analysis on AlphaFold3 structures')
    parser.add_argument('--input', '-i', default='../alphafold_outputs',
                       help='Input directory with AlphaFold3 outputs')
    parser.add_argument('--output', '-o', default='../results/pocket_analysis',
                       help='Output directory for fpocket results')
    parser.add_argument('--csv-output', default='../results/per_structure/fpocket_results.csv',
                       help='Output CSV file path')
    parser.add_argument('--workers', '-w', type=int, default=8,
                       help='Number of parallel workers')
    parser.add_argument('--seed', '-s', default=None,
                       help='Process only specific seed (e.g., seed-1_sample-0)')
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all structure/seed combinations
    tasks = []
    for structure_dir in sorted(input_dir.iterdir()):
        if not structure_dir.is_dir():
            continue

        for seed_dir in sorted(structure_dir.iterdir()):
            if not seed_dir.is_dir() or not seed_dir.name.startswith('seed-'):
                continue

            if args.seed and seed_dir.name != args.seed:
                continue

            if (seed_dir / 'model.cif').exists():
                tasks.append((str(structure_dir), str(seed_dir), str(output_dir)))

    print(f"Found {len(tasks)} structures to analyze")

    # Process in parallel
    all_results = []

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(analyze_structure, task): task for task in tasks}

        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
            result = future.result()
            if result:
                all_results.extend(result)

    # Save results
    if all_results:
        df = pd.DataFrame(all_results)

        # Ensure output directory exists
        Path(args.csv_output).parent.mkdir(parents=True, exist_ok=True)

        df.to_csv(args.csv_output, index=False)
        print(f"\nSaved {len(df)} pocket records to {args.csv_output}")

        # Print summary
        print("\n=== Summary ===")
        print(f"Total pockets detected: {len(df)}")
        print(f"Structures analyzed: {df['structure'].nunique()}")

        if 'mapped_active_site' in df.columns:
            mapped = df[df['mapped_active_site'].notna()]
            print(f"\nPockets mapped to active sites: {len(mapped)}")
            print(mapped['mapped_active_site'].value_counts())
    else:
        print("No results generated. Check that fpocket is installed and structures exist.")


if __name__ == '__main__':
    main()
