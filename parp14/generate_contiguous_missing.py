#!/usr/bin/env python3
"""
Generate AF3 inputs for the 18 missing CONTIGUOUS domain constructs.

These are all N-terminal-extending constructs (RRM1/2/3 + core domains)
that haven't been predicted yet. Uses the same JSON format as
generate_all_missing.py but targets only the contiguous subsequences.

Usage:
    python generate_contiguous_missing.py
    python generate_contiguous_missing.py --dry-run
"""

import os
import json
import csv
import glob
from argparse import ArgumentParser

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FASTA_FILE = os.path.join(SCRIPT_DIR, "input/PARP14.fasta")
DOMAIN_FILE = os.path.join(SCRIPT_DIR, "input/domain_boundaries.csv")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "alphafold_outputs")
JSON_DIR = os.path.join(SCRIPT_DIR, "alphafold_inputs_missing")

UNITS = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a', 'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']


def main():
    parser = ArgumentParser()
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    # Read sequence
    with open(FASTA_FILE) as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines[1:])

    # Read domain boundaries
    domain_bounds = {}
    with open(DOMAIN_FILE) as f:
        for row in csv.DictReader(f):
            domain_bounds[row['Domain'].strip()] = (int(row['Start']), int(row['End']))

    # Enumerate all contiguous subsequences
    all_contiguous = []
    for i in range(len(UNITS)):
        for j in range(i, len(UNITS)):
            units = UNITS[i:j+1]
            name = '_'.join(units)
            all_contiguous.append((name, units))

    # Check which are missing
    missing = []
    for name, units in all_contiguous:
        name_lower = name.lower()
        # Check exact dir or timestamped
        found = False
        exact = os.path.join(OUTPUT_DIR, name_lower)
        if glob.glob(os.path.join(exact, 'seed-*/model.cif')):
            found = True
        else:
            for d in glob.glob(exact + '_*'):
                if glob.glob(os.path.join(d, 'seed-*/model.cif')):
                    found = True
                    break
        if not found:
            # Build sequence
            combo_seq = ''
            for domain in units:
                start, end = domain_bounds[domain]
                combo_seq += sequence[start - 1:end]
            missing.append((name, units, combo_seq))

    print(f"Total contiguous constructs: {len(all_contiguous)}")
    print(f"Already completed: {len(all_contiguous) - len(missing)}")
    print(f"Missing: {len(missing)}")
    print()

    for name, units, seq in missing:
        n_units = len(units)
        print(f"  {n_units:>2d} domains  {name:<65s}  {len(seq)} aa")

    if args.dry_run:
        print("\n[DRY RUN] No files written.")
        return

    # Generate JSON inputs
    os.makedirs(JSON_DIR, exist_ok=True)
    written = 0
    for name, units, combo_seq in missing:
        data = {
            "name": name,
            "sequences": [
                {
                    "protein": {
                        "id": ["A"],
                        "sequence": combo_seq
                    }
                }
            ],
            "modelSeeds": [1, 2, 3, 4, 5],
            "dialect": "alphafold3",
            "version": 1
        }
        out_path = os.path.join(JSON_DIR, f"{name}.json")
        with open(out_path, 'w') as f:
            json.dump(data, f, indent=2)
        written += 1

    print(f"\nWrote {written} JSON files to {JSON_DIR}/")
    print(f"\nTo run AF3:")
    print(f"  cd {SCRIPT_DIR}")
    print(f"  bash run_af3_missing.sh              # full pipeline (MSA + inference)")
    print(f"  bash run_af3_missing.sh --skip-msa   # inference only (if MSA done)")
    print(f"  bash run_af3_missing.sh --dry-run    # preview")


if __name__ == '__main__':
    main()
