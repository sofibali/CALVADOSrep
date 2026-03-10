#!/usr/bin/env python3
"""
Generate the missing AF3 JSON input files to complete the full 2^11 - 1 = 2047 library.

Domain units (11):
  RRM1, RRM2, RRM3, KH1-KH6, KH7a, MD1L1, MD2, MD3, KHb-KH8, WWE, ART

Existing: 1535 of 2047 (plus 15 MD1-without-L1 inputs that are ignored)
Missing: 512 combinations to generate
"""
import os
import json
import csv
from itertools import combinations

# ── Configuration ──
FASTA_FILE = "/home/sbali/CALVADOS/parp14/input/PARP14.fasta"
DOMAIN_FILE = "/home/sbali/CALVADOS/parp14/input/domain_boundaries.csv"
EXISTING_DIR = "/home/sbali/CALVADOS/parp14/alphafold_inputs"
OUTPUT_DIR = "/home/sbali/CALVADOS/parp14/alphafold_inputs_missing"
MISSING_LIST = "/home/sbali/CALVADOS/parp14/missing_predictions.txt"

DOMAINS_ORDERED = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a', 'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']

# ── Read full-length sequence ──
with open(FASTA_FILE) as f:
    lines = f.readlines()
sequence = ''.join(line.strip() for line in lines[1:])
print(f"Full-length sequence: {len(sequence)} aa")

# ── Read domain boundaries ──
domain_bounds = {}
with open(DOMAIN_FILE) as f:
    reader = csv.DictReader(f)
    for row in reader:
        domain_bounds[row['Domain'].strip()] = (int(row['Start']), int(row['End']))

print(f"Domain boundaries: {len(domain_bounds)}")
for d in DOMAINS_ORDERED:
    s, e = domain_bounds[d]
    print(f"  {d}: {s}-{e} ({e - s + 1} aa)")

# ── Find missing combinations ──
# All 2^11 - 1 target subsets
all_targets = set()
for r in range(1, len(DOMAINS_ORDERED) + 1):
    for c in combinations(DOMAINS_ORDERED, r):
        all_targets.add(c)

# Existing inputs (case-insensitive match)
existing_lower = set()
for f in os.listdir(EXISTING_DIR):
    if f.endswith('.json'):
        existing_lower.add(f.replace('.json', '').lower())

# Find missing
missing = []
for combo in sorted(all_targets):
    name = '_'.join(combo)
    if name.lower() not in existing_lower:
        missing.append(combo)

print(f"\nTotal target: {len(all_targets)}")
print(f"Already exist: {len(all_targets) - len(missing)}")
print(f"Missing: {len(missing)}")

# ── Generate missing JSON files ──
os.makedirs(OUTPUT_DIR, exist_ok=True)

for combo in missing:
    name = '_'.join(combo)

    # Extract and concatenate domain sequences
    combo_seq = ''
    for domain in combo:
        start, end = domain_bounds[domain]
        combo_seq += sequence[start - 1:end]

    # Create AF3 JSON
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

    json_path = os.path.join(OUTPUT_DIR, f"{name}.json")
    with open(json_path, 'w') as f:
        json.dump(data, f, indent=2)

print(f"\nGenerated {len(missing)} JSON files in {OUTPUT_DIR}/")

# ── Write missing predictions list (for the run script) ──
with open(MISSING_LIST, 'w') as f:
    for combo in missing:
        f.write('_'.join(combo) + '\n')

print(f"Missing predictions list: {MISSING_LIST}")

# ── Summary by combination size ──
from collections import Counter
size_counts = Counter(len(c) for c in missing)
print(f"\nMissing by combination size:")
for size in sorted(size_counts):
    print(f"  {size} domains: {size_counts[size]}")

# ── Sequence length distribution ──
lengths = []
for combo in missing:
    seq = ''
    for domain in combo:
        s, e = domain_bounds[domain]
        seq += sequence[s - 1:e]
    lengths.append(len(seq))

print(f"\nSequence length range: {min(lengths)}-{max(lengths)} aa")
print(f"Mean: {sum(lengths)/len(lengths):.0f} aa")
