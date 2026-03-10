#!/usr/bin/env python3
"""
Generate ALL missing AF3 inputs and a complete missing predictions list.

Target: full 2^11 - 1 = 2047 domain combinations using MD1L1 (not MD1).
Checks which already have completed model.cif outputs and generates
JSON inputs + a run list for everything still needed.
"""
import os
import json
import csv
from itertools import combinations
from collections import Counter

FASTA_FILE = "/home/sbali/CALVADOS/parp14/input/PARP14.fasta"
DOMAIN_FILE = "/home/sbali/CALVADOS/parp14/input/domain_boundaries.csv"
EXISTING_INPUT_DIR = "/home/sbali/CALVADOS/parp14/alphafold_inputs"
OUTPUT_DIR = "/home/sbali/CALVADOS/parp14/alphafold_outputs"
MISSING_JSON_DIR = "/home/sbali/CALVADOS/parp14/alphafold_inputs_missing"
MISSING_LIST = "/home/sbali/CALVADOS/parp14/missing_predictions.txt"

DOMAINS = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a', 'MD1L1', 'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']

# ── Read sequence and domain boundaries ──
with open(FASTA_FILE) as f:
    lines = f.readlines()
sequence = ''.join(line.strip() for line in lines[1:])

domain_bounds = {}
with open(DOMAIN_FILE) as f:
    for row in csv.DictReader(f):
        domain_bounds[row['Domain'].strip()] = (int(row['Start']), int(row['End']))

# ── Full target set ──
target = {}
for r in range(1, len(DOMAINS) + 1):
    for c in combinations(DOMAINS, r):
        name = '_'.join(c)
        target[name.lower()] = c

print(f"Target: {len(target)} combinations")

# ── Find completed outputs (have model.cif) ──
completed = set()
for d in os.listdir(OUTPUT_DIR):
    full = os.path.join(OUTPUT_DIR, d)
    if not os.path.isdir(full) or d in ('logs', 'test_run', 'missing_l1'):
        continue
    for sub in os.listdir(full):
        subpath = os.path.join(full, sub)
        if os.path.isdir(subpath) and os.path.exists(os.path.join(subpath, 'model.cif')):
            completed.add(d.lower())
            break

print(f"Completed: {len(completed)}")

# ── Find existing input JSONs ──
existing_inputs = set()
for f in os.listdir(EXISTING_INPUT_DIR):
    if f.endswith('.json'):
        existing_inputs.add(f.replace('.json', '').lower())

# ── Determine what's missing ──
needed = {k: v for k, v in target.items() if k not in completed}
print(f"Still needed: {len(needed)}")

# ── Generate missing JSON inputs (only those without existing input) ──
os.makedirs(MISSING_JSON_DIR, exist_ok=True)

# Clear old files in missing dir
for f in os.listdir(MISSING_JSON_DIR):
    if f.endswith('.json'):
        os.remove(os.path.join(MISSING_JSON_DIR, f))

new_jsons = 0
reuse_jsons = 0
for name_lower, combo in sorted(needed.items()):
    name = '_'.join(combo)

    if name_lower in existing_inputs:
        # Copy from existing inputs dir
        reuse_jsons += 1
        src = None
        for f in os.listdir(EXISTING_INPUT_DIR):
            if f.replace('.json', '').lower() == name_lower:
                src = os.path.join(EXISTING_INPUT_DIR, f)
                break
        if src:
            with open(src) as f:
                data = json.load(f)
            # Ensure 5 seeds
            if data.get('modelSeeds') != [1, 2, 3, 4, 5]:
                data['modelSeeds'] = [1, 2, 3, 4, 5]
            with open(os.path.join(MISSING_JSON_DIR, f"{name}.json"), 'w') as f:
                json.dump(data, f, indent=2)
        continue

    # Generate new JSON
    new_jsons += 1
    combo_seq = ''
    for domain in combo:
        start, end = domain_bounds[domain]
        combo_seq += sequence[start - 1:end]

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
    with open(os.path.join(MISSING_JSON_DIR, f"{name}.json"), 'w') as f:
        json.dump(data, f, indent=2)

print(f"\nJSON files in {MISSING_JSON_DIR}/:")
print(f"  Copied from existing inputs: {reuse_jsons}")
print(f"  Newly generated: {new_jsons}")
print(f"  Total: {reuse_jsons + new_jsons}")

# ── Write missing predictions list ──
with open(MISSING_LIST, 'w') as f:
    for name_lower in sorted(needed.keys()):
        combo = needed[name_lower]
        f.write('_'.join(combo) + '\n')

print(f"\nMissing list: {MISSING_LIST} ({len(needed)} entries)")

# ── Summary ──
size_counts = Counter(len(c) for c in needed.values())
print(f"\nMissing by combination size:")
for size in sorted(size_counts):
    print(f"  {size} domains: {size_counts[size]}")

lengths = []
for combo in needed.values():
    seq = ''
    for domain in combo:
        s, e = domain_bounds[domain]
        seq += sequence[s - 1:e]
    lengths.append(len(seq))

print(f"\nSequence lengths: {min(lengths)}-{max(lengths)} aa (mean {sum(lengths)/len(lengths):.0f})")
print(f"\nEstimated models to generate: {len(needed)} x 25 = {len(needed) * 25}")
