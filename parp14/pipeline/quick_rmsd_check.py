#!/usr/bin/env python3
"""
Quick RMSD calculation for AlphaFold3 structures
"""
import numpy as np
from Bio import PDB
from pathlib import Path
import sys

def calculate_rmsd(coords1, coords2):
    """Calculate RMSD between two sets of coordinates after alignment"""
    # Center both structures
    coords1_centered = coords1 - coords1.mean(axis=0)
    coords2_centered = coords2 - coords2.mean(axis=0)
    
    # Calculate RMSD
    diff = coords1_centered - coords2_centered
    rmsd = np.sqrt((diff ** 2).sum() / len(coords1))
    return rmsd

def get_ca_coords(structure):
    """Extract CA coordinates from structure"""
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
    return np.array(coords)

def main():
    structure_dir = sys.argv[1] if len(sys.argv) > 1 else "/home/sbali/parp14/alphafold_outputs/kh1_kh2_kh3_kh4_kh5_kh6_art"
    
    parser = PDB.MMCIFParser(QUIET=True)
    
    # Get all structure files
    cif_files = sorted(Path(structure_dir).glob("seed-*/model.cif"))[:10]  # First 10 structures
    
    if not cif_files:
        print(f"No CIF files found in {structure_dir}")
        return
    
    print(f"Analyzing structures from: {structure_dir}")
    print(f"Found {len(cif_files)} structures")
    print("\nRMSD Matrix (Å) - CA atoms only:")
    print("=" * 80)
    
    # Load all structures
    structures = []
    names = []
    for cif_file in cif_files:
        name = cif_file.parent.name
        names.append(name)
        structure = parser.get_structure(name, str(cif_file))
        structures.append(get_ca_coords(structure))
    
    # Calculate RMSD matrix
    print(f"\n{'':20}", end="")
    for name in names[:6]:  # Show first 6
        print(f"{name:15}", end="")
    print()
    
    for i, (name1, coords1) in enumerate(zip(names[:6], structures[:6])):
        print(f"{name1:20}", end="")
        for j, coords2 in enumerate(structures[:6]):
            if len(coords1) == len(coords2):
                rmsd = calculate_rmsd(coords1, coords2)
                print(f"{rmsd:15.3f}", end="")
            else:
                print(f"{'N/A':15}", end="")
        print()
    
    # Summary statistics
    print("\n" + "=" * 80)
    print("Summary Statistics:")
    print("=" * 80)
    
    # Compare seed-1 samples
    seed1_structures = [(n, s) for n, s in zip(names, structures) if n.startswith('seed-1')]
    if len(seed1_structures) >= 2:
        print("\nRMSD between seed-1 samples:")
        for i in range(min(5, len(seed1_structures))):
            for j in range(i+1, min(5, len(seed1_structures))):
                name1, coords1 = seed1_structures[i]
                name2, coords2 = seed1_structures[j]
                if len(coords1) == len(coords2):
                    rmsd = calculate_rmsd(coords1, coords2)
                    print(f"  {name1} vs {name2}: {rmsd:.3f} Å")
    
    # Compare different seeds
    print("\nRMSD between different seeds (sample-0):")
    sample0_structures = [(n, s) for n, s in zip(names, structures) if n.endswith('sample-0')]
    for i in range(min(3, len(sample0_structures))):
        for j in range(i+1, min(3, len(sample0_structures))):
            name1, coords1 = sample0_structures[i]
            name2, coords2 = sample0_structures[j]
            if len(coords1) == len(coords2):
                rmsd = calculate_rmsd(coords1, coords2)
                print(f"  {name1} vs {name2}: {rmsd:.3f} Å")

if __name__ == "__main__":
    main()
