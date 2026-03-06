#!/usr/bin/env python
"""
extract_fasta_from_pdb.py - Extract FASTA sequences from PDB files

Usage:
    # Single PDB file
    python extract_fasta_from_pdb.py --pdb input/protein.pdb --output sequences.fasta

    # Multiple PDB files
    python extract_fasta_from_pdb.py --pdb_dir input/structures/ --output sequences.fasta

    # Specify custom name
    python extract_fasta_from_pdb.py --pdb input/protein.pdb --name MyProtein --output sequences.fasta

This uses CALVADOS built-in functions from calvados.sequence module.
"""

import os
import sys
import glob
from argparse import ArgumentParser

# Add CALVADOS to path if needed
try:
    from calvados.sequence import seq_from_pdb, write_fasta, record_from_seq
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from calvados.sequence import seq_from_pdb, write_fasta, record_from_seq


def extract_single_pdb(pdb_file, name=None):
    """
    Extract sequence from a single PDB file.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file
    name : str
        Sequence name. If None, uses PDB filename without extension.

    Returns
    -------
    record : Bio.SeqRecord
        FASTA record
    """

    if name is None:
        name = os.path.splitext(os.path.basename(pdb_file))[0]

    seq, n_termini, c_termini = seq_from_pdb(pdb_file)

    record = record_from_seq(seq, name)

    print(f"  {name}: {len(seq)} residues")

    return record, seq


def extract_from_directory(pdb_dir, pattern='*.pdb'):
    """
    Extract sequences from all PDB files in a directory.

    Parameters
    ----------
    pdb_dir : str
        Directory containing PDB files
    pattern : str
        Glob pattern for PDB files

    Returns
    -------
    records : list
        List of Bio.SeqRecord objects
    """

    pdb_files = sorted(glob.glob(f'{pdb_dir}/{pattern}'))

    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in {pdb_dir} with pattern {pattern}")

    records = []
    for pdb_file in pdb_files:
        try:
            record, _ = extract_single_pdb(pdb_file)
            records.append(record)
        except Exception as e:
            print(f"  Error processing {pdb_file}: {e}")

    return records


def main(pdb_file=None, pdb_dir=None, name=None, output_file='sequences.fasta',
         append=False):
    """
    Main extraction routine.

    Parameters
    ----------
    pdb_file : str
        Single PDB file to process
    pdb_dir : str
        Directory of PDB files to process
    name : str
        Custom name for single PDB (ignored for directory mode)
    output_file : str
        Output FASTA file
    append : bool
        Append to existing FASTA file
    """

    records = []

    if pdb_file:
        print(f"Extracting sequence from: {pdb_file}")
        record, seq = extract_single_pdb(pdb_file, name)
        records.append(record)

    elif pdb_dir:
        print(f"Extracting sequences from directory: {pdb_dir}")
        records = extract_from_directory(pdb_dir)

    else:
        raise ValueError("Must specify either --pdb or --pdb_dir")

    if not records:
        print("No sequences extracted!")
        return

    # Write FASTA file
    if append and os.path.exists(output_file):
        print(f"\nAppending to: {output_file}")
    else:
        print(f"\nWriting to: {output_file}")

    write_fasta(records, output_file)

    print(f"Wrote {len(records)} sequence(s) to {output_file}")

    # Print summary
    print("\n--- Sequences ---")
    for rec in records:
        print(f">{rec.id}")
        seq_str = str(rec.seq)
        # Print first 60 chars + ... if longer
        if len(seq_str) > 60:
            print(f"{seq_str[:60]}...")
        else:
            print(seq_str)


def print_sequence_stats(pdb_file):
    """Print detailed sequence statistics."""

    seq, n_termini, c_termini = seq_from_pdb(pdb_file)

    # Count residue types
    from collections import Counter
    counts = Counter(seq)

    # Charged residues
    positive = counts.get('K', 0) + counts.get('R', 0)
    negative = counts.get('D', 0) + counts.get('E', 0)
    net_charge = positive - negative

    # Aromatic residues
    aromatic = counts.get('F', 0) + counts.get('Y', 0) + counts.get('W', 0)

    # Hydrophobic
    hydrophobic = sum(counts.get(aa, 0) for aa in 'AILMFWV')

    print(f"\n--- Sequence Statistics ---")
    print(f"Length: {len(seq)}")
    print(f"N-termini: {n_termini}")
    print(f"C-termini: {c_termini}")
    print(f"Positive residues (K, R): {positive}")
    print(f"Negative residues (D, E): {negative}")
    print(f"Net charge: {net_charge:+d}")
    print(f"Aromatic (F, Y, W): {aromatic} ({100*aromatic/len(seq):.1f}%)")
    print(f"Hydrophobic (AILMFWV): {hydrophobic} ({100*hydrophobic/len(seq):.1f}%)")

    print(f"\nResidue composition:")
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
        count = counts.get(aa, 0)
        if count > 0:
            print(f"  {aa}: {count} ({100*count/len(seq):.1f}%)")


if __name__ == "__main__":
    parser = ArgumentParser(description="Extract FASTA from PDB files")
    parser.add_argument('--pdb', default=None, type=str,
                        help='Single PDB file')
    parser.add_argument('--pdb_dir', default=None, type=str,
                        help='Directory of PDB files')
    parser.add_argument('--name', default=None, type=str,
                        help='Custom name for sequence (single PDB mode)')
    parser.add_argument('--output', default='sequences.fasta', type=str,
                        help='Output FASTA file')
    parser.add_argument('--append', action='store_true',
                        help='Append to existing FASTA file')
    parser.add_argument('--stats', action='store_true',
                        help='Print sequence statistics')

    args = parser.parse_args()

    if not args.pdb and not args.pdb_dir:
        parser.print_help()
        print("\nError: Must specify either --pdb or --pdb_dir")
        sys.exit(1)

    main(
        pdb_file=args.pdb,
        pdb_dir=args.pdb_dir,
        name=args.name,
        output_file=args.output,
        append=args.append
    )

    if args.stats and args.pdb:
        print_sequence_stats(args.pdb)
