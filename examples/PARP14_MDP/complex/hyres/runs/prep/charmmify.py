#!/usr/bin/env python3
"""All-atom PDB -> heavy-atom PDB with a CHARMM segment id.

    python charmmify.py in.pdb out.pdb SEGID

Rebuilt 2026-10-07; the original was lost with a scratchpad.

VERIFICATION. Re-running the full convert.sh on binding/p9_dtx3l/input/
parp9.pdb reproduces the tracked parp9_hyres.pdb exactly in particle count
(5,346) and topology, and all 3,416 backbone N/CA/C/O atoms land at
identical coordinates (0.000 A). Sidechain beads differ by 0.37 A on
average, 0.79 A at worst, because HyRes places each bead at a centroid
computed over the all-atom sidechain INCLUDING its hydrogens, and
addH.py's OpenMM Modeller.addHydrogens does not place rotatable hydrogens
reproducibly between runs -- the same OpenMM 8.3.1 gives different answers.
So the pipeline is reproducible in topology but not bit-for-bit in
coordinates, which is why the *_hyres.pdb files are committed: they, not a
re-run, are the authoritative inputs to the published HyRes trajectories.

Three jobs, all of them things psfgen/at2hyres downstream require:

1. Drop every hydrogen EXCEPT the backbone amide HN, which HyRes keeps
   explicitly and at2hyres will not rebuild. The element column (77-78) is
   used when present and the atom name only as a fallback.
2. Write SEGID into columns 73-76, which CHARMM-style tooling keys on and
   which the PDB written by OpenMM leaves blank.
3. Renumber atom serials from 1, since dropping hydrogens leaves gaps.

The B-factor column is zeroed. The archived *_cg.pdb files in this
directory show real B-factors from an earlier generation of the pipeline,
but nothing downstream reads them: pdbfix_res writes 0.0 and the final
*_hyres.pdb carries 0.00.
"""
import sys


def is_dropped_hydrogen(line):
    """True for every hydrogen EXCEPT the backbone amide HN.

    HyRes keeps the amide hydrogen explicitly -- it drives the hydrogen-bond
    term -- so at2hyres expects to find it and only maps atoms that already
    exist. addH.py builds it (OpenMM names it H) and convert.sh renames it to
    HN just before this step; dropping it here would leave the final structure
    815 particles short of its own PSF, one per non-proline residue.
    """
    name = line[12:16].strip()
    if name == 'HN':
        return False
    el = line[76:78].strip()
    if el:
        return el.upper() == 'H'
    return name.lstrip('0123456789').startswith('H')


def main(src, dst, segid):
    if len(segid) > 4:
        raise SystemExit(f'segid {segid!r} is longer than the 4 columns PDB allows')
    n = nh = 0
    with open(src) as f, open(dst, 'w') as g:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            if is_dropped_hydrogen(line):
                continue
            n += 1
            nh += line[12:16].strip() == 'HN'
            row = line.rstrip('\n').ljust(80)
            # The element column (77-78) is deliberately NOT carried through.
            # pdbfix_res.py reads the LAST whitespace token of the line as the
            # segment id, so leaving the element there makes it read "N"/"C"/"O"
            # and collapse all 854 residues into one. Dropping it leaves SEGID
            # last, which is what that script actually wants.
            row = (row[:6] + f'{n:5d}' + row[11:54]
                   + f'{1.00:6.2f}{0.00:6.2f}' + ' ' * 6
                   + f'{segid:<4s}')
            g.write(row.rstrip() + '\n')
        g.write('END\n')
    print(f'  {dst}: {n} atoms ({nh} amide HN kept), segid {segid}')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
