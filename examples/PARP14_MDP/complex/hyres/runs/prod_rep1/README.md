# HyRes PARP9-DTX3L, one replicate

Conditions matched to the CALVADOS `p9_dtx3l` runs: 40 nm cubic box
(25.9 uM per chain), 293 K, 150 mM salt, separated start (chains 17 nm apart),
500 ns, 50 ps/frame. **No restraints of any kind** -- HyRes has no working
folded-protein restraint module (restraints/folded_restriants.py is an
unfinished stub with a syntax error and is not wired into the run scripts), so
this is a genuinely unrestrained test.

Two questions from one trajectory:
1. FOLD STABILITY (the gate). HyRes is validated on ~100 *disordered* proteins;
   nothing establishes it holds a folded macrodomain or KH fold. Per-domain
   RMSD/Rg vs the start decides whether anything else here is interpretable.
   If the domains unfold, the reported HyRes result -- PARP14-PARP9 in contact
   99.98% of the time -- is explained trivially, because unfolded chains with
   exposed hydrophobic cores stick to everything.
2. CROSSLINK AGREEMENT. Score inter-chain contacts against the three BS3
   crosslinks of Ashok et al. 2022 (doi 10.1042/BCJ20210722) using the same
   3.0 nm CA-CA criterion and buried/exposed filter as the CALVADOS side.

## Build notes (four upstream issues had to be worked around)
1. `pdbfix_res.py:28` takes the LAST whitespace token as the segment id. A
   standard PDB's last token is the element symbol, so it saw a new segment at
   every N->C->O and collapsed 854 residues to 1. Fixed by writing a SEGID into
   columns 73-76 of the input and blanking the element columns.
2. HyRes keeps the backbone amide H explicitly, and `at2hyres` only MAPS atoms
   that exist -- it does not build them. Hydrogens added with OpenMM Modeller
   against CHARMM36. This is the only thing CHARMM-GUI was needed for, so
   registration is NOT actually required, contrary to the README.
3. `at2hyres.py:74` expects the amide hydrogen named `HN` (CHARMM); OpenMM
   writes `H` (PDB). Renamed. PSF and PDB then match exactly at 9,968.
4. The repo run script sets velocities and then calls a bare `minimizeEnergy()`,
   which NaNs immediately on this system. Minimising with
   `maxIterations=10000` BEFORE setting velocities fixes it, and 4 fs is then
   stable (295 K, equilibrated).

Throughput measured here: ~3,450 ns/day, i.e. ~3.5 h for the 500 ns replicate.
