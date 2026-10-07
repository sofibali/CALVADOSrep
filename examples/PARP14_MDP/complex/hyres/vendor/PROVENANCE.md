# Vendored from HyRes_GPU

Source: https://github.com/lslumass/HyRes_GPU
Commit: ff02ac14b5da7e1d0dd820079741d938003e54bf
Vendored: 2026-10-07

## Why vendored rather than a submodule

Only two files of the upstream repo's 53 are used by this project
(`at2hyres/at2hyres.py` and `at2hyres/pdbfix_res.py`, both called by
`hyres/runs/prep/convert.sh`). Of upstream's 46 MB, 44 MB is
`examples/phase_separation/` demo data this project never touches, and a
submodule would leave a plain `git clone` with an empty directory.

Pinning the exact bytes also matters because this project works *around*
bugs in these two scripts rather than patching them, so the behaviour
being compensated for has to stay fixed. See
`hyres/runs/prod_rep1/README.md`; in brief:

- `pdbfix_res.py` reads the last whitespace token of a PDB line as a
  segment id, which on a standard PDB is the element symbol, collapsing
  854 residues to 1.
- `at2hyres.py` does not build the backbone amide H that the PSF expects,
  and wants it named `HN`, not `H`.

Files here are unmodified. To work with full upstream instead:

    git clone https://github.com/lslumass/HyRes_GPU.git \
        examples/PARP14_MDP/complex/hyres/HyRes_GPU
    cd examples/PARP14_MDP/complex/hyres/HyRes_GPU && git checkout ff02ac14b5da7e1d0dd820079741d938003e54bf

## Licence

Upstream ships no LICENSE file (checked at the pinned commit), so these
files carry no explicit grant. They are kept here only to reproduce this
project's own results. Before redistributing them or using them beyond
that, ask the upstream authors for terms.
