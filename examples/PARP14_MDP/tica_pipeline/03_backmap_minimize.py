#!/usr/bin/env python
"""
STAGE 3 - Back-map the CG representative states (stage 2) to full atom and remove
the clashes that back-mapping introduces.

Thin wrapper around ../fullatom_minimize_states.py (same flags), kept here so the
three pipeline stages live together. Two-tier clash removal:
  1. restrained amber14 energy minimization (CA restraints preserve the arrangement),
  2. rigid-body declash + re-minimize for deep inter-domain interpenetration.

Flags:
  --states-dir DIR     directory of state_*.pdb from stage 2
  --reference  PDB     all-atom reference with side chains
                         md_full -> md_full/input/ref_allatom.pdb
                         fl_optimized -> input/parp14.pdb
  --fdomains   YAML    domain ranges (md_full/input/domains.yaml or input/domains.yaml)
  --k FLOAT            CA restraint constant (default 1000 kJ/mol/nm^2)
  --passes INT         minimization retries with weaker restraint (default 2)

Example:
  python 03_backmap_minimize.py --states-dir ../states/md_full_pose \
      --reference ../md_full/input/ref_allatom.pdb --fdomains ../md_full/input/domains.yaml
"""
import os, sys, runpy

_TARGET = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'fullatom_minimize_states.py')

if __name__ == '__main__':
    # Re-dispatch to the canonical implementation with the same CLI args.
    sys.argv[0] = _TARGET
    runpy.run_path(_TARGET, run_name='__main__')
