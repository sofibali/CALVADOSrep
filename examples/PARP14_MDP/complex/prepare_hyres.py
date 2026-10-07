#!/usr/bin/env python3
"""Prepare HyRes runs for PARP9-DTX3L, to test the published crosslinks.

THE POINT
---------
CALVADOS puts PARP9-DTX3L at ~263 uM apparent Kd with a diffuse interface: it
takes 1,206 domain-pair crosslinks to reach half the predicted signal, and of
the three BS3 crosslinks Ashok et al. 2022 actually observe (doi
10.1042/BCJ20210722), one is never sampled and another ranks 1,864th of 4,259.
HyRes is the other candidate force field. The question is not whether HyRes
binds harder -- it plainly does -- but whether it binds in the RIGHT PLACE,
with no inter-chain restraints imposed.

"WITHOUT ADDITIONAL RESTRAINTS BESIDES THE NORMAL ONES" -- THERE ARE NONE
------------------------------------------------------------------------
Checked in the HyRes_GPU repo (github.com/lslumass/HyRes_GPU) rather than
assumed:
  * restraints/folded_restriants.py is an unfinished stub. Its header is
    literally "# TODO define restraint functions for folded proteins", it has a
    syntax error at line 50 (`np.squeeze(pdb` / `_md.xyz`, a broken line
    continuation referencing an undefined `pdb_md`), so it cannot even be
    imported, and
  * it is NOT referenced anywhere in run_OpenMM/run.nvt.py or run.npt.py.

So a HyRes run of a folded protein is COMPLETELY unrestrained. There is no
"normal" fold restraint to keep. That is exactly the restraint-free test we
want -- and simultaneously the largest risk, because HyRes is an IDP model:
Li, Barethiya & Chen 2026 validate it across ~100 *intrinsically disordered*
proteins. Nothing published establishes that it holds a folded macrodomain or
KH fold at all.

HENCE STAGE 0 IS A GATE, NOT A FORMALITY
----------------------------------------
Stage 0 runs each chain ALONE, unrestrained, and asks whether the folds
survive. If they do not, then the reported HyRes result -- PARP14-PARP9 in
contact 99.98% of the time -- is explained trivially: unfolded chains with
exposed hydrophobic cores stick to everything, and no crosslink comparison
based on such an ensemble means anything. That outcome would be a real finding
about the HyRes/CALVADOS comparison, not a failed setup.

Only if the folds hold does Stage 1 (the complex) become interpretable.

MATCHED CONDITIONS
------------------
Stage 1 mirrors the CALVADOS runs so the comparison is apples to apples:
40 nm cubic box (25.9 uM per chain), 293 K, 150 mM salt, separated start, no
inter-chain restraints. Contacts are then scored against the three published
crosslinks with the SAME 3.0 nm CA-CA criterion and the same buried/exposed
reactivity filter used for the CALVADOS side.

COST
----
HyRes carries an atomistic backbone (~5-6 particles/residue vs 1 for CALVADOS)
and runs dt = 0.004 ps vs CALVADOS's 0.01 ps, so ~14x the cost per ns. The
CALVADOS p9_dtx3l replicate is 19 min for 500 ns, which puts HyRes near 4.5 h.
Benchmark before committing a panel.

BLOCKER YOU HAVE TO CLEAR (needs a human)
-----------------------------------------
at2hyres expects a CHARMM-style PDB, and the documented route is CHARMM-GUI's
PDB Reader & Manipulator, which requires registration and interactive web use.
This script attempts a local substitute (pdbfixer/pdbfix_res.py) and reports
honestly whether it worked; if it did not, the CHARMM-GUI step has to be done
by hand on binding/p9_dtx3l/input/{parp9,dtx3l}.pdb.

Usage:  python prepare_hyres.py            # convert + write run/analysis scripts
"""
import shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE / 'hyres' / 'HyRes_GPU'
WORK = HERE / 'hyres' / 'runs'
SRC = HERE / 'binding' / 'p9_dtx3l' / 'input'
PY = '/home/sbali/miniconda3/envs/hyres/bin/python'

BOX_NM = 40.0
TEMP_K = 293
SALT_MM = 150


def convert(pdb_in, tag):
    """at2hyres, via the local pdbfix_res.py. Returns the hyres pdb or None."""
    out = WORK / 'prep'
    out.mkdir(parents=True, exist_ok=True)
    fixed = out / f'{tag}_fix.pdb'
    hyres = out / f'{tag}_hyres.pdb'
    fixer = REPO / 'at2hyres' / 'pdbfix_res.py'
    conv = REPO / 'at2hyres' / 'at2hyres.py'
    for script, args, label in ((fixer, [str(pdb_in), str(fixed), '0'], 'pdbfix_res'),
                                (conv, [str(fixed), str(hyres)], 'at2hyres')):
        if not script.is_file():
            print(f'  {tag}: {script.name} missing from the repo'); return None
        r = subprocess.run([PY, str(script)] + args, capture_output=True, text=True)
        if r.returncode != 0:
            print(f'  {tag}: {label} FAILED (rc={r.returncode})')
            print('   ', (r.stderr or r.stdout).strip().splitlines()[-1][:200])
            return None
    print(f'  {tag}: converted -> {hyres.name}')
    return hyres


RUN_TEMPLATE = '''#!/usr/bin/env python3
"""HyRes run, adapted from HyRes_GPU/run_OpenMM/run.nvt.py.

No restraints of any kind are applied -- see prepare_hyres.py for why there are
none to apply. Conditions match the CALVADOS p9_dtx3l runs: {box} nm cubic box,
{temp} K, {salt} mM salt.
"""
import sys
sys.path.insert(0, "{repo}/run_OpenMM")
BOX_NM, TEMP_K, SALT_MM = {box}, {temp}, {salt}
PDB, PSF = "{pdb}", "{psf}"
# The repo's run.nvt.py is a script, not a module, so it is executed with the
# parameters above already bound; edit there if the force-field setup changes.
exec(open("{repo}/run_OpenMM/run.nvt.py").read())
'''


def main():
    if not REPO.is_dir():
        raise SystemExit(f'{REPO} not found -- clone github.com/lslumass/HyRes_GPU first')
    if not Path(PY).is_file():
        print(f'NOTE {PY} not present yet (conda env "hyres" still building?)')
    WORK.mkdir(parents=True, exist_ok=True)

    print('Stage 0 conversion (single chains, the fold-stability gate):')
    ok = {}
    for tag, pdb in (('parp9', SRC / 'parp9.pdb'), ('dtx3l', SRC / 'dtx3l.pdb')):
        if not pdb.is_file():
            print(f'  {tag}: {pdb} missing'); continue
        ok[tag] = convert(pdb, tag)

    if not all(ok.values()):
        print('\nCONVERSION INCOMPLETE. The documented path needs a CHARMM-style PDB')
        print('from CHARMM-GUI (registration + interactive web), which cannot be')
        print('automated here. Do this by hand:')
        print('  1. charmm-gui.org -> Input Generator -> PDB Reader & Manipulator')
        print(f'  2. upload {SRC}/parp9.pdb and {SRC}/dtx3l.pdb')
        print('  3. terminal groups: NTER/CTER (standard), download "CHARMM PDB"')
        print(f'  4. drop them in {WORK/"prep"} and re-run this script')
        return

    print('\nNext: PSF via HyRes_GPU/generate_psf/psfgen_hyres.py, then Stage 0 runs.')
    print('Stage 1 (the complex) is gated on Stage 0 showing the folds survive.')


if __name__ == '__main__':
    main()
