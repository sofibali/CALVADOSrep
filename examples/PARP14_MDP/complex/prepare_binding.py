#!/usr/bin/env python3
"""
Prepare SMALL-BOX association simulations for the PARP14 / PARP9 / DTX3L pairs.

WHY: the existing complexes/ runs were built with BOX_MARGIN = 80 nm, giving
~94-110 nm boxes. Every one of them dissociates -- the chains start from the AF3
docked pose, come apart within 10-60 ns, and never find each other again, because
at ~1 molecule per (100 nm)^3 (1.7 uM) the re-encounter time is far beyond 100 ns.
Those runs therefore measure unbinding only.

THIS SCRIPT instead puts the chains in the smallest box that does not distort
them, starts them APART (topol 'grid', no docked pose, so no AF3 bias), and runs
long enough for many binding/unbinding cycles. The observable becomes the
equilibrium bound fraction, which is what distinguishes the ternary complex from
the individual homo- and heterodimers.

BOX SIZING. A chain must not interact with its own periodic image:
    L  >  (chain span)  +  cutoff_yu       [cutoff_yu = 4.0 nm]
Spans measured over the existing production trajectories (99th percentile):
    parp14  32.0 nm  -> L >= 36     parp9  20.8 nm -> L >= 25
    dtx3l   26.8 nm  -> L >= 31
Default L = 40 nm gives margin on the largest chain (parp14) and puts every
system at the SAME concentration, which is required to compare bound fractions
across systems. --box overrides; --box 50 is the box-size sensitivity control.

Concentration per chain = 1 / (L^3 * N_A):   40 nm -> 26 uM,  50 nm -> 13 uM.

PAIR COUNTING. A 1+1 heterodimer box and a 2x homodimer box each contain exactly
one potential pair, so their bound fractions are directly comparable. The ternary
box contains three pairs, each at the same per-partner concentration.

RESTRAINTS. Harmonic domain restraints (k=700) only -- these preserve the folds
(FNC 0.62-0.67). The Go-KH scheme is deliberately NOT offered here: a difference
map against its unrestrained control showed it manufactures a parp14:MD3 -
dtx3l:linker contact (+25 contacts/frame) that does not otherwise exist.

Usage:
    python prepare_binding.py --all                      # 7 systems, 40 nm, 20 reps
    python prepare_binding.py --all --box 50 --tag box50 # sensitivity control
    python prepare_binding.py --name ternary --components parp14:1 parp9:1 dtx3l:1
    bash binding/<name>/run.sh 8                         # 8 replicates in parallel
"""
import os, argparse, shutil, yaml, json, importlib.util
import numpy as np
from pathlib import Path

CWD = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('pc', str(CWD / 'prepare_complex.py'))
pc = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(pc)

OUT_DIR   = CWD / 'binding'
BOX_NM    = 40.0      # cubic edge; see BOX SIZING above
N_REP     = 20        # independent replicates (different RNG seeds)
PROD_NS   = 500.0     # ns per replicate
DT_PS     = 0.01      # stock CALVADOS timestep; his output_0.01dt_* runs completed at this
FRAME_PS  = 50.0      # trajectory frame spacing
CUTOFF_YU = 4.0       # nm, must match config cutoff_yu
N_A       = 6.02214076e23

SET_LAUNCHER = """#!/bin/bash
# Launch this system's rep-* replicates across GPUs (round-robin).
#   bash run.sh [PER_GPU]          GPU_LIST="0 1" bash run.sh 2
set -e
CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHON_EXE="${PYTHON_EXE:-/home/vle/.conda/envs/calvados/bin/python}"
read -r -a GPUS <<< "${GPU_LIST:-0 1 2 3 4 5 6 7}"
NGPU=${#GPUS[@]}; PER_GPU=${1:-3}; PARALLEL=$(( NGPU * PER_GPU ))
cd "$CWD"
i=0
for d in $(ls -d rep-* | sort -V); do echo "${GPUS[$((i % NGPU))]} $d"; i=$((i + 1)); done \
| xargs -P "$PARALLEL" -L1 bash -c '
    gpu="$0"; d="$1"
    cd "$d" || exit 1
    if compgen -G "*.dcd" > /dev/null; then echo "  $d: has trajectory, skipping"; exit 0; fi
    echo "  starting $d on GPU $gpu"
    CUDA_VISIBLE_DEVICES="$gpu" "$PYTHON_EXE" run.py > run.log 2>&1 || echo "  $d: FAILED"
    echo "  $d: done"
'
"""

MASTER_LAUNCHER = """#!/bin/bash
# Run every prepared binding set, one after another.
#   bash run_all.sh [PER_GPU]
set -e
CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for s in "$CWD"/*/run.sh; do
    echo "=== $(basename $(dirname $s)) ==="
    bash "$s" "${1:-3}"
done
"""

# AF3 complex prediction backing each set's DOCKED arm.
# None -> no AF3 complex model exists, so that set is separated-start only.
AF3_COMPLEX = {
    'ternary':    ('parp14_parp9_dtx3l', ['parp14', 'parp9', 'dtx3l']),
    'p14_p9':     ('parp14_parp9',       ['parp14', 'parp9']),
    'p14_dtx3l':  ('parp14_dtx3l',       ['parp14', 'dtx3l']),
    'p9_dtx3l':   ('parp9_dtx3l',        ['parp9',  'dtx3l']),
    'dtx3l_homo': ('dtx3l_dtx3l',        ['dtx3l',  'dtx3l']),
    'p14_homo':   None,
    'p9_homo':    None,
}

SETS = [
    ('ternary',    ['parp14:1', 'parp9:1', 'dtx3l:1']),
    ('p14_p9',     ['parp14:1', 'parp9:1']),
    ('p14_dtx3l',  ['parp14:1', 'dtx3l:1']),
    ('p9_dtx3l',   ['parp9:1',  'dtx3l:1']),
    ('p14_homo',   ['parp14:2']),
    ('p9_homo',    ['parp9:2']),
    ('dtx3l_homo', ['dtx3l:2']),
]

def parse_components(items):
    out = []
    for it in items:
        p, _, n = it.partition(':')
        out.append((p, int(n) if n else 1))
    return out

def monomer_span_nm(pdb):
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure('m', str(pdb))
    ca = np.array([a.coord for a in s[0].get_atoms() if a.name == 'CA'])
    return float((ca.max(0) - ca.min(0)).max()) / 10.0

def conc_uM(box_nm, n):
    vol_L = (box_nm * 1e-7) ** 3 * 1e-3            # nm -> cm; cm^3 -> L
    return n / (N_A * vol_L) * 1e6

def build_set(name, comp_items, box_nm, nreps, prod_ns, dt_ps, tag):
    comps = parse_components(comp_items)
    set_dir = OUT_DIR / (f'{name}_{tag}' if tag else name)
    shared = set_dir / 'input'; shared.mkdir(parents=True, exist_ok=True)
    with open(pc.DEFAULT_CONFIG) as f:
        base_config = yaml.safe_load(f)

    domains_yaml, meta, spans = {}, [], []
    for protein, nmol in comps:
        base = pc.base_protein(protein)
        job = pc.find_job_dir(base)
        best = pc.best_sample_per_seed(job / 'ranking_scores.csv')
        seed, (sample, score) = max(best.items(), key=lambda kv: kv[1][1])
        pc.cif_to_chain_pdbs(job / f'seed-{seed}_sample-{sample}' / 'model.cif',
                             [protein], shared)
        domains_yaml[protein] = pc.build_domain_list(protein)
        spans.append(monomer_span_nm(shared / f'{protein}.pdb'))
        meta.append({'protein': protein, 'nmol': nmol,
                     'monomer_model': str((job / f'seed-{seed}_sample-{sample}' / 'model.cif'))})
    with open(shared / 'domains.yaml', 'w') as f:
        yaml.dump(domains_yaml, f, default_flow_style=True)
    shutil.copy2(pc.FL_RESIDUES, shared / 'residues_CALVADOS3.csv')

    # split-KH custom restraints, written for every copy of every chain
    cres = shared / 'custom_restraints.txt'
    if cres.exists(): cres.unlink()
    n_custom = 0
    for protein, nmol in comps:
        base = pc.base_protein(protein)
        if base not in pc.INTERDOMAIN:
            continue
        ra, rb = pc.INTERDOMAIN[base]
        pairs = pc.generate_interdomain_pairs(shared / f'{protein}.pdb', ra, rb)
        with open(cres, 'a') as f:
            for copy in range(1, nmol + 1):
                for i, j, d_nm in sorted(pairs):
                    f.write(f'{protein} {copy} {i} | {protein} {copy} {j} | '
                            f'{d_nm:.3f} {pc.K_CUSTOM:.1f}\n')
        n_custom += len(pairs) * nmol
    has_custom = n_custom > 0

    steps = int(prod_ns * 1000 / dt_ps)
    wfreq = int(FRAME_PS / dt_ps)
    n_chain = sum(n for _, n in comps)
    max_span = max(spans)
    if box_nm <= max_span + CUTOFF_YU:
        raise SystemExit(f'{name}: box {box_nm} nm too small for span {max_span:.1f} nm '
                         f'(need > span + {CUTOFF_YU})')

    for rep in range(1, nreps + 1):
        sim = set_dir / f'rep-{rep}'; idir = sim / 'input'; idir.mkdir(parents=True, exist_ok=True)
        shared_files = ['domains.yaml', 'residues_CALVADOS3.csv'] + \
                       (['custom_restraints.txt'] if has_custom else []) + \
                       [f'{p}.pdb' for p, _ in comps]
        for fname in shared_files:
            link, src = idir / fname, shared / fname
            if not link.exists() and src.exists():
                link.symlink_to(os.path.relpath(src, idir))
        config = dict(base_config)
        config.update({
            'sysname': name, 'box': [box_nm, box_nm, box_nm],
            'temp': pc.TEMP, 'ionic': pc.IONIC, 'pH': pc.PH,
            'topol': 'grid',                 # chains placed APART, no docked pose
            'slab_eq': False, 'box_eq': False, 'bilayer_eq': False,
            'restart': 'checkpoint', 'frestart': 'restart.chk',
            'dt': dt_ps, 'steps': steps, 'wfreq': wfreq, 'steps_eq': 1000, 'runtime': 0,
            'cutoff_lj': 2.0, 'cutoff_yu': CUTOFF_YU,
            'platform': pc.PLATFORM, 'threads': pc.THREADS,
            'verbose': True, 'random_number_seed': rep * 1000,
        })
        if has_custom:
            config.update({'custom_restraints': True,
                           'custom_restraint_type': 'harmonic',
                           'fcustom_restraints': 'input/custom_restraints.txt'})
        with open(sim / 'config.yaml', 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        components_dict = {
            'defaults': {
                'molecule_type': 'protein', 'nmol': 1, 'charge_termini': 'both',
                'alpha': 0, 'ffasta': 'fastabib.fasta', 'kb': 8033.0,
                'ext_restraint': True, 'restraint': True, 'cutoff_restr': pc.CUTOFF_RESTR,
                'pdb_folder': 'input', 'restraint_type': 'harmonic',
                'k_harmonic': pc.K_HARMONIC, 'fdomains': 'input/domains.yaml',
                'k_go': pc.K_GO, 'use_com': True, 'periodic': False, 'colabfold': 1,
                'bfac_shift': 0.8, 'bfac_width': 50.0, 'pae_shift': 0.3, 'pae_width': 15.0,
                'rna_kb1': 8033.0, 'rna_kb2': 8033.0, 'rna_ka': 7.24, 'rna_pa': 3.14,
                'rna_nb_sigma': 0.4, 'rna_nb_scale': 15, 'rna_nb_cutoff': 0.6,
                'n_ends': 1, 'ptm_name': 'example_ptm', 'ptm_locations': [],
                'fresidues': 'input/residues_CALVADOS3.csv',
            },
            'system': {p: {'nmol': n} for p, n in comps},
        }
        with open(sim / 'components.yaml', 'w') as f:
            yaml.dump(components_dict, f, default_flow_style=False, sort_keys=False)
        pc.write_run_py(sim)
    lf = set_dir / 'run.sh'
    lf.write_text(SET_LAUNCHER)
    lf.chmod(0o755)
    c = conc_uM(box_nm, 1)
    with open(set_dir / 'binding_summary.json', 'w') as f:
        json.dump({'name': name, 'box_nm': box_nm, 'n_chains': n_chain,
                   'conc_uM_per_chain': round(c, 1), 'prod_ns': prod_ns, 'dt_ps': dt_ps,
                   'nreps': nreps, 'topol': 'grid', 'start': 'separated',
                   'restraints': f'harmonic k={pc.K_HARMONIC}',
                   'n_custom_restraint_pairs': n_custom, 'components': meta}, f, indent=2)
    print(f'  {set_dir.name:20s} {"+".join(f"{p}x{n}" for p,n in comps):28s} '
          f'{n_chain} chains  box {box_nm:.0f} nm  {c:.0f} uM/chain  '
          f'{nreps}x{prod_ns:.0f} ns  (max span {max_span:.1f} nm)')

def build_docked(name, box_nm, nreps, prod_ns, dt_ps, tag):
    """DOCKED arm: same box, same length, but chains start from the AF3 complex
    pose (restart='pdb'). Paired with the separated arm this brackets the bound
    fraction: if the two arms agree the result is thermodynamic, if they disagree
    it is hysteresis -- exactly the test the HyRes docked/separated pair failed."""
    entry = AF3_COMPLEX.get(name)
    if entry is None:
        print(f'  {name:20s} SKIPPED - no AF3 complex model, separated-start only')
        return
    job_name, proteins = entry
    comp_names = pc.uniquify(proteins)
    set_dir = OUT_DIR / (f'{name}_docked_{tag}' if tag else f'{name}_docked')
    shared = set_dir / 'input'; shared.mkdir(parents=True, exist_ok=True)
    with open(pc.DEFAULT_CONFIG) as f:
        base_config = yaml.safe_load(f)
    job = pc.find_job_dir(job_name)
    best = pc.best_sample_per_seed(job / 'ranking_scores.csv')
    models = [job / f'seed-{sd}_sample-{sm}' / 'model.cif'
              for sd, (sm, _) in sorted(best.items())]

    domains_yaml = {c: pc.build_domain_list(c) for c in comp_names}
    with open(shared / 'domains.yaml', 'w') as f:
        yaml.dump(domains_yaml, f, default_flow_style=True)
    shutil.copy2(pc.FL_RESIDUES, shared / 'residues_CALVADOS3.csv')

    steps = int(prod_ns * 1000 / dt_ps); wfreq = int(FRAME_PS / dt_ps)
    spans = []
    for rep in range(1, nreps + 1):
        sim = set_dir / f'rep-{rep}'; idir = sim / 'input'; idir.mkdir(parents=True, exist_ok=True)
        cif = models[(rep - 1) % len(models)]          # spread AF3 poses over replicates
        pc.cif_to_chain_pdbs(cif, comp_names, idir)
        edge = pc.write_system_pdb(cif, comp_names, idir / 'system.pdb', box_abs=box_nm)
        spans.append(edge)
        meta_c = pc.write_custom_restraints(idir, comp_names)
        has_custom = (idir / 'custom_restraints.txt').exists()
        for fname in ['domains.yaml', 'residues_CALVADOS3.csv']:
            link, src = idir / fname, shared / fname
            if not link.exists() and src.exists():
                link.symlink_to(os.path.relpath(src, idir))
        config = dict(base_config)
        config.update({
            'sysname': name, 'box': [box_nm, box_nm, box_nm],
            'temp': pc.TEMP, 'ionic': pc.IONIC, 'pH': pc.PH,
            'topol': 'grid', 'restart': 'pdb', 'frestart': 'input/system.pdb',
            'slab_eq': False, 'box_eq': False, 'bilayer_eq': False,
            'dt': dt_ps, 'steps': steps, 'wfreq': wfreq, 'steps_eq': 1000, 'runtime': 0,
            'cutoff_lj': 2.0, 'cutoff_yu': CUTOFF_YU,
            'platform': pc.PLATFORM, 'threads': pc.THREADS,
            'verbose': True, 'random_number_seed': rep * 1000,
        })
        if has_custom:
            config.update({'custom_restraints': True,
                           'custom_restraint_type': 'harmonic',
                           'fcustom_restraints': 'input/custom_restraints.txt'})
        with open(sim / 'config.yaml', 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        components_dict = {
            'defaults': {
                'molecule_type': 'protein', 'nmol': 1, 'charge_termini': 'both',
                'alpha': 0, 'ffasta': 'fastabib.fasta', 'kb': 8033.0,
                'ext_restraint': True, 'restraint': True, 'cutoff_restr': pc.CUTOFF_RESTR,
                'pdb_folder': 'input', 'restraint_type': 'harmonic',
                'k_harmonic': pc.K_HARMONIC, 'fdomains': 'input/domains.yaml',
                'k_go': pc.K_GO, 'use_com': True, 'periodic': False, 'colabfold': 1,
                'bfac_shift': 0.8, 'bfac_width': 50.0, 'pae_shift': 0.3, 'pae_width': 15.0,
                'rna_kb1': 8033.0, 'rna_kb2': 8033.0, 'rna_ka': 7.24, 'rna_pa': 3.14,
                'rna_nb_sigma': 0.4, 'rna_nb_scale': 15, 'rna_nb_cutoff': 0.6,
                'n_ends': 1, 'ptm_name': 'example_ptm', 'ptm_locations': [],
                'fresidues': 'input/residues_CALVADOS3.csv',
            },
            'system': {c: {} for c in comp_names},
        }
        with open(sim / 'components.yaml', 'w') as f:
            yaml.dump(components_dict, f, default_flow_style=False, sort_keys=False)
        pc.write_run_py(sim)
    lf = set_dir / 'run.sh'; lf.write_text(SET_LAUNCHER); lf.chmod(0o755)
    with open(set_dir / 'binding_summary.json', 'w') as f:
        json.dump({'name': name, 'arm': 'docked', 'af3_job': str(job), 'box_nm': box_nm,
                   'n_chains': len(comp_names), 'conc_uM_per_chain': round(conc_uM(box_nm,1),1),
                   'prod_ns': prod_ns, 'dt_ps': dt_ps, 'nreps': nreps,
                   'start': 'af3_docked', 'components': comp_names,
                   'af3_models': [str(m) for m in models]}, f, indent=2)
    print(f'  {set_dir.name:26s} {"+".join(comp_names):28s} {len(comp_names)} chains  '
          f'box {box_nm:.0f} nm  {len(models)} AF3 poses over {nreps} reps')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--all', action='store_true', help='prepare all 7 systems')
    ap.add_argument('--name'); ap.add_argument('--components', nargs='+')
    ap.add_argument('--box', type=float, default=BOX_NM)
    ap.add_argument('--nreps', type=int, default=N_REP)
    ap.add_argument('--ns', type=float, default=PROD_NS)
    ap.add_argument('--dt', type=float, default=DT_PS)
    ap.add_argument('--tag', default='')
    ap.add_argument('--start', choices=['separated','docked','both'],
                    default='separated',
                    help="'separated' (default) = grid start, chains apart; "
                         "'docked' = start from the AF3 complex pose in the same "
                         "small box; 'both' = prepare both arms for comparison.")
    a = ap.parse_args()
    OUT_DIR.mkdir(exist_ok=True)
    sets = SETS if a.all else [(a.name, a.components)]
    print(f'box {a.box:.0f} nm -> {conc_uM(a.box,1):.0f} uM per chain; '
          f'{a.nreps} reps x {a.ns:.0f} ns at dt={a.dt} ps')
    if a.start in ('separated', 'both'):
        print('-- separated-start arm (chains placed apart, must associate)')
        for name, comps in sets:
            build_set(name, comps, a.box, a.nreps, a.ns, a.dt, a.tag)
    if a.start in ('docked', 'both'):
        print('-- docked-start arm (AF3 complex pose, may dissociate)')
        for name, _ in sets:
            build_docked(name, a.box, a.nreps, a.ns, a.dt, a.tag)
    mf = OUT_DIR / 'run_all.sh'
    mf.write_text(MASTER_LAUNCHER); mf.chmod(0o755)
    print(f'\nprepared -> {OUT_DIR}')
    print(f'launch one set: bash {OUT_DIR.name}/<set>/run.sh 3')
    print(f'launch all:     bash {OUT_DIR.name}/run_all.sh 3')

if __name__ == '__main__':
    main()
