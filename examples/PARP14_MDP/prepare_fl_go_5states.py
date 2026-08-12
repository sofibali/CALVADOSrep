#!/usr/bin/env python
"""
Set up 5 long full-length PARP14 runs (2 us each -> 10 us total) seeded from the
5 MOST DISTINCT conformations among the 25 fl_optimized replicates, with the KH7a-KHb
"split domain" restraints switched from custom HARMONIC (k=350) to Go contacts at the
CALVADOS standard force (k_go = 15).

WHY FRESH RUNS (not checkpoint-continue): switching the custom restraint from harmonic
to Go changes the OpenMM System's forces, so an existing checkpoint cannot be continued
(loadCheckpoint requires an identical System). We therefore start fresh runs whose
INITIAL COORDINATES are each selected replicate's current conformation (its
checkpoint.pdb, CG, Angstrom) via restart='pdb', while the restraint REFERENCE stays
the folded AF2 structure (input/parp14.pdb) so domain + Go native distances are native.

STATE SELECTION: fl_optimized's 25 reps all started from the SAME AF2 structure
(input/parp14.pdb) and diverged only via random seed over 100 ns -- they were NOT
started from per-replicate AF3 rank-0 structures. So "5 most distinct states" is
computed by farthest-point sampling on each replicate's CURRENT inter-domain
arrangement (11-unit COM-COM distance fingerprint of its checkpoint.pdb).

RESTRAINTS:
  - domain (intra-domain) restraints: unchanged -- harmonic, k_harmonic = 700
  - KH7a-KHb custom (inter-domain) restraints: Go, standard k_go = 15
    (custom_restraint_type: go; native distances kept from the AF2-derived pairs)

Output: fl_go/state-{1..5}_from_{srcrep}/ , 2 us each. Does NOT launch.

  python prepare_fl_go_5states.py            # 5 most distinct, 2 us each
  python prepare_fl_go_5states.py --n 5 --us-per-state 2
  bash fl_go/run_all.sh parallel             # launch when ready (heavy, multi-day)
"""
import os, sys, glob, argparse, shutil, warnings
warnings.simplefilter('ignore')
import numpy as np, yaml
import MDAnalysis as mda

CWD = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(CWD, 'fl_optimized')
OUT = os.path.join(CWD, 'fl_go')
STEPS_PER_NS = 100000                 # dt = 0.01 ps
K_GO = 15.0                           # CALVADOS standard Go force constant (kJ/mol)

# TICA/PCA state selection reuses the clustering engine (in tica_pipeline/)
sys.path.insert(0, os.path.join(CWD, 'tica_pipeline'))

# FL domain units (full-length numbering) for the distinctness fingerprint
UNITS = {'RRM1': (1,145),'RRM2':(146,224),'RRM3':(225,314),'KH1-KH6':(315,737),
         'KH7a':(738,789),'MD1L1':(790,1004),'MD2':(1004,1193),'MD3':(1207,1388),
         'KHb-KH8':(1389,1533),'WWE':(1534,1602),'ART':(1603,1801)}


def fingerprint(pdb):
    """11-unit COM-COM distance vector (nm) of a rep's conformation."""
    u = mda.Universe(pdb)
    coms = []
    for lo, hi in UNITS.values():
        ag = u.select_atoms(f'name CA and resid {lo}:{hi}')
        coms.append(ag.center_of_mass() / 10.0)
    coms = np.array(coms)
    d = []
    for i in range(len(coms)):
        for j in range(i + 1, len(coms)):
            d.append(np.linalg.norm(coms[i] - coms[j]))
    return np.array(d)


def farthest_point(F, n):
    """Indices of n most mutually-distinct rows (max-min), seeded at the row
    farthest from the mean."""
    Z = (F - F.mean(0)) / (F.std(0) + 1e-9)
    start = int(np.argmax(np.linalg.norm(Z - Z.mean(0), axis=1)))
    chosen = [start]
    dmin = np.linalg.norm(Z - Z[start], axis=1)
    while len(chosen) < n:
        nxt = int(np.argmax(dmin))
        chosen.append(nxt)
        dmin = np.minimum(dmin, np.linalg.norm(Z - Z[nxt], axis=1))
    return chosen


def select_landscape(method, features, ca_stride, tica_lag, n, workers):
    """Pick n REPRESENTATIVE states from the pooled fl_optimized landscape (TICA or
    PCA) using farthest-point spread + Voronoi tiling. Returns a list of
    (seed, sample, frame) for the n representative frames (any frame from any rep)."""
    from cluster_states import (collect_features, compute_tica, find_spread_frames)
    ALL = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a', 'MD1L1', 'MD2', 'MD3',
           'KHb-KH8', 'WWE', 'ART']
    print(f"Featurizing fl_optimized (pooled) feature='{features}' stride={ca_stride} ...")
    d = collect_features('fl_optimized', ALL, False, workers,
                         features_mode=features, ca_stride=ca_stride)
    X = d['features']; meta = d['metadata']; rb = d['rep_boundaries']
    print(f"  {X.shape[0]} frames x {X.shape[1]} features across {len(rb)} replicates")
    if method == 'tica':
        print(f"  TICA lag {tica_lag} frames ({tica_lag*0.01:.1f} ns) ...")
        Y, _, _ = compute_tica(X, rb, lag=tica_lag, n_components=10)
    else:  # pca
        print("  PCA (SVD on standardized features) ...")
        Z = (X - X.mean(0)) / (X.std(0) + 1e-9)
        _, _, Vt = np.linalg.svd(Z, full_matrices=False)
        Y = Z @ Vt[:10].T
    idx, _ = find_spread_frames(Y[:, :2], n)
    return [tuple(int(v) for v in meta[i]) for i in idx]


def write_run_py(d):
    open(os.path.join(d, 'run.py'), 'w').write(
        'from calvados import sim\nfrom argparse import ArgumentParser\n\n'
        'if __name__ == "__main__":\n'
        '    p = ArgumentParser()\n'
        "    p.add_argument('--path', nargs='?', default='.', const='.', type=str)\n"
        "    p.add_argument('--config', nargs='?', default='config.yaml', const='config.yaml', type=str)\n"
        "    p.add_argument('--components', nargs='?', default='components.yaml', const='components.yaml', type=str)\n"
        '    a = p.parse_args()\n'
        '    sim.run(path=a.path, fconfig=a.config, fcomponents=a.components)\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=5, help='Number of states (default 5)')
    ap.add_argument('--us-per-state', type=float, default=2.0, help='us per state (default 2)')
    ap.add_argument('--select', choices=['tica', 'pca', 'distinct'], default='tica',
                    help="How to pick the n starting states: 'tica'/'pca' = representative "
                         "states from the pooled fl_optimized landscape (spread selection); "
                         "'distinct' = farthest-point on per-rep checkpoint conformations.")
    ap.add_argument('--features', default='ca',
                    help="Feature for tica/pca selection (ca | pose | interface_ca | com; default ca)")
    ap.add_argument('--ca-stride', type=int, default=25)
    ap.add_argument('--tica-lag', type=int, default=200, help='TICA lag in frames (default 200 = 2 ns)')
    ap.add_argument('--threads', type=int, default=16,
                    help='OpenMM CPU threads PER run (default 16; 5 runs x 16 = 80 cores)')
    ap.add_argument('--workers', type=int, default=16, help='parallel workers for featurization')
    a = ap.parse_args()
    steps = int(round(a.us_per_state * 1000 * STEPS_PER_NS))   # 2 us -> 2e8

    # 1. + 2. pick the n starting states -> each yields a way to write restart.pdb
    #    states[k] = (label, kind, payload) ; kind 'frame'->(seed,sample,frame), 'ckpt'->path
    reps = sorted(glob.glob(os.path.join(SRC, 'seed-*_sample-*')))
    reps = [r for r in reps if os.path.isfile(os.path.join(r, 'checkpoint.pdb'))]
    if len(reps) < a.n:
        print(f"Only {len(reps)} reps with checkpoint.pdb (< {a.n})"); return
    states = []
    if a.select in ('tica', 'pca'):
        frames = select_landscape(a.select, a.features, a.ca_stride, a.tica_lag, a.n, a.workers)
        print(f"\n{a.n} representative states ({a.select.upper()} landscape, spread):")
        for k, (s, sm, fr) in enumerate(frames, 1):
            print(f"  state-{k}: seed-{s}_sample-{sm} frame {fr} ({fr*0.01:.1f} ns)")
            states.append((f"state-{k}_{a.select}_seed-{s}_sample-{sm}_fr{fr}", 'frame', (s, sm, fr)))
    else:
        print(f"Fingerprinting {len(reps)} replicates (11-unit COM distances)...")
        F = np.array([fingerprint(os.path.join(r, 'checkpoint.pdb')) for r in reps])
        sel = farthest_point(F, a.n)
        print(f"\n{a.n} most distinct states (farthest-point on checkpoints):")
        for k, i in enumerate(sel, 1):
            print(f"  state-{k}: {os.path.basename(reps[i])}")
            states.append((f"state-{k}_from_{os.path.basename(reps[i])}", 'ckpt',
                           os.path.join(reps[i], 'checkpoint.pdb')))

    # 3. shared input: AF2 reference + domains + residues + Go custom restraints
    shared = os.path.join(OUT, 'input'); os.makedirs(shared, exist_ok=True)
    for f in ['parp14.pdb', 'domains.yaml', 'residues_CALVADOS3.csv']:
        s = os.path.join(SRC, 'input', f)
        if os.path.isfile(s) and not os.path.exists(os.path.join(shared, f)):
            shutil.copy2(os.path.realpath(s), os.path.join(shared, f))
    # Go custom restraints: same native distances, k -> standard k_go
    go_lines = []
    for ln in open(os.path.join(SRC, 'input', 'custom_restraints.txt')):
        parts = ln.split('|')
        r = parts[2].split()[0]
        go_lines.append(f"{parts[0].strip()} | {parts[1].strip()} | {r} {K_GO}\n")
    open(os.path.join(shared, 'custom_restraints_go.txt'), 'w').writelines(go_lines)
    print(f"\nWrote {len(go_lines)} KH7a-KHb Go restraints (k_go={K_GO}) -> "
          f"{shared}/custom_restraints_go.txt")

    # 4. one fresh 2-us sim per selected state
    from cluster_states import extract_pdb_for_frame
    for k, (name, kind, payload) in enumerate(states, 1):
        d = os.path.join(OUT, name); idir = os.path.join(d, 'input')
        os.makedirs(idir, exist_ok=True)
        for f in ['parp14.pdb', 'domains.yaml', 'residues_CALVADOS3.csv', 'custom_restraints_go.txt']:
            dest = os.path.join(idir, f)
            if not os.path.exists(dest):
                os.symlink(os.path.join(shared, f), dest)
        # initial coordinates = this state's conformation (CG, Angstrom)
        restart = os.path.join(d, 'restart.pdb')
        if kind == 'ckpt':
            shutil.copy2(payload, restart)
        else:  # 'frame' -> extract representative frame from its source replicate
            s, sm, fr = payload
            if not extract_pdb_for_frame('fl_optimized', s, sm, fr, restart):
                print(f"  WARN: could not extract frame for {name}; skipping"); continue

        # start from the known-good fl_optimized config (all default keys present),
        # then override only what changes for a fresh Go-restraint 2-us run
        config = yaml.safe_load(open(os.path.join(SRC, 'seed-1_sample-0', 'config.yaml')))
        config.update({
            'sysname': 'parp14', 'box': [300, 300, 300],
            'steps': steps, 'wfreq': 1000, 'runtime': 0, 'threads': a.threads,
            'restart': 'pdb', 'frestart': 'restart.pdb',      # fresh, seeded from this state
            'custom_restraints': True,
            'custom_restraint_type': 'go',                    # KH7a-KHb -> Go
            'fcustom_restraints': 'input/custom_restraints_go.txt',
            'random_number_seed': 1000 + k,
        })
        yaml.dump(config, open(os.path.join(d, 'config.yaml'), 'w'), default_flow_style=False)

        components = {'defaults': {
            'molecule_type': 'protein', 'nmol': 1, 'charge_termini': 'both', 'alpha': 0,
            'ffasta': 'fastabib.fasta', 'kb': 8033.0,
            'ext_restraint': True, 'restraint': True, 'cutoff_restr': 0.9,
            'pdb_folder': str(idir),
            'restraint_type': 'harmonic', 'k_harmonic': 700.0,   # domain restraints unchanged
            'fdomains': str(os.path.join(idir, 'domains.yaml')),
            'k_go': K_GO, 'use_com': True, 'periodic': False, 'colabfold': 1,
            'bfac_shift': 0.8, 'bfac_width': 50.0, 'pae_shift': 0.3, 'pae_width': 15.0,
            'rna_kb1': 8033.0, 'rna_kb2': 8033.0, 'rna_ka': 7.24, 'rna_pa': 3.14,
            'rna_nb_sigma': 0.4, 'rna_nb_scale': 15, 'rna_nb_cutoff': 0.6,
            'n_ends': 1, 'ptm_name': 'example_ptm', 'ptm_locations': [],
            'fresidues': str(os.path.join(idir, 'residues_CALVADOS3.csv')),
        }, 'system': {'parp14': {}}}
        yaml.dump(components, open(os.path.join(d, 'components.yaml'), 'w'), default_flow_style=False)
        write_run_py(d)
        print(f"  prepared {name}  ({a.us_per_state:g} us, {steps} steps)")

    # 5. launcher
    run_sh = os.path.join(OUT, 'run_all.sh')
    open(run_sh, 'w').write(
        '#!/bin/bash\n# Launch the 5 fl_go 2-us states. Arg: "parallel" or "serial" (default).\n'
        'set -e\nHERE="$(cd "$(dirname "$0")" && pwd)"\nMODE="${1:-serial}"\n'
        'PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python\n'
        'for d in "$HERE"/state-*/; do\n'
        '  [ -f "$d/run.py" ] || continue\n'
        '  echo "RUN $d"\n'
        '  if [ "$MODE" = parallel ]; then ( cd "$d" && $PY run.py ) & else ( cd "$d" && $PY run.py ); fi\n'
        'done\n[ "$MODE" = parallel ] && wait\necho "fl_go done."\n')
    os.chmod(run_sh, 0o755)
    print(f"\nLauncher: {run_sh}  ->  bash fl_go/run_all.sh parallel")
    print(f"Total: {a.n} x {a.us_per_state:g} us = {a.n * a.us_per_state:g} us | "
          f"{a.threads} threads/run x {a.n} runs = {a.threads * a.n} cores in parallel.")
    print("Not launched (heavy multi-day job). Tune parallelism with --threads.")


if __name__ == '__main__':
    main()
