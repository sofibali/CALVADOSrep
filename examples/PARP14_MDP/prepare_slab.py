#!/usr/bin/env python
"""
Prepare PARP14 slab (direct-coexistence) simulations for the puncta / phase-
separation arm of the project.

WHY
---
Every PARP14 simulation to date is SINGLE-CHAIN. Self-association is by
definition a multi-chain property, so nothing currently on disk can predict the
diffuse-vs-punctate sort-seq screen. Slab direct-coexistence gives the actual
quantitative observable -- the saturation concentration c_sat -- via
`calvados.analysis.SlabAnalysis`.

TWO ARMS
--------
  homotypic   N copies of one PARP14 construct.
              Baseline: does this construct self-associate at all?
  rna         the same, plus polyU RNA.
              Tests the RNA-multivalency leg of the puncta hypothesis.

An ADPr-substrate arm is deliberately NOT built here -- see LIMITATIONS.

WHERE THIS RUNS
---------------
NOT on pollux (no GPU: nvidia-smi fails, OpenMM sees only Reference/CPU). The
measured CPU throughput from this project's own md_full runs is ~5.8e5
bead-steps/s, which puts a single 100-chain slab at roughly TWO YEARS. These
inputs target the SLURM GPU server that shares this filesystem (same host the
BindCraft campaign targets; see bindcraft_md/GPU_SERVER_SETUP.md).

Because no GPU is reachable from here, the per-construct cost below is an
ESTIMATE from a CPU-derived scaling argument, not a measurement. Run
`--benchmark` FIRST on the GPU; it writes a real ns/day number that
`--report` then uses to cost the panel honestly.

USAGE
-----
    python prepare_slab.py --benchmark          # tiny calibration run, do this first
    python prepare_slab.py --report             # panel + cost table, writes nothing
    python prepare_slab.py --arm homotypic      # prepare the homotypic arm
    python prepare_slab.py --arm rna            # prepare the +RNA arm
    python prepare_slab.py --arm both --panel core
    # then, on the GPU server:
    sbatch --array=0-5 slab/submit_slab.slurm slab/homotypic

LIMITATIONS (verified in source, not assumed)
---------------------------------------------
* RNA is CALVADOS's generic 2-bead-per-nucleotide model with NO base identity
  (components.py:394-396); polyU is the only meaningful choice. Its bead
  parameters (RBC/RNA rows) were fit alongside CALVADOS2, while the PARP14
  protein parameters here are CALVADOS3 -- mixing the two is a real, unresolved
  approximation and any RNA-arm result inherits it.
* NO ADPr-substrate arm. `PTMProtein` has no example and no test in this repo;
  it sets c_termini to the last PTM bead rather than the protein C-terminus
  (components.py:737), and it never reads from PDB, so it is incompatible with
  `restraint: True` -- which every multi-domain PARP14 construct requires.
  Additionally no ADP-ribose bead parameters exist anywhere and would have to
  be derived. The viable route is a short UNRESTRAINED substrate peptide
  carrying the ADPr beads, which sidesteps the restraint bug; that is a
  separate piece of work, not a config change.
* calvados/analysis.py:772 returns before the "NOT CONVERGED" check at 774-776,
  so a failed tanh fit never warns and c_sat can come back silently wrong.
  `analyze_slab.py` re-implements that check.
"""
import os
import shutil
import csv as _csv
from argparse import ArgumentParser

import yaml

from calvados.cfg import Config, Components

CWD = os.path.dirname(os.path.abspath(__file__))
SLAB_ROOT = os.path.join(CWD, 'slab')

TEMP = 293.15      # K, matches the orig_examples slab protocol
IONIC = 0.15       # M
PH = 7.5
PLATFORM = 'CUDA'

# RNA bead parameters, from examples/orig_examples/slab_mixed/input/residues_C2RNA.csv
# (and single_dsRNA for the structured 's' bead). Columns are read by NAME
# (components.py:28 does read_csv(...).set_index('one')), so column ORDER in the
# merged file does not matter -- only that every name is present.
RNA_BEADS = [
    # one, three, MW,    lambdas, sigmas, q,  bondlength
    ('p', 'RBC', 194.1, 0.00,   0.6954, -1, 0.59),   # phosphate / backbone
    ('r', 'RNA', 126.3, 1.18,   0.6238,  0, 0.54),   # unstructured base
    ('s', 'SRN', 126.3, 0.13,   0.6238,  0, 0.54),   # structured base
]

# Published RNA force-field settings (single_RNA/prepare.py, slab_mixed/prepare.py)
RNA_OPTS = dict(rna_kb1=1400.0, rna_kb2=2200.0, rna_ka=4.20, rna_pa=3.14,
                rna_nb_sigma=0.4, rna_nb_scale=15, rna_nb_cutoff=2.0)

RNA_NAME = 'polyU40'
RNA_NT = 40                 # nucleotides; 40 matches slab_mixed's polyU40
RNA_PER_PROTEIN = 0.30      # molar ratio, ~ slab_mixed's 60 RNA : 200 protein

# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------
# Chosen to span the two valences the multivalency hypothesis turns on, using
# ONLY constructs that already have a prepared structure + domains.yaml on disk
# (so no new AF3 inference is needed -- that also requires a GPU).
#
#   rna_val   number of individual RNA-binding domains (RRMs + KHs)
#   adpr_val  number of ADPr-reader modules (MD2, MD3, WWE)
#
# 'core' is the minimum panel that still contains a +/-ART matched pair at each
# end of the RNA-valence range; 'full' adds the intermediate points.
PANEL = [
    # key,               tier,   note
    ('md2_md3',          'core', 'RNA-null baseline, no eraser, no writer'),
    ('md_full',          'full', 'RNA-null + eraser'),
    ('md3_wwe_full',     'full', 'low RNA valence, small'),
    ('mka_wwe_full',     'core', '-ART  of matched pair A'),
    ('mka_full',         'core', '+ART  of matched pair A'),
    ('core_full_go',     'full', 'mid RNA valence'),
    ('kh1_wwe_full',     'full', '-ART  of matched pair B'),
    ('kh1_art_full',     'full', '+ART  of matched pair B'),
    ('fl_wwe_full_go',   'core', '-ART  of matched pair C'),
    ('fl',               'core', '+ART  of matched pair C, full length'),
]

# Cost model. TARGET_BEADS keeps every run in the same cost bracket by trading
# chain count against construct size, so one construct does not dominate the
# queue. MIN_CHAINS is a floor: below ~30 chains a slab is too thin for the
# tanh interface fit in SlabAnalysis to be trustworthy.
TARGET_BEADS = 50_000
MIN_CHAINS, MAX_CHAINS = 30, 100
N_STEPS = 200_000_000        # 2 us at dt = 0.01 ps
WFREQ = 100_000
STEPS_EQ = 5_000_000         # slab-centering equilibration

# CPU throughput measured from this project's own md_full runs:
#   599 beads * 9e7 steps / (26 h) = 5.76e5 bead-steps/s
CPU_BEAD_STEPS_PER_S = 5.76e5
# Conservative GPU assumption ONLY used until --benchmark provides a real number.
ASSUMED_GPU_SPEEDUP = 50


def load_registry():
    import sys
    sys.path.insert(0, CWD)
    import sim_registry as reg
    return reg


def _first_existing(*paths):
    for p in paths:
        if p and os.path.isfile(p):
            return p
    return paths[0] if paths else None


def construct_info(reg, key):
    """Residue count, sysname and input paths for a prepared construct.

    The sets do not all store their structure in the same place: most keep
    `{set}/input/{sysname}.pdb`, but md_full only has `ref_allatom.pdb` there
    and keeps the real CG structure under a replicate's own input dir. Look in
    both rather than special-casing, since new sets may follow either layout.
    """
    v = reg.SETS[key]
    sysname = v['sysname']
    inp = os.path.join(CWD, key, 'input')
    rep_inp = os.path.join(CWD, key, 'seed-1_sample-0', 'input')
    pdb = _first_existing(os.path.join(inp, f'{sysname}.pdb'),
                          os.path.join(rep_inp, f'{sysname}.pdb'))
    dom = _first_existing(os.path.join(inp, 'domains.yaml'),
                          os.path.join(rep_inp, 'domains.yaml'))
    res = _first_existing(os.path.join(inp, 'residues_CALVADOS3.csv'),
                          os.path.join(rep_inp, 'residues_CALVADOS3.csv'))
    nres = None
    if os.path.isfile(pdb):
        rs = set()
        with open(pdb) as fh:
            for line in fh:
                if line.startswith(('ATOM', 'HETATM')):
                    rs.add(int(line[22:26]))
        nres = len(rs)
    return dict(key=key, sysname=sysname, pdb=pdb, domains=dom, residues=res,
                nres=nres, ok=all(os.path.isfile(p) for p in (pdb, dom, res)))


def n_chains_for(nres):
    n = round(TARGET_BEADS / nres)
    return int(max(MIN_CHAINS, min(MAX_CHAINS, n)))


def box_for(nres, nchain):
    """Slab box [Lx, Lx, Lz] in nm.

    Lx must comfortably exceed twice the chain's radius of gyration or a chain
    interacts with its own periodic image across the slab. Rg is estimated as
    ~0.25*sqrt(N) nm, which for these compact multi-domain constructs lands in
    the 4-7 nm range and matches the boxes the orig_examples slabs use for
    comparable sizes. Lz is held at 13x Lx, the aspect ratio of slab_MDP
    ([20,20,270]), which leaves enough dilute phase for the tanh fit.
    """
    rg = 0.25 * nres ** 0.5
    lx = max(20.0, 4.0 * rg)
    lx = float(round(lx / 5.0) * 5)      # round to 5 nm
    return [lx, lx, round(13.0 * lx)]


def merged_residues_csv(src_c3, dst):
    """CALVADOS3 amino acids + the RNA bead rows, written as one file.

    The PARP14 sets use residues_CALVADOS3.csv, which has no RNA bead types;
    slab_mixed's residues_C2RNA.csv has them but carries CALVADOS2 protein
    lambdas. Keep the C3 protein rows (what every PARP14 run so far used) and
    append only the RNA rows.
    """
    with open(src_c3) as fh:
        rows = list(_csv.DictReader(fh))
    cols = list(rows[0].keys())
    have = {r['one'] for r in rows}
    for one, three, mw, lam, sig, q, bl in RNA_BEADS:
        if one in have:
            continue
        rows.append({'one': one, 'three': three, 'MW': mw, 'lambdas': lam,
                     'sigmas': sig, 'q': q, 'bondlength': bl})
    with open(dst, 'w', newline='') as fh:
        w = _csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return dst


def write_rna_fasta(path):
    """Unstructured polyU. CALVADOS's RNA model has no base identity, so the
    letter is the generic unstructured-base bead 'r' (components.py:389-392),
    exactly as slab_mixed/input/mix.fasta does."""
    with open(path, 'w') as fh:
        fh.write(f'>{RNA_NAME}\n{"r" * RNA_NT}\n')
    return path


RUN_PY = '''import os
from calvados import sim

path = os.path.dirname(os.path.abspath(__file__))
sim.run(path=path, fconfig='config.yaml', fcomponents='components.yaml')
'''


def prepare_one(reg, key, arm, steps, outroot, note=''):
    info = construct_info(reg, key)
    if not info['ok']:
        missing = [p for p in (info['pdb'], info['domains'], info['residues'])
                   if not os.path.isfile(p)]
        return None, f"missing input: {', '.join(os.path.relpath(m, CWD) for m in missing)}"

    nres = info['nres']
    nchain = n_chains_for(nres)
    box = box_for(nres, nchain)
    n_rna = int(round(nchain * RNA_PER_PROTEIN)) if arm == 'rna' else 0

    d = os.path.join(outroot, key)
    inp = os.path.join(d, 'input')
    os.makedirs(inp, exist_ok=True)
    shutil.copy(info['pdb'], os.path.join(inp, os.path.basename(info['pdb'])))
    # Copy the domains file and point the config at the COPY, so a slab run is
    # self-contained: pointing fdomains back at the source set's input/ would
    # silently break if that set is moved, renamed, or re-prepared.
    local_domains = os.path.join(inp, 'domains.yaml')
    shutil.copy(info['domains'], local_domains)

    if arm == 'rna':
        fres = merged_residues_csv(info['residues'],
                                   os.path.join(inp, 'residues_CALVADOS3_RNA.csv'))
        ffasta = write_rna_fasta(os.path.join(inp, 'rna.fasta'))
    else:
        fres = os.path.join(inp, 'residues_CALVADOS3.csv')
        shutil.copy(info['residues'], fres)
        ffasta = None

    cfg = Config(
        sysname=info['sysname'],
        box=box,
        temp=TEMP, ionic=IONIC, pH=PH,
        topol='slab',
        slab_eq=True,
        k_eq=0.02,
        steps_eq=STEPS_EQ,
        slab_width=max(20, int(box[0])),
        wfreq=WFREQ,
        steps=steps,
        runtime=0,
        platform=PLATFORM,
        restart='checkpoint',
        frestart='restart.chk',
        verbose=True,
        random_number_seed=12345,
    )
    cfg.write(d, name='config.yaml')

    comp_kwargs = dict(fresidues=fres, pdb_folder=inp, fdomains=local_domains)
    if arm == 'rna':
        comp_kwargs['ffasta'] = ffasta
        comp_kwargs.update(RNA_OPTS)
    comps = Components(**comp_kwargs)
    comps.add(name=info['sysname'], molecule_type='protein', nmol=nchain,
              restraint=True, charge_termini='both')
    if arm == 'rna':
        comps.add(name=RNA_NAME, molecule_type='rna', nmol=n_rna)
    comps.write(d, name='components.yaml')

    with open(os.path.join(d, 'run.py'), 'w') as fh:
        fh.write(RUN_PY)

    beads = nchain * nres + n_rna * RNA_NT * 2     # RNA is 2 beads / nucleotide
    meta = dict(construct=key, arm=arm, sysname=info['sysname'], nres=nres,
                nchain=nchain, n_rna=n_rna, box=box, steps=steps, beads=beads,
                note=note)
    with open(os.path.join(d, 'slab_meta.yaml'), 'w') as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)
    return meta, None


SLURM = '''#!/bin/bash
#SBATCH --job-name=parp14_slab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/slab_%x_%A_%a.log
#
# One PARP14 slab coexistence run per array task.
#
#   sbatch --array=0-5 submit_slab.slurm slab/homotypic
#
# ARM_DIR (positional) holds one subdirectory per construct; the array index
# selects one. Partition/qos/gres mirror bindcraft_job.slurm, which is the only
# verified-working GPU submission in this repo -- adjust if the GPU server uses
# different names.
#
# CALVADOS restarts from restart.chk run `steps` ADDITIONAL steps and append to
# the DCD, so re-submitting a timed-out task continues rather than restarting.

set -u
ARM_DIR=${1:?usage: sbatch --array=0-N submit_slab.slurm <arm_dir>}
CAL_ENV=/home/sbali/miniconda3/envs/CALVADOS

mapfile -t DIRS < <(find "$ARM_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
D=${DIRS[$SLURM_ARRAY_TASK_ID]}
echo "host=$(hostname)  task=$SLURM_ARRAY_TASK_ID  dir=$D"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

cd "$D" || exit 1
srun "$CAL_ENV/bin/python" -u run.py
'''

BENCH_NOTE = '''# GPU benchmark for the PARP14 slab campaign.
#
# pollux has no GPU, so the panel cost in prepare_slab.py --report is an
# ESTIMATE extrapolated from CPU throughput. Run this once on the GPU server to
# replace it with a measurement:
#
#     sbatch --array=0-0 submit_slab.slurm slab/benchmark
#
# It runs the smallest construct for a short fixed number of steps. Divide
# steps*beads by the wall time in the log to get bead-steps/s, then re-run
#     python prepare_slab.py --report --gpu-rate <bead_steps_per_s>
# to cost the real panel.
'''


def cmd_report(reg, args):
    rate = args.gpu_rate or CPU_BEAD_STEPS_PER_S * ASSUMED_GPU_SPEEDUP
    measured = 'MEASURED' if args.gpu_rate else f'ASSUMED ({ASSUMED_GPU_SPEEDUP}x CPU)'
    keys = [k for k, tier, _ in PANEL if args.panel == 'full' or tier == 'core']
    print(f"\nPARP14 slab panel -- throughput {rate:,.0f} bead-steps/s  [{measured}]")
    print(f"steps per run: {args.steps:,}\n")
    print(f"{'construct':18s}{'res':>6}{'chains':>7}  {'box (nm)':<18}{'beads':>9}"
          f"{'homotypic':>11}{'+RNA':>9}")
    print('-' * 82)
    tot = 0.0
    for k in keys:
        info = construct_info(reg, k)
        if not info['ok']:
            print(f"{k:18s}  -- inputs missing, skipped")
            continue
        n = info['nres']
        nc = n_chains_for(n)
        box = box_for(n, nc)
        b_h = nc * n
        b_r = b_h + int(round(nc * RNA_PER_PROTEIN)) * RNA_NT * 2
        d_h = b_h * args.steps / rate / 86400
        d_r = b_r * args.steps / rate / 86400
        tot += d_h + d_r
        bx = f"{box[0]:.0f}x{box[1]:.0f}x{box[2]:.0f}"
        print(f"{k:18s}{n:>6}{nc:>7}  {bx:<18}{b_h:>9,}"
              f"{d_h:>10.1f}d{d_r:>8.1f}d")
    print('-' * 82)
    print(f"{'TOTAL (both arms, sequential)':<57}{tot:>10.1f} GPU-days")
    print(f"\nOn CPU this would be ~{rate / CPU_BEAD_STEPS_PER_S:.0f}x longer "
          f"({tot * rate / CPU_BEAD_STEPS_PER_S:,.0f} days) -- not viable, which is "
          f"why these target the GPU server.")
    if not args.gpu_rate:
        print("\nThe GPU rate above is ASSUMED, not measured. Run "
              "`--benchmark` first, then re-run --report --gpu-rate <value>.")


def main():
    ap = ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--arm', choices=('homotypic', 'rna', 'both'), default=None)
    ap.add_argument('--panel', choices=('core', 'full'), default='core')
    ap.add_argument('--steps', type=int, default=N_STEPS)
    ap.add_argument('--report', action='store_true',
                    help='print the panel + cost table and exit, writing nothing')
    ap.add_argument('--benchmark', action='store_true',
                    help='prepare a short calibration run of the smallest construct')
    ap.add_argument('--gpu-rate', type=float, default=None,
                    help='measured GPU throughput in bead-steps/s, from --benchmark')
    args = ap.parse_args()

    reg = load_registry()

    if args.report or (not args.arm and not args.benchmark):
        cmd_report(reg, args)
        return

    os.makedirs(SLAB_ROOT, exist_ok=True)
    os.makedirs(os.path.join(SLAB_ROOT, 'logs'), exist_ok=True)
    with open(os.path.join(SLAB_ROOT, 'submit_slab.slurm'), 'w') as fh:
        fh.write(SLURM)

    if args.benchmark:
        root = os.path.join(SLAB_ROOT, 'benchmark')
        os.makedirs(root, exist_ok=True)
        with open(os.path.join(root, 'README'), 'w') as fh:
            fh.write(BENCH_NOTE)
        smallest = min((construct_info(reg, k) for k, _, _ in PANEL
                        if construct_info(reg, k)['ok']), key=lambda i: i['nres'])
        meta, err = prepare_one(reg, smallest['key'], 'homotypic',
                                steps=1_000_000, outroot=root,
                                note='GPU throughput calibration')
        print(f"benchmark: {smallest['key']} "
              f"({meta['nchain']} chains, {meta['beads']:,} beads, 1e6 steps)"
              if not err else f"benchmark FAILED: {err}")
        print(f"  -> {os.path.relpath(root, CWD)}")
        print(f"  sbatch --array=0-0 slab/submit_slab.slurm slab/benchmark")
        return

    arms = ('homotypic', 'rna') if args.arm == 'both' else (args.arm,)
    keys = [(k, note) for k, tier, note in PANEL
            if args.panel == 'full' or tier == 'core']
    for arm in arms:
        root = os.path.join(SLAB_ROOT, arm)
        os.makedirs(root, exist_ok=True)
        print(f"\n=== arm: {arm} ===")
        made = 0
        for k, note in keys:
            meta, err = prepare_one(reg, k, arm, args.steps, root, note)
            if err:
                print(f"  {k:18s} SKIPPED -- {err}")
                continue
            made += 1
            print(f"  {k:18s} {meta['nchain']:3d} chains"
                  + (f" + {meta['n_rna']:2d} RNA" if arm == 'rna' else "        ")
                  + f"  box {meta['box']}  {meta['beads']:,} beads")
        print(f"  prepared {made} construct(s) -> {os.path.relpath(root, CWD)}")
        print(f"  sbatch --array=0-{made-1} slab/submit_slab.slurm slab/{arm}")


if __name__ == '__main__':
    main()
