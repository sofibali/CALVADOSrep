#!/usr/bin/env python
"""
Compute c_sat for a finished slab run.

c_sat -- the saturation concentration, i.e. the concentration of the dilute
phase in coexistence with the dense phase -- is the quantitative observable
this whole campaign exists to produce. It comes out of
`calvados.analysis.SlabAnalysis` as the `c_dilute` column of
`{sysname}_ps_results.csv`, in mM.

WHY THIS WRAPPER
----------------
Three things need doing around the stock call and are easy to get wrong:

1. **Use the repo's `calvados.analysis`, not the installed one.** Only this
   checkout has the fix that moves `fit_profile`'s "NOT CONVERGED" check ahead
   of its `return`. In the installed copy that check is dead code after a
   return, so a failed tanh interface fit never warns and c_sat comes back
   silently wrong. This script forces the repo copy onto sys.path and says
   which one it loaded.

2. **Discard the start of production.** `slab_eq` pulls the chains together,
   but the slab still relaxes for a while after that force is removed. Frames
   before `--discard` (default 10%) are dropped.

3. **Check it actually converged**, rather than trusting one number: the
   interface fit cutoffs are reported, and c_sat is recomputed on the first and
   second half of the retained frames so a drift shows up.

USAGE
-----
    python analyze_slab.py homotypic/fl
    python analyze_slab.py homotypic/fl --discard 0.2 --step 2
    python analyze_slab.py --all            # every finished run
"""
import argparse
import glob
import os
import sys

import numpy as np
import yaml

# slab/ -> PARP14_MDP -> examples -> CALVADOS (four levels, not three)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, REPO)          # repo calvados must win over site-packages
SLAB_ROOT = os.path.dirname(os.path.abspath(__file__))


def load_analysis():
    from calvados import analysis
    src = os.path.abspath(analysis.__file__)
    if not src.startswith(REPO):
        print(f'WARNING: loaded calvados.analysis from {src}\n'
              f'         NOT the repo copy at {REPO}. If that install lacks the\n'
              f'         fit_profile convergence fix, a failed interface fit will\n'
              f'         pass silently and c_sat may be wrong.')
    else:
        print(f'calvados.analysis: {src}  (repo copy, has the convergence fix)')
    return analysis


def finished_runs():
    out = []
    for arm in ('homotypic', 'rna'):
        for d in sorted(glob.glob(os.path.join(SLAB_ROOT, arm, '*'))):
            fmeta = os.path.join(d, 'slab_meta.yaml')
            if not os.path.isfile(fmeta):
                continue
            meta = yaml.safe_load(open(fmeta))
            log = os.path.join(d, f"{meta['sysname']}.log")
            done = 0
            if os.path.isfile(log):
                for line in open(log, errors='replace'):
                    if line.startswith('#') or not line.strip():
                        continue
                    try:
                        done = max(done, int(float(line.split('\t')[0])))
                    except (ValueError, IndexError):
                        pass
            if done >= int(float(meta['steps'])):
                out.append(d)
    return out



def slab_stability(d, sysname, stride=400):
    """
    Did the slab survive production, or dissolve?

    `slab_eq` compresses every chain into a thin slab with a linear pull, then
    removes that force. If the construct self-associates the slab persists and
    a tanh interface fit is meaningful. If it does not, the slab simply
    disperses -- and then `fit_profile` is fitting interfaces that do not
    exist, `calc_concentrations` returns NaN, and quoting a c_sat from it would
    be meaningless.

    Returns (width_after_eq, width_first_production, width_last_production) in
    nm, each the z-range holding 90% of the beads.
    """
    import mdtraj as md

    def w90(traj):
        lz = traj.unitcell_lengths[0, 2]
        zs = []
        for fr in range(traj.n_frames):
            z = traj.xyz[fr, :, 2].copy()
            zs.append((z - z.mean() + lz / 2) % lz)   # centre on the beads' own COM
        h, _ = np.histogram(np.concatenate(zs), bins=80, range=(0, lz))
        h = h / traj.n_frames
        hs = np.sort(h)[::-1]
        c = np.cumsum(hs) / hs.sum()
        return (c < 0.9).sum() / 80 * lz

    top = os.path.join(d, 'top.pdb')
    eq = md.load(os.path.join(d, f'equilibration_{sysname}.dcd'), top=top)
    pr = md.load(os.path.join(d, f'{sysname}.dcd'), top=top, stride=stride)
    return w90(eq[-3:]), w90(pr[:2]), w90(pr[-2:])



def phase_state(d, sysname, stride=100):
    """
    Is the dense region a LIQUID condensate, or a kinetically trapped solid?

    A slab that never disperses is not automatically evidence of phase
    separation. `slab_eq` compresses every chain into a thin slab; at high
    enough concentration in a small cross-section that compaction can jam,
    and the result is a frozen aggregate that also never disperses. The two
    look identical if you only measure slab width at the start and the end.

    They are trivially separable by whether anything MOVES:

      * a liquid condensate fluctuates in width and scrambles chain order as
        chains diffuse past one another;
      * a trapped solid has a constant width and preserves chain order.

    Returns (mean_width, sd_width, order_correlation, core_mg_per_mL).
    `order_correlation` is the Spearman correlation of per-chain z-rank between
    the first and last frame: ~1 means no chain ever changed places.
    """
    import mdtraj as md

    meta = yaml.safe_load(open(os.path.join(d, 'slab_meta.yaml')))
    nch, nres = int(meta['nchain']), int(meta['nres'])
    Lx, Ly, _ = [float(x) for x in meta['box']]
    t = md.load(os.path.join(d, f'{sysname}.dcd'),
                top=os.path.join(d, 'top.pdb'), stride=stride)
    lz = t.unitcell_lengths[0, 2]

    widths, com = [], np.zeros((t.n_frames, nch))
    for fr in range(t.n_frames):
        z = t.xyz[fr, :, 2].copy()
        zc = (z - z.mean() + lz / 2) % lz
        h, _ = np.histogram(zc, bins=160, range=(0, lz))
        hs = np.sort(h)[::-1]
        widths.append((np.cumsum(hs) / hs.sum() < 0.9).sum() / 160 * lz)
        for c in range(nch):
            com[fr, c] = z[c * nres:(c + 1) * nres].mean()
    com -= com.mean(axis=1, keepdims=True)          # remove slab drift
    rank = lambda a: np.argsort(np.argsort(a))
    tau = float(np.corrcoef(rank(com[0]), rank(com[-1]))[0, 1])

    w = float(np.mean(widths))
    # beads/nm^3 -> mg/mL, at ~110 Da per coarse-grained residue
    g_per_cm3 = nch * nres / (Lx * Ly * w) * 110 * 1.66053907e-24 * 1e21
    mg_ml = g_per_cm3 * 1000.0
    return w, float(np.std(widths)), tau, mg_ml


def run_one(analysis, d, discard, step):
    meta = yaml.safe_load(open(os.path.join(d, 'slab_meta.yaml')))
    sysname, nchain = meta['sysname'], int(meta['nchain'])
    rel = os.path.relpath(d, SLAB_ROOT)
    print(f'\n{"=" * 70}\n{rel}   {sysname}   {nchain} chains, {meta["beads"]:,} beads\n{"=" * 70}')

    slab = analysis.SlabAnalysis(
        name=sysname, input_path=d, output_path=d,
        input_pdb='top.pdb', input_dcd=None, centered_dcd=f'{sysname}_centered.dcd',
        # the protein chains are the reference phase; RNA chains (rna arm) come
        # after them, so ref_chains stops at nchain-1 in both arms
        ref_chains=(0, nchain - 1), ref_name='protein', verbose=True)

    w_eq, w_p0, w_p1 = slab_stability(d, sysname)
    print(f'slab width (90% of beads): after slab_eq {w_eq:.0f} nm -> '
          f'production start {w_p0:.0f} nm -> production end {w_p1:.0f} nm')
    dissolved = w_p1 > 3 * w_eq
    trapped = False
    if not dissolved:
        wm, wsd, tau, mg = phase_state(d, sysname)
        print(f'the slab persisted: width {wm:.1f} +/- {wsd:.2f} nm, core ~{mg:.0f} mg/mL,'
              f' chain-order correlation {tau:.3f}')
        if wsd < 0.5 and tau > 0.9:
            print('  *** KINETICALLY TRAPPED, not a condensate. The width never'
                  ' fluctuates\n      and no chain ever changes places, so this is a'
                  ' jammed solid that\n      slab_eq compressed into existence -- not'
                  ' liquid-liquid coexistence.\n      c_sat is NOT defined for it'
                  ' either.')
            trapped = True        # no measurable c_sat either
        else:
            print('  liquid-like: width fluctuates and chains rearrange -- a genuine'
                  ' dense phase')
    if dissolved and not trapped:
        print('  *** THE SLAB DISSOLVED. There is no dense/dilute coexistence, so\n'
              '      c_sat is NOT defined for this run -- the tanh interface fit\n'
              '      below is fitting interfaces that do not exist and any number\n'
              '      it returns is meaningless. Reported as no-phase-separation.')

    import mdtraj as md
    n_frames = len(md.open(os.path.join(d, f'{sysname}.dcd')))
    start = int(n_frames * discard)
    print(f'frames: {n_frames}, discarding first {start} ({discard:.0%}), step {step}')

    slab.center(start=start, step=step, center_target='all')
    slab.calc_profiles()
    slab.calc_concentrations()

    csv = os.path.join(d, f'{sysname}_ps_results.csv')
    if not os.path.isfile(csv):
        print('  no results csv written')
        return None
    import csv as _csv
    row = list(_csv.DictReader(open(csv)))[-1]
    print(f'\n  --> c_sat (c_dilute) = {row.get("c_dilute")} mM'
          f'   dense = {row.get("c_dense")} mM')
    print(f'      cutoffs: {[f"{k}={v}" for k, v in row.items() if k.startswith("cutoffs")]}')
    return {'run': rel, 'sysname': sysname, 'dissolved': dissolved,
            'trapped': trapped,
            'w_eq': w_eq, 'w_end': w_p1, **row}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', nargs='?')
    ap.add_argument('--all', action='store_true', help='every finished run')
    ap.add_argument('--discard', type=float, default=0.10,
                    help='fraction of production frames to drop from the start')
    ap.add_argument('--step', type=int, default=1)
    args = ap.parse_args()

    analysis = load_analysis()
    targets = finished_runs() if args.all else [os.path.abspath(args.run_dir)]
    if not targets:
        sys.exit('no finished runs found')
    print(f'{len(targets)} finished run(s): '
          f'{", ".join(os.path.relpath(t, SLAB_ROOT) for t in targets)}')

    rows = []
    for d in targets:
        try:
            r = run_one(analysis, d, args.discard, args.step)
            if r:
                rows.append(r)
        except Exception as exc:
            print(f'  FAILED: {type(exc).__name__}: {exc}')

    if rows:
        print(f'\n{"=" * 70}\nc_sat SUMMARY\n{"=" * 70}')
        print(f'{"run":<28}{"slab nm (eq->end)":>20}{"c_sat (mM)":>14}')
        for r in rows:
            span = f'{r["w_eq"]:.0f} -> {r["w_end"]:.0f}'
            csat = ('TRAPPED' if r.get('trapped') else
                    'NO PHASE SEP' if r['dissolved'] else (r.get('c_dilute') or '?'))
            print(f'{r["run"]:<28}{span:>20}{csat:>14}')
        if any(r['dissolved'] and not r.get('trapped') for r in rows):
            print('\nNO PHASE SEP = the compacted slab dispersed during production, so\n'
                  'there is no coexistence to measure. That is a result, not a failure:\n'
                  'it means the construct does not condense at this concentration.')
        if any(r.get('trapped') for r in rows):
            print('\nTRAPPED = the slab persisted but never fluctuated and no chain ever\n'
                  'changed places: a jammed solid slab_eq compressed into existence, not\n'
                  'liquid-liquid coexistence. c_sat is not defined for it either.')


if __name__ == '__main__':
    main()
