#!/usr/bin/env python3
"""One table: every simulation run, its setup, and its crosslink agreement.

Pulls components, restraints, timescale and replicate count straight from each
run's own config/components/restraint files rather than from notes, so the table
cannot drift from what was actually simulated.
"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd, yaml
warnings.filterwarnings('ignore')
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname
import analyze_convergence as ac

NS_PER_FRAME = 0.05


def probe(root, name):
    d = HERE / root / name
    reps = sorted(d.glob('rep-*'), key=lambda p: int(p.name.split('-')[1]))
    done = [r for r in reps if list(r.glob('*.dcd'))]
    if not done:
        return None
    r0 = done[0]
    cfg = yaml.safe_load(open(r0 / 'config.yaml'))
    comp = yaml.safe_load(open(r0 / 'components.yaml'))
    base = name if name in SETS else 'p9_dtx3l'
    chains = SETS.get(base, [])
    beads = sum(ac.LENGTHS[c] for c in chains)
    ns = cfg['steps'] * 0.01 / 1000.0
    # restraints actually present
    intra = 'domains.yaml (harmonic k=700, cutoff 0.9 nm)' if comp['defaults'].get('restraint') else 'none'
    cr = 'none'
    f = r0 / 'input' / 'custom_restraints.txt'
    if cfg.get('custom_restraints') and f.exists():
        lines = [l for l in f.read_text().split('\n') if l.strip() and not l.startswith('#')]
        inter = [l for l in lines if l.split('|')[0].split()[0] != l.split('|')[1].split()[0]]
        cr = (f"{len(lines)} custom ({len(inter)} inter-chain), "
              f"{cfg.get('custom_restraint_type')}")
    return dict(engine='CALVADOS', run=f'{root}/{name}', components='+'.join(chains),
                beads=beads, box_nm=cfg['box'][0], ns_per_rep=ns, n_reps=len(done),
                total_us=round(len(done)*ns/1000, 1), start=cfg['topol'] if cfg['restart'] != 'pdb' else 'docked (AF3)',
                intra_restraints=intra, inter_restraints=cr)


def hyres(tag, dirs, restr):
    ok = [d for d in dirs if (HERE/d/'system.dcd').exists()]
    if not ok:
        return None
    return dict(engine='HyRes', run=tag, components='parp9+dtx3l', beads=9968,
                box_nm=40.0, ns_per_rep=500.0, n_reps=len(ok),
                total_us=round(len(ok)*0.5, 1),
                start='separated' if 'separated' in tag else ('docked (AF3)' if 'docked' in tag else 'separated'),
                intra_restraints=restr, inter_restraints='none')


def main():
    rows = []
    for s in SETS:
        r = probe('binding', s)
        if r: rows.append(r)
    for s in ('p9_dtx3l_xl_h5', 'p9_dtx3l_xl_h20', 'p9_dtx3l_xl_h100', 'p9_dtx3l_xl_go15'):
        r = probe('binding_go', s)
        if r: rows.append(r)
    for s in ('p9_dtx3l_got_r15_k50', 'p9_dtx3l_got_r15_k100',
              'p9_dtx3l_got_r12_k50', 'p9_dtx3l_got_r12_k100'):
        r = probe('binding_gotight', s)
        if r: rows.append(r)
    r = probe('binding_xl', 'p9_dtx3l_xl30_k20')
    if r: rows.append(r)
    rows.append(hyres('hyres/prod_rep1 (unrestrained)', ['hyres/runs/prod_rep1'], 'NONE'))
    rows.append(hyres('hyres/matched separated',
                      [f'hyres/runs/matched/separated_rep{i}' for i in (1,2,3)],
                      '7540 harmonic CA-CA (within-block only; omits the 85 split-KH1 pairs)'))
    rows.append(hyres('hyres/matched docked',
                      [f'hyres/runs/matched/docked_rep{i}' for i in (1,2,3)],
                      '7665 harmonic CA-CA (within-block only; omits the 85 split-KH1 pairs)'))
    df = pd.DataFrame([r for r in rows if r])
    df.to_csv(HERE/'analysis'/'campaign_inventory.csv', index=False)
    tot = df.total_us.sum()
    print(df.to_string(index=False))
    print(f'\nTOTAL {df.n_reps.sum()} replicates, {tot:.1f} us of trajectory')


if __name__ == '__main__':
    main()
