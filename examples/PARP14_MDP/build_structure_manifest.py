#!/usr/bin/env python3
"""One row per simulation run across the whole campaign: what exists, and where.

campaign_inventory.py answers this for the 24 complex runs. Nothing answered it
for the single-chain tree, the fragments or the slab arm, and the counts quoted
in different places disagree four ways -- which is the reason for this file.

Each row is read from the run's own config.yaml, so it cannot drift from what
was simulated. Frame counts come from the DCD header (NSET at byte 8), not by
loading the trajectory, so the whole campaign walks in seconds.

TWO COLUMNS THAT MATTER MORE THAN THE REST

  status        Three states, not two. A directory can hold a trajectory; or
                hold none but carry the dated snapshot PDBs and analysis
                figures of a run whose trajectory was later deleted; or be
                prepared and never started. md/, mka/ and md3art/ are the
                second kind -- 75 replicates with 37 analysis figures each and
                not one .dcd on disk. Their conclusions stand; they cannot be
                re-analysed without re-running. A two-state flag would file
                them with the never-started directories, which is wrong.

  custom_restraints
                The optimised restraint scheme is NOT simply "single-chain" vs
                "slab": md/, core/ and mka/ run with custom_restraints false
                while fl_optimized runs with it true. Flagged per run so the
                index can separate them on the truth rather than on the folder.

Representative-state bookkeeping (n_state_pdb, has_pse) describes THIS disk.
*.pse is gitignored, so those columns do not survive a clone -- the state PDBs
and the CSVs do.

Usage:  python build_structure_manifest.py
        -> analysis/structure_manifest.csv
"""
import hashlib, struct, sys
from datetime import datetime
from pathlib import Path
import pandas as pd, yaml

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'; OUT.mkdir(exist_ok=True)

# (label, glob of run dirs, how to name the construct from the path)
ROOTS = [
    ('single-chain', 'fl/seed-*_sample-*',            lambda p: p.parts[-2]),
    ('single-chain', 'fl_optimized/seed-*_sample-*',  lambda p: p.parts[-2]),
    ('single-chain', 'md/seed-*_sample-*',            lambda p: p.parts[-2]),
    ('single-chain', 'mka/seed-*_sample-*',           lambda p: p.parts[-2]),
    ('single-chain', 'core/seed-*_sample-*',          lambda p: p.parts[-2]),
    ('single-chain', 'norrm/seed-*_sample-*',         lambda p: p.parts[-2]),
    ('single-chain', 'noart/seed-*_sample-*',         lambda p: p.parts[-2]),
    ('single-chain', 'md3art/seed-*_sample-*',        lambda p: p.parts[-2]),
    ('fragment',     'fragments/*/seed-*_sample-*',   lambda p: p.parts[-2]),
    ('slab',         'slab/*/*',                      lambda p: p.name),
    ('complex',      'complex/binding/*/rep-*',       lambda p: p.parts[-2]),
    ('complex',      'complex/binding_go/*/rep-*',    lambda p: p.parts[-2]),
    ('complex',      'complex/binding_gotight/*/rep-*', lambda p: p.parts[-2]),
    ('complex',      'complex/binding_xl/*/rep-*',    lambda p: p.parts[-2]),
]


def dcd_frames(path):
    """NSET from the DCD header, without reading the coordinates."""
    try:
        with open(path, 'rb') as f:
            head = f.read(12)
            if len(head) < 12 or head[4:8] != b'CORD':
                return None
            return int(struct.unpack('<i', head[8:12])[0])
    except OSError:
        return None


def topology_for(d):
    for cand in ('top.pdb', 'checkpoint.pdb'):
        if (d / cand).is_file():
            return d / cand
    pdbs = sorted((d / 'input').glob('*.pdb')) if (d / 'input').is_dir() else []
    return pdbs[0] if pdbs else None


def scan():
    rows = []
    for arm, pattern, construct_of in ROOTS:
        for d in sorted(HERE.glob(pattern)):
            cfg_p = d / 'config.yaml'
            if not (d.is_dir() and cfg_p.is_file()):
                continue
            try:
                cfg = yaml.safe_load(open(cfg_p)) or {}
            except Exception:
                cfg = {}
            sysname = cfg.get('sysname') or ''
            dcds = sorted(p for p in d.glob('*.dcd')
                          if not p.name.startswith('equilibration'))
            traj = next((p for p in dcds if p.stem == sysname), dcds[0] if dcds else None)
            comp_p = d / 'components.yaml'
            comps = ''
            if comp_p.is_file():
                try:
                    comps = '+'.join((yaml.safe_load(open(comp_p)) or {})
                                     .get('system', {}).keys())
                except Exception:
                    pass
            top = topology_for(d)
            # dated snapshot PDBs are written during a run, so they are
            # evidence it executed even when the trajectory is gone
            snaps = len(list(d.glob('*_2026*.pdb')))
            status = ('has_traj' if traj is not None else
                      'traj_removed' if (snaps or (d / 'checkpoint.pdb').is_file())
                      else 'prepared_only')
            rows.append(dict(
                arm=arm, construct=construct_of(d), run=str(d.relative_to(HERE)),
                engine='CALVADOS', components=comps,
                box_nm=(cfg.get('box') or [None])[0],
                steps=cfg.get('steps'), wfreq=cfg.get('wfreq'),
                temp=cfg.get('temp'), topol=cfg.get('topol'),
                custom_restraints=cfg.get('custom_restraints'),
                custom_restraint_type=cfg.get('custom_restraint_type'),
                has_traj=traj is not None, status=status,
                n_snapshot_pdb=snaps,
                n_frames=dcd_frames(traj) if traj else None,
                traj=str(traj.relative_to(HERE)) if traj else '',
                topology=str(top.relative_to(HERE)) if top else '',
                has_checkpoint=(d / 'checkpoint.pdb').is_file(),
                config_md5=hashlib.md5(cfg_p.read_bytes()).hexdigest()[:10],
                modified=datetime.fromtimestamp(
                    (traj or cfg_p).stat().st_mtime).strftime('%Y-%m-%d')))
    return pd.DataFrame(rows)


def hyres_rows():
    """HyRes runs carry no config.yaml, so they are walked explicitly."""
    base = HERE / 'complex' / 'hyres' / 'runs'
    out = []
    for d in sorted(list(base.glob('matched/*_rep*')) + [base / 'prod_rep1']):
        if not d.is_dir():
            continue
        traj = d / 'system.dcd'
        top = next((d / n for n in ('complex_start.pdb', 'start.pdb')
                    if (d / n).is_file()), None)
        out.append(dict(
            arm='complex', construct='p9_dtx3l (HyRes)',
            run=str(d.relative_to(HERE)), engine='HyRes',
            components='parp9+dtx3l', box_nm=40.0, steps=None, wfreq=None,
            temp=293.0, topol='', custom_restraints=False,
            custom_restraint_type='harmonic (within-block only)',
            has_traj=traj.is_file(),
            status='has_traj' if traj.is_file() else 'prepared_only',
            n_snapshot_pdb=0,
            n_frames=dcd_frames(traj) if traj.is_file() else None,
            traj=str(traj.relative_to(HERE)) if traj.is_file() else '',
            topology=str(top.relative_to(HERE)) if top else '',
            has_checkpoint=(d / 'system.chk').is_file(),
            config_md5='', modified=datetime.fromtimestamp(
                (traj if traj.is_file() else d).stat().st_mtime).strftime('%Y-%m-%d')))
    return pd.DataFrame(out)


def add_structure_bookkeeping(df):
    """Per CONSTRUCT: do representative states exist, and from what."""
    cx_states = HERE / 'complex' / 'analysis' / 'states'
    cx_pse = HERE / 'complex' / 'analysis' / 'pymol'
    rep_csv = HERE / 'complex' / 'analysis' / 'representative_states.csv'
    known = {}
    if rep_csv.is_file():
        r = pd.read_csv(rep_csv)
        for s, g in r.groupby('set'):
            known[s] = dict(n_state_pdb=len(g), state_method=g.reduce.iloc[0],
                            state_source='representative_states.csv')
    n_pdb, meth, src, pse = [], [], [], []
    for _, row in df.iterrows():
        k = known.get(row.construct, {})
        n_pdb.append(k.get('n_state_pdb', 0))
        meth.append(k.get('state_method', ''))
        src.append(k.get('state_source', ''))
        pse.append(bool(list(cx_pse.glob(f'{row.construct}*.pse'))) if cx_pse.is_dir() else False)
    df['n_state_pdb'] = n_pdb
    df['state_method'] = meth
    df['state_source'] = src
    df['has_pse'] = pse
    return df


def main():
    df = add_structure_bookkeeping(
        pd.concat([scan(), hyres_rows()], ignore_index=True))
    df.to_csv(OUT / 'structure_manifest.csv', index=False)
    print(f'{len(df)} run directories -> {OUT / "structure_manifest.csv"}\n')
    g = pd.crosstab(df.arm, df.status)
    g['runs'] = g.sum(axis=1)
    print(g.to_string())
    print(f'\nTOTAL {len(df)} dirs across {df.construct.nunique()} constructs')
    for k, v in df.status.value_counts().items():
        print(f'   {k:14s} {v}')
    rs = df[df.n_state_pdb > 0].construct.nunique()
    print(f'constructs with representative states: {rs} of {df.construct.nunique()}')


if __name__ == '__main__':
    main()
