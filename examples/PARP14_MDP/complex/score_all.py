#!/usr/bin/env python3
"""Crosslink agreement for every run that contains both PARP9 and DTX3L.

WHICH CROSSLINKS ARE BEING COMPARED
-----------------------------------
The 51 inter-protein BS3 crosslinks of Ashok et al. 2022 (Biochem J 479:289,
doi 10.1042/BCJ20210722), deduplicated across constructs to 36 unique PARP9-
DTX3L residue pairs. Three nested subsets are scored, because they carry very
different evidential weight:

  ALL (36 unique)   every inter-protein crosslink reported, pooled over the
                    three constructs (FL PARP9 + D3 / + D3RD / + FL DTX3L)
  REPRODUCIBLE (21) seen in >=2 of the 4 BS3 replicates
  ALL-4 (8 rows)    seen in every replicate -- the most trustworthy restraints

Sets without both chains (the homodimers, p14_dtx3l, p14_p9) cannot be scored:
there is no published PARP9-DTX3L crosslink to test in them. The ternary set
can, since it contains both.
"""
import subprocess, sys, re
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent
PY = '/home/sbali/miniconda3/envs/CALVADOS/bin/python'

JOBS = [
    ('CALVADOS unrestrained (separated)', ['--root','binding','--set','p9_dtx3l']),
    ('CALVADOS unrestrained (docked)',    ['--root','binding','--set','p9_dtx3l_docked']),
    ('CALVADOS ternary (separated)',      ['--root','binding','--set','ternary']),
    ('CALVADOS ternary (docked)',         ['--root','binding','--set','ternary_docked']),
    ('CALVADOS +3 XL harmonic k=5',       ['--root','binding_go','--set','p9_dtx3l_xl_h5']),
    ('CALVADOS +3 XL harmonic k=20',      ['--root','binding_go','--set','p9_dtx3l_xl_h20']),
    ('CALVADOS +3 XL harmonic k=100',     ['--root','binding_go','--set','p9_dtx3l_xl_h100']),
    ('CALVADOS +3 XL Go k=15',            ['--root','binding_go','--set','p9_dtx3l_xl_go15']),
    ('CALVADOS +3 XL Go tight r=1.5 k=50',  ['--root','binding_gotight','--set','p9_dtx3l_got_r15_k50']),
    ('CALVADOS +3 XL Go tight r=1.5 k=100', ['--root','binding_gotight','--set','p9_dtx3l_got_r15_k100']),
    ('CALVADOS +3 XL Go tight r=1.2 k=50',  ['--root','binding_gotight','--set','p9_dtx3l_got_r12_k50']),
    ('CALVADOS +3 XL Go tight r=1.2 k=100', ['--root','binding_gotight','--set','p9_dtx3l_got_r12_k100']),
    ('CALVADOS +21 XL IMProv-tiered',      ['--root','binding_xl','--set','p9_dtx3l_xl30_k20']),
    ('HyRes unrestrained (separated)',    ['--set','p9_dtx3l','--hyres','hyres/runs/prod_rep1']),
    ('HyRes matched (separated)',         ['--set','p9_dtx3l','--hyres'] +
                                          [f'hyres/runs/matched/separated_rep{i}' for i in (1,2,3)]),
    ('HyRes matched (docked)',            ['--set','p9_dtx3l','--hyres'] +
                                          [f'hyres/runs/matched/docked_rep{i}' for i in (1,2,3)]),
]

rows = []
for label, args in JOBS:
    r = subprocess.run([PY, str(HERE/'rescore_crosslinks.py')] + args + ['--stride','20'],
                       capture_output=True, text=True, cwd=HERE)
    csv = HERE/'analysis'/'crosslink_rescore.csv'
    if r.returncode != 0 or not csv.exists():
        print(f'  SKIP {label}: {(r.stderr or "").strip().splitlines()[-1][:80] if r.stderr else "no output"}')
        continue
    # read the per-crosslink table rather than parsing stdout -- the summary
    # block's ragged "all 4 replicates" label breaks whitespace splitting
    d = pd.read_csv(csv)
    d = d[d.model == d.model.iloc[-1]]
    for key, sub in (('ALL (36)', d),
                     ('reproducible (21)', d[d.repro]),
                     ('all-4-replicate (8)', d[d.n_samples == 4])):
        if not len(sub):
            continue
        rows.append(dict(model=label, subset=key, n=len(sub),
                         per_frame_30A=round(sub.pct_30A.mean(), 2),
                         ever_30A=round(100*(sub.pct_30A > 0).mean(), 1),
                         pct1_30A=round(100*(sub.pct_30A >= 1).mean(), 1)))
    print(f'  scored {label}', flush=True)

df = pd.DataFrame(rows)
df.to_csv(HERE/'analysis'/'crosslink_agreement_all.csv', index=False)
piv = df.pivot_table(index='model', columns='subset', values='per_frame_30A')
order = [j[0] for j in JOBS]
piv = piv.reindex([o for o in order if o in piv.index])
print('\n=== per-frame % within 30 A ===')
print(piv.round(2).to_string())
print(f'\nsaved {HERE/"analysis"/"crosslink_agreement_all.csv"}')
