#!/usr/bin/env python3
"""One row per simulation run: what it is, and every score it carries.

Joins the three result tables that are each keyed differently:
    analysis/campaign_inventory.csv      run   = 'binding/p9_dtx3l'
    analysis/crosslink_agreement_all.csv model = 'CALVADOS unrestrained (separated)'
    analysis/domain_contact_maps.csv     model = 'binding/p9_dtx3l'

Not every run is crosslink-scorable. The published BS3 set has inter-protein
crosslinks for PARP9-DTX3L only, so a homodimer or a PARP14 pair has nothing to
be scored against -- that is a property of the experiment, not a gap in the
analysis, and those rows say so rather than showing a blank.

Writes analysis/campaign_scorecard.csv and prints the markdown table.
"""
from pathlib import Path
import re
import pandas as pd, numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'

# inventory run -> label in crosslink_agreement_all.csv
XL_LABEL = {
    'binding/p9_dtx3l':                     'CALVADOS unrestrained (separated)',
    'binding/p9_dtx3l_docked':              'CALVADOS unrestrained (docked)',
    'binding/ternary':                      'CALVADOS ternary (separated)',
    'binding/ternary_docked':               'CALVADOS ternary (docked)',
    'binding_go/p9_dtx3l_xl_h5':            'CALVADOS +3 XL harmonic k=5',
    'binding_go/p9_dtx3l_xl_h20':           'CALVADOS +3 XL harmonic k=20',
    'binding_go/p9_dtx3l_xl_h100':          'CALVADOS +3 XL harmonic k=100',
    'binding_go/p9_dtx3l_xl_go15':          'CALVADOS +3 XL Go k=15',
    'binding_gotight/p9_dtx3l_got_r15_k50': 'CALVADOS +3 XL Go tight r=1.5 k=50',
    'binding_gotight/p9_dtx3l_got_r15_k100':'CALVADOS +3 XL Go tight r=1.5 k=100',
    'binding_gotight/p9_dtx3l_got_r12_k50': 'CALVADOS +3 XL Go tight r=1.2 k=50',
    'binding_gotight/p9_dtx3l_got_r12_k100':'CALVADOS +3 XL Go tight r=1.2 k=100',
    'binding_xl/p9_dtx3l_xl30_k20':         'CALVADOS +21 XL IMProv-tiered',
    'hyres/prod_rep1 (unrestrained)':       'HyRes unrestrained (separated)',
    'hyres/matched separated':              'HyRes matched (separated)',
    'hyres/matched docked':                 'HyRes matched (docked)',
}
# inventory run -> label in domain_contact_maps.csv (identical except for HyRes)
DC_LABEL = {'hyres/prod_rep1 (unrestrained)': 'HyRes unrestrained',
            'hyres/matched separated':        'HyRes matched separated',
            'hyres/matched docked':           'HyRes matched docked'}


def inter_restraints(v):
    """'88 custom (3 inter-chain), harmonic' -> '3 crosslinks, harmonic'.

    The 'custom' count includes the intra-chain restraints that hold split
    domains together (PARP9's split KH1 and the like); only the inter-chain
    ones distinguish the runs, and a set with 0 of them is unrestrained between
    chains however many custom restraints it carries.
    """
    v = str(v)
    m = re.search(r'\((\d+) inter-chain\)(?:,\s*(\w+))?', v)
    if not m:
        return 'none' if v in ('none', 'nan') else v
    n, kind = int(m.group(1)), (m.group(2) or '')
    if n == 0:
        return 'none'
    return f'{n} crosslink{"s" if n > 1 else ""}' + (f', {kind}' if kind else '')


def main():
    inv = pd.read_csv(OUT / 'campaign_inventory.csv')
    agr = pd.read_csv(OUT / 'crosslink_agreement_all.csv')
    dc = pd.read_csv(OUT / 'domain_contact_maps.csv')

    rows = []
    for _, r in inv.iterrows():
        run = r['run']
        d = dict(engine=r.engine, run=run, components=r.components,
                 beads=r.beads, n_reps=r.n_reps, total_us=r.total_us,
                 start={'grid': 'separated',
                        'docked (AF3)': 'docked'}.get(r.start, r.start),
                 inter_restraints=inter_restraints(r.inter_restraints))

        # --- crosslink score -------------------------------------------------
        lab = XL_LABEL.get(run)
        has_both = ('parp9' in r.components) and ('dtx3l' in r.components)
        if lab and lab in set(agr.model):
            g = agr[agr.model == lab].set_index('subset')
            d.update(xl_scored='yes',
                     xl_all36=g.loc['ALL (36)', 'per_frame_30A'],
                     xl_repro21=g.loc['reproducible (21)', 'per_frame_30A'],
                     xl_all4=g.loc['all-4-replicate (8)', 'per_frame_30A'],
                     xl_ever=g.loc['reproducible (21)', 'ever_30A'],
                     xl_pct1=g.loc['reproducible (21)', 'pct1_30A'])
        else:
            d.update(xl_scored='not scorable - no PARP9+DTX3L pair'
                     if not has_both else 'not scored',
                     xl_all36=np.nan, xl_repro21=np.nan, xl_all4=np.nan,
                     xl_ever=np.nan, xl_pct1=np.nan)

        # --- domain contact map ----------------------------------------------
        dlab = DC_LABEL.get(run, run)
        g = dc[(dc.model == dlab) & (dc.kind == 'inter')]
        if len(g):
            # the PARP9-DTX3L block where it exists, else the largest chain pair
            sub = g[(g.chain_i == 'parp9') & (g.chain_j == 'dtx3l')]
            if not len(sub):
                key = g.groupby(['chain_i', 'chain_j']).size().idxmax()
                sub = g[(g.chain_i == key[0]) & (g.chain_j == key[1])]
            top = sub.nlargest(1, 'pct_contact_bound').iloc[0]
            tot = sub.pct_contact_bound.sum()
            # Interface BREADTH, not the top pair's share of contact. Share is
            # ambiguous: HyRes docked scores 7% not because it is diffuse but
            # because ~14 pairs are near-saturated, so the top one is a small
            # slice of a large total. Counting pairs above a threshold separates
            # "broad interface" from "diffuse nothing" the way share cannot.
            n50 = int((sub.pct_contact_bound > 50).sum())
            n20 = int((sub.pct_contact_bound > 20).sum())
            d.update(dmap='yes', n_pairs_over50=n50, n_pairs_over20=n20,
                     interface=('never associated'
                                if not sub.n_bound_frames.iloc[0] else
                                'none' if top.pct_contact_bound < 20 else
                                'single-pair lock' if n20 <= 1 else
                                'broad' if n50 >= 6 else 'narrow'),
                     assoc_pct=round(100 * sub.n_bound_frames.iloc[0]
                                     / sub.n_frames.iloc[0], 1),
                     top_pair=f'{top.chain_i} {top.domain_i} ↔ {top.chain_j} {top.domain_j}',
                     top_pair_pct=top.pct_contact_bound,
                     top_pair_nm=top.mean_min_bound_nm,
                     top_pair_share=round(100 * top.pct_contact_bound / tot, 1)
                                    if tot else np.nan)
        else:
            d.update(dmap='no', n_pairs_over50=np.nan, n_pairs_over20=np.nan,
                     interface='', assoc_pct=np.nan, top_pair='',
                     top_pair_pct=np.nan, top_pair_nm=np.nan,
                     top_pair_share=np.nan)
        rows.append(d)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'campaign_scorecard.csv', index=False)
    print(f'{len(df)} runs, {df.n_reps.sum()} replicates, '
          f'{df.total_us.sum():.1f} us')
    print(f"crosslink-scored: {(df.xl_scored == 'yes').sum()}  "
          f"domain-mapped: {(df.dmap == 'yes').sum()}")
    print('\nwrote', OUT / 'campaign_scorecard.csv')
    return df


if __name__ == '__main__':
    main()
