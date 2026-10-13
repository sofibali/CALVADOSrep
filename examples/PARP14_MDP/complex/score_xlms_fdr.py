#!/usr/bin/env python3
"""Target-decoy quality control on the 2026-08-18 DMTMM crosslinking screen.

The six `*_pairs_acid_K.xls` files under
`docs/20260818_SB31_P14-7_XLMS/` are tab-separated search output (no header),
three replicates each of two samples: SB31 and P14-7 (+ADPr). Each row is one
crosslink PSM. Columns 3 and 4 name the two proteins the two peptides came
from, and a `reverse_` prefix marks a decoy sequence, so every PSM falls into
one of three classes:

    TT  both peptides target   -- candidate real crosslink
    TD  one target, one decoy
    DD  both decoy

Under Walzthoeni et al. 2012 (Nat Methods 9:901) the false fraction of TT is
estimated as

    FDR = (TD - DD) / TT

and a search in which identifications are purely random produces the classes
in a 1:2:1 ratio. The distance from that ratio is the other half of the
readout: a screen whose TT:TD:DD approaches 1:2:1 has no signal to rank.

Column 27 is the search score, sorted descending within each file; this script
re-sorts rather than trusting that. The decoy content of the top-scoring PSMs
says whether the score separates anything, and the replicate-overlap
enrichment says whether reproducibility does.

Usage:  python score_xlms_fdr.py
        -> analysis/xlms_screen_fdr.csv
"""
from pathlib import Path
from collections import Counter
import pandas as pd

HERE = Path(__file__).resolve().parent
XLDIR = HERE.parent / 'docs' / '20260818_SB31_P14-7_XLMS'
OUT = HERE / 'analysis'; OUT.mkdir(exist_ok=True)

# 1-based column positions in the search output, as documented above
C_PROT1, C_PROT2, C_POS1, C_POS2, C_SCORE, C_SITES = 3, 4, 8, 9, 27, 29
TOP_N = 20          # how deep to look for decoys among the best-scoring PSMs

SAMPLES = {'SB31': 'SB31', 'P14_ADPr': 'P14-7 (+ADPr)'}


def read_psms(path):
    """One row per crosslink PSM, with its target/decoy class."""
    rows = []
    for line in path.read_text().splitlines():
        f = line.split('\t')
        if len(f) < C_SITES:
            continue
        p1, p2 = f[C_PROT1 - 1], f[C_PROT2 - 1]
        d1, d2 = p1.startswith('reverse_'), p2.startswith('reverse_')
        try:
            score = float(f[C_SCORE - 1])
        except ValueError:
            continue
        rows.append(dict(
            cls='DD' if (d1 and d2) else ('TT' if not (d1 or d2) else 'TD'),
            score=score, site=f[C_SITES - 1],
            # unordered residue-pair key, for replicate overlap
            pair=tuple(sorted((f[C_POS1 - 1], f[C_POS2 - 1]))),
        ))
    return pd.DataFrame(rows)


def fdr(c):
    """Walzthoeni 2012 target-decoy FDR. None when there are no TT at all."""
    return (c['TD'] - c['DD']) / c['TT'] if c['TT'] else None


def main():
    recs, per_rep = [], {}
    for key, label in SAMPLES.items():
        files = sorted(XLDIR.glob(f'*_{key}_ZL_*_pairs_acid_K.xls'))
        if not files:
            raise SystemExit(f'no files matched sample {key} in {XLDIR}')
        pooled = []
        for path in files:
            df = read_psms(path)
            pooled.append(df)
            c = Counter(df.cls)
            recs.append(dict(sample=label, replicate=path.name, scope='replicate',
                             n_psm=len(df), **_stats(df, c)))
        allf = pd.concat(pooled, ignore_index=True)
        c = Counter(allf.cls)
        recs.append(dict(sample=label, replicate='pooled', scope='pooled',
                         n_psm=len(allf), **_stats(allf, c)))

        # reproducibility: does appearing in all 3 replicates enrich for targets?
        seen = [set(d[d.cls == 'TT'].pair) | set(d[d.cls != 'TT'].pair) for d in pooled]
        in_all = set.intersection(*seen)
        sub = allf[allf.pair.isin(in_all)]
        base = (allf.cls == 'TT').mean()
        rep = (sub.cls == 'TT').mean() if len(sub) else float('nan')
        per_rep[label] = dict(n_pairs_in_all_3=len(in_all),
                              tt_frac_all=round(base, 4),
                              tt_frac_in_all_3=round(rep, 4),
                              enrichment=round(rep / base, 2) if base else None)

    tab = pd.DataFrame(recs)
    tab.to_csv(OUT / 'xlms_screen_fdr.csv', index=False)
    print(tab.to_string(index=False))
    print()
    for label, d in per_rep.items():
        print(f'{label}: present in all 3 replicates -> {d["n_pairs_in_all_3"]} residue pairs; '
              f'target fraction {d["tt_frac_all"]:.3f} overall vs {d["tt_frac_in_all_3"]:.3f} '
              f'among reproducible ({d["enrichment"]}x)')
    print(f'\nwrote {OUT / "xlms_screen_fdr.csv"}')


def _stats(df, c):
    tt, td, dd = c['TT'], c['TD'], c['DD']
    top = df.sort_values('score', ascending=False).head(TOP_N)
    f = fdr(c)
    return dict(
        TT=tt, TD=td, DD=dd,
        fdr=round(f, 4) if f is not None else None,
        # random matching gives TT:TD:DD = 1:2:1; these are TD and DD per TT
        td_per_tt=round(td / tt, 3) if tt else None,
        dd_per_tt=round(dd / tt, 3) if tt else None,
        decoys_in_top20=int((top.cls != 'TT').sum()),
    )


if __name__ == '__main__':
    main()
