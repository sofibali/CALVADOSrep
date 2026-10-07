#!/usr/bin/env python3
"""Extract every crosslink from the Ashok 2022 supplement into a flat CSV.

Rebuilt 2026-10-07; the original was lost with a scratchpad. Verified to
reproduce data/published_xlinks_v2.csv byte for byte.

    python data/extract_v2.py bcj-2021-0722_supp.xlsx data/published_xlinks_v2.csv

LOCATE BLOCKS BY HEADER, NEVER BY COLUMN LETTER. Each of the three Summary
sheets holds three side-by-side blocks -- intra-PARP9, inter-protein,
intra-DTX3L -- and they do NOT sit at the same columns in every sheet. In
"Summary FL PARP9 + D3" the inter-protein block's title sits in column K
while its Peptides column is J; the other two blocks have title and Peptides
in the same column. An earlier pass hard-coded the letter and silently got
nothing from two of the three sheets, reporting 3 inter-protein crosslinks
instead of 51. Blocks are therefore found by scanning the header row for
"Peptides", and each block runs to the next one.

Per-replicate detection is an "x" in the Sample 1-4 columns. The authors'
own ">=2/4" column is AUTHORITATIVE for the reproducibility flag and is not
recomputed from the x marks: the two disagree on 4 of 300 crosslinks, and
the published analysis follows the authors. The script prints those 4 so
they stay visible. One of them is inter-protein and so reaches the scoring:
FL PARP9 + D3, 567-363, marked reproducible on one x mark.
"""
import csv, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readxlsx import rows

SHEETS = ['Summary FL PARP9 + D3',
          'Summary FL PARP9 + D3RD',
          'Summary FL PARP9 + FL DTX3L']


def col_index(letter):
    n = 0
    for ch in letter:
        n = n * 26 + (ord(ch) - 64)
    return n


def blocks_of(sheet):
    """[(block name, peptides col, [sample cols], repro col)] for one sheet."""
    title_row, header_row = sheet[0], None
    for r in sheet[:6]:
        if 'Peptides' in r.values():
            header_row = r
            break
    if header_row is None:
        raise SystemExit('no header row containing "Peptides"')

    pep = sorted((c for c, v in header_row.items() if v == 'Peptides'),
                 key=col_index)
    out = []
    for i, c in enumerate(pep):
        lo = col_index(c)
        hi = col_index(pep[i + 1]) if i + 1 < len(pep) else 10**6
        inblock = lambda cc: lo <= col_index(cc) < hi
        samples = sorted((cc for cc, v in header_row.items()
                          if inblock(cc) and v.startswith('Sample')),
                         key=col_index)
        repro = [cc for cc, v in header_row.items()
                 if inblock(cc) and v.replace(' ', '').upper() == 'YES/NO']
        # the block's title sits at or just after its Peptides column
        name = next((v for cc, v in sorted(title_row.items(), key=lambda kv: col_index(kv[0]))
                     if inblock(cc)), None)
        if name is None or not samples:
            continue
        out.append((name, c, samples, repro[0] if repro else None))
    return out


def main(xlsx, out_csv):
    recs, disagree = [], []
    for name in SHEETS:
        sheet = list(rows(xlsx, name))
        construct = name.replace('Summary ', '')
        for block, pep_c, sample_cs, repro_c in blocks_of(sheet):
            for r in sheet:
                v = r.get(pep_c, '')
                if '-' not in v:
                    continue
                a, _, b = v.partition('-')
                if not (a.strip().isdigit() and b.strip().isdigit()):
                    continue          # skips the header and any stray text
                n = sum(1 for cc in sample_cs if r.get(cc, '').lower() == 'x')
                flagged = bool(repro_c and r.get(repro_c, ''))
                if flagged != (n >= 2):
                    disagree.append((construct, block, v, n, flagged))
                recs.append(dict(construct=construct, block=block,
                                 res_a=int(a), res_b=int(b),
                                 n_samples=n, repro_2of4=flagged))

    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, ['construct', 'block', 'res_a', 'res_b',
                               'n_samples', 'repro_2of4'])
        w.writeheader()
        for r in recs:
            w.writerow({**r, 'repro_2of4': str(r['repro_2of4'])})

    if disagree:
        print(f'  NOTE: author >=2/4 flag disagrees with the x marks on '
              f'{len(disagree)} of {len(recs)}; the flag wins:')
        for c, b, v, n, f in disagree:
            print(f'    {c} / {b} / {v}: {n} x marks, flag {f}')
    inter = [r for r in recs if r['block'] == 'P9-Dtx3L']
    print(f'{len(recs)} crosslinks -> {out_csv}')
    print(f'  inter-protein {len(inter)}, reproducible '
          f'{sum(r["repro_2of4"] for r in inter)}, '
          f'unique pairs {len({(r["res_a"], r["res_b"]) for r in inter})}')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'bcj-2021-0722_supp.xlsx',
         sys.argv[2] if len(sys.argv) > 2 else 'data/published_xlinks_v2.csv')
