#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Which predicted crosslinks are CONSERVED across constructs?

Q10 asks which crosslinks DIFFER between constructs -- those discriminate
between architectures and tell you which constructs to compare. This script asks
the complement, which is what you want when picking what to actually test:

    which crosslinks form in (almost) every construct that could form them?

A crosslink that survives across the whole panel is a robust structural
restraint: it does not depend on which truncation happens to express well, so it
is the safest thing to take to the bench. One that appears in a single construct
is either a real architectural signature (Q10's subject) or noise.

CONSERVATION
------------
For each unique lysine pair, over the constructs that contain BOTH its domains:

    conservation = (constructs where it is reachable) / (constructs that could form it)

The denominator matters. A pair absent because the domain is not in the
construct is not evidence against it, so it must not count against conservation.

Usage:
    python consolidate_crosslinks.py
    python consolidate_crosslinks.py --min-capable 4 --min-persistence 0.3
"""

import os
import sys
import csv
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from collections import defaultdict

import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sim_registry as reg

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(CWD, 'data')
sys.path.insert(0, CWD)
from _fig_layout import get_fig_dir as _get_fig_dir

UNIT_OF_DOMAIN = {
    'RRM1': 'rrm1', 'RRM2': 'rrm2', 'RRM3': 'rrm3', 'KH1-6': 'kh1-kh6',
    'KH7a': 'kh7a', 'MD1': 'md1l1', 'MD2': 'md2', 'MD3': 'md3',
    'KHb-KH8': 'khb-kh8', 'WWE': 'wwe', 'ART': 'art',
}
ROLE = {'MD1': 'eraser', 'MD2': 'reader', 'MD3': 'reader', 'WWE': 'reader',
        'ART': 'writer'}


def domain_of(resid):
    for name, (lo, hi) in reg.FL_DOMAINS.items():
        if lo <= resid <= hi:
            return name
    return None


def load_all():
    """set -> {(i,j): (persistence, sasd, dom_i, dom_j, reachable)}"""
    out = {}
    for f in sorted(glob.glob(os.path.join(DATA_PATH, '*_lys_sasd.csv'))):
        s = os.path.basename(f).replace('_lys_sasd.csv', '')
        if s not in reg.SETS:
            continue
        with open(f) as fh:
            head = fh.readline().rstrip('\n').split(',')
            first = fh.readline().rstrip('\n').split(',')
        if first and first != [''] and len(head) != len(first):
            print(f"  WARN: {s} has a malformed CSV (header {len(head)} vs row "
                  f"{len(first)}) -- skipped")
            continue
        d = {}
        with open(f) as fh:
            for r in csv.DictReader(fh):
                i, j = int(r['resid_i']), int(r['resid_j'])
                di = r.get('domain_i') or domain_of(i)
                dj = r.get('domain_j') or domain_of(j)
                d[(i, j)] = (float(r['euclid_persistence']),
                             float(r['sasd_min_nm']) if r['sasd_min_nm'] != 'inf' else np.inf,
                             di, dj, r['reachable'] == '1')
        out[s] = d
    return out


def can_form(set_key, di, dj):
    u = set(reg.SETS[set_key]['units'])
    return (UNIT_OF_DOMAIN.get(di) in u) and (UNIT_OF_DOMAIN.get(dj) in u)


def main():
    ap = ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--min-capable', type=int, default=3,
                    help='Only report pairs that at least N constructs could '
                         'form (default 3). Below that "conservation" is not '
                         'meaningful.')
    ap.add_argument('--min-persistence', type=float, default=0.25)
    args = ap.parse_args()

    data = load_all()
    if not data:
        sys.exit("ERROR: no *_lys_sasd.csv found. Run run_crosslinks_all.sh")
    sets = sorted(data)
    print(f"  {len(sets)} constructs: {', '.join(sets)}\n")

    # per-set summary
    print(f"{'construct':20s}{'scored':>8}{'reachable':>11}{'buried FP':>11}"
          f"{'inter-dom':>11}")
    print('-' * 62)
    tot = defaultdict(int)
    for s in sets:
        d = data[s]
        n = len(d)
        r = sum(1 for v in d.values() if v[4])
        inter = sum(1 for v in d.values() if v[4] and v[2] and v[3] and v[2] != v[3])
        tot['n'] += n; tot['r'] += r; tot['i'] += inter
        print(f"{s:20s}{n:>8}{r:>11}{n - r:>11}{inter:>11}")
    print('-' * 62)
    print(f"{'TOTAL':20s}{tot['n']:>8}{tot['r']:>11}{tot['n']-tot['r']:>11}"
          f"{tot['i']:>11}")
    print(f"\n  SASD removed {tot['n']-tot['r']:,} of {tot['n']:,} "
          f"({100*(tot['n']-tot['r'])/tot['n']:.0f}%) Euclidean-passing pairs as "
          f"buried.\n")

    # conservation
    pair_dom = {}
    forms = defaultdict(set)
    persist = defaultdict(list)
    sasd = defaultdict(list)
    for s, d in data.items():
        for (i, j), (p, sd, di, dj, ok) in d.items():
            if not (di and dj):
                continue
            pair_dom[(i, j)] = (di, dj)
            if ok and p >= args.min_persistence:
                forms[(i, j)].add(s)
                persist[(i, j)].append(p)
                sasd[(i, j)].append(sd)

    rows = []
    for pair, (di, dj) in pair_dom.items():
        capable = [s for s in sets if can_form(s, di, dj)]
        if len(capable) < args.min_capable:
            continue
        f = forms.get(pair, set())
        rows.append(dict(
            resid_i=pair[0], resid_j=pair[1], domain_i=di, domain_j=dj,
            inter_domain=int(di != dj),
            n_capable=len(capable), n_forms=len(f),
            conservation=len(f) / len(capable),
            mean_persistence=float(np.mean(persist[pair])) if pair in persist else 0.0,
            min_sasd=float(np.min(sasd[pair])) if pair in sasd else float('inf'),
            forms_in=';'.join(sorted(f)),
        ))
    rows.sort(key=lambda r: (-r['conservation'], -r['mean_persistence']))

    inter = [r for r in rows if r['inter_domain']]
    print("=" * 84)
    print("MOST CONSERVED INTER-DOMAIN CROSSLINKS")
    print("   (formable in every construct that contains both domains)")
    print("=" * 84)
    print(f"{'pair':16s}{'domains':22s}{'capable':>8}{'forms':>7}{'cons':>7}"
          f"{'persist':>9}{'SASD':>7}")
    print('-' * 84)
    for r in inter[:15]:
        pair = f"K{r['resid_i']}-K{r['resid_j']}"
        doms = f"{r['domain_i']}-{r['domain_j']}"
        print(f"{pair:16s}{doms:22s}{r['n_capable']:>8}{r['n_forms']:>7}"
              f"{r['conservation']:>7.2f}{r['mean_persistence']:>9.2f}"
              f"{r['min_sasd']:>7.2f}")

    out = os.path.join(DATA_PATH, 'crosslink_conservation.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\n  Saved: {os.path.relpath(out, CWD)}  ({len(rows)} pairs)")

    make_figures(rows, inter, data, sets)


def make_figures(rows, inter, data, sets):
    figdir = str(_get_fig_dir('07_lysine_contacts', subname='conservation'))

    # 1. conservation distribution, inter vs intra
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.3))
    intra = [r['conservation'] for r in rows if not r['inter_domain']]
    inter_c = [r['conservation'] for r in inter]
    bins = np.linspace(0, 1, 21)
    a1.hist([inter_c, intra], bins=bins, stacked=False,
            label=[f'inter-domain (n={len(inter_c)})', f'intra-domain (n={len(intra)})'],
            color=['#4363d8', '#bbbbbb'], edgecolor='black', linewidth=0.4)
    a1.set_xlabel('conservation  (constructs forming / constructs capable)')
    a1.set_ylabel('crosslinks')
    a1.set_title('How reproducible is each predicted crosslink?', fontsize=11)
    a1.legend(frameon=False, fontsize=8.5)

    # 2. per-construct SASD filtering rate
    ss = sorted(data)
    kept = [100 * sum(1 for v in data[s].values() if v[4]) / max(1, len(data[s]))
            for s in ss]
    a2.barh(range(len(ss)), kept, color='#2ca02c', edgecolor='black', linewidth=0.4)
    a2.set_yticks(range(len(ss))); a2.set_yticklabels(ss, fontsize=7.5)
    a2.invert_yaxis()
    a2.axvline(np.mean(kept), ls='--', color='#c0392b', lw=1)
    a2.set_xlabel('% of Euclidean-passing pairs with a solvent path')
    a2.set_title(f'SASD survival rate (mean {np.mean(kept):.0f}%)', fontsize=11)
    a2.grid(axis='x', alpha=0.25)
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(figdir, f'crosslink_conservation.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: crosslink_conservation.png/.svg")
    print(f"  Figures -> {os.path.relpath(figdir, CWD)}")


if __name__ == '__main__':
    main()
