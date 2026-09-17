#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Which constructs should we crosslink to map eraser / reader / writer exposure?

THE QUESTION
------------
Different PARP14 constructs arrange their domains differently. If those
differences produce *different feasible crosslinks*, then XL-MS on a well-chosen
panel of constructs reads out the architecture -- and specifically the exposure
of the eraser (MD1), the readers (MD2 / MD3 / WWE) and the writer (ART).

This script answers three sub-questions, in order:

  1. Which INTER-DOMAIN crosslinks are feasible in each construct?
     Intra-domain crosslinks report on a fold that does not change between
     constructs; only inter-domain ones report on architecture.

  2. Which of those DISCRIMINATE between constructs?
     A crosslink feasible in every construct that contains both its domains is
     uninformative for construct comparison. The useful ones are feasible in
     some and not others.

  3. Does the crosslink signature actually TRACK site exposure?
     This is the part that decides whether the whole idea works, and it is
     tested rather than assumed: correlate each site's inter-domain crosslink
     count against its solid-angle accessibility from Q3/Q4. A correlation means
     crosslink pattern is a usable proxy for exposure; no correlation means
     crosslinks report on contacts but NOT on how open the pocket is, and the
     two must be measured separately.

INPUTS (produced by analyze_lys_contacts.py --sasd, e.g. via run_crosslinks_all.sh)
  data/{set}_lys_sasd.csv        SASD-filtered crosslinks, FL numbering
  data/enzyme_dominance_saa.csv  per-site accessibility (Q3/Q4)

Usage:
    python predict_crosslink_discrimination.py
    python predict_crosslink_discrimination.py --min-persistence 0.5
"""

import os
import sys
import csv
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from argparse import ArgumentParser

import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sim_registry as reg

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(CWD, 'data')
sys.path.insert(0, CWD)
from _fig_layout import get_fig_dir as _get_fig_dir

# Functional role of each structured domain. The whole point of the panel is to
# separate these three, so crosslinks BETWEEN roles are worth more than within.
ROLE = {
    'MD1': 'eraser',
    'MD2': 'reader', 'MD3': 'reader', 'WWE': 'reader',
    'ART': 'writer',
    'RRM1': 'scaffold', 'RRM2': 'scaffold', 'RRM3': 'scaffold',
    'KH1-6': 'scaffold', 'KH7a': 'scaffold', 'KHb-KH8': 'scaffold',
}
ROLE_COLOR = {'eraser': '#e6194b', 'reader': '#3cb44b',
              'writer': '#f58231', 'scaffold': '#7f8c8d', 'mixed': '#4363d8'}

# Which catalytic site each domain carries (for the exposure correlation).
DOMAIN_SITE = {'MD1': 'MD1', 'MD2': 'MD2', 'MD3': 'MD3', 'WWE': 'WWE', 'ART': 'ART'}


def domain_of(resid):
    """FL domain containing this residue, or None if it sits in a linker."""
    for name, (lo, hi) in reg.FL_DOMAINS.items():
        if lo <= resid <= hi:
            return name
    return None


def load_crosslinks(min_persistence):
    """set_key -> {(dom_i, dom_j): [(resid_i, resid_j, persistence, sasd), ...]}

    Only reachable, inter-domain pairs are kept.
    """
    out, skipped, malformed = {}, [], []
    for key in reg.SETS:
        path = os.path.join(DATA_PATH, f'{key}_lys_sasd.csv')
        if not os.path.isfile(path):
            skipped.append(key)
            continue
        # Refuse a misaligned CSV rather than silently reading shifted columns.
        # A header/row width mismatch once made every construct report zero
        # inter-domain crosslinks, which looks exactly like a real negative
        # result -- fail loud instead.
        with open(path) as fh:
            head = fh.readline().rstrip('\n').split(',')
            first = fh.readline().rstrip('\n').split(',')
        if first and first != [''] and len(head) != len(first):
            malformed.append((key, len(head), len(first)))
            continue
        missing = {'domain_i', 'domain_j', 'inter_domain'} - set(head)
        if missing:
            malformed.append((key, f"missing {','.join(sorted(missing))}", ''))
            continue
        per_pair = {}
        with open(path) as fh:
            for r in csv.DictReader(fh):
                if r['reachable'] != '1':
                    continue
                p = float(r['euclid_persistence'])
                if p < min_persistence:
                    continue
                ri, rj = int(r['resid_i']), int(r['resid_j'])
                di, dj = domain_of(ri), domain_of(rj)
                if di is None or dj is None or di == dj:
                    continue      # intra-domain or linker: not architectural
                dp = tuple(sorted((di, dj)))
                sasd = float(r['sasd_min_nm']) if r['sasd_min_nm'] != 'inf' else np.inf
                per_pair.setdefault(dp, []).append((ri, rj, p, sasd))
        out[key] = per_pair
    if malformed:
        msg = '\n'.join(f"    {k}: header {a} vs row {b}" if b != '' else f"    {k}: {a}"
                         for k, a, b in malformed)
        sys.exit(f"ERROR: {len(malformed)} malformed *_lys_sasd.csv file(s):\n{msg}\n"
                 f"These predate the domain-column fix. Regenerate with:\n"
                 f"    python analyze_lys_contacts.py --set <name> --sasd --sasd-top 500")
    return out, skipped


def load_exposure():
    """set_key -> site -> mean SAA, from the Q3/Q4 table."""
    path = os.path.join(DATA_PATH, 'enzyme_dominance_saa.csv')
    if not os.path.isfile(path):
        return {}
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out[r['set']] = {s: float(r[f'{s}_mean'])
                             for s in reg.SITE_NAMES if int(r[f'{s}_n']) > 0}
    return out


def pair_role(dp):
    ra, rb = ROLE.get(dp[0], 'scaffold'), ROLE.get(dp[1], 'scaffold')
    if ra == rb:
        return ra
    if 'scaffold' in (ra, rb):
        return rb if ra == 'scaffold' else ra
    return 'mixed'


def units_present(set_key):
    return set(reg.SETS[set_key]['units'])


UNIT_OF_DOMAIN = {
    'RRM1': 'rrm1', 'RRM2': 'rrm2', 'RRM3': 'rrm3', 'KH1-6': 'kh1-kh6',
    'KH7a': 'kh7a', 'MD1': 'md1l1', 'MD2': 'md2', 'MD3': 'md3',
    'KHb-KH8': 'khb-kh8', 'WWE': 'wwe', 'ART': 'art',
}


def construct_can_form(set_key, dp):
    """True if this construct contains BOTH domains of the pair.

    Essential for the discrimination logic: a crosslink absent because the
    domain is not in the construct is not evidence about architecture.
    """
    u = units_present(set_key)
    return all(UNIT_OF_DOMAIN.get(d) in u for d in dp)


def main():
    ap = ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--min-persistence', type=float, default=0.25,
                    help='Minimum Euclidean persistence for a crosslink to count '
                         '(default 0.25). A crosslink seen in a quarter of frames '
                         'is plausibly capturable; below that the prediction is '
                         'too marginal to spend bench time on.')
    args = ap.parse_args()

    xl, skipped = load_crosslinks(args.min_persistence)
    expo = load_exposure()
    if not xl:
        sys.exit("ERROR: no {set}_lys_sasd.csv found. Run:\n"
                 "    bash run_crosslinks_all.sh")

    sets = sorted(xl)
    print(f"  Loaded {len(sets)} construct(s) with crosslink data"
          + (f"; {len(skipped)} without: {', '.join(skipped[:6])}"
             + (' ...' if len(skipped) > 6 else '') if skipped else ''))
    print(f"  Persistence floor: {args.min_persistence}\n")

    # ---- 1. inter-domain crosslinks per construct -------------------------
    all_pairs = sorted({dp for v in xl.values() for dp in v})
    print("=" * 96)
    print("1. INTER-DOMAIN CROSSLINKS PER CONSTRUCT")
    print("=" * 96)
    print(f"{'construct':20s}{'inter-dom XLs':>14}{'domain pairs':>14}"
          f"{'role-crossing':>15}   roles touched")
    print("-" * 96)
    rows = []
    for k in sets:
        n_xl = sum(len(v) for v in xl[k].values())
        pairs = list(xl[k])
        crossing = [dp for dp in pairs if pair_role(dp) == 'mixed']
        roles = sorted({pair_role(dp) for dp in pairs})
        rows.append(dict(set=k, n_xl=n_xl, n_pairs=len(pairs),
                         n_cross=len(crossing), roles=roles))
        print(f"{k:20s}{n_xl:>14d}{len(pairs):>14d}{len(crossing):>15d}   "
              f"{', '.join(roles) if roles else '-'}")

    # ---- 2. discriminating crosslinks -------------------------------------
    # A domain pair discriminates if, among constructs that COULD form it,
    # some do and some do not.
    print("\n" + "=" * 96)
    print("2. DISCRIMINATING DOMAIN PAIRS")
    print("   (feasible in some constructs but not others, among those containing both domains)")
    print("=" * 96)
    disc = []
    for dp in all_pairs:
        capable = [k for k in sets if construct_can_form(k, dp)]
        if len(capable) < 2:
            continue
        have = [k for k in capable if dp in xl[k]]
        lack = [k for k in capable if dp not in xl[k]]
        if have and lack:
            disc.append(dict(pair=dp, role=pair_role(dp), have=have, lack=lack,
                             n_capable=len(capable)))
    disc.sort(key=lambda d: (d['role'] != 'mixed', -min(len(d['have']), len(d['lack']))))

    if disc:
        print(f"{'domain pair':24s}{'role':10s}{'forms in':>10}{'absent in':>11}   "
              f"discriminates")
        print("-" * 96)
        for d in disc:
            print(f"{d['pair'][0] + '-' + d['pair'][1]:24s}{d['role']:10s}"
                  f"{len(d['have']):>10d}{len(d['lack']):>11d}   "
                  f"{', '.join(d['have'][:3])}{' ...' if len(d['have']) > 3 else ''}"
                  f"  vs  {', '.join(d['lack'][:3])}{' ...' if len(d['lack']) > 3 else ''}")
    else:
        print("  None. Every domain pair behaves the same way in every construct")
        print("  that contains it -- crosslinking the panel would not distinguish them.")

    # ---- 3. does the signature track exposure? ----------------------------
    print("\n" + "=" * 96)
    print("3. DOES THE CROSSLINK SIGNATURE TRACK SITE EXPOSURE?")
    print("=" * 96)
    if not expo:
        print("  data/enzyme_dominance_saa.csv missing -- run predict_enzyme_dominance.py")
    else:
        print(f"{'site':6s}{'n constructs':>14}{'Pearson r':>12}{'   reading'}")
        print("-" * 96)
        for site in reg.SITE_NAMES:
            dom = next((d for d, s in DOMAIN_SITE.items() if s == site), None)
            xs, ys = [], []
            for k in sets:
                if k not in expo or site not in expo[k]:
                    continue
                n = sum(len(v) for dp, v in xl[k].items() if dom in dp)
                xs.append(n)
                ys.append(expo[k][site])
            if len(xs) >= 4 and np.std(xs) > 0 and np.std(ys) > 0:
                r = float(np.corrcoef(xs, ys)[0, 1])
                if abs(r) >= 0.6:
                    rd = 'strong - XL count is a usable proxy for exposure'
                elif abs(r) >= 0.35:
                    rd = 'weak - suggestive, not a proxy on its own'
                else:
                    rd = 'none - XLs report contacts, NOT pocket openness'
                print(f"{site:6s}{len(xs):>14d}{r:>12.3f}   {rd}")
            else:
                print(f"{site:6s}{len(xs):>14d}{'-':>12}   too few constructs / no variance")

    # ---- 4. recommended panel ---------------------------------------------
    print("\n" + "=" * 96)
    print("4. RECOMMENDED CROSSLINKING PANEL")
    print("=" * 96)
    # Greedy: repeatedly add the construct that resolves the most
    # still-unresolved discriminating pairs.
    if disc:
        unresolved = {d['pair'] for d in disc}
        chosen, guard = [], 0
        while unresolved and guard < 12:
            guard += 1
            best, best_gain = None, 0
            for k in sets:
                if k in chosen:
                    continue
                gain = sum(1 for dp in unresolved
                           if construct_can_form(k, dp) and dp in xl[k])
                if gain > best_gain:
                    best, best_gain = k, gain
            if not best:
                break
            chosen.append(best)
            unresolved -= {dp for dp in list(unresolved)
                           if construct_can_form(best, dp) and dp in xl[best]}
        print(f"  Greedy cover of the {len(disc)} discriminating domain pairs:\n")
        for i, k in enumerate(chosen, 1):
            got = [dp for dp in (d['pair'] for d in disc)
                   if construct_can_form(k, dp) and dp in xl[k]]
            print(f"  {i}. {k:20s} contributes {len(got):2d} pair(s): "
                  f"{', '.join(a + '-' + b for a, b in got[:4])}"
                  f"{' ...' if len(got) > 4 else ''}")
        if unresolved:
            print(f"\n  Not covered by any construct: "
                  f"{', '.join(a + '-' + b for a, b in sorted(unresolved))}")
    else:
        print("  No discriminating pairs -- no panel would separate these constructs.")

    # ---- outputs -----------------------------------------------------------
    out_csv = os.path.join(DATA_PATH, 'crosslink_discrimination.csv')
    with open(out_csv, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['domain_pair', 'role', 'n_constructs_capable',
                    'forms_in', 'absent_in', 'discriminating'])
        for dp in all_pairs:
            capable = [k for k in sets if construct_can_form(k, dp)]
            have = [k for k in capable if dp in xl[k]]
            lack = [k for k in capable if dp not in xl[k]]
            w.writerow([f'{dp[0]}-{dp[1]}', pair_role(dp), len(capable),
                        ';'.join(have), ';'.join(lack),
                        int(bool(have and lack))])
    print(f"\n  Saved: {os.path.relpath(out_csv, CWD)}")

    # shopping list of the actual residue pairs to watch
    out_xl = os.path.join(DATA_PATH, 'crosslink_shopping_list.csv')
    with open(out_xl, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['set', 'domain_pair', 'role', 'resid_i', 'resid_j',
                    'persistence', 'sasd_min_nm'])
        for k in sets:
            for dp, lst in sorted(xl[k].items()):
                for ri, rj, p, s in sorted(lst, key=lambda t: -t[2]):
                    w.writerow([k, f'{dp[0]}-{dp[1]}', pair_role(dp), ri, rj,
                                f'{p:.4f}', f'{s:.3f}'])
    print(f"  Saved: {os.path.relpath(out_xl, CWD)}")

    figdir = str(_get_fig_dir('07_lysine_contacts', subname='discrimination'))
    make_figure(sets, xl, all_pairs, figdir)
    print(f"  Figures -> {os.path.relpath(figdir, CWD)}")


def make_figure(sets, xl, all_pairs, figdir):
    """Construct x domain-pair matrix of inter-domain crosslink counts."""
    if not all_pairs:
        return
    mat = np.full((len(sets), len(all_pairs)), np.nan)
    for i, k in enumerate(sets):
        for j, dp in enumerate(all_pairs):
            if not construct_can_form(k, dp):
                continue          # domain absent -> stays NaN (grey), not zero
            mat[i, j] = len(xl[k].get(dp, []))

    fig, ax = plt.subplots(figsize=(0.46 * len(all_pairs) + 5,
                                    0.34 * len(sets) + 3.2))
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('#d8d8d8')
    im = ax.imshow(np.ma.masked_invalid(mat), aspect='auto', cmap=cmap)
    ax.set_xticks(range(len(all_pairs)))
    ax.set_xticklabels([f'{a}-{b}' for a, b in all_pairs], rotation=90, fontsize=8)
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels(sets, fontsize=8)
    for j, dp in enumerate(all_pairs):
        ax.get_xticklabels()[j].set_color(ROLE_COLOR[pair_role(dp)])
    for i in range(len(sets)):
        for j in range(len(all_pairs)):
            if np.isnan(mat[i, j]):
                continue
            ax.text(j, i, int(mat[i, j]), ha='center', va='center', fontsize=7,
                    color='white' if mat[i, j] < np.nanmax(mat) * 0.55 else 'black')
    ax.set_title('Feasible inter-domain crosslinks per construct\n'
                 'grey = domain not in construct;  0 = domains present but never '
                 'crosslinkable\n'
                 'label colour = role  (red eraser, green reader, orange writer, '
                 'blue role-crossing)', fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.7, label='n crosslinks')
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(figdir, f'crosslink_discrimination.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
