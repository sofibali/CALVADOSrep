#!/usr/bin/env python3
"""Inter-chain contact maps resolved by DOMAIN PAIR.

Two quantities per (domain of chain A) x (domain of chain B) pair:

    pct_contact    % of frames with min C-alpha--C-alpha distance < 1.0 nm
    mean_min_nm    mean over frames of that minimum distance

They answer different questions and both are needed. pct_contact is an
occupancy: how much of the ensemble has these two domains touching. mean_min_nm
is a geometry: how close they get on average, whether or not they ever touch. A
pair can sit at 2 nm for the whole run -- never in contact, but clearly the
closest approach of the two chains -- and only the second column shows it.

DOMAIN BOUNDARIES are the restrained blocks of <root>/<set>/input/domains.yaml,
the same definition figure_face_separation.py uses: those are what the
simulation actually holds rigid. Residues outside every block are linker and are
excluded, so a flexible tail brushing the partner does not register as a domain
contact.

DISTANCES are minimum-image throughout (the chains are not kept in the same
periodic image -- see pbc_center.py), and are computed between restrained-block
beads only.

Both engines are supported. CALVADOS reads the bead trajectory directly; HyRes
selects name CA, PARP9 being the first 854 and DTX3L the rest, matching
rescore_crosslinks.dists_hyres.

Usage:
    python figure_domain_contacts.py                        # every binding set
    python figure_domain_contacts.py --root binding_go
    python figure_domain_contacts.py --sets p9_dtx3l --target-frames 600
    python figure_domain_contacts.py --hyres hyres/runs/matched/docked_rep1
"""
import sys, warnings
from pathlib import Path
from collections import defaultdict
from argparse import ArgumentParser
import numpy as np, yaml, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'
FIG = OUT / 'figures'; FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology
import analyze_convergence as ac
from figure_face_separation import DOMAIN_NAMES, load_domains, NS_PER_RAW_FRAME, SKIP_FRAMES

CONTACT_NM = 1.0                  # same cutoff as analyze_binding / figure_face_separation
DEFAULT_TARGET_FRAMES = 300
HYRES_BOX_A = 400.0               # 40 nm, as built
INK = {'primary': '#0B0B0B', 'secondary': '#52514E', 'grid': '#E1E0D9'}


def domain_index(blocks, resids, base, chain):
    """[(name, global bead indices)] for each restrained block of one chain."""
    names = DOMAIN_NAMES.get(chain, [f'D{i+1}' for i in range(len(blocks))])
    out = []
    for i, (lo, hi) in enumerate(blocks):
        sel = np.where((resids >= lo) & (resids <= hi))[0]
        if len(sel) < 4:
            continue
        out.append((names[i] if i < len(names) else f'D{i+1}', sel + base))
    return out


def _accumulate(P, boxv, doms, acc):
    """Per-frame min distance for every domain pair, appended to acc.

    Both INTER-chain pairs (i < j) and INTRA-chain pairs (i == j, two different
    domains of the same chain). The intra-chain maps answer the companion
    question to the interface one: how the domains of a single protein are
    arranged relative to each other, and whether binding rearranges them.

    The full chain-vs-chain distance matrix is computed once per chain pair and
    then sliced, rather than calling distance_array per domain pair: one call on
    ~700x600 beads is far cheaper than 48 calls on the blocks.
    """
    # --- intra-chain: domain i vs domain j of the SAME chain -----------------
    for ci, di in doms:
        idx = np.concatenate([ix for _, ix in di])
        D = distance_array(P[idx], P[idx], box=boxv) / 10.0
        offs, o = [], 0
        for nm, ix in di:
            offs.append((nm, o, o + len(ix))); o += len(ix)
        for a in range(len(offs)):
            for b in range(a + 1, len(offs)):
                na, la, ha = offs[a]; nb, lb, hb = offs[b]
                acc[('intra', ci, na, ci, nb)].append(float(D[la:ha, lb:hb].min()))

    # --- inter-chain ---------------------------------------------------------
    for i in range(len(doms)):
        for j in range(i + 1, len(doms)):
            (ci, di), (cj, dj) = doms[i], doms[j]
            idx_i = np.concatenate([ix for _, ix in di])
            idx_j = np.concatenate([ix for _, ix in dj])
            D = distance_array(P[idx_i], P[idx_j], box=boxv) / 10.0
            # whole-chain-pair minimum: the frame's association state, recorded
            # alongside the per-domain minima so the two stay index-aligned
            acc[('__chainmin__', ci, cj, '', '')].append(float(D.min()))
            oi = 0
            for ni, ix in di:
                oj = 0
                for nj, jx in dj:
                    sub = D[oi:oi + len(ix), oj:oj + len(jx)]
                    acc[('inter', ci, ni, cj, nj)].append(float(sub.min()))
                    oj += len(jx)
                oi += len(ix)


def analyse_calvados(s, root, target_frames, per_rep=False):
    chains = SETS[s if s in SETS else 'p9_dtx3l']
    blocks_by_chain = load_domains(s if (HERE / 'binding' / s).is_dir() else 'p9_dtx3l',
                                   chains)
    offs, o = [], 0
    for c in chains:
        offs.append((c, o, o + ac.LENGTHS[c])); o += ac.LENGTHS[c]

    acc = defaultdict(list)
    step_set = None
    nrep = nfr_tot = 0
    per_rep_out = []          # [(rep name, acc, nframes)] when per_rep
    for rep in sorted((HERE / root / s).glob('rep-*'),
                      key=lambda p: int(p.name.split('-')[1])):
        dcd = rep / f'{sysname(s if s in SETS else "p9_dtx3l")}.dcd'
        top = topology(rep)
        if top is None or not dcd.exists():
            continue
        cfg = yaml.safe_load(open(rep / 'config.yaml'))
        b = np.array(cfg['box'], float)
        boxv = np.array([b[0]*10, b[1]*10, b[2]*10, 90., 90., 90.], np.float32)
        u = mda.Universe(str(top), str(dcd))
        if u.atoms.n_atoms != o:
            continue
        usable = max(0, len(u.trajectory) - SKIP_FRAMES)
        if usable < 10:
            continue
        if step_set is None:                       # one stride for the whole set
            step_set = max(1, usable // target_frames) if target_frames else 1

        resids = u.atoms.resids
        doms = [(c, domain_index(blocks_by_chain[k], resids[a:b2], a, c))
                for k, (c, a, b2) in enumerate(offs)]
        racc, nfr = (defaultdict(list) if per_rep else acc), 0
        for ts in u.trajectory[SKIP_FRAMES::step_set]:
            P = u.atoms.positions.astype(np.float32)
            _accumulate(P, boxv, doms, acc if not per_rep else racc)
            nfr += 1; nfr_tot += 1
        if per_rep:
            per_rep_out.append((rep.name, racc, nfr))
        nrep += 1
    if per_rep:
        return per_rep_out
    return (acc, nrep, nfr_tot) if nrep else (None, 0, 0)


def analyse_hyres(run_dirs, target_frames):
    blocks = load_domains('p9_dtx3l', ['parp9', 'dtx3l'])
    boxv = np.array([HYRES_BOX_A]*3 + [90., 90., 90.], np.float32)
    acc = defaultdict(list)
    nrep = nfr_tot = 0
    for rd in run_dirs:
        run = Path(rd) if Path(rd).is_absolute() else HERE / rd
        top = run / 'complex_start.pdb'
        if not top.exists():
            top = run / 'start.pdb'
        dcd = run / 'system.dcd'
        if not top.exists() or not dcd.exists():
            print(f'  skip {rd}: missing topology or dcd')
            continue
        u = mda.Universe(str(top), str(dcd))
        ca = u.select_atoms('name CA')
        n9 = ac.LENGTHS['parp9']
        if ca.n_atoms != n9 + ac.LENGTHS['dtx3l']:
            print(f'  skip {rd}: {ca.n_atoms} CA, expected '
                  f'{n9 + ac.LENGTHS["dtx3l"]}')
            continue
        # HyRes PDBs renumber per chain, so index by position, not by resid
        rid9 = np.arange(1, n9 + 1)
        ridD = np.arange(1, ac.LENGTHS['dtx3l'] + 1)
        doms = [('parp9', domain_index(blocks[0], rid9, 0, 'parp9')),
                ('dtx3l', domain_index(blocks[1], ridD, n9, 'dtx3l'))]
        step = max(1, len(u.trajectory) // target_frames) if target_frames else 1
        for ts in u.trajectory[::step]:
            _accumulate(ca.positions.astype(np.float32), boxv, doms, acc)
            nfr_tot += 1
        nrep += 1
    return (acc, nrep, nfr_tot) if nrep else (None, 0, 0)


def tabulate(acc, label, nrep, nfr):
    """Per-domain-pair statistics, both unconditional and conditioned on contact.

    The unconditional mean of the minimum distance is dominated by the frames in
    which the two chains are simply apart -- for a weak binder that is most of
    the run, and the number comes back at box scale (~16 nm) no matter which
    domains are involved. The *_bound columns restrict to frames where the two
    chains touch anywhere, which is the question the map is actually asked:
    given that the complex is formed, which domains are at the interface and how
    close do they sit.
    """
    assoc = {(ci, cj): np.asarray(v) < CONTACT_NM
             for (tag, ci, cj, _, _), v in acc.items() if tag == '__chainmin__'}
    rows = []
    for key, v in acc.items():
        kind, ci, ni, cj, nj = key
        if kind == '__chainmin__':
            continue
        a = np.asarray(v)
        # An intra-chain pair is always "associated" -- it is one molecule, so
        # the conditional and unconditional columns coincide by construction.
        m = (assoc.get((ci, cj)) if kind == 'inter'
             else np.ones(len(a), bool))
        nb = int(m.sum()) if m is not None else 0
        ab = a[m] if nb else np.array([])
        rows.append(dict(
            model=label, kind=kind, chain_i=ci, domain_i=ni,
            chain_j=cj, domain_j=nj,
            n_rep=nrep, n_frames=len(a), n_bound_frames=nb,
            pct_contact=round(100.0 * float((a < CONTACT_NM).mean()), 3),
            pct_contact_bound=round(100.0 * float((ab < CONTACT_NM).mean()), 2)
                              if nb else np.nan,
            mean_min_nm=round(float(a.mean()), 2),
            mean_min_bound_nm=round(float(ab.mean()), 2) if nb else np.nan,
            median_min_bound_nm=round(float(np.median(ab)), 2) if nb else np.nan,
            closest_nm=round(float(a.min()), 2)))
    return pd.DataFrame(rows)


def heatmaps(df, label, stem):
    """One figure per chain pair: occupancy beside closest approach."""
    df = df[df.kind == 'inter']
    pairs = df[['chain_i', 'chain_j']].drop_duplicates().values.tolist()
    n = len(pairs)
    fig, axes = plt.subplots(n, 3, figsize=(19.5, 4.6 * n), squeeze=False)
    for r, (ci, cj) in enumerate(pairs):
        sub = df[(df.chain_i == ci) & (df.chain_j == cj)]
        rows = [d for d in DOMAIN_NAMES.get(ci, []) if d in set(sub.domain_i)]
        cols = [d for d in DOMAIN_NAMES.get(cj, []) if d in set(sub.domain_j)]
        piv = lambda col: (sub.pivot(index='domain_i', columns='domain_j',
                                     values=col).reindex(index=rows, columns=cols)
                           .to_numpy(float))
        for c, (col, cmap, title, unit) in enumerate([
                ('pct_contact', 'Purples', 'contact occupancy',
                 '% of ALL frames < 1.0 nm'),
                ('pct_contact_bound', 'Purples', 'interface composition',
                 '% of ASSOCIATED frames < 1.0 nm'),
                ('mean_min_bound_nm', 'Greens_r', 'closest approach',
                 'mean min Cα–Cα over associated frames (nm)')]):
            ax = axes[r, c]
            M = piv(col)
            im = ax.imshow(M, cmap=cmap, aspect='auto',
                           vmin=0 if col.startswith('pct_') else None)
            ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=45,
                                                               ha='right', fontsize=9)
            ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=9)
            ax.set_xlabel(cj.upper(), fontsize=9.5, color=INK['secondary'])
            ax.set_ylabel(ci.upper(), fontsize=9.5, color=INK['secondary'])
            # direct labels: the grid is small enough that a colourbar alone
            # would make readers estimate values off a ramp
            lo, hi = np.nanmin(M), np.nanmax(M)
            for y in range(M.shape[0]):
                for x in range(M.shape[1]):
                    v = M[y, x]
                    if np.isnan(v):
                        continue
                    dark = (v - lo) / (hi - lo + 1e-9) > .55
                    if col.startswith('mean_'):
                        dark = not dark
                    ax.text(x, y, f'{v:.2f}' if col == 'pct_contact' else f'{v:.1f}',
                            ha='center', va='center', fontsize=7.5,
                            color='white' if dark else INK['primary'])
            ax.set_title(f'{title}\n{unit}', loc='left', fontsize=10.5,
                         color=INK['primary'], fontweight='semibold', pad=8)
            fig.colorbar(im, ax=ax, fraction=.035, pad=.02)
            ax.tick_params(colors=INK['secondary'], length=0)
    fig.suptitle(f'Domain-resolved inter-chain contacts — {label}',
                 x=0.006, y=1.0, ha='left', fontsize=13.5,
                 color=INK['primary'], fontweight='semibold')
    fig.tight_layout(rect=[0, 0, 1, .97])
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'{stem}.{ext}', dpi=200, bbox_inches='tight',
                    facecolor='white')
    plt.close(fig)
    print('  wrote', FIG / f'{stem}.png')


def main():
    ap = ArgumentParser()
    ap.add_argument('--root', default='binding')
    ap.add_argument('--sets', nargs='*', default=None)
    ap.add_argument('--hyres', nargs='*', default=None)
    ap.add_argument('--label', default=None)
    ap.add_argument('--target-frames', type=int, default=DEFAULT_TARGET_FRAMES)
    ap.add_argument('--out', default='domain_contact_maps.csv')
    ap.add_argument('--per-rep', action='store_true',
                    help='one model row per replicate, for spread and convergence')
    ap.add_argument('--compare', action='store_true',
                    help='build the cross-model figure from the existing CSV')
    a = ap.parse_args()

    if a.compare:
        compare(a.out); return

    allrows = []
    if a.hyres:
        label = a.label or f'HyRes ({len(a.hyres)} rep)'
        print(f'== {label}')
        acc, nrep, nfr = analyse_hyres(a.hyres, a.target_frames)
        if acc:
            df = tabulate(acc, label, nrep, nfr)
            allrows.append(df)
            print(f'   {nrep} rep, {nfr} frames, {len(df)} domain pairs')
            heatmaps(df, label, f'18_domain_contacts_{(a.label or "hyres").replace(" ", "_")}')
    else:
        sets = a.sets or [s for s in SETS if (HERE / a.root / s).is_dir()]
        if a.sets is None and a.root != 'binding':
            sets = sorted(p.name for p in (HERE / a.root).iterdir() if p.is_dir())
        for s in sets:
            if not (HERE / a.root / s).is_dir():
                continue
            label = a.label or f'{a.root}/{s}'
            print(f'== {label}')
            if a.per_rep:
                for rname, racc, nfr in analyse_calvados(s, a.root,
                                                         a.target_frames, True):
                    lab = f'{label} {rname}'
                    allrows.append(tabulate(racc, lab, 1, nfr))
                    print(f'   {rname}: {nfr} frames')
                continue
            acc, nrep, nfr = analyse_calvados(s, a.root, a.target_frames)
            if not acc:
                print('   no usable replicates'); continue
            df = tabulate(acc, label, nrep, nfr)
            allrows.append(df)
            print(f'   {nrep} rep, {nfr} frames, {len(df)} domain pairs')
            heatmaps(df, label, f'18_domain_contacts_{a.root}_{s}')

    if not allrows:
        print('nothing analysed'); return
    out = pd.concat(allrows, ignore_index=True)
    dest = OUT / a.out
    if dest.exists():
        prev = pd.read_csv(dest)
        prev = prev[~prev.model.isin(set(out.model))]
        out = pd.concat([prev, out], ignore_index=True)
    out.to_csv(dest, index=False)
    print('\nwrote', dest, f'({len(out)} rows)')




# --------------------------------------------------------------------------
# Cross-model comparison, run separately once every model is in the CSV:
#     python figure_domain_contacts.py --compare
# --------------------------------------------------------------------------
COMPARE = [('binding/p9_dtx3l',                        'CALVADOS\nunrestrained, separated'),
           ('binding/p9_dtx3l_docked',                 'CALVADOS\nunrestrained, docked'),
           ('binding_go/p9_dtx3l_xl_h20',              'CALVADOS\n+3 XL harmonic k=20'),
           ('HyRes matched separated',                 'HyRes matched\nseparated'),
           ('HyRes matched docked',                    'HyRes matched\ndocked')]


def compare(csv='domain_contact_maps.csv'):
    d = pd.read_csv(OUT / csv)
    d = d[(d.kind == 'inter') & (d.chain_i == 'parp9') & (d.chain_j == 'dtx3l')]
    models = [(m, lab) for m, lab in COMPARE if m in set(d.model)]
    rows = DOMAIN_NAMES['parp9']; cols = DOMAIN_NAMES['dtx3l']
    n = len(models)

    fig, axes = plt.subplots(2, n, figsize=(3.55 * n + 1.4, 8.4), squeeze=False)
    specs = [('pct_contact_bound', 'Purples', 0, 100,
              'Where the two chains actually touch, by domain\n'
              'contact occupancy — % of associated frames within 1.0 nm'),
             ('mean_min_bound_nm', 'Greens_r', 0, 10,
              'closest approach — mean min Cα–Cα over associated frames (nm)')]
    for r, (col, cmap, vmin, vmax, rowtitle) in enumerate(specs):
        for c, (m, lab) in enumerate(models):
            sub = d[d.model == m]
            M = (sub.pivot(index='domain_i', columns='domain_j', values=col)
                 .reindex(index=rows, columns=cols).to_numpy(float))
            ax = axes[r, c]
            im = ax.imshow(M, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
            ax.set_xticks(range(len(cols)))
            ax.set_xticklabels(cols if r == 1 else [], rotation=45, ha='right',
                               fontsize=8.5)
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels(rows if c == 0 else [], fontsize=8.5)
            ax.tick_params(colors=INK['secondary'], length=0)
            for y in range(M.shape[0]):
                for x in range(M.shape[1]):
                    v = M[y, x]
                    if np.isnan(v):
                        continue
                    frac = (v - vmin) / (vmax - vmin)
                    dark = frac > .55 if r == 0 else frac < .45
                    ax.text(x, y, f'{v:.0f}' if r == 0 else f'{v:.1f}',
                            ha='center', va='center', fontsize=6.8,
                            color='white' if dark else INK['primary'])
            if r == 0:
                a = sub.n_bound_frames.iloc[0] / sub.n_frames.iloc[0] * 100
                ax.set_title(f'{lab}\nassociated in {a:.0f}% of frames',
                             fontsize=9.5, color=INK['primary'],
                             fontweight='semibold', pad=6)
            if c == 0:
                ax.set_ylabel('PARP9', fontsize=9.5, color=INK['secondary'])
        fig.colorbar(im, ax=axes[r, :].tolist(), fraction=.018, pad=.012)
        axes[r, 0].text(0, 1.33 if r == 0 else 1.045, rowtitle,
                        transform=axes[r, 0].transAxes, ha='left', va='bottom',
                        fontsize=11, color=INK['primary'], fontweight='semibold')
    axes[1, 0].set_xlabel('DTX3L', fontsize=9.5, color=INK['secondary'])
    for ext in ('png', 'svg'):
        fig.savefig(FIG / f'19_domain_contacts_compare.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print('wrote', FIG / '19_domain_contacts_compare.png')


if __name__ == '__main__':
    main()
