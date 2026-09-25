#!/usr/bin/env python
"""
Figures for the CALVADOS campaign summary page.

Every figure follows the PARP14 Figures design system: locked domain colours,
Arial, flat fills, 1 pt lines, no top/right spines, no gradients.

For each question the page pairs two figures:
  * a MEASUREMENT panel -- a schematic of what is actually being computed
  * a RESULT panel -- the measured numbers, as a graph

and for each simulation type, a rendering of the simulation itself from real
coordinates.

    python make_campaign_figures.py --out <dir>
"""
import argparse
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

# --- design system tokens ---------------------------------------------------
INK, INK2, MUTED = '#0b0b0b', '#52514e', '#898781'
GRID, BASE, SURF = '#e1e0d9', '#c3c2b7', '#ffffff'
DOM = {'rrm1': '#93700f', 'rrm2': '#b88c13', 'rrm3': '#dda816',
       'kh1-kh6': '#117c60', 'kh7a': '#16a07b', 'khb-kh8': '#1cc497',
       'md1': '#0072b2', 'md2': '#6a3fa0', 'md3': '#b31973',
       'wwe': '#f584cf', 'art': '#c8481a'}
SITE = {'MD1': '#0072b2', 'MD2': '#6a3fa0', 'MD3': '#b31973',
        'WWE': '#f584cf', 'ART': '#c8481a'}
FL, GOOD, CRIT, WARN = '#000000', '#0ca30c', '#d03b3b', '#fab219'
SEQ = ['#cde2fb', '#86b6ef', '#3987e5', '#2a78d6', '#256abf', '#184f95', '#0d366b']

# Single-chain sets excluded from the results.
#
# 'fl' is the only set still on the untrimmed input/domains.yaml AND without the
# KH7a-KHb custom restraints, i.e. the pre-optimisation full-length run.
# fl_optimized supersedes it: per-domain trims from the restraint_tests sweep
# plus custom_restraints: true. Keeping both would show the same construct twice
# under two different restraint schemes.
EXCLUDE_SETS = {'fl'}
FL_REF = 'fl_optimized'          # the full-length reference in every figure

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 10, 'axes.titleweight': 'bold',
    'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 8.5, 'legend.frameon': False,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': INK, 'axes.linewidth': 0.8,
    'xtick.color': INK, 'ytick.color': INK, 'text.color': INK,
    'axes.labelcolor': INK, 'grid.color': GRID, 'grid.linewidth': 0.5,
    'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
})


def save(fig, out, name):
    p = os.path.join(out, name)
    fig.savefig(p, dpi=170, bbox_inches='tight')
    plt.close(fig)
    print(f'  {name}')


# ---------------------------------------------------------------- measurement
def fig_measure_accessibility(out):
    """What solid-angle accessibility actually computes."""
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    rng = np.random.default_rng(3)
    # protein body: a blob of beads, with a pocket notch on the right
    ang = rng.uniform(0, 2 * np.pi, 320)
    rad = 1.0 * np.sqrt(rng.uniform(0, 1, 320))
    bx, by = rad * np.cos(ang), rad * np.sin(ang)
    keep = ~((bx > 0.45) & (np.abs(by) < 0.38))
    ax.scatter(bx[keep], by[keep], s=10, c=MUTED, lw=0, zorder=1)
    px, py = 0.60, 0.0
    ax.add_patch(Circle((px, py), 0.14, fc=SITE['MD1'], ec='none', zorder=3))
    ax.text(px, py - 0.30, 'pocket COM', ha='center', fontsize=8, color=INK)
    # rays: open to the right, blocked to the left
    for a in np.linspace(-np.pi, np.pi, 28, endpoint=False):
        dx, dy = np.cos(a), np.sin(a)
        blocked = dx < 0.15
        ax.annotate('', xy=(px + 1.05 * dx, py + 1.05 * dy), xytext=(px, py),
                    arrowprops=dict(arrowstyle='-', lw=0.6,
                                    color=(CRIT if blocked else GOOD),
                                    alpha=0.85), zorder=2)
    ax.plot([], [], color=GOOD, lw=1, label='escapes to 5 nm')
    ax.plot([], [], color=CRIT, lw=1, label='blocked by own chain')
    ax.legend(loc='upper left', bbox_to_anchor=(-0.02, 1.06))
    ax.set_xlim(-1.4, 1.9); ax.set_ylim(-1.35, 1.35); ax.set_aspect('equal')
    ax.axis('off')
    ax.text(0.5, -1.28, 'SAA  =  escaping rays / 200', ha='center',
            fontsize=9, color=INK, style='italic')
    save(fig, out, 'measure_accessibility.png')


def fig_measure_dominance(out):
    """The writer-vs-eraser dominance score."""
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    for x, lab, col, frac in ((-0.75, 'MD1\neraser', SITE['MD1'], 0.35),
                              (0.75, 'ART\nwriter', SITE['ART'], 0.72)):
        ax.add_patch(Circle((x, 0.35), 0.34, fc=col, ec='none'))
        ax.text(x, 0.35, lab, ha='center', va='center', fontsize=8.5,
                color='white', fontweight='bold')
        ax.add_patch(Rectangle((x - 0.30, -0.62), 0.60, 0.42,
                               fc='#f0efec', ec='none'))
        ax.add_patch(Rectangle((x - 0.30, -0.62), 0.60 * frac, 0.42,
                               fc=col, ec='none'))
        ax.text(x, -0.78, f'SAA {frac:.2f}', ha='center', fontsize=8.5, color=INK)
    ax.annotate('', xy=(0.36, 0.35), xytext=(-0.36, 0.35),
                arrowprops=dict(arrowstyle='<->', lw=1, color=INK2))
    ax.text(0, 1.02, 'D  =  (SAA$_{ART}$ − SAA$_{MD1}$) / (SAA$_{ART}$ + SAA$_{MD1}$)',
            ha='center', fontsize=9.5, color=INK)
    ax.text(-0.75, -1.08, 'D < 0  eraser-leaning', ha='center', fontsize=8, color=INK2)
    ax.text(0.75, -1.08, 'D > 0  writer-leaning', ha='center', fontsize=8, color=INK2)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.3, 1.25); ax.set_aspect('equal'); ax.axis('off')
    save(fig, out, 'measure_dominance.png')


def fig_measure_interdomain(out):
    """Inter-site COM-COM distance."""
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    pts = {'MD1': (-0.9, 0.5), 'MD2': (0.0, 0.85), 'MD3': (0.85, 0.35),
           'ART': (0.35, -0.75)}
    for a in ('MD1', 'MD2', 'MD3'):
        for b in ('MD2', 'MD3', 'ART'):
            if a == b:
                continue
            ax.plot(*zip(pts[a], pts[b]), lw=0.6, color=MUTED, zorder=1, ls=(0, (3, 2)))
    ax.annotate('', xy=pts['ART'], xytext=pts['MD1'],
                arrowprops=dict(arrowstyle='<->', lw=1.2, color=INK))
    mx = (pts['MD1'][0] + pts['ART'][0]) / 2, (pts['MD1'][1] + pts['ART'][1]) / 2
    ax.text(mx[0] - 0.28, mx[1], 'd', fontsize=11, style='italic', color=INK)
    for n, (x, y) in pts.items():
        ax.add_patch(Circle((x, y), 0.26, fc=SITE[n], ec='none', zorder=3))
        ax.text(x, y, n, ha='center', va='center', fontsize=8, color='white',
                fontweight='bold', zorder=4)
    ax.text(0, -1.22, 'distance between active-site centres of mass,\npooled over 25 replicates',
            ha='center', fontsize=8.5, color=INK2)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.45, 1.3); ax.set_aspect('equal'); ax.axis('off')
    save(fig, out, 'measure_interdomain.png')


def fig_measure_slab(out):
    """What the slab run measures: a density profile along z."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 2.7),
                                 gridspec_kw={'width_ratios': [1, 1.25]})
    rng = np.random.default_rng(11)
    z = np.concatenate([rng.normal(0, 0.10, 700), rng.uniform(-1, 1, 40)])
    x = rng.uniform(-0.35, 0.35, len(z))
    a1.scatter(x, z, s=3, c=SEQ[4], lw=0)
    a1.axhline(0.22, color=INK, lw=0.8, ls=(0, (4, 2)))
    a1.axhline(-0.22, color=INK, lw=0.8, ls=(0, (4, 2)))
    a1.text(0.40, 0.0, 'dense', fontsize=8.5, color=INK, rotation=90, va='center')
    a1.text(0.40, 0.62, 'dilute', fontsize=8.5, color=INK2, rotation=90, va='center')
    a1.set_xlim(-0.55, 0.62); a1.set_ylim(-1, 1)
    a1.set_ylabel('z (box long axis)'); a1.set_xticks([])
    a1.set_yticks([]); a1.spines['left'].set_visible(False); a1.spines['bottom'].set_visible(False)
    zz = np.linspace(-1, 1, 400)
    prof = 1 / (1 + np.exp((np.abs(zz) - 0.22) / 0.035)) * 0.96 + 0.02
    a2.plot(prof, zz, color=SEQ[5], lw=1.2)
    a2.axvline(0.02, color=MUTED, lw=0.8, ls=(0, (3, 2)))
    a2.annotate('c$_{sat}$', xy=(0.02, 0.72), xytext=(0.30, 0.80), fontsize=10,
                color=INK, arrowprops=dict(arrowstyle='->', lw=0.8, color=INK2))
    a2.annotate('c$_{dense}$', xy=(0.96, 0.0), xytext=(0.52, 0.22), fontsize=10,
                color=INK, arrowprops=dict(arrowstyle='->', lw=0.8, color=INK2))
    a2.set_xlabel('concentration'); a2.set_ylabel('z')
    a2.set_yticks([]); a2.set_xticks([])
    fig.text(0.5, -0.06, 'c$_{sat}$ is the dilute-phase plateau — it only exists if the slab survives',
             ha='center', fontsize=8.5, color=INK2)
    fig.tight_layout()
    save(fig, out, 'measure_slab.png')


# -------------------------------------------------------------------- results
def _load(npz):
    return np.load(os.path.join(ROOT, 'data', npz), allow_pickle=True)


def fig_result_accessibility(out):
    d = _load('accessibility_stats.npz')
    sets = ['fl_optimized', 'core', 'norrm', 'noart', 'md_full', 'mka_full',
            'core_full_go']
    sets = [s for s in sets if s not in EXCLUDE_SETS
            and any(k.startswith(s + '_') for k in d.keys())]
    sites = ['MD1', 'MD2', 'MD3', 'WWE', 'ART']
    fig, ax = plt.subplots(figsize=(6.4, 3.1))
    w = 0.15
    for i, st in enumerate(sites):
        xs, ys, es = [], [], []
        for j, s in enumerate(sets):
            k = f'{s}_{st}_saa'
            if k in d:
                v = np.asarray(d[k], dtype=float)
                v = v[np.isfinite(v)]
                if len(v):
                    xs.append(j + (i - 2) * w); ys.append(v.mean()); es.append(v.std())
        ax.bar(xs, ys, width=w * 0.92, color=SITE[st], label=st, lw=0)
        ax.errorbar(xs, ys, yerr=es, fmt='none', ecolor=INK, elinewidth=0.7, capsize=1.6)
    ax.set_xticks(range(len(sets)))
    ax.set_xticklabels([s.replace('_', '\n') for s in sets])
    ax.set_ylabel('SAA  (fraction of rays escaping)')
    ax.grid(axis='y', lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(ncol=5, loc='upper center', bbox_to_anchor=(0.5, 1.16))
    save(fig, out, 'result_accessibility.png')


def fig_result_dominance(out):
    d = _load('accessibility_stats.npz')
    rows = []
    for k in d.keys():
        if not k.endswith('_ART_saa'):
            continue
        s = k[:-len('_ART_saa')]
        if s in EXCLUDE_SETS:
            continue
        km = f'{s}_MD1_saa'
        if km not in d:
            continue
        a = np.asarray(d[k], dtype=float); m = np.asarray(d[km], dtype=float)
        ok = np.isfinite(a) & np.isfinite(m) & ((a + m) > 0)
        if ok.sum() < 2:
            continue
        D = (a[ok] - m[ok]) / (a[ok] + m[ok])
        rows.append((s, D.mean(), D.std() / max(np.sqrt(ok.sum()), 1)))
    rows.sort(key=lambda r: r[1])
    fig, ax = plt.subplots(figsize=(6.4, max(2.6, 0.30 * len(rows) + 1.0)))
    y = np.arange(len(rows))
    cols = [SITE['ART'] if r[1] > 0 else SITE['MD1'] for r in rows]
    ax.barh(y, [r[1] for r in rows], xerr=[r[2] for r in rows], color=cols, lw=0,
            error_kw=dict(ecolor=INK, elinewidth=0.7, capsize=1.6), height=0.66)
    ax.axvline(0, color=BASE, lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xlabel('D   ←  eraser-leaning        writer-leaning  →')
    ax.grid(axis='x', lw=0.5); ax.set_axisbelow(True)
    save(fig, out, 'result_dominance.png')


def fig_result_interdomain(out):
    d = _load('active_site_stats.npz')
    pairs = ['MD1-MD2', 'MD2-MD3', 'MD1-MD3', 'MD1-ART', 'MD3-ART', 'MD2-ART']
    sets = ['fl_optimized', 'md', 'md_full', 'core', 'norrm', 'core_full_go']
    sets = [s for s in sets if s not in EXCLUDE_SETS
            and any(k.startswith(s + '_dist_') for k in d.keys())]
    fig, ax = plt.subplots(figsize=(6.4, 3.1))
    for i, s in enumerate(sets):
        xs, ys = [], []
        for j, p in enumerate(pairs):
            k = f'{s}_dist_{p}'
            if k in d:
                v = np.asarray(d[k], dtype=float); v = v[np.isfinite(v)]
                if len(v):
                    xs.append(j); ys.append(v.mean())
        ax.plot(xs, ys, marker='o', ms=4, lw=1, label=s,
                color=FL if s == FL_REF else SEQ[min(i + 1, 6)])
    ax.set_xticks(range(len(pairs))); ax.set_xticklabels(pairs, rotation=20, ha='right')
    ax.set_ylabel('active-site COM–COM distance (nm)')
    ax.grid(axis='y', lw=0.5); ax.set_axisbelow(True)
    ax.legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.18))
    save(fig, out, 'result_interdomain.png')


# ------------------------------------------------------- simulation renderings
def _domain_colors(nres, units):
    import sim_registry as reg
    c = np.array([MUTED] * nres, dtype=object)
    fl2c = reg.build_fl_to_construct_map(units, None)
    for u in units:
        a, b = reg.DOMAIN_UNITS[u]
        for fl in range(a, b + 1):
            i = fl2c(fl)
            if i and 1 <= i <= nres:
                c[i - 1] = DOM.get(u, MUTED)
    return c


def fig_sim_monomer(out):
    import mdtraj as md
    import sim_registry as reg
    d = reg.get_sim_dir(FL_REF, 1, 0)
    dcd = [f for f in os.listdir(d) if f.endswith('.dcd')
           and not f.startswith(('backup', 'equilibration'))][0]
    t = md.load(os.path.join(d, dcd), top=os.path.join(d, 'top.pdb'), stride=2000)
    x = t.xyz[-1]
    x = x - x.mean(0)
    cols = _domain_colors(t.n_atoms, list(reg.DOMAIN_UNITS))
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    ax.plot(x[:, 0], x[:, 1], lw=0.35, color=MUTED, zorder=1)
    ax.scatter(x[:, 0], x[:, 1], s=7, c=list(cols), lw=0, zorder=2)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title('Full-length PARP14, optimised restraints (1801 beads)', fontsize=9.5)
    hs = [plt.Line2D([], [], marker='o', ls='', ms=5, color=DOM[u],
                     label=u.upper()) for u in ('rrm1', 'kh1-kh6', 'md1', 'md2',
                                                'md3', 'wwe', 'art')]
    ax.legend(handles=hs, ncol=4, loc='lower center', bbox_to_anchor=(0.5, -0.13),
              fontsize=7.5, handletextpad=0.3, columnspacing=0.9)
    save(fig, out, 'sim_monomer.png')


def _slab_frame(ax, d, sysname, title):
    """
    The LAST production frame, drawn against the full box height.

    Loading with a stride >= n_frames silently returns frame 0, i.e. the start
    of production, where every run still looks like the slab slab_eq made. Index
    the final frame explicitly instead -- that is the whole point of the panel.
    """
    import mdtraj as md
    path = os.path.join(d, f'{sysname}.dcd')
    n = len(md.open(path))
    t = md.load_frame(path, n - 1, top=os.path.join(d, 'top.pdb'))
    x = t.xyz[0]
    lz = t.unitcell_lengths[0, 2]
    z = (x[:, 2] - np.median(x[:, 2]) + lz / 2) % lz
    # small opaque marks: the dispersed tail is sparse and vanishes if alpha-blended
    ax.scatter(x[:, 0] - x[:, 0].min(), z, s=0.5, c=SEQ[4], lw=0)
    lo, hi = np.percentile(z, [5, 95])
    ax.axhspan(lo, hi, color=SEQ[0], zorder=0)
    ax.set_ylim(0, lz); ax.set_xlabel('x (nm)'); ax.set_ylabel('z (nm)')
    ax.set_title(title, fontsize=9.5)
    ax.text(0.97, 0.97, f'90% within\n{hi - lo:.0f} nm', transform=ax.transAxes,
            ha='right', va='top', fontsize=8, color=INK)
    return t


def fig_sim_slabs(out):
    import sim_registry as reg
    slab = os.path.join(ROOT, 'slab')
    import yaml
    panels = [('homotypic/fl', 'Disperses — no condensate'),
              ('rna/fl', '+RNA — disperses more slowly'),
              ('homotypic/md_full', 'Arrested — jammed, not liquid')]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.4))
    for ax, (rel, title) in zip(axes, panels):
        d = os.path.join(slab, rel)
        m = yaml.safe_load(open(os.path.join(d, 'slab_meta.yaml')))
        _slab_frame(ax, d, m['sysname'], title)
    fig.tight_layout()
    save(fig, out, 'sim_slabs.png')


def fig_result_slab(out):
    """Dispersal of the compacted slab over 2 us -- the actual slab result."""
    import mdtraj as md
    import yaml
    slab = os.path.join(ROOT, 'slab')
    runs = [('homotypic/fl', 'fl', FL, '-'),
            ('rna/fl', 'fl +RNA', FL, '--'),
            ('homotypic/fl_wwe_full_go', 'fl −ART', '#623d37', '-'),
            ('rna/fl_wwe_full_go', 'fl −ART +RNA', '#623d37', '--'),
            ('homotypic/core_full_go', 'core', '#697e96', '-'),
            ('homotypic/md_full', 'md_full (406 mM)', SITE['MD1'], '-')]
    fig, ax = plt.subplots(figsize=(6.4, 3.3))
    for rel, lab, col, ls in runs:
        d = os.path.join(slab, rel)
        f = os.path.join(d, 'slab_meta.yaml')
        if not os.path.isfile(f):
            continue
        m = yaml.safe_load(open(f))
        p = os.path.join(d, f"{m['sysname']}.dcd")
        if not os.path.isfile(p):
            continue
        t = md.load(p, top=os.path.join(d, 'top.pdb'), stride=200)
        lz = t.unitcell_lengths[0, 2]
        w = []
        for fr in range(t.n_frames):
            z = t.xyz[fr, :, 2].copy()
            z = (z - z.mean() + lz / 2) % lz
            h, _ = np.histogram(z, bins=100, range=(0, lz))
            hs = np.sort(h)[::-1]
            w.append((np.cumsum(hs) / hs.sum() < 0.9).sum() / 100 * lz)
        ts = np.linspace(0, 2.0, len(w))
        ax.plot(ts, w, lw=1.1, color=col, ls=ls, label=lab)
    ax.set_xlabel('production time (µs)')
    ax.set_ylabel('slab width, 90% of beads (nm)')
    ax.grid(lw=0.5); ax.set_axisbelow(True)
    ax.legend(ncol=2, loc='upper left')
    ax.text(0.98, 0.30, 'flat = arrested\nrising = dispersing', transform=ax.transAxes,
            ha='right', va='top', fontsize=8, color=INK2)
    save(fig, out, 'result_slab.png')



def fig_result_ladder(out):
    """Concentration ladder: the three regimes, side by side."""
    import mdtraj as md
    import yaml
    slab = os.path.join(ROOT, 'slab')
    runs = [('homotypic/fl', 'FL  108 mM', FL, '-'),
            ('ladder/fl_c2x', 'FL  216 mM', '#697e96', '-'),
            ('ladder/fl_c4x', 'FL  431 mM', '#2f2f6a', '-'),
            ('ladder/md_full_lowc', 'MD1–MD3  108 mM', SITE['MD1'], '-'),
            ('homotypic/md_full', 'MD1–MD3  406 mM', SITE['MD3'], '--')]
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for rel, lab, col, ls in runs:
        d = os.path.join(slab, rel)
        f = os.path.join(d, 'slab_meta.yaml')
        if not os.path.isfile(f):
            continue
        m = yaml.safe_load(open(f))
        pth = os.path.join(d, f"{m['sysname']}.dcd")
        if not os.path.isfile(pth):
            continue
        t = md.load(pth, top=os.path.join(d, 'top.pdb'), stride=100)
        lz = t.unitcell_lengths[0, 2]
        w = []
        for fr in range(t.n_frames):
            z = t.xyz[fr, :, 2].copy()
            z = (z - z.mean() + lz / 2) % lz
            h, _ = np.histogram(z, bins=200, range=(0, lz))
            hs = np.sort(h)[::-1]
            w.append((np.cumsum(hs) / hs.sum() < 0.9).sum() / 200 * lz)
        ax.plot(np.linspace(0, 2, len(w)), w, lw=1.2, color=col, ls=ls, label=lab)
    ax.set_xlabel('production time (µs)')
    ax.set_ylabel('slab width, 90% of beads (nm)')
    ax.grid(lw=0.5); ax.set_axisbelow(True)
    ax.legend(loc='upper left', ncol=2)
    ax.annotate('full length disperses at every\nconcentration tested',
                xy=(1.6, 150), xytext=(0.75, 185), fontsize=8.5, color=INK2,
                arrowprops=dict(arrowstyle='->', lw=0.8, color=MUTED))
    ax.annotate('macrodomains hold —\nstable dense phase', xy=(1.6, 8),
                xytext=(0.95, 55), fontsize=8.5, color=SITE['MD1'],
                arrowprops=dict(arrowstyle='->', lw=0.8, color=SITE['MD1']))
    save(fig, out, 'result_ladder.png')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    print('writing figures:')
    for fn in (fig_measure_accessibility, fig_measure_dominance,
               fig_measure_interdomain, fig_measure_slab,
               fig_result_accessibility, fig_result_dominance,
               fig_result_interdomain, fig_sim_monomer, fig_sim_slabs,
               fig_result_slab, fig_result_ladder):
        try:
            fn(a.out)
        except Exception as exc:
            print(f'  FAILED {fn.__name__}: {type(exc).__name__}: {exc}')


if __name__ == '__main__':
    main()
