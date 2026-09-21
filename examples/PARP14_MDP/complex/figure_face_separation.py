#!/usr/bin/env python3
"""Where, and how often, the chains separate -- resolved by domain and by face.

Method follows sim_analysis/figure_sasa_faces.py and figure_face_contacts.py:

DOMAIN BOUNDARIES come from the RESTRAINED residue blocks (binding/<set>/input/
domains.yaml), not from full unit ranges. Those blocks are what the simulation
actually holds rigid; a flexible linker brushing past the partner chain is not a
domain interface. Residues outside every block are linker and are excluded.

FACES are the three-way split along the domain's own axis, normalised by the
domain's 95th-percentile half-extent so the rim band means the same thing for
every domain:

    |norm| <= 0.35   rim      equatorial band; the sign of a near-zero
                              projection carries no geometric information
     norm  >  0.35   active / front
     norm  < -0.35   back

  * Domains WITH a catalytic site (PARP14 MD1, MD2, MD3, WWE, ART, from
    parp14/input/active_sites.yaml) use the domain-COM -> catalytic-COM axis,
    so the poles are 'active' and 'back'.
  * Domains WITHOUT one -- every PARP9 and DTX3L domain, and PARP14's RRM/KH
    units -- have no catalytic reference, so the axis is the domain's first
    principal component and the poles are the geometric 'front' and 'back'.
    PARP9's macrodomains are ADP-ribose binders by homology, but no residue-level
    site is defined for them in this project, and inventing one would put a
    label on the plot that no data supports.

Faces are computed once from the first frame: the domains are harmonically
restrained, so their internal geometry is fixed for the run.

CONTACTS are inter-CHAIN CA-CA within 1.0 nm under the minimum image convention,
the same cutoff figure_md_distances.py and analyze_binding.py use.

ENRICHMENT, not raw counts. The rim holds ~35% of a domain's residues so it wins
on size alone. Every face comparison is therefore
    enrichment(face) = (contacts on face / contacts on domain)
                     / (residues in face / residues in domain)
1.0 = exactly what the face's size predicts; >1 = preferred.

Usage:
    python figure_face_separation.py                 # every set with trajectories
    python figure_face_separation.py p9_dtx3l
    python figure_face_separation.py --target-frames 500
"""
import os, sys, warnings
from pathlib import Path
from collections import defaultdict
from argparse import ArgumentParser
import numpy as np, yaml, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
warnings.filterwarnings('ignore')

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'binding'
OUT  = HERE / 'analysis'
FIG  = OUT / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(HERE))
from analyze_binding import SETS, sysname, topology

ACTIVE_SITES = {
    'parp14': Path('/home/sbali/CALVADOS/parp14/input/active_sites.yaml'),
    # PARP9 + DTX3L sites derived for this campaign; see that file's header for
    # the PDB entries and the structure/homology/mutagenesis provenance of each.
    'parp9':  HERE / 'active_sites_parp9_dtx3l.yaml',
    'dtx3l':  HERE / 'active_sites_parp9_dtx3l.yaml',
}
CONTACT_NM  = 1.0
RIM_FRACTION = 0.35
NS_PER_RAW_FRAME = 0.05   # wfreq 5000 x 0.01 ps
SKIP_NS = 25.0            # equilibration discarded
SKIP_FRAMES = int(SKIP_NS / NS_PER_RAW_FRAME)   # RAW dcd frames, not analysed samples
DEFAULT_TARGET_FRAMES = 2000   # ~0.25 ns/sample; episodes here are ~1-2 ns

# canonical domain names, in restraint-block order, per chain
DOMAIN_NAMES = {
 'parp14': ['RRM1','RRM2','RRM3','KH1','KH2','KH3','KH4','KH5','KH6','1/2KH7a',
            'MD1','MD2','MD3','1/2KH7b','KH8','WWE','ART'],
 'parp9':  ['1/2KH1a','MD1','MD2','1/2KH1b','KH2','ART'],
 'dtx3l':  ['RRM','KH1','KH2','KH3','KH4','KH5','RING','DTC'],
}
FACE_COLOR = {'active':'#27ae60','front':'#2e86c1','rim':'#f0b323','back':'#7f8c8d'}


def load_active_sites():
    """chain -> {domain name: pocket resids}, in each chain's own numbering.

    Prefers pocket_residues over catalytic_residues: the face axis is the
    domain-COM -> site-COM vector, and a pocket centroid locates the ligand far
    better than three scattered catalytic side chains.
    """
    out = {}
    for chain, path in ACTIVE_SITES.items():
        if not path.is_file():
            continue
        d = yaml.safe_load(open(path))
        sub = d.get(chain, d)            # parp9/dtx3l file is nested by chain
        out[chain] = {k: (v.get('pocket_residues') or v.get('catalytic_residues'))
                      for k, v in sub.items()
                      if isinstance(v, dict)
                      and (v.get('pocket_residues') or v.get('catalytic_residues'))}
    return out


def load_domains(s, chains):
    """Per chain INSTANCE -> [(lo, hi), ...] restrained blocks.

    domains.yaml is keyed by component name, which is not always the chain name:
    the homodimer-docked sets were built as two separate components, so
    dtx3l_homo_docked is keyed dtx3l_1/dtx3l_2 while SETS says ['dtx3l','dtx3l'].
    Resolve per instance, falling back to <chain>_<n> and then to any key with
    the chain name as prefix.
    """
    raw = yaml.safe_load(open(ROOT / s / 'input' / 'domains.yaml'))
    out = []
    for i, c in enumerate(chains):
        if c in raw:
            out.append(raw[c]); continue
        k = f'{c}_{i + 1}'
        if k in raw:
            out.append(raw[k]); continue
        cand = [v for kk, v in raw.items() if kk.startswith(c)]
        if not cand:
            raise KeyError(f'{s}: no domains.yaml entry for chain {c} '
                           f'(keys: {sorted(raw)})')
        out.append(cand[min(i, len(cand) - 1)])
    return out


def classify_faces(coords, resids, blocks, chain, sites):
    """resid -> (domain_name, face). coords in nm, one frame, this chain only."""
    names = DOMAIN_NAMES.get(chain, [f'D{i+1}' for i in range(len(blocks))])
    out = {}
    for i, (lo, hi) in enumerate(blocks):
        dname = names[i] if i < len(names) else f'D{i+1}'
        sel = (resids >= lo) & (resids <= hi)
        if sel.sum() < 4:
            continue
        X = coords[sel]; rid = resids[sel]
        com = X.mean(axis=0)

        cat = sites.get(chain, {}).get(dname)
        cat = [r for r in (cat or []) if lo <= r <= hi]
        if cat:
            m = np.isin(rid, cat)
            if m.sum():
                axis = X[m].mean(axis=0) - com
                pole_hi, pole_lo = 'active', 'back'
            else:
                cat = None
        if not cat:
            # No catalytic reference -> first principal component. The SVD sign is
            # arbitrary, so it must be fixed deterministically.
            #
            # It used to be fixed to point away from the PARENT CHAIN's centre of
            # mass ('outward'). That was wrong: the chain COM is not restrained, it
            # reorients continuously as the linkers flex, so the same residue got a
            # different label at different frames -- measured at 15-19% of residues
            # relabelled between frames, with whole domains (PARP9 1/2KH1a, DTX3L
            # RRM/KH1-3) flipping front<->back ~70% of the time. Those columns
            # carried no information.
            #
            # The sign is now anchored to the domain's OWN sequence: the axis points
            # from the N-terminal half's COM to the C-terminal half's COM. Both
            # endpoints sit inside the same harmonically restrained block, so the
            # sign is rigid for the whole run and identical across replicates.
            #
            # NOTE: 'front'/'back' is therefore a STABLE GEOMETRIC label, not an
            # outward/inward one. Do not read solvent exposure into it.
            U, S, Vt = np.linalg.svd(X - com, full_matrices=False)
            axis = Vt[0]
            mid = len(X) // 2
            seq_dir = X[mid:].mean(axis=0) - X[:mid].mean(axis=0)
            if np.dot(axis, seq_dir) < 0:
                axis = -axis
            pole_hi, pole_lo = 'front', 'back'

        n = np.linalg.norm(axis)
        if n == 0:
            continue
        axis = axis / n
        proj = (X - com) @ axis
        half_extent = float(np.percentile(np.abs(proj), 95)) or 1.0
        # Exposure tag. A buried residue cannot make an inter-chain contact, so
        # counting it in the enrichment denominator dilutes whichever face holds
        # more buried residues. The rim is systematically the most buried band
        # (measured intra-chain coordination: rim 18.0, active 17.5, back 15.6,
        # front 13.9), which alone pushed every rim row below 1.0. Exposure here
        # is intra-chain CA coordination within CONTACT_NM; a residue counts as
        # exposed if it is below its own domain's median.
        nb = (np.linalg.norm(coords[:, None, :] - X[None, :, :], axis=2)
              < CONTACT_NM).sum(axis=0)
        med = float(np.median(nb))
        for r, pr, cn in zip(rid, proj / half_extent, nb):
            out[int(r)] = (dname,
                           pole_hi if pr > RIM_FRACTION else
                           pole_lo if pr < -RIM_FRACTION else 'rim',
                           bool(cn <= med))
    return out


def analyse(s, target_frames=DEFAULT_TARGET_FRAMES):
    chains = SETS[s]
    blocks_by_chain = load_domains(s, chains)
    sites = load_active_sites()
    sysn = sysname(s)

    # chain bead offsets, from analyze_convergence's canonical lengths
    import analyze_convergence as ac
    offs, o = [], 0
    for idx, c in enumerate(chains):
        offs.append((idx, c, o, o + ac.LENGTHS[c])); o += ac.LENGTHS[c]

    step_set = None
    per_res  = defaultdict(float)                 # (cidx, resid) -> contacts/frame
    dom_series = defaultdict(list)                # (cidx, domain) -> per-frame bound flag
    face_map = {}
    nrep = nframes_tot = 0

    for rep in sorted((ROOT / s).glob('rep-*'), key=lambda p: int(p.name.split('-')[1])):
        dcd = rep / f'{sysn}.dcd'
        top = topology(rep)
        if top is None or not dcd.exists():
            continue
        cfg = yaml.safe_load(open(rep / 'config.yaml'))
        box = np.array(cfg['box'], float)
        boxv = np.array([box[0]*10, box[1]*10, box[2]*10, 90., 90., 90.], np.float32)
        u = mda.Universe(str(top), str(dcd))
        if u.atoms.n_atoms != o:
            continue

        usable = max(0, len(u.trajectory) - SKIP_FRAMES)
        if usable < 10:
            continue
        # One stride for the whole set. It used to be recomputed per replicate and
        # only the last value survived into the return dict, so a set containing
        # one short replicate converted pooled episode lengths with the wrong
        # scalar (up to 4x off) and pooled samples of unequal time weight.
        if step_set is None:
            step_set = max(1, usable // target_frames) if target_frames else 1
        step = step_set
        ns_per_sample = cfg['wfreq'] * NS_PER_RAW_FRAME / 5000.0 * step

        if not face_map:                           # faces from the first frame
            u.trajectory[0]
            P0 = u.atoms.positions / 10.0
            resids = u.atoms.resids
            for cidx, c, a, b in offs:
                face_map[cidx] = classify_faces(P0[a:b], resids[a:b],
                                                blocks_by_chain[cidx], c, sites)

        rep_series = defaultdict(list)
        nfr = 0
        for ts in u.trajectory[SKIP_FRAMES::step]:
            P = u.atoms.positions.astype(np.float32)
            touched = set()
            for i in range(len(offs)):
                for j in range(i + 1, len(offs)):
                    _, _, ai, bi = offs[i]; _, _, aj, bj = offs[j]
                    pairs = capped_distance(P[ai:bi], P[aj:bj], CONTACT_NM * 10.0,
                                            box=boxv, return_distances=False)
                    if not len(pairs):
                        continue
                    for k, cidx, base in ((0, offs[i][0], ai), (1, offs[j][0], aj)):
                        rid = u.atoms.resids[base + pairs[:, k]]
                        for r in rid:
                            per_res[(cidx, int(r))] += 1
                            fm = face_map.get(cidx, {}).get(int(r))
                            if fm:
                                touched.add((cidx, fm[0]))
            for cidx in face_map:
                for dname in {v[0] for v in face_map[cidx].values()}:
                    rep_series[(cidx, dname)].append(1 if (cidx, dname) in touched else 0)
            nfr += 1
        for k, v in rep_series.items():
            dom_series[k].append(np.array(v, dtype=np.int8))
        nrep += 1; nframes_tot += nfr

    if not nrep:
        return None
    for k in per_res:
        per_res[k] /= nframes_tot
    return dict(set=s, chains=chains, offs=offs, per_res=per_res,
                dom_series=dom_series, face_map=face_map, nrep=nrep,
                ns_per_sample=ns_per_sample, step=step)


def episodes(b):
    d = np.diff(np.concatenate(([0], np.asarray(b), [0])))
    return np.where(d == -1)[0] - np.where(d == 1)[0]


def tabulate(res, ns_per_sample):
    rows = []
    for cidx, c, a, b in res['offs']:
        fm = res['face_map'].get(cidx, {})
        doms = sorted({v[0] for v in fm.values()},
                      key=lambda n: DOMAIN_NAMES.get(c, []).index(n)
                      if n in DOMAIN_NAMES.get(c, []) else 99)
        for dname in doms:
            resids_d = [r for r, v in fm.items() if v[0] == dname]
            ctot = sum(res['per_res'].get((cidx, r), 0.0) for r in resids_d)
            series = res['dom_series'].get((cidx, dname), [])
            allb = np.concatenate(series) if series else np.array([0])
            eps = np.concatenate([episodes(x) for x in series]) if series else np.array([])
            exposed_d = [r for r in resids_d if fm[r][2]]
            for face in ('active', 'front', 'rim', 'back'):
                rf = [r for r in resids_d if fm[r][1] == face]
                if not rf:
                    continue
                cf = sum(res['per_res'].get((cidx, r), 0.0) for r in rf)
                # Primary metric: exposed residues only in BOTH numerator and
                # denominator -- a buried residue cannot contact the partner chain,
                # so leaving it in the denominator measures burial, not preference.
                ef = [r for r in rf if fm[r][2]]
                cfe = sum(res['per_res'].get((cidx, r), 0.0) for r in ef)
                cte = sum(res['per_res'].get((cidx, r), 0.0) for r in exposed_d)
                enr = ((cfe / cte) / (len(ef) / len(exposed_d))) \
                    if (cte > 0 and ef and exposed_d) else np.nan
                # Kept for comparison: the uncorrected, burial-confounded version.
                enr_raw = ((cf / ctot) / (len(rf) / len(resids_d))) if ctot > 0 else np.nan
                rows.append(dict(
                    set=res['set'], chain=f'{c}#{cidx+1}', domain=dname, face=face,
                    n_res=len(rf), n_res_exposed=len(ef), contacts_per_frame=round(cf, 4),
                    frac_domain_contacts=round(cf / ctot, 4) if ctot > 0 else np.nan,
                    enrichment=round(enr, 3) if enr == enr else np.nan,
                    enrichment_uncorrected=round(enr_raw, 3) if enr_raw == enr_raw else np.nan,
                    domain_contact_freq=round(float(allb.mean()), 4),
                    n_separations=int(len(eps)),
                    tau_mean_ns=round(float(eps.mean() * ns_per_sample), 2) if len(eps) else 0.0))
    return pd.DataFrame(rows)


def plots(df, res, ns_per_sample):
    s = res['set']
    dom = (df.groupby(['chain', 'domain'], sort=False)
             .agg(freq=('domain_contact_freq', 'first'),
                  seps=('n_separations', 'first'),
                  tau=('tau_mean_ns', 'first'),
                  contacts=('contacts_per_frame', 'sum')).reset_index())
    dom = dom[dom.contacts > 0] if (dom.contacts > 0).any() else dom
    lab = dom.chain.str.replace('#', ' #') + ' · ' + dom.domain

    # ---- Fig 1: where, and how often it lets go ----
    fig, axes = plt.subplots(1, 3, figsize=(13, max(3.2, 0.26 * len(dom) + 1.4)), sharey=True)
    y = np.arange(len(dom))
    axes[0].barh(y, dom.freq, color='#2e86c1')
    axes[0].set_xlabel('fraction of frames in inter-chain contact')
    axes[0].set_title('WHERE the chains touch', fontsize=9, loc='left')
    axes[1].barh(y, dom.seps, color='#c0392b')
    axes[1].set_xlabel('separation events (contact → no contact)')
    axes[1].set_title('HOW OFTEN it lets go', fontsize=9, loc='left')
    axes[2].barh(y, dom.tau, color='#27ae60')
    axes[2].set_xlabel(f'mean contact lifetime (ns)')
    axes[2].set_title('HOW LONG it holds', fontsize=9, loc='left')
    axes[0].set_yticks(y); axes[0].set_yticklabels(lab, fontsize=6); axes[0].invert_yaxis()
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle(f'{s} — inter-chain contact by restrained domain  '
                 f'({res["nrep"]} replicates, {CONTACT_NM} nm CA-CA)',
                 fontsize=10, x=0.01, ha='left')
    fig.tight_layout(); fig.savefig(FIG / f'07_{s}_domain_separation.png', dpi=150)
    plt.close(fig)

    # ---- Fig 2: face enrichment ----
    d2 = df[df.contacts_per_frame > 0].copy()
    if len(d2):
        d2['lab'] = d2.chain.str.replace('#', ' #') + ' · ' + d2.domain
        order = [l for l in lab if l in set(d2.lab)]
        faces = [f for f in ('active', 'front', 'rim', 'back') if f in set(d2.face)]
        fig, ax = plt.subplots(figsize=(8.5, max(3.0, 0.30 * len(order) + 1.4)))
        h = 0.8 / len(faces)
        for k, face in enumerate(faces):
            sub = d2[d2.face == face].set_index('lab').reindex(order)
            ax.barh(np.arange(len(order)) + k * h - 0.4 + h / 2, sub.enrichment.fillna(0),
                    height=h, label=face, color=FACE_COLOR[face])
        ax.axvline(1.0, color='k', ls='--', lw=0.9)
        ax.set_yticks(np.arange(len(order))); ax.set_yticklabels(order, fontsize=6)
        ax.invert_yaxis(); ax.set_xlabel('contact enrichment (1.0 = as expected from face size)')
        ax.legend(fontsize=7, frameon=False, ncol=len(faces))
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_title(f'{s} — which face of each domain carries the inter-chain contacts\n'
                     "poles are active/back where a catalytic site is defined, else front/back (PC1)",
                     fontsize=9, loc='left')
        fig.tight_layout(); fig.savefig(FIG / f'08_{s}_face_enrichment.png', dpi=150)
        plt.close(fig)

    # ---- Fig 3: contact/separation timeline, replicate 1 ----
    keys = [(c, d) for c, d in zip(dom.chain, dom.domain)]
    mats = []
    for chain, dname in keys:
        cidx = int(chain.split('#')[1]) - 1
        ser = res['dom_series'].get((cidx, dname), [])
        mats.append(ser[0] if ser else np.zeros(1, np.int8))
    if mats:
        n = min(len(m) for m in mats)
        M = np.array([m[:n] for m in mats])
        fig, ax = plt.subplots(figsize=(11, max(2.6, 0.22 * len(keys) + 1.2)))
        ax.imshow(M, aspect='auto', cmap='Blues', interpolation='nearest',
                  extent=[SKIP_FRAMES * NS_PER_RAW_FRAME,
                          SKIP_FRAMES * NS_PER_RAW_FRAME + n * ns_per_sample,
                          len(keys) - 0.5, -0.5])
        ax.set_yticks(np.arange(len(keys))); ax.set_yticklabels(lab, fontsize=6)
        ax.set_xlabel('time (ns)')
        ax.set_title(f'{s} — replicate 1: when each domain is in inter-chain contact '
                     '(dark = touching). Gaps are separations.', fontsize=9, loc='left')
        fig.tight_layout(); fig.savefig(FIG / f'09_{s}_separation_timeline.png', dpi=150)
        plt.close(fig)


def main():
    ap = ArgumentParser()
    ap.add_argument('sets', nargs='*')
    ap.add_argument('--target-frames', type=int, default=DEFAULT_TARGET_FRAMES)
    a = ap.parse_args()
    todo = a.sets or [s for s in SETS if (ROOT / s).is_dir()
                      and list((ROOT / s).glob(f'rep-*/{sysname(s)}.dcd'))]
    allrows = []
    for s in todo:
        print(f'== {s}', flush=True)
        res = analyse(s, a.target_frames)
        if res is None:
            print('   no usable replicates'); continue
        ns_per_sample = res['ns_per_sample']
        df = tabulate(res, ns_per_sample)
        allrows.append(df)
        plots(df, res, ns_per_sample)
        print(f'   {res["nrep"]} reps, {ns_per_sample:.2f} ns/sample, '
              f'{len(df)} domain-face rows', flush=True)
    if allrows:
        out = pd.concat(allrows, ignore_index=True)
        # merge with earlier sets so the table grows as the campaign finishes
        f = OUT / 'domain_face_contacts.csv'
        if f.exists():
            try:
                prev = pd.read_csv(f)
                out = pd.concat([prev[~prev.set.isin(out.set.unique())], out],
                                ignore_index=True)
            except Exception as e:
                print(f'  ! could not merge existing csv: {e}')
        out.to_csv(f, index=False)
        print(f'\nsaved {f}  ({len(out)} rows, {out.set.nunique()} sets)')


if __name__ == '__main__':
    main()
