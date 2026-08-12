#!/usr/bin/env python
"""
AF3 ports of CALVADOS analysis modules.

Subcommands:
  sasa       — atomistic relative SASA per active-site residue (FREE: reuses per_combo residues parquet)
  wcn        — weighted coordination number around each active site COM
  saa        — Solid-Angle Accessibility (Fibonacci ray-cast), matches analyze_accessibility.py
  sites      — inter-active-site COM distances, radial position (vs chain COM / Rg),
               per-site contact profile (counts from each atS domain)
  topology   — per-construct Rg, Ree, 11x11 inter-atS-domain COM distance matrix,
               and Fraction Native Contacts (FNC) vs FL_AF3 (Cα<10Å, |i-j|>3)
  all        — runs sasa, wcn, saa, sites, topology in order

Inputs:
  parp14/analysis/structure_analysis/per_combo/<combo>__residues.parquet
  parp14/analysis/structure_analysis/aggregated/all_dompairs.parquet
  parp14/input/active_sites.yaml

Outputs (under parp14/analysis/structure_analysis/aggregated/):
  site_sasa.parquet         (combo, model, site, n_cat, mean_cat_sasa_rel, mean_pocket_sasa_rel, n_cat_exposed)
  site_wcn.parquet          (combo, model, site, wcn_5A, wcn_8A, wcn_12A)
  site_saa.parquet          (combo, model, site, accessible_frac, max_cone_deg, shell_density)
  site_geometry.parquet     (combo, model, site_i, site_j, distance_A) — inter-site distances
  site_radial.parquet       (combo, model, site, dist_to_chainCOM_A, rg_A, ratio)
  site_contacts.parquet     (combo, model, site, atS_domain, n_contacts_sc8_proxy_via_ca10)
  combo_topology.parquet    (combo, model, n_res, rg_A, ree_A)
  combo_dmap.parquet        (combo, model, dom_i, dom_j, com_distance_A)
  combo_fnc.parquet         (combo, model, fnc_ca10, n_native_applicable, n_preserved)

Each subcommand is resumable per-combo via in-place parquet appending where useful;
subcommand outputs are full overwrites. Run one subcommand at a time as you like.
"""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path('/home/sbali/CALVADOS/parp14')
PC = ROOT / 'analysis' / 'structure_analysis' / 'per_combo'
A = ROOT / 'analysis' / 'structure_analysis' / 'aggregated'
A.mkdir(parents=True, exist_ok=True)
ACTIVE_SITES_YAML = ROOT / 'input' / 'active_sites.yaml'

ATS_ORDER = ['RRM1', 'RRM2', 'RRM3', 'KH_N', 'KH_Ca',
             'MD1', 'MD2', 'MD3', 'KH_b', 'WWE', 'ART']

# Match analyze_accessibility.py CG params (CG used nm; here we work in Å, ×10)
PROBE_RADIUS_A = 5.0     # = 0.5 nm
MAX_DIST_A = 50.0        # = 5.0 nm
N_RAYS = 200
SEQ_SEP = 10             # exclude blockers within ±SEQ_SEP residues of any catalytic residue
SHELL_INNER_A = 20.0     # 2.0 nm
SHELL_OUTER_A = 50.0     # 5.0 nm
EXPOSED_RSASA = 0.20     # consistent with analyze_structures.py


# ============================================================================
# Active site loading
# ============================================================================

def load_active_sites():
    cfg = yaml.safe_load(open(ACTIVE_SITES_YAML))
    sites = {}
    for site, info in cfg.items():
        if site == 'annotations':
            continue
        if not isinstance(info, dict) or 'catalytic_residues' not in info:
            continue
        sites[site] = {
            'catalytic': sorted(set(info.get('catalytic_residues', []))),
            'pocket': sorted(set(info.get('pocket_residues', []))),
        }
    return sites


# ============================================================================
# Geometry helpers
# ============================================================================

def fibonacci_sphere(n=N_RAYS):
    golden = (1 + np.sqrt(5)) / 2
    idx = np.arange(n)
    theta = np.arccos(1 - 2 * (idx + 0.5) / n)
    phi = 2 * np.pi * idx / golden
    return np.column_stack([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta),
    ])


_RAY_DIRS = fibonacci_sphere()  # cached


def ray_accessibility(site_com, blocker_pos, probe=PROBE_RADIUS_A, max_dist=MAX_DIST_A):
    if blocker_pos.shape[0] == 0:
        return 1.0, 180.0
    vecs = blocker_pos - site_com  # (M,3)
    proj = _RAY_DIRS @ vecs.T       # (R, M)
    ahead = (proj > 0) & (proj < max_dist)
    # perpendicular distance from each blocker to each ray
    # |v - (v·d)d| = sqrt(|v|^2 - (v·d)^2)
    v_norm2 = np.sum(vecs ** 2, axis=1)             # (M,)
    perp2 = v_norm2[None, :] - proj ** 2            # (R, M)
    np.maximum(perp2, 0, out=perp2)
    blocked = (perp2 < probe ** 2) & ahead          # (R, M)
    ray_blocked = blocked.any(axis=1)
    accessible_frac = float(1 - ray_blocked.sum() / _RAY_DIRS.shape[0])
    # Cone angle: largest cone (around any axis) of unblocked rays. Approx via
    # mean direction of unblocked rays then half-angle to nearest blocked direction.
    if ray_blocked.all():
        cone_deg = 0.0
    elif not ray_blocked.any():
        cone_deg = 180.0
    else:
        free = _RAY_DIRS[~ray_blocked]
        # mean axis
        axis = free.mean(0)
        norm = np.linalg.norm(axis)
        if norm < 1e-6:
            cone_deg = 90.0
        else:
            axis /= norm
            # half-angle to nearest blocked ray
            blocked_dirs = _RAY_DIRS[ray_blocked]
            cos_to_axis = blocked_dirs @ axis
            half_angle_rad = np.arccos(np.clip(cos_to_axis.max(), -1, 1))
            cone_deg = float(np.degrees(half_angle_rad))
    return accessible_frac, cone_deg


def shell_density(site_com, all_pos, r_in=SHELL_INNER_A, r_out=SHELL_OUTER_A):
    d = np.linalg.norm(all_pos - site_com, axis=1)
    n_in = int(((d >= r_in) & (d < r_out)).sum())
    return n_in


def kabsch_rmsd(P, Q):
    if P.shape != Q.shape or P.shape[0] < 3:
        return np.nan
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((Pc @ R.T - Qc) ** 2, axis=1))))


# ============================================================================
# Worker: per-combo SASA / WCN / SAA / sites / topology
# ============================================================================

def _site_present(catalytic_set, fl_set):
    return set(catalytic_set).issubset(fl_set)


def process_combo(parquet_path, sites_def, fl_native_contacts=None,
                  modules=('sasa', 'wcn', 'saa', 'sites', 'topology')):
    """Returns dict module->list[dict rows]. fl_native_contacts: set of (res_i_fl, res_j_fl) tuples."""
    combo = parquet_path.name.replace('__residues.parquet', '')
    df = pd.read_parquet(parquet_path)
    fl_set = set(df.res_fl.unique())
    out = {m: [] for m in modules}

    sites_in_combo = {s: d for s, d in sites_def.items() if _site_present(d['catalytic'], fl_set)}

    for model, g in df.groupby('model'):
        g = g.sort_values('res_construct')
        ca = g[['ca_x', 'ca_y', 'ca_z']].values.astype(np.float32)
        res_fl_arr = g['res_fl'].values.astype(np.int32)
        atS_arr = g['atS_domain'].values
        sasa_rel_arr = g['sasa_rel'].values.astype(np.float32)
        # Map FL residue -> local construct index
        fl_to_idx = {int(r): i for i, r in enumerate(res_fl_arr)}

        # ---------------- SASA ----------------
        if 'sasa' in modules:
            for site, sd in sites_in_combo.items():
                cat_idx = [fl_to_idx[r] for r in sd['catalytic'] if r in fl_to_idx]
                pkt_idx = [fl_to_idx[r] for r in sd['pocket'] if r in fl_to_idx]
                cat_sasa = sasa_rel_arr[cat_idx]
                pkt_sasa = sasa_rel_arr[pkt_idx]
                out['sasa'].append({
                    'combo': combo, 'model': model, 'site': site,
                    'n_cat': len(cat_idx), 'n_pocket': len(pkt_idx),
                    'mean_cat_sasa_rel': float(np.nanmean(cat_sasa)) if len(cat_sasa) else np.nan,
                    'mean_pocket_sasa_rel': float(np.nanmean(pkt_sasa)) if len(pkt_sasa) else np.nan,
                    'n_cat_exposed': int(np.nansum(cat_sasa > EXPOSED_RSASA)),
                    'frac_cat_exposed': float(np.nanmean(cat_sasa > EXPOSED_RSASA)) if len(cat_sasa) else np.nan,
                })

        # Pre-compute COMs and chain-COM for downstream metrics
        chain_com = ca.mean(0)
        rg = float(np.sqrt(np.mean(np.sum((ca - chain_com) ** 2, axis=1))))
        ree = float(np.linalg.norm(ca[-1] - ca[0]))

        # Site COMs (catalytic) and exclusion masks
        site_info = {}
        for site, sd in sites_in_combo.items():
            cat_idx = np.array([fl_to_idx[r] for r in sd['catalytic'] if r in fl_to_idx], dtype=np.int32)
            site_info[site] = {
                'cat_idx': cat_idx,
                'com': ca[cat_idx].mean(0) if len(cat_idx) else None,
                'cat_fl': set(sd['catalytic']),
            }

        # ---------------- WCN ----------------
        if 'wcn' in modules:
            for site, info in site_info.items():
                if info['com'] is None:
                    continue
                # exclude residues within SEQ_SEP of any catalytic residue (in FL numbering)
                excl = np.zeros(len(ca), dtype=bool)
                for c in info['cat_fl']:
                    near = (res_fl_arr >= c - SEQ_SEP) & (res_fl_arr <= c + SEQ_SEP)
                    excl |= near
                d = np.linalg.norm(ca[~excl] - info['com'], axis=1)
                d = d[d > 1e-3]
                row = {'combo': combo, 'model': model, 'site': site}
                for cut in (5.0, 8.0, 12.0):
                    in_r = d < cut
                    # weighted coord: sum 1/r within cutoff; raw count too
                    row[f'wcn_{int(cut)}A'] = float(np.sum(1.0 / d[in_r])) if in_r.any() else 0.0
                    row[f'n_in_{int(cut)}A'] = int(in_r.sum())
                out['wcn'].append(row)

        # ---------------- SAA ----------------
        if 'saa' in modules:
            for site, info in site_info.items():
                if info['com'] is None:
                    continue
                excl = np.zeros(len(ca), dtype=bool)
                for c in info['cat_fl']:
                    near = (res_fl_arr >= c - SEQ_SEP) & (res_fl_arr <= c + SEQ_SEP)
                    excl |= near
                blockers = ca[~excl]
                acc, cone = ray_accessibility(info['com'], blockers)
                shell_n = shell_density(info['com'], ca)  # all residues
                out['saa'].append({
                    'combo': combo, 'model': model, 'site': site,
                    'accessible_frac': acc,
                    'max_cone_deg': cone,
                    'shell_n': shell_n,
                })

        # ---------------- sites: distances, radial, contact profile ----------------
        if 'sites' in modules:
            site_list = sorted(site_info.keys())
            # Inter-site distances
            for i, si in enumerate(site_list):
                if site_info[si]['com'] is None: continue
                for sj in site_list[i + 1:]:
                    if site_info[sj]['com'] is None: continue
                    d = float(np.linalg.norm(site_info[si]['com'] - site_info[sj]['com']))
                    out['sites'].append({
                        'combo': combo, 'model': model,
                        'kind': 'pair_dist', 'site_i': si, 'site_j': sj,
                        'value': d,
                    })
            # Radial position
            for site, info in site_info.items():
                if info['com'] is None: continue
                d_to_com = float(np.linalg.norm(info['com'] - chain_com))
                out['sites'].append({
                    'combo': combo, 'model': model,
                    'kind': 'radial', 'site_i': site, 'site_j': '',
                    'value': d_to_com,
                    'rg': rg, 'ratio': d_to_com / rg if rg > 0 else np.nan,
                })
            # Per-site contact profile: count residues from each atS domain within Cα 10Å
            #   of any catalytic residue of the site (exclude residues within SEQ_SEP)
            for site, info in site_info.items():
                if len(info['cat_idx']) == 0: continue
                cat_pos = ca[info['cat_idx']]  # (n_cat, 3)
                excl = np.zeros(len(ca), dtype=bool)
                for c in info['cat_fl']:
                    near = (res_fl_arr >= c - SEQ_SEP) & (res_fl_arr <= c + SEQ_SEP)
                    excl |= near
                pool_ca = ca[~excl]
                pool_atS = atS_arr[~excl]
                # min distance from each pool residue to any catalytic residue
                d_min = np.linalg.norm(pool_ca[:, None, :] - cat_pos[None, :, :], axis=2).min(axis=1)
                in_contact = d_min < 10.0
                # group by atS domain
                for d_atS in ATS_ORDER:
                    n = int(np.sum(in_contact & (pool_atS == d_atS)))
                    if n > 0:
                        out['sites'].append({
                            'combo': combo, 'model': model,
                            'kind': 'contact_profile',
                            'site_i': site, 'site_j': d_atS,
                            'value': n,
                        })

        # ---------------- topology ----------------
        if 'topology' in modules:
            # per-domain COM distances (11x11 fragment for those domains present)
            doms_present = sorted(set(atS_arr.tolist()) - {'none'})
            dom_coms = {}
            for d in doms_present:
                m = atS_arr == d
                dom_coms[d] = ca[m].mean(0)
            # store per pair
            for i, di in enumerate(doms_present):
                for dj in doms_present[i + 1:]:
                    out['topology'].append({
                        'combo': combo, 'model': model,
                        'kind': 'dmap', 'a': di, 'b': dj,
                        'value': float(np.linalg.norm(dom_coms[di] - dom_coms[dj])),
                    })
            # rg, ree
            out['topology'].append({
                'combo': combo, 'model': model,
                'kind': 'rg', 'a': '', 'b': '',
                'value': rg, 'n_res': int(len(ca)),
            })
            out['topology'].append({
                'combo': combo, 'model': model,
                'kind': 'ree', 'a': '', 'b': '',
                'value': ree,
            })
            # FNC
            if fl_native_contacts is not None:
                # compute current Cα contacts at <10 Å, |i-j|>3
                D = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=2)
                seq_diff = np.abs(res_fl_arr[:, None] - res_fl_arr[None, :])
                pair_mask = (D < 10.0) & (seq_diff > 3) & np.triu(np.ones_like(D, dtype=bool), 1)
                ii, jj = np.where(pair_mask)
                cur_pairs = set(zip(res_fl_arr[ii].tolist(), res_fl_arr[jj].tolist()))
                # restrict native to pairs whose residues are present here
                applicable = {(a, b) for (a, b) in fl_native_contacts
                              if a in fl_set and b in fl_set}
                preserved = applicable & cur_pairs
                out['topology'].append({
                    'combo': combo, 'model': model,
                    'kind': 'fnc', 'a': '', 'b': '',
                    'value': float(len(preserved) / len(applicable)) if applicable else np.nan,
                    'n_native_applicable': len(applicable),
                    'n_preserved': len(preserved),
                })
    return out


# ============================================================================
# FL native contacts (cached)
# ============================================================================

def compute_fl_native_contacts():
    """Cα<10Å, |fl_i - fl_j|>3, in FL combo seed-1_sample-0. Cached as JSON."""
    cache = A / 'fl_native_contacts_ca10.json'
    if cache.exists():
        data = json.load(open(cache))
        return set(tuple(p) for p in data['pairs'])
    fl_files = sorted(PC.glob('rrm1_rrm2_rrm3_kh1-kh6_kh7a_md1_md2_md3_khb-kh8_wwe_art*__residues.parquet'))
    if not fl_files:
        # fall back to glob with different name
        fl_files = sorted(PC.glob('*md1_md2_md3*wwe_art__residues.parquet'))
    if not fl_files:
        raise RuntimeError('FL combo residues parquet not found')
    fp = fl_files[0]
    df = pd.read_parquet(fp)
    sub = df[df.model == 'seed-1_sample-0'].sort_values('res_construct')
    if sub.empty:
        sub = df[df.model == df.model.unique()[0]].sort_values('res_construct')
    ca = sub[['ca_x', 'ca_y', 'ca_z']].values
    res_fl = sub['res_fl'].values
    D = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=2)
    seq_diff = np.abs(res_fl[:, None] - res_fl[None, :])
    mask = (D < 10.0) & (seq_diff > 3) & np.triu(np.ones_like(D, dtype=bool), 1)
    ii, jj = np.where(mask)
    pairs = list(zip(res_fl[ii].tolist(), res_fl[jj].tolist()))
    json.dump({'source_file': str(fp), 'n_pairs': len(pairs), 'pairs': pairs}, open(cache, 'w'))
    print(f'  cached {len(pairs)} FL native contacts -> {cache}')
    return set(pairs)


# ============================================================================
# CLI dispatcher
# ============================================================================

def run(modules, workers=16):
    sites_def = load_active_sites()
    print(f'Active sites: {list(sites_def.keys())}')
    fl_native = None
    if 'topology' in modules:
        fl_native = compute_fl_native_contacts()
        print(f'  FL native contacts: {len(fl_native)} pairs')

    files = sorted(PC.glob('*__residues.parquet'))
    print(f'Combos: {len(files)} | modules: {modules} | workers: {workers}')

    bins = {m: [] for m in modules}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(process_combo, p, sites_def, fl_native, modules): p for p in files}
        n = 0
        for fut in as_completed(futs):
            try:
                r = fut.result()
                for m in modules:
                    bins[m].extend(r.get(m, []))
            except Exception as e:
                print(f'  err {futs[fut].name}: {e}')
            n += 1
            if n % 200 == 0:
                print(f'  [{n}/{len(files)}]')

    # Save
    if 'sasa' in modules and bins['sasa']:
        pd.DataFrame(bins['sasa']).to_parquet(A / 'site_sasa.parquet', index=False)
        print(f'wrote site_sasa.parquet ({len(bins["sasa"])} rows)')
    if 'wcn' in modules and bins['wcn']:
        pd.DataFrame(bins['wcn']).to_parquet(A / 'site_wcn.parquet', index=False)
        print(f'wrote site_wcn.parquet ({len(bins["wcn"])} rows)')
    if 'saa' in modules and bins['saa']:
        pd.DataFrame(bins['saa']).to_parquet(A / 'site_saa.parquet', index=False)
        print(f'wrote site_saa.parquet ({len(bins["saa"])} rows)')
    if 'sites' in modules and bins['sites']:
        df = pd.DataFrame(bins['sites'])
        df[df.kind == 'pair_dist'].to_parquet(A / 'site_geometry.parquet', index=False)
        df[df.kind == 'radial'].to_parquet(A / 'site_radial.parquet', index=False)
        df[df.kind == 'contact_profile'].to_parquet(A / 'site_contacts.parquet', index=False)
        print(f'wrote site_geometry / site_radial / site_contacts')
    if 'topology' in modules and bins['topology']:
        df = pd.DataFrame(bins['topology'])
        df[df.kind.isin(['rg', 'ree'])].to_parquet(A / 'combo_topology.parquet', index=False)
        df[df.kind == 'dmap'].to_parquet(A / 'combo_dmap.parquet', index=False)
        df[df.kind == 'fnc'].to_parquet(A / 'combo_fnc.parquet', index=False)
        print(f'wrote combo_topology / combo_dmap / combo_fnc')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['sasa', 'wcn', 'saa', 'sites', 'topology', 'all'])
    ap.add_argument('--workers', type=int, default=16)
    args = ap.parse_args()
    if args.cmd == 'all':
        modules = ('sasa', 'wcn', 'saa', 'sites', 'topology')
    else:
        modules = (args.cmd,)
    run(modules, workers=args.workers)
