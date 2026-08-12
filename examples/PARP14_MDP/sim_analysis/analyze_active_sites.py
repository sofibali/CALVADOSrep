#!/usr/bin/env python3
"""
Active site analysis for PARP14 CALVADOS simulations.

Analyzes how active sites move and change in exposure during MD:
  1. Active site RMSF (per-residue fluctuation)
  2. Contact number time series (exposure proxy: fewer contacts = more exposed)
  3. Pairwise inter-site distances (MD1-MD2-MD3-ART active site COM distances)
  4. Active site burial relative to protein surface (distance-to-COM / Rg)
  5. Per-site contact profiles (which domains contact each active site)

Set/construct handling, active-site numbering, domain boundaries and trajectory
paths all come from the shared `sim_registry` module, so this works for ANY named
set, fragment, or new simulation folder with no edits — point it at a folder with
--sim-folder.

Usage:
    python analyze_active_sites.py --set fl                 # named set, all 25 reps (ensemble)
    python analyze_active_sites.py --set md --seed 1 --sample 0   # single replicate, full plots
    python analyze_active_sites.py --sim-folder mynewconstruct    # any new simulation folder
    python analyze_active_sites.py --sim-folder path --units md1l1 md2 md3
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
import os
import warnings
from argparse import ArgumentParser

import os as _os_boot, sys as _sys_boot  # __ROOTBOOT__ (script lives in a subfolder; root=parent)
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__))))
import sim_registry as reg

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(CWD, 'data')
os.makedirs(DATA_PATH, exist_ok=True)

SITE_NAMES = reg.SITE_NAMES
SITE_COLORS = reg.SITE_COLORS


def get_active_sites(set_key):
    """Active sites with residue numbers in the set's construct numbering."""
    return reg.get_active_sites_for_set(set_key)


def get_domains(set_key):
    """Domain boundaries (construct numbering) for the set's contact profile."""
    return reg.get_construct_domains_for_set(set_key)


# ============================================================
# Analysis functions
# ============================================================

def compute_contact_number(u, site_ag, all_ag, cutoff_nm=1.0, seq_sep=5):
    """
    Compute contact number for active site residues over trajectory.

    Contact number = count of non-local residues (|i-j| > seq_sep) within cutoff.
    Lower contact number = more exposed. Higher = more buried.

    Returns: (n_frames,) array of mean contact number per site residue.
    """
    cutoff_A = cutoff_nm * 10.0  # nm -> Angstrom
    site_resids = site_ag.resids
    all_resids = all_ag.resids

    contact_numbers = np.zeros(len(u.trajectory))
    for t, ts in enumerate(u.trajectory):
        site_pos = site_ag.positions  # (n_site, 3) Angstrom
        all_pos = all_ag.positions    # (n_all, 3) Angstrom

        n_contacts = 0
        for i, (spos, sresid) in enumerate(zip(site_pos, site_resids)):
            dists = np.linalg.norm(all_pos - spos, axis=1)
            # Non-local contacts within cutoff
            mask = (dists < cutoff_A) & (np.abs(all_resids - sresid) > seq_sep)
            n_contacts += np.sum(mask)

        # Mean contact number per site residue
        contact_numbers[t] = n_contacts / len(site_resids)

    return contact_numbers


def compute_site_com(u, site_ag):
    """Compute COM trajectory for an active site atom group. Returns (n_frames, 3) in nm."""
    coms = np.zeros((len(u.trajectory), 3))
    for t, ts in enumerate(u.trajectory):
        coms[t] = site_ag.center_of_mass() / 10.0  # Angstrom -> nm
    return coms


def compute_rmsf_site(u, site_ag):
    """Compute per-residue RMSF for active site residues. Returns RMSF in nm."""
    positions = np.zeros((len(u.trajectory), len(site_ag), 3))
    for t, ts in enumerate(u.trajectory):
        positions[t] = site_ag.positions / 10.0  # nm

    mean_pos = positions.mean(axis=0)
    deviations = positions - mean_pos
    rmsf = np.sqrt(np.mean(np.sum(deviations**2, axis=2), axis=0))
    return rmsf


def compute_domain_contacts_per_site(u, site_ag, domains):
    """
    For each active site, compute fraction of frames where it contacts each domain.
    Contact = any site residue within cutoff of any domain residue.
    Returns dict: {domain_name: contact_fraction}
    """
    cutoff_A = 10.0  # 1.0 nm in Angstrom
    site_resids = set(site_ag.resids)

    domain_ags = {}
    for dname, (s, e) in domains.items():
        ag = u.select_atoms(f'resid {s}:{e}')
        # Exclude site residues from domain to avoid self-contacts
        non_site = [a for a in ag if a.resid not in site_resids]
        if non_site:
            domain_ags[dname] = u.select_atoms(
                f'resid {s}:{e}') if len(non_site) == len(ag) else mda.AtomGroup(non_site)
        else:
            domain_ags[dname] = None

    contact_fracs = {}
    for dname, dag in domain_ags.items():
        if dag is None or len(dag) == 0:
            contact_fracs[dname] = 0.0
            continue

        n_contact_frames = 0
        for ts in u.trajectory:
            site_pos = site_ag.positions
            dom_pos = dag.positions
            # Check minimum distance between any site-domain residue pair
            for sp in site_pos:
                dists = np.linalg.norm(dom_pos - sp, axis=1)
                if np.min(dists) < cutoff_A:
                    n_contact_frames += 1
                    break

        contact_fracs[dname] = n_contact_frames / len(u.trajectory)

    return contact_fracs


# ============================================================
# Main analysis
# ============================================================

def analyze_simulation(pdb_path, dcd_path, set_key, label, data_prefix,
                       skip_frames=0):
    """Run active site analysis on a single simulation."""
    print(f"\n{'='*60}")
    print(f"Active Site Analysis: {label}")
    print(f"{'='*60}")

    # Load trajectory
    print(f"  Loading: {pdb_path}")
    print(f"           {dcd_path}")
    u = mda.Universe(pdb_path, dcd_path, in_memory=True)
    n_total = len(u.trajectory)
    print(f"  Total frames: {n_total}")

    if skip_frames > 0 and skip_frames < n_total:
        # Slice trajectory to skip equilibration
        # MDAnalysis doesn't easily slice in_memory, so we track frame range
        print(f"  Skipping first {skip_frames} frames (equilibration)")

    active_sites = get_active_sites(set_key)
    domains = get_domains(set_key)
    # Only the active sites actually present in this construct.
    site_names = [s for s in SITE_NAMES if s in active_sites]
    all_ag = u.select_atoms('all')

    # --- 1. Active site RMSF ---
    print("\n--- Active site RMSF ---")
    site_rmsfs = {}
    for sname in site_names:
        cat_resids = active_sites[sname]['catalytic']
        pocket_resids = active_sites[sname]['pocket']

        cat_sel = ' or '.join([f'resid {r}' for r in cat_resids])
        pocket_sel = ' or '.join([f'resid {r}' for r in pocket_resids])

        cat_ag = u.select_atoms(cat_sel)
        pocket_ag = u.select_atoms(pocket_sel)

        if len(cat_ag) == 0:
            print(f"  {sname}: no catalytic residues found, skipping")
            continue

        rmsf_cat = compute_rmsf_site(u, cat_ag)
        rmsf_pocket = compute_rmsf_site(u, pocket_ag)

        site_rmsfs[sname] = {
            'catalytic_resids': cat_resids,
            'catalytic_rmsf': rmsf_cat,
            'pocket_resids': pocket_resids,
            'pocket_rmsf': rmsf_pocket,
        }
        print(f"  {sname} catalytic RMSF: {np.mean(rmsf_cat):.3f} +/- {np.std(rmsf_cat):.3f} nm")
        print(f"  {sname} pocket RMSF:    {np.mean(rmsf_pocket):.3f} +/- {np.std(rmsf_pocket):.3f} nm")

    # --- 2. Contact number time series (exposure proxy) ---
    print("\n--- Active site contact number (exposure proxy) ---")
    site_contacts = {}
    for sname in site_names:
        pocket_resids = active_sites[sname]['pocket']
        pocket_sel = ' or '.join([f'resid {r}' for r in pocket_resids])
        pocket_ag = u.select_atoms(pocket_sel)

        if len(pocket_ag) == 0:
            continue

        cn = compute_contact_number(u, pocket_ag, all_ag, cutoff_nm=1.0, seq_sep=5)
        site_contacts[sname] = cn
        print(f"  {sname}: mean contact number = {np.mean(cn):.2f} +/- {np.std(cn):.2f}")

    # --- 3. Pairwise inter-site distances ---
    print("\n--- Inter-active-site COM distances ---")
    site_coms = {}
    for sname in site_names:
        cat_resids = active_sites[sname]['catalytic']
        cat_sel = ' or '.join([f'resid {r}' for r in cat_resids])
        cat_ag = u.select_atoms(cat_sel)
        if len(cat_ag) > 0:
            site_coms[sname] = compute_site_com(u, cat_ag)

    available_sites = [s for s in SITE_NAMES if s in site_coms]
    n_sites = len(available_sites)
    inter_site_dists = {}
    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            si, sj = available_sites[i], available_sites[j]
            d = np.linalg.norm(site_coms[si] - site_coms[sj], axis=1)
            pair = f'{si}-{sj}'
            inter_site_dists[pair] = d
            print(f"  {pair}: {np.mean(d):.2f} +/- {np.std(d):.2f} nm")

    # --- 4. Relative exposure (distance to protein COM / Rg) ---
    print("\n--- Active site relative position (dist_to_COM / Rg) ---")
    # Compute protein COM and Rg per frame
    prot_com = np.zeros((len(u.trajectory), 3))
    prot_rg = np.zeros(len(u.trajectory))
    for t, ts in enumerate(u.trajectory):
        pos = all_ag.positions / 10.0  # nm
        com = np.mean(pos, axis=0)
        prot_com[t] = com
        prot_rg[t] = np.sqrt(np.mean(np.sum((pos - com)**2, axis=1)))

    site_relative_pos = {}
    for sname in available_sites:
        dist_to_com = np.linalg.norm(site_coms[sname] - prot_com, axis=1)
        relative = dist_to_com / prot_rg
        site_relative_pos[sname] = relative
        print(f"  {sname}: relative position = {np.mean(relative):.3f} +/- {np.std(relative):.3f}")
        print(f"         (>1 = outside Rg shell = surface-exposed, <1 = buried)")

    # --- 5. Domain contact profiles ---
    print("\n--- Domain contact profiles per active site ---")
    site_domain_contacts = {}
    for sname in site_names:
        pocket_resids = active_sites[sname]['pocket']
        pocket_sel = ' or '.join([f'resid {r}' for r in pocket_resids])
        pocket_ag = u.select_atoms(pocket_sel)

        if len(pocket_ag) == 0:
            continue

        contacts = compute_domain_contacts_per_site(u, pocket_ag, domains)
        site_domain_contacts[sname] = contacts
        top_contacts = sorted(contacts.items(), key=lambda x: -x[1])[:5]
        print(f"  {sname} top contacts: {', '.join(f'{d}={v:.2f}' for d, v in top_contacts)}")

    # ============================================================
    # Save data
    # ============================================================
    print("\n--- Saving data ---")

    # Save all results as npz
    save_dict = {
        'prot_rg': prot_rg,
    }
    for sname in available_sites:
        if sname in site_contacts:
            save_dict[f'{sname}_contact_number'] = site_contacts[sname]
        save_dict[f'{sname}_com'] = site_coms[sname]
        save_dict[f'{sname}_relative_pos'] = site_relative_pos[sname]

    for pair, d in inter_site_dists.items():
        save_dict[f'dist_{pair}'] = d

    np.savez(os.path.join(DATA_PATH, f'{data_prefix}_active_sites.npz'), **save_dict)
    print(f"  Saved: {data_prefix}_active_sites.npz")

    # ============================================================
    # Plots
    # ============================================================

    # --- Plot 1: Active site RMSF bar chart ---
    fig, ax = plt.subplots(figsize=(10, 5))
    x_pos = 0
    xticks, xlabels = [], []
    for sname in SITE_NAMES:
        if sname not in site_rmsfs:
            continue
        data = site_rmsfs[sname]
        resids = data['pocket_resids']
        rmsf = data['pocket_rmsf']
        cat_set = set(data['catalytic_resids'])

        colors = [SITE_COLORS[sname] if r in cat_set else 'lightgray' for r in resids]
        positions = np.arange(x_pos, x_pos + len(resids))
        ax.bar(positions, rmsf, color=colors, width=0.8, edgecolor='none')

        mid = x_pos + len(resids) / 2
        xticks.append(mid)
        xlabels.append(sname)
        x_pos += len(resids) + 2  # gap between sites

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=12)
    ax.set_ylabel('RMSF (nm)', fontsize=12)
    ax.set_title(f'Active Site Residue Fluctuations — {label}', fontsize=13)
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='gray', label='Pocket residue')]
    for sname in SITE_NAMES:
        legend_elements.append(Patch(facecolor=SITE_COLORS[sname], label=f'{sname} catalytic'))
    ax.legend(handles=legend_elements, fontsize=9, loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(DATA_PATH, f'{data_prefix}_active_site_rmsf.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {data_prefix}_active_site_rmsf.png")

    # --- Plot 2: Contact number time series (exposure) ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for idx, sname in enumerate(SITE_NAMES):
        ax = axes.flat[idx]
        if sname in site_contacts:
            cn = site_contacts[sname]
            frames = np.arange(len(cn))
            ax.plot(frames, cn, lw=0.5, color=SITE_COLORS[sname], alpha=0.7)
            # Running average
            window = max(1, len(cn) // 50)
            if window > 1:
                cn_smooth = np.convolve(cn, np.ones(window)/window, mode='valid')
                ax.plot(np.arange(window-1, len(cn)), cn_smooth,
                        lw=2, color='black', label=f'running avg (w={window})')
            ax.axhline(np.mean(cn), color='red', ls='--', lw=1,
                       label=f'mean={np.mean(cn):.1f}')
            ax.legend(fontsize=8)
        ax.set_title(f'{sname} Active Site', fontsize=11)
        ax.set_ylabel('Contact number')
    for ax in axes[1]:
        ax.set_xlabel('Frame')
    fig.suptitle(f'Active Site Exposure (contact number) — {label}\n'
                 'Lower = more exposed, Higher = more buried', fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(DATA_PATH, f'{data_prefix}_active_site_exposure.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {data_prefix}_active_site_exposure.png")

    # --- Plot 3: Inter-site distance time series ---
    n_pairs = len(inter_site_dists)
    if n_pairs > 0:
        fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
        for idx, (pair, d) in enumerate(inter_site_dists.items()):
            if idx >= 6:
                break
            ax = axes.flat[idx]
            frames = np.arange(len(d))
            ax.plot(frames, d, lw=0.5, color='C0', alpha=0.7)
            ax.axhline(np.mean(d), color='red', ls='--', lw=1,
                       label=f'mean={np.mean(d):.1f} nm')
            ax.set_title(pair, fontsize=11)
            ax.set_ylabel('Distance (nm)')
            ax.legend(fontsize=8)
        # Hide unused subplots
        for idx in range(n_pairs, 6):
            axes.flat[idx].set_visible(False)
        for ax in axes[1]:
            if ax.get_visible():
                ax.set_xlabel('Frame')
        fig.suptitle(f'Inter-Active-Site Distances — {label}', fontsize=13, y=1.01)
        fig.tight_layout()
        fig.savefig(os.path.join(DATA_PATH, f'{data_prefix}_intersite_distances.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {data_prefix}_intersite_distances.png")

    # --- Plot 4: Relative position (exposure vs burial) ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for sname in available_sites:
        rel = site_relative_pos[sname]
        frames = np.arange(len(rel))
        ax.plot(frames, rel, lw=0.5, alpha=0.5, color=SITE_COLORS[sname])
        # Running average
        window = max(1, len(rel) // 50)
        if window > 1:
            rel_smooth = np.convolve(rel, np.ones(window)/window, mode='valid')
            ax.plot(np.arange(window-1, len(rel)), rel_smooth,
                    lw=2, color=SITE_COLORS[sname], label=f'{sname} ({np.mean(rel):.2f})')
        else:
            ax.plot([], [], color=SITE_COLORS[sname], lw=2,
                    label=f'{sname} ({np.mean(rel):.2f})')
    ax.axhline(1.0, color='gray', ls=':', lw=1, label='Rg shell')
    ax.set_xlabel('Frame')
    ax.set_ylabel('Distance to COM / Rg')
    ax.set_title(f'Active Site Radial Position — {label}\n'
                 '>1 = beyond Rg (surface), <1 = within Rg (core)', fontsize=12)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(DATA_PATH, f'{data_prefix}_active_site_radial.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {data_prefix}_active_site_radial.png")

    # --- Plot 5: Domain contact profile heatmap ---
    if site_domain_contacts:
        domain_names = list(domains.keys())
        contact_matrix = np.zeros((len(SITE_NAMES), len(domain_names)))
        for i, sname in enumerate(SITE_NAMES):
            if sname in site_domain_contacts:
                for j, dname in enumerate(domain_names):
                    contact_matrix[i, j] = site_domain_contacts[sname].get(dname, 0.0)

        fig, ax = plt.subplots(figsize=(12, 4))
        im = ax.imshow(contact_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)
        ax.set_yticks(range(len(SITE_NAMES)))
        ax.set_yticklabels(SITE_NAMES, fontsize=11)
        ax.set_xticks(range(len(domain_names)))
        ax.set_xticklabels(domain_names, rotation=45, ha='right', fontsize=10)
        plt.colorbar(im, ax=ax, label='Contact fraction', shrink=0.8)
        # Annotate
        for i in range(len(SITE_NAMES)):
            for j in range(len(domain_names)):
                v = contact_matrix[i, j]
                if v > 0.01:
                    ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=8,
                            color='white' if v > 0.5 else 'black')
        ax.set_title(f'Active Site Domain Contacts — {label}', fontsize=13)
        fig.tight_layout()
        fig.savefig(os.path.join(DATA_PATH, f'{data_prefix}_active_site_domain_contacts.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {data_prefix}_active_site_domain_contacts.png")

    # --- Summary ---
    print(f"\n{'='*60}")
    print("ACTIVE SITE SUMMARY")
    print(f"{'='*60}")
    print(f"{'Site':<6} {'RMSF_cat (nm)':<15} {'Contact#':<15} {'Radial pos':<15}")
    print("-" * 51)
    for sname in SITE_NAMES:
        rmsf_str = f"{np.mean(site_rmsfs[sname]['catalytic_rmsf']):.3f}" if sname in site_rmsfs else "N/A"
        cn_str = f"{np.mean(site_contacts[sname]):.1f}" if sname in site_contacts else "N/A"
        rp_str = f"{np.mean(site_relative_pos[sname]):.3f}" if sname in site_relative_pos else "N/A"
        print(f"{sname:<6} {rmsf_str:<15} {cn_str:<15} {rp_str:<15}")
    print(f"{'='*60}")


def analyze_ensemble(set_key, seeds=reg.SEEDS, samples=reg.SAMPLES, skip_frames=50):
    """Analyze all replicates of a set/folder and compute ensemble statistics."""
    label = reg.SETS[set_key]['label']
    print("\n" + "=" * 70)
    print(f"ENSEMBLE ACTIVE SITE ANALYSIS — {label}")
    print("=" * 70)

    active_sites = get_active_sites(set_key)
    site_names = [s for s in SITE_NAMES if s in active_sites]
    all_contacts = {s: [] for s in site_names}
    all_relative_pos = {s: [] for s in site_names}
    all_inter_dists = {}

    n_analyzed = 0
    for seed in seeds:
        for sample in samples:
            pdb, dcd = reg.get_sim_paths(set_key, seed, sample)
            if not os.path.isfile(pdb) or not os.path.isfile(dcd):
                continue

            print(f"\n  Processing: seed-{seed}_sample-{sample}")
            u = mda.Universe(pdb, dcd, in_memory=True)
            all_ag = u.select_atoms('all')

            # Contact numbers
            for sname in site_names:
                pocket_resids = active_sites[sname]['pocket']
                pocket_sel = ' or '.join([f'resid {r}' for r in pocket_resids])
                pocket_ag = u.select_atoms(pocket_sel)
                if len(pocket_ag) > 0:
                    cn = compute_contact_number(u, pocket_ag, all_ag, cutoff_nm=1.0, seq_sep=5)
                    all_contacts[sname].append(np.mean(cn))

            # Site COMs and relative positions
            site_coms = {}
            prot_com = np.zeros((len(u.trajectory), 3))
            prot_rg = np.zeros(len(u.trajectory))
            for t, ts in enumerate(u.trajectory):
                pos = all_ag.positions / 10.0
                com = np.mean(pos, axis=0)
                prot_com[t] = com
                prot_rg[t] = np.sqrt(np.mean(np.sum((pos - com)**2, axis=1)))

            for sname in site_names:
                cat_resids = active_sites[sname]['catalytic']
                cat_sel = ' or '.join([f'resid {r}' for r in cat_resids])
                cat_ag = u.select_atoms(cat_sel)
                if len(cat_ag) > 0:
                    site_com = compute_site_com(u, cat_ag)
                    site_coms[sname] = site_com
                    dist_to_com = np.linalg.norm(site_com - prot_com, axis=1)
                    relative = dist_to_com / prot_rg
                    all_relative_pos[sname].append(np.mean(relative))

            # Inter-site distances
            available = [s for s in site_names if s in site_coms]
            for i in range(len(available)):
                for j in range(i + 1, len(available)):
                    pair = f'{available[i]}-{available[j]}'
                    d = np.linalg.norm(site_coms[available[i]] - site_coms[available[j]], axis=1)
                    if pair not in all_inter_dists:
                        all_inter_dists[pair] = []
                    all_inter_dists[pair].append(np.mean(d))

            n_analyzed += 1

    if n_analyzed == 0:
        print(f"  No completed simulations found for set '{set_key}'.")
        return

    # --- Ensemble summary plot ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: Contact number distribution
    ax = axes[0]
    for sname in site_names:
        if all_contacts[sname]:
            ax.hist(all_contacts[sname], bins=10, alpha=0.6, color=SITE_COLORS[sname],
                    label=f'{sname} ({np.mean(all_contacts[sname]):.1f})', edgecolor='black')
    ax.set_xlabel('Mean contact number')
    ax.set_ylabel('Count (replicates)')
    ax.set_title('Active Site Burial')
    ax.legend(fontsize=9)

    # Panel 2: Radial position distribution
    ax = axes[1]
    for sname in site_names:
        if all_relative_pos[sname]:
            ax.hist(all_relative_pos[sname], bins=10, alpha=0.6, color=SITE_COLORS[sname],
                    label=f'{sname} ({np.mean(all_relative_pos[sname]):.2f})', edgecolor='black')
    ax.axvline(1.0, color='gray', ls=':', lw=1)
    ax.set_xlabel('Distance to COM / Rg')
    ax.set_ylabel('Count (replicates)')
    ax.set_title('Radial Position (>1 = surface)')
    ax.legend(fontsize=9)

    # Panel 3: Inter-site distances
    ax = axes[2]
    pairs = sorted(all_inter_dists.keys())
    means = [np.mean(all_inter_dists[p]) for p in pairs]
    stds = [np.std(all_inter_dists[p]) for p in pairs]
    ax.bar(range(len(pairs)), means, yerr=stds, color='steelblue',
           edgecolor='black', capsize=3)
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels(pairs, rotation=30, ha='right')
    ax.set_ylabel('Distance (nm)')
    ax.set_title('Inter-Site Distances')

    fig.suptitle(f'{label} — Ensemble ({n_analyzed} replicates)', fontsize=14, y=1.02)
    fig.tight_layout()
    out = os.path.join(DATA_PATH, f'{set_key}_ensemble_active_sites.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved: {os.path.basename(out)}")

    # Print ensemble summary
    print(f"\n{'='*60}")
    print(f"ENSEMBLE SUMMARY ({n_analyzed} replicates)")
    print(f"{'='*60}")
    print(f"{'Site':<6} {'Contact# (mean+/-std)':<25} {'Radial pos (mean+/-std)':<25}")
    print("-" * 56)
    for sname in site_names:
        cn = all_contacts[sname]
        rp = all_relative_pos[sname]
        cn_str = f"{np.mean(cn):.1f} +/- {np.std(cn):.1f}" if cn else "N/A"
        rp_str = f"{np.mean(rp):.3f} +/- {np.std(rp):.3f}" if rp else "N/A"
        print(f"{sname:<6} {cn_str:<25} {rp_str:<25}")
    print(f"{'='*60}")


# ============================================================
# CLI
# ============================================================

def resolve_set_key(args):
    """Resolve which registered set to analyze from --set / --sim-folder."""
    if args.sim_folder:
        # Register the (first) external folder and return its key.
        keys = [reg.register_external_folder(f, units=args.units)
                for f in args.sim_folder]
        if len(keys) > 1:
            print(f"  Note: multiple --sim-folder given; using first ('{keys[0]}'). "
                  f"Run again per folder for the rest.")
        return keys[0]
    s = args.set
    if s in reg.SETS:
        return s
    # bare fragment name?
    reg.discover_fragments()
    if f'frag_{s}' in reg.SETS:
        return f'frag_{s}'
    if s in reg.SETS:
        return s
    raise SystemExit(f"unknown set '{s}'. Known: {', '.join(sorted(reg.SETS))} "
                     f"(or pass --sim-folder PATH)")


def main():
    parser = ArgumentParser(description='PARP14 active site analysis')
    parser.add_argument('--set', default='fl',
                        help='Named set (fl/md/core/mka/norrm/noart/md3art/'
                             'fl_optimized) or fragment name. Default: fl')
    reg.add_sim_folder_args(parser)
    parser.add_argument('--seed', type=int, default=None,
                        help='Specific seed (1-5) for single-replicate analysis')
    parser.add_argument('--sample', type=int, default=None,
                        help='Specific sample (0-4) for single-replicate analysis')
    parser.add_argument('--ensemble', action='store_true',
                        help='Ensemble analysis across all replicates (default '
                             'when no --seed/--sample given)')
    parser.add_argument('--skip-frames', type=int, default=50,
                        help='Initial frames to skip (equilibration, default=50 = 0.5 ns)')
    args = parser.parse_args()

    set_key = resolve_set_key(args)
    label = reg.SETS[set_key]['label']

    if args.seed is not None and args.sample is not None:
        # Single replicate — full per-replicate plots
        pdb, dcd = reg.get_sim_paths(set_key, args.seed, args.sample)
        if not (os.path.isfile(pdb) and os.path.isfile(dcd)):
            raise SystemExit(f"missing trajectory: {dcd}")
        analyze_simulation(pdb, dcd, set_key,
                           label=f'{label} (seed-{args.seed}_sample-{args.sample})',
                           data_prefix=f'{set_key}_seed-{args.seed}_sample-{args.sample}',
                           skip_frames=args.skip_frames)
    else:
        # Ensemble across all replicates (the common case)
        analyze_ensemble(set_key, skip_frames=args.skip_frames)


if __name__ == '__main__':
    main()
