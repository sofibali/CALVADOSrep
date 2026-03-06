#!/usr/bin/env python3
"""
Analyze existing PARP14 0.5 ns CALVADOS simulation.

Computes per-domain Rg, RMSD, FNC, RMSF, interdomain COM distances,
and interdomain contact maps.

Usage:
    python analyze_parp14.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import MDAnalysis as mda
from calvados.analysis import calc_rg, calc_rmsd, calc_fnc, calc_cmap, cmap_traj, calc_dmap
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================

CWD = os.path.dirname(os.path.abspath(__file__))
SIM_PATH = os.path.join(CWD, 'parp14')
DATA_PATH = os.path.join(CWD, 'data')
os.makedirs(DATA_PATH, exist_ok=True)

PDB = os.path.join(SIM_PATH, 'top.pdb')
DCD = os.path.join(SIM_PATH, 'parp14.dcd')

# Domain definitions (full-length PARP14, 1-indexed residue ranges)
DOMAINS = {
    'RRM1':    (6, 88),
    'RRM2':    (150, 223),
    'RRM3':    (227, 301),
    'KH1':	(315,	396),
    'KH2':	(397,	465),
    'KH3':	(466,	530),
    'KH4':	(531,	596),
    'KH5':	(597,	669),
    'KH6':	(670,	737),
    'KH7a':	(738,	776),
    'MD1':     (791, 978),
    'MD2':     (1003, 1190),
    'MD3':     (1216, 1387),
    'KHb':      (1425,	1453),
    'KH8':      (1454,	1533),
    'WWE':     (1534, 1601),
    'ART':     (1605, 1801),
}

DOMAIN_NAMES = list(DOMAINS.keys())
N_DOMAINS = len(DOMAIN_NAMES)

# ============================================================
# Helper functions
# ============================================================

def resid_selection(start, end):
    """MDAnalysis selection string for 1-indexed residue range."""
    return f'resid {start}:{end}'


def compute_domain_com(u, domain_selections, start=None, stop=None, step=None):
    """Compute COM trajectory for each domain."""
    n_frames = len(u.trajectory[start:stop:step])
    com_traj = np.zeros((n_frames, len(domain_selections), 3))
    for t, ts in enumerate(u.trajectory[start:stop:step]):
        for d, ag in enumerate(domain_selections):
            com_traj[t, d] = ag.center_of_mass() / 10.0  # Angstrom -> nm
    return com_traj


# ============================================================
# Main analysis
# ============================================================

def main():
    print("=" * 60)
    print("PARP14 Simulation Analysis (0.5 ns)")
    print("=" * 60)

    # Load universe
    print(f"\nLoading: {PDB}")
    print(f"         {DCD}")
    u = mda.Universe(PDB, DCD, in_memory=True)
    print(f"  Frames: {len(u.trajectory)}")
    print(f"  Atoms:  {len(u.atoms)}")

    # Reference universe for RMSD/FNC
    uref = mda.Universe(PDB)

    # Create atom groups for each domain
    domain_ags = {}
    for name, (s, e) in DOMAINS.items():
        ag = u.select_atoms(resid_selection(s, e))
        domain_ags[name] = ag
        print(f"  {name}: resid {s}-{e}, {len(ag)} atoms")

    # ----------------------------------------------------------
    # 1. Per-domain Rg
    # ----------------------------------------------------------
    print("\n--- Per-domain Rg ---")
    domain_rgs = {}
    for name, ag in domain_ags.items():
        rg = calc_rg(u, ag)
        domain_rgs[name] = rg
        print(f"  {name}: Rg = {np.mean(rg):.3f} +/- {np.std(rg):.3f} nm")
    np.savez(os.path.join(DATA_PATH, 'domain_rgs.npz'), **domain_rgs)

    # Plot Rg time series
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharex=True)
    for idx, name in enumerate(DOMAIN_NAMES):
        ax = axes.flat[idx]
        rg = domain_rgs[name]
        frames = np.arange(len(rg))
        ax.plot(frames, rg, lw=0.5, color='C0')
        ax.axhline(np.mean(rg), color='red', ls='--', lw=1, label=f'mean={np.mean(rg):.2f}')
        ax.set_title(name, fontsize=11)
        ax.set_ylabel('Rg (nm)')
        ax.legend(fontsize=8)
    for ax in axes[1]:
        ax.set_xlabel('Frame')
    fig.suptitle('Per-domain Radius of Gyration', fontsize=14, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(DATA_PATH, 'domain_rg_timeseries.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: domain_rg_timeseries.png")

    # ----------------------------------------------------------
    # 2. Per-domain RMSD and RMSF
    # ----------------------------------------------------------
    print("\n--- Per-domain RMSD and RMSF ---")
    domain_rmsd_ref = {}
    domain_rmsd_mean = {}
    domain_rmsf = {}
    for name, (s, e) in DOMAINS.items():
        sel = resid_selection(s, e)
        # Reload universe for each domain (calc_rmsd modifies trajectory in memory)
        u_dom = mda.Universe(PDB, DCD, in_memory=True)
        uref_dom = mda.Universe(PDB)
        rmsd_ref, rmsd_mean, rmsf = calc_rmsd(u_dom, uref_dom, select=sel)
        domain_rmsd_ref[name] = rmsd_ref
        domain_rmsd_mean[name] = rmsd_mean
        domain_rmsf[name] = rmsf
        print(f"  {name}: RMSD_ref = {np.mean(rmsd_ref[2]):.3f} A, RMSD_mean = {np.mean(rmsd_mean[2]):.3f} A")

    # Plot RMSD time series
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharex=True)
    for idx, name in enumerate(DOMAIN_NAMES):
        ax = axes.flat[idx]
        rmsd = domain_rmsd_ref[name]
        ax.plot(rmsd[0], rmsd[2] / 10., lw=0.5, color='C0', label='vs ref')
        rmsd_m = domain_rmsd_mean[name]
        ax.plot(rmsd_m[0], rmsd_m[2] / 10., lw=0.5, color='C1', label='vs mean')
        ax.set_title(name, fontsize=11)
        ax.set_ylabel('RMSD (nm)')
        ax.legend(fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel('Frame')
    fig.suptitle('Per-domain RMSD', fontsize=14, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(DATA_PATH, 'domain_rmsd_timeseries.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: domain_rmsd_timeseries.png")

    # Plot per-residue RMSF (all domains on one plot)
    fig, ax = plt.subplots(figsize=(16, 5))
    for name, (s, e) in DOMAINS.items():
        rmsf = domain_rmsf[name] / 10.  # Angstrom -> nm
        resids = np.arange(s, e + 1)
        ax.plot(resids, rmsf, lw=1, label=name)
    # Add domain boundary lines
    for name, (s, e) in DOMAINS.items():
        ax.axvspan(s, e, alpha=0.1)
    ax.set_xlabel('Residue')
    ax.set_ylabel('RMSF (nm)')
    ax.set_title('Per-residue RMSF (within-domain alignment)')
    ax.legend(fontsize=8, ncol=4)
    fig.tight_layout()
    fig.savefig(os.path.join(DATA_PATH, 'domain_rmsf.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: domain_rmsf.png")

    # ----------------------------------------------------------
    # 3. Per-domain Fraction of Native Contacts (FNC)
    # ----------------------------------------------------------
    print("\n--- Per-domain FNC ---")
    domain_fnc = {}
    for name, (s, e) in DOMAINS.items():
        sel = resid_selection(s, e)
        u_fnc = mda.Universe(PDB, DCD, in_memory=True)
        uref_fnc = mda.Universe(PDB)
        fnc = calc_fnc(u_fnc, uref_fnc, sel, cutoff=1.0)
        domain_fnc[name] = fnc
        print(f"  {name}: FNC = {np.mean(fnc):.3f} +/- {np.std(fnc):.3f}")
    np.savez(os.path.join(DATA_PATH, 'domain_fnc.npz'), **domain_fnc)

    # Plot FNC time series
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharex=True)
    for idx, name in enumerate(DOMAIN_NAMES):
        ax = axes.flat[idx]
        fnc = domain_fnc[name]
        frames = np.arange(len(fnc))
        ax.plot(frames, fnc, lw=0.5, color='C0')
        ax.axhline(np.mean(fnc), color='red', ls='--', lw=1, label=f'mean={np.mean(fnc):.2f}')
        ax.set_title(name, fontsize=11)
        ax.set_ylabel('FNC')
        ax.set_ylim(0, 1.1)
        ax.legend(fontsize=8)
    for ax in axes[1]:
        ax.set_xlabel('Frame')
    fig.suptitle('Per-domain Fraction of Native Contacts', fontsize=14, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(DATA_PATH, 'domain_fnc_timeseries.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: domain_fnc_timeseries.png")

    # ----------------------------------------------------------
    # 4. Interdomain COM distances
    # ----------------------------------------------------------
    print("\n--- Interdomain COM distances ---")
    u_com = mda.Universe(PDB, DCD, in_memory=True)
    domain_ags_com = [u_com.select_atoms(resid_selection(s, e)) for s, e in DOMAINS.values()]
    com_traj = compute_domain_com(u_com, domain_ags_com)

    # Mean pairwise distance matrix
    mean_dist = np.zeros((N_DOMAINS, N_DOMAINS))
    std_dist = np.zeros((N_DOMAINS, N_DOMAINS))
    for i in range(N_DOMAINS):
        for j in range(N_DOMAINS):
            dists = np.linalg.norm(com_traj[:, i, :] - com_traj[:, j, :], axis=1)
            mean_dist[i, j] = np.mean(dists)
            std_dist[i, j] = np.std(dists)

    np.save(os.path.join(DATA_PATH, 'interdomain_com_dist.npy'), mean_dist)
    np.save(os.path.join(DATA_PATH, 'interdomain_com_dist_std.npy'), std_dist)

    # Plot interdomain COM distance heatmap
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(mean_dist, cmap='viridis', origin='lower')
    ax.set_xticks(range(N_DOMAINS))
    ax.set_xticklabels(DOMAIN_NAMES, rotation=45, ha='right')
    ax.set_yticks(range(N_DOMAINS))
    ax.set_yticklabels(DOMAIN_NAMES)
    cbar = plt.colorbar(im, ax=ax, label='COM distance (nm)')
    # Annotate with values
    for i in range(N_DOMAINS):
        for j in range(N_DOMAINS):
            ax.text(j, i, f'{mean_dist[i,j]:.1f}', ha='center', va='center', fontsize=7,
                    color='white' if mean_dist[i,j] > np.median(mean_dist) else 'black')
    ax.set_title('Mean Interdomain COM Distance (nm)')
    fig.tight_layout()
    fig.savefig(os.path.join(DATA_PATH, 'interdomain_com_distance.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: interdomain_com_distance.png")

    # ----------------------------------------------------------
    # 5. Interdomain contact frequency
    # ----------------------------------------------------------
    print("\n--- Interdomain contact map ---")
    u_cmap = mda.Universe(PDB, DCD, in_memory=True)
    domain_ags_cmap = {name: u_cmap.select_atoms(resid_selection(s, e)) for name, (s, e) in DOMAINS.items()}

    contact_freq = np.zeros((N_DOMAINS, N_DOMAINS))
    for i, name_i in enumerate(DOMAIN_NAMES):
        for j, name_j in enumerate(DOMAIN_NAMES):
            if j < i:
                contact_freq[i, j] = contact_freq[j, i]
                continue
            cmap = cmap_traj(u_cmap, domain_ags_cmap[name_i], domain_ags_cmap[name_j], cutoff=1.0)
            # Average contact frequency: mean over all residue pairs
            contact_freq[i, j] = np.mean(cmap)

    np.save(os.path.join(DATA_PATH, 'interdomain_contact_freq.npy'), contact_freq)

    # Plot interdomain contact frequency heatmap
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(contact_freq, cmap='hot_r', origin='lower', vmin=0)
    ax.set_xticks(range(N_DOMAINS))
    ax.set_xticklabels(DOMAIN_NAMES, rotation=45, ha='right')
    ax.set_yticks(range(N_DOMAINS))
    ax.set_yticklabels(DOMAIN_NAMES)
    cbar = plt.colorbar(im, ax=ax, label='Mean contact frequency')
    for i in range(N_DOMAINS):
        for j in range(N_DOMAINS):
            ax.text(j, i, f'{contact_freq[i,j]:.3f}', ha='center', va='center', fontsize=7,
                    color='white' if contact_freq[i,j] > 0.5 * np.max(contact_freq) else 'black')
    ax.set_title('Interdomain Contact Frequency')
    fig.tight_layout()
    fig.savefig(os.path.join(DATA_PATH, 'interdomain_contact_freq.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: interdomain_contact_freq.png")

    # ----------------------------------------------------------
    # 6. Global properties
    # ----------------------------------------------------------
    print("\n--- Global properties ---")
    from calvados.analysis import calc_ete
    u_glob = mda.Universe(PDB, DCD, in_memory=True)
    ag_all = u_glob.select_atoms('all')
    rg_global = calc_rg(u_glob, ag_all)
    etes, ete_m, ete_sem = calc_ete(u_glob, ag_all)
    print(f"  Global Rg: {np.mean(rg_global):.3f} +/- {np.std(rg_global):.3f} nm")
    print(f"  Global Ree: {ete_m:.3f} +/- {ete_sem:.3f} nm")

    # ----------------------------------------------------------
    # 7. Summary table
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"{'Domain':<10} {'Rg (nm)':<15} {'FNC':<15} {'RMSD_ref (nm)':<15}")
    print("-" * 55)
    for name in DOMAIN_NAMES:
        rg_m = np.mean(domain_rgs[name])
        rg_s = np.std(domain_rgs[name])
        fnc_m = np.mean(domain_fnc[name])
        fnc_s = np.std(domain_fnc[name])
        rmsd_m = np.mean(domain_rmsd_ref[name][2]) / 10.
        print(f"{name:<10} {rg_m:.3f}+/-{rg_s:.3f}  {fnc_m:.3f}+/-{fnc_s:.3f}  {rmsd_m:.3f}")
    print("-" * 55)
    print(f"Global Rg: {np.mean(rg_global):.3f} +/- {np.std(rg_global):.3f} nm")
    print(f"Global Ree: {ete_m:.3f} +/- {ete_sem:.3f} nm")
    print(f"\nAll outputs saved to: {DATA_PATH}")
    print("=" * 60)


if __name__ == '__main__':
    main()
