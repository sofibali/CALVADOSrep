#!/usr/bin/env python3
"""
Compare restraint test results per domain and generate:
  1. Per-domain RMSD-to-crystal/AF2 bar plots across all tests
  2. Per-domain Rg ratio + cmap correlation summary figure
  3. PyMOL session with KH domains from each test as separate states

Usage:
    python compare_restraints.py
"""
import os
import json
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

CWD = Path(__file__).resolve().parent
TEST_DIR = CWD / 'restraint_tests'
PLOT_DIR = TEST_DIR / 'comparison_plots'
AF2_PDB = CWD / 'input' / 'parp14.pdb'

TESTS = [
    'baseline_current', 'tight_uniform_700', 'tight_uniform_350',
    'medium_uniform_700', 'medium_uniform_350',
    'loose_uniform_700', 'loose_uniform_350',
]
TEST_LABELS = {
    'baseline_current': 'Baseline\n(no KH)',
    'tight_uniform_700': 'Tight\nk=700',
    'tight_uniform_350': 'Tight\nk=350',
    'medium_uniform_700': 'Medium\nk=700',
    'medium_uniform_350': 'Medium\nk=350',
    'loose_uniform_700': 'Loose\nk=700',
    'loose_uniform_350': 'Loose\nk=350',
}

ANALYSIS_DOMAINS = {
    'RRM1': (1, 145), 'RRM2': (146, 224), 'RRM3': (225, 314),
    'KH1': (315, 384), 'KH2': (385, 454), 'KH3': (455, 520),
    'KH4': (521, 593), 'KH5': (594, 665), 'KH6': (666, 737),
    'KH7a': (738, 789),
    'MD1L1': (790, 1004), 'MD2': (1004, 1193), 'MD3': (1207, 1388),
    'KHb': (1389, 1461), 'KH8': (1462, 1533),
    'WWE': (1534, 1602), 'ART': (1603, 1801),
}

DOMAIN_ORDER = ['RRM1','RRM2','RRM3','KH1','KH2','KH3','KH4','KH5','KH6',
                'KH7a','MD1L1','MD2','MD3','KHb','KH8','WWE','ART']

KH_DOMAINS = ['KH1','KH2','KH3','KH4','KH5','KH6','KH7a','KHb','KH8']

DOMAIN_COLORS = {
    'RRM1': '#1f77b4', 'RRM2': '#2ca02c', 'RRM3': '#aec7e8',
    'KH1': '#ff7f0e', 'KH2': '#ff9f40', 'KH3': '#ffbf70',
    'KH4': '#ffd700', 'KH5': '#e6c200', 'KH6': '#ccac00',
    'KH7a': '#d62728', 'MD1L1': '#17becf', 'MD2': '#9467bd',
    'MD3': '#8c564b', 'KHb': '#e377c2', 'KH8': '#f7b6d2',
    'WWE': '#7f7f7f', 'ART': '#bcbd22',
}

# Crystal structures for RMSD reference
XTAL_DIR = CWD / 'input' / 'xtal_refs'
XTAL_REFS = {
    'MD1L1': str(XTAL_DIR / '3q6z.pdb'),   # 3Q6Z: FL 792-978
    'MD2':   str(XTAL_DIR / '3vfq.pdb'),   # 3VFQ: FL 1005-1191 (MD2 portion)
    'WWE':   str(XTAL_DIR / '3goy.pdb'),   # 3GOY: FL 1532-1720 (WWE portion)
    'ART':   str(XTAL_DIR / '3goy.pdb'),   # 3GOY: FL 1532-1720 (ART portion)
}
XTAL_RANGES = {
    'MD1L1': (792, 978),
    'MD2':   (1005, 1191),
    'WWE':   (1534, 1602),
    'ART':   (1603, 1720),
}


def load_test_traj(test_name):
    """Load trajectory and return CA-only mdtraj trajectory."""
    test_dir = TEST_DIR / test_name
    dcd_files = list(test_dir.glob('*.dcd'))
    top_file = test_dir / 'top.pdb'
    if not dcd_files or not top_file.exists():
        return None
    traj = md.load(str(dcd_files[0]), top=str(top_file))
    ca = traj.topology.select('name CA')
    return traj.atom_slice(ca)


def load_xtal_ca(pdb_path, resid_range):
    """Load CA coords from crystal PDB for a residue range."""
    traj = md.load(str(pdb_path))
    ca = traj.topology.select('name CA')
    ca_traj = traj.atom_slice(ca)
    # Get residue indices within range
    resids = [r.resSeq for r in ca_traj.topology.residues]
    mask = [(resid_range[0] <= r <= resid_range[1]) for r in resids]
    idx = [i for i, m in enumerate(mask) if m]
    if not idx:
        return None, []
    return ca_traj.xyz[0, idx, :], [resids[i] for i in idx]


def compute_domain_rmsd(sim_traj, ref_coords, domain_fl_range, ref_resid_range):
    """Compute per-frame RMSD of a domain in simulation vs reference structure.

    Aligns domain residues and computes RMSD. Handles mismatched lengths
    by finding the overlapping residue range.
    """
    # Simulation domain indices (0-based)
    sim_idx = list(range(domain_fl_range[0] - 1, domain_fl_range[1]))

    # Find overlap between simulation domain and reference
    ref_start, ref_end = ref_resid_range
    overlap_start = max(domain_fl_range[0], ref_start)
    overlap_end = min(domain_fl_range[1], ref_end)

    if overlap_end <= overlap_start:
        return None

    # Indices within simulation array
    sim_overlap = list(range(overlap_start - 1, overlap_end))
    # Indices within reference array
    ref_overlap = list(range(overlap_start - ref_start, overlap_end - ref_start))

    n = min(len(sim_overlap), len(ref_overlap), ref_coords.shape[0])
    if n < 5:
        return None

    sim_overlap = sim_overlap[:n]
    ref_overlap = ref_overlap[:n]

    ref_xyz = ref_coords[ref_overlap]

    # Per-frame RMSD with Kabsch alignment
    rmsds = []
    for frame in range(len(sim_traj)):
        sim_xyz = sim_traj.xyz[frame, sim_overlap, :]

        # Center both
        sim_c = sim_xyz - sim_xyz.mean(axis=0)
        ref_c = ref_xyz - ref_xyz.mean(axis=0)

        # Kabsch (SVD)
        H = sim_c.T @ ref_c
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        sign_matrix = np.diag([1, 1, np.sign(d)])
        R = Vt.T @ sign_matrix @ U.T
        sim_aligned = sim_c @ R

        rmsd = np.sqrt(np.mean(np.sum((sim_aligned - ref_c) ** 2, axis=1)))
        rmsds.append(rmsd)

    return np.array(rmsds)


def main():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    # Load AF2 reference
    af2_traj = md.load(str(AF2_PDB))
    af2_ca = af2_traj.atom_slice(af2_traj.topology.select('name CA'))
    af2_xyz = af2_ca.xyz[0]  # (1801, 3) in nm

    # Load all test trajectories
    print("Loading trajectories...")
    trajs = {}
    for t in TESTS:
        traj = load_test_traj(t)
        if traj is not None:
            trajs[t] = traj
            print(f"  {t}: {traj.n_frames} frames, {traj.n_residues} residues")

    # Load crystal references
    print("\nLoading crystal references...")
    xtal_coords = {}
    for dname, pdb_path in XTAL_REFS.items():
        coords, resids = load_xtal_ca(pdb_path, XTAL_RANGES[dname])
        if coords is not None:
            xtal_coords[dname] = coords
            print(f"  {dname}: {len(resids)} CA from {os.path.basename(pdb_path)} "
                  f"({resids[0]}-{resids[-1]})")

    # Load analysis results
    print("\nLoading analysis results...")
    results = {}
    for t in TESTS:
        f = TEST_DIR / t / 'analysis_results.json'
        if f.exists():
            with open(f) as fh:
                results[t] = json.load(fh)

    # ════════════════════════════════════════════════════════════
    # Compute per-domain RMSD to crystal/AF2 for each test
    # ════════════════════════════════════════════════════════════
    print("\nComputing per-domain RMSD...")

    rmsd_data = {}  # {domain: {test: array_of_rmsds}}
    for dname in DOMAIN_ORDER:
        ds, de = ANALYSIS_DOMAINS[dname]
        rmsd_data[dname] = {}

        # Choose reference: crystal if available, else AF2
        if dname in xtal_coords:
            ref_coords = xtal_coords[dname]
            ref_range = XTAL_RANGES[dname]
            ref_label = 'xtal'
        else:
            ref_coords = af2_xyz[ds-1:de]
            ref_range = (ds, de)
            ref_label = 'AF2'

        for t, traj in trajs.items():
            rmsds = compute_domain_rmsd(traj, ref_coords, (ds, de), ref_range)
            if rmsds is not None:
                rmsd_data[dname][t] = rmsds

        if rmsd_data[dname]:
            means = {t: v.mean() for t, v in rmsd_data[dname].items()}
            best = min(means, key=means.get)
            print(f"  {dname} (vs {ref_label}): best={best} "
                  f"({means[best]:.3f} nm)")

    # ════════════════════════════════════════════════════════════
    # Plot 1: Per-domain RMSD bar chart
    # ════════════════════════════════════════════════════════════
    print("\nGenerating RMSD comparison plot...")
    fig, ax = plt.subplots(figsize=(22, 8))

    x = np.arange(len(DOMAIN_ORDER))
    width = 0.11
    n_tests = len(TESTS)

    test_colors = plt.cm.Set2(np.linspace(0, 1, n_tests))

    for i, t in enumerate(TESTS):
        means = []
        stds = []
        for dname in DOMAIN_ORDER:
            if t in rmsd_data.get(dname, {}):
                rmsds = rmsd_data[dname][t]
                means.append(rmsds.mean())
                stds.append(rmsds.std())
            else:
                means.append(0)
                stds.append(0)

        offset = (i - n_tests/2 + 0.5) * width
        label = TEST_LABELS.get(t, t).replace('\n', ' ')
        ax.bar(x + offset, means, width, yerr=stds, label=label,
               color=test_colors[i], alpha=0.85, capsize=2)

    # Mark crystal-referenced domains
    for j, dname in enumerate(DOMAIN_ORDER):
        if dname in xtal_coords:
            ax.text(j, -0.02, '★', ha='center', fontsize=10, color='red')

    ax.set_xticks(x)
    ax.set_xticklabels(DOMAIN_ORDER, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('RMSD to Reference (nm)', fontsize=12)
    ax.set_title('Per-Domain RMSD Across Restraint Configurations\n'
                 '(★ = crystal structure reference, others = AF2)',
                 fontsize=14)
    ax.legend(fontsize=7, ncol=4, loc='upper left')
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    fig.savefig(PLOT_DIR / 'rmsd_per_domain.png', dpi=200, bbox_inches='tight')
    fig.savefig(PLOT_DIR / 'rmsd_per_domain.pdf', bbox_inches='tight')
    plt.close()

    # ════════════════════════════════════════════════════════════
    # Plot 2: Rg ratio + CMap heatmap
    # ════════════════════════════════════════════════════════════
    print("Generating Rg/CMap heatmap...")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 8), sharex=True)

    rg_matrix = np.full((len(TESTS), len(DOMAIN_ORDER)), np.nan)
    cmap_matrix = np.full((len(TESTS), len(DOMAIN_ORDER)), np.nan)

    for i, t in enumerate(TESTS):
        if t not in results:
            continue
        for j, dname in enumerate(DOMAIN_ORDER):
            dm = results[t]['domain_metrics'].get(dname, {})
            rg_matrix[i, j] = dm.get('rg_ratio', np.nan)
            cmap_matrix[i, j] = dm.get('cmap_corr', np.nan)

    im1 = ax1.imshow(rg_matrix, cmap='RdBu_r', vmin=0.5, vmax=2.0, aspect='auto')
    ax1.set_yticks(range(len(TESTS)))
    ax1.set_yticklabels([TEST_LABELS.get(t, t).replace('\n', ' ') for t in TESTS],
                        fontsize=8)
    plt.colorbar(im1, ax=ax1, label='Rg Ratio (sim/ref)', shrink=0.8)
    ax1.set_title('Rg Ratio (1.0 = ideal)', fontsize=12)

    # Annotate cells
    for i in range(len(TESTS)):
        for j in range(len(DOMAIN_ORDER)):
            val = rg_matrix[i, j]
            if not np.isnan(val):
                color = 'white' if abs(val - 1.0) > 0.4 else 'black'
                ax1.text(j, i, f'{val:.2f}', ha='center', va='center',
                        fontsize=5, color=color)

    im2 = ax2.imshow(cmap_matrix, cmap='YlGn', vmin=0.3, vmax=1.0, aspect='auto')
    ax2.set_xticks(range(len(DOMAIN_ORDER)))
    ax2.set_xticklabels(DOMAIN_ORDER, rotation=45, ha='right', fontsize=10)
    ax2.set_yticks(range(len(TESTS)))
    ax2.set_yticklabels([TEST_LABELS.get(t, t).replace('\n', ' ') for t in TESTS],
                        fontsize=8)
    plt.colorbar(im2, ax=ax2, label='Contact Map Correlation', shrink=0.8)
    ax2.set_title('Contact Map Correlation with Reference (higher = better)', fontsize=12)

    for i in range(len(TESTS)):
        for j in range(len(DOMAIN_ORDER)):
            val = cmap_matrix[i, j]
            if not np.isnan(val):
                color = 'white' if val < 0.5 else 'black'
                ax2.text(j, i, f'{val:.2f}', ha='center', va='center',
                        fontsize=5, color=color)

    plt.tight_layout()
    fig.savefig(PLOT_DIR / 'rg_cmap_heatmap.png', dpi=200, bbox_inches='tight')
    fig.savefig(PLOT_DIR / 'rg_cmap_heatmap.pdf', bbox_inches='tight')
    plt.close()

    # ════════════════════════════════════════════════════════════
    # Plot 3: KH domains focus — RMSD boxplots
    # ════════════════════════════════════════════════════════════
    print("Generating KH domain RMSD boxplots...")
    n_kh = len(KH_DOMAINS)
    fig, axes = plt.subplots(3, 3, figsize=(18, 14))
    axes = axes.flatten()

    for idx, dname in enumerate(KH_DOMAINS):
        ax = axes[idx]
        data = []
        labels = []
        colors = []

        for t in TESTS:
            if t in rmsd_data.get(dname, {}):
                rmsds = rmsd_data[dname][t]
                data.append(rmsds * 10)  # nm -> Å for readability
                label = TEST_LABELS.get(t, t).replace('\n', ' ')
                labels.append(label)

        if data:
            bp = ax.boxplot(data, labels=labels, showfliers=False,
                           patch_artist=True)
            test_cols = plt.cm.Set2(np.linspace(0, 1, len(TESTS)))
            for patch, color in zip(bp['boxes'], test_cols[:len(data)]):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

        ds, de = ANALYSIS_DOMAINS[dname]
        ax.set_title(f'{dname} (FL {ds}-{de})', fontsize=11, fontweight='bold')
        ax.set_ylabel('RMSD (Å)', fontsize=9)
        ax.tick_params(axis='x', rotation=45, labelsize=6)

    plt.suptitle('KH Domain RMSD to AF2 Reference per Restraint Configuration',
                 fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(PLOT_DIR / 'kh_rmsd_boxplots.png', dpi=200, bbox_inches='tight')
    fig.savefig(PLOT_DIR / 'kh_rmsd_boxplots.pdf', bbox_inches='tight')
    plt.close()

    # ════════════════════════════════════════════════════════════
    # Plot 4: Best restraint recommendation per domain
    # ════════════════════════════════════════════════════════════
    print("\nGenerating recommendation plot...")

    fig, ax = plt.subplots(figsize=(18, 6))

    best_per_domain = {}
    for j, dname in enumerate(DOMAIN_ORDER):
        # Score = 0.5*(1 - |rg_ratio - 1|) + 0.5*cmap_corr
        best_score = -1
        best_test = None
        for t in TESTS:
            if t not in results:
                continue
            dm = results[t]['domain_metrics'].get(dname, {})
            rg_r = dm.get('rg_ratio', 2.0)
            cmap = dm.get('cmap_corr', 0)
            score = 0.5 * (1 - min(abs(rg_r - 1.0), 1.0)) + 0.5 * cmap
            if score > best_score:
                best_score = score
                best_test = t
        best_per_domain[dname] = (best_test, best_score)

        color = DOMAIN_COLORS.get(dname, '#666')
        ax.bar(j, best_score, color=color, alpha=0.8, edgecolor='black', linewidth=0.5)
        ax.text(j, best_score + 0.01,
                best_test.replace('uniform_', 'u').replace('baseline_', 'bl_'),
                ha='center', va='bottom', fontsize=6, rotation=45)

    ax.set_xticks(range(len(DOMAIN_ORDER)))
    ax.set_xticklabels(DOMAIN_ORDER, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Domain Score (Rg match + CMap fidelity)', fontsize=12)
    ax.set_title('Best Restraint Configuration Per Domain', fontsize=14)
    ax.set_ylim(0, 1.1)

    plt.tight_layout()
    fig.savefig(PLOT_DIR / 'best_per_domain.png', dpi=200, bbox_inches='tight')
    fig.savefig(PLOT_DIR / 'best_per_domain.pdf', bbox_inches='tight')
    plt.close()

    # Print recommendation table
    print("\n" + "=" * 80)
    print("RECOMMENDED RESTRAINT PER DOMAIN")
    print("=" * 80)
    print(f"{'Domain':<8} {'Best Config':<25} {'Score':>6} {'Rg Ratio':>9} {'CMap':>6}")
    print("-" * 60)
    for dname in DOMAIN_ORDER:
        best_test, score = best_per_domain[dname]
        dm = results[best_test]['domain_metrics'].get(dname, {})
        print(f"{dname:<8} {best_test:<25} {score:>6.3f} "
              f"{dm.get('rg_ratio', 0):>9.3f} {dm.get('cmap_corr', 0):>6.3f}")

    # ════════════════════════════════════════════════════════════
    # PyMOL session: KH domains from each test as states
    # ════════════════════════════════════════════════════════════
    print("\nGenerating PyMOL session...")
    generate_pymol_session(trajs, af2_xyz)

    print(f"\nAll plots saved to {PLOT_DIR}/")
    print("Done!")


def generate_pymol_session(trajs, af2_xyz):
    """Generate PyMOL PSE with KH domains from each test as separate states."""
    try:
        import pymol
        from pymol import cmd
    except ImportError:
        print("  PyMOL not available — writing PDB files instead")
        write_kh_pdbs(trajs, af2_xyz)
        return

    pymol.finish_launching(['pymol', '-cq'])
    cmd.reinitialize()

    # Color scheme for tests
    test_rgb = {
        'baseline_current': [0.6, 0.6, 0.6],
        'tight_uniform_700': [0.9, 0.2, 0.2],
        'tight_uniform_350': [1.0, 0.5, 0.5],
        'medium_uniform_700': [0.2, 0.5, 0.9],
        'medium_uniform_350': [0.5, 0.7, 1.0],
        'loose_uniform_700': [0.2, 0.8, 0.3],
        'loose_uniform_350': [0.5, 0.9, 0.6],
    }

    # For each KH domain, create an object with test states
    for dname in KH_DOMAINS:
        ds, de = ANALYSIS_DOMAINS[dname]
        n_res = de - ds + 1
        obj_name = f'KH_{dname}'

        state_idx = 1
        for t_name, traj in trajs.items():
            if traj is None:
                continue
            # Use the middle frame (most equilibrated)
            mid = traj.n_frames // 2
            # Extract domain CA coords
            domain_idx = list(range(ds - 1, min(de, traj.n_residues)))
            if len(domain_idx) != n_res:
                continue

            coords = traj.xyz[mid, domain_idx, :] * 10  # nm -> Å

            # Create pseudoatom-based structure
            if state_idx == 1:
                # First state: create the object
                for i, (x, y, z) in enumerate(coords):
                    cmd.pseudoatom(obj_name, pos=[float(x), float(y), float(z)],
                                   resv=ds + i, chain='A', name='CA',
                                   state=state_idx)
            else:
                for i, (x, y, z) in enumerate(coords):
                    cmd.pseudoatom(obj_name, pos=[float(x), float(y), float(z)],
                                   resv=ds + i, chain='A', name='CA',
                                   state=state_idx)

            # Name the state
            cmd.set('state', state_idx, obj_name)
            state_idx += 1

        # Style
        cmd.show('spheres', obj_name)
        cmd.set('sphere_scale', 0.3, obj_name)

    # Save
    pse_path = str(PLOT_DIR / 'kh_restraint_comparison.pse')
    cmd.save(pse_path)
    cmd.quit()
    print(f"  PyMOL session saved: {pse_path}")


def write_kh_pdbs(trajs, af2_xyz):
    """Write KH domain PDBs for manual PyMOL loading (fallback)."""
    pdb_dir = PLOT_DIR / 'kh_pdbs'
    pdb_dir.mkdir(exist_ok=True)

    for dname in KH_DOMAINS:
        ds, de = ANALYSIS_DOMAINS[dname]

        for t_name, traj in trajs.items():
            if traj is None:
                continue

            # Use middle frame
            mid = traj.n_frames // 2
            domain_idx = list(range(ds - 1, min(de, traj.n_residues)))
            coords = traj.xyz[mid, domain_idx, :] * 10  # nm -> Å

            pdb_path = pdb_dir / f'{dname}_{t_name}.pdb'
            with open(pdb_path, 'w') as f:
                for i, (x, y, z) in enumerate(coords):
                    f.write(f"ATOM  {i+1:>5}  CA  ALA A{ds+i:>4}    "
                            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00\n")
                f.write("END\n")

        # Also write AF2 reference
        domain_idx = list(range(ds - 1, de))
        coords = af2_xyz[domain_idx] * 10  # nm -> Å
        pdb_path = pdb_dir / f'{dname}_AF2_reference.pdb'
        with open(pdb_path, 'w') as f:
            for i, (x, y, z) in enumerate(coords):
                f.write(f"ATOM  {i+1:>5}  CA  ALA A{ds+i:>4}    "
                        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00\n")
            f.write("END\n")

    print(f"  KH domain PDBs written to {pdb_dir}/")

    # Write a PyMOL loading script
    script_path = PLOT_DIR / 'load_kh_comparison.pml'
    with open(script_path, 'w') as f:
        f.write("# Load KH domain comparison PDBs\n")
        f.write("# Run: pymol load_kh_comparison.pml\n\n")

        test_colors = {
            'baseline_current': 'gray60',
            'tight_uniform_700': 'red',
            'tight_uniform_350': 'salmon',
            'medium_uniform_700': 'blue',
            'medium_uniform_350': 'lightblue',
            'loose_uniform_700': 'green',
            'loose_uniform_350': 'palegreen',
            'AF2_reference': 'yellow',
        }

        for dname in KH_DOMAINS:
            f.write(f"\n# {dname}\n")
            # Load AF2 reference first
            f.write(f"load kh_pdbs/{dname}_AF2_reference.pdb, {dname}_AF2\n")
            f.write(f"color yellow, {dname}_AF2\n")

            for t_name in TESTS:
                pdb_file = f'kh_pdbs/{dname}_{t_name}.pdb'
                obj = f'{dname}_{t_name.replace("uniform_", "u")}'
                color = test_colors.get(t_name, 'white')
                f.write(f"load {pdb_file}, {obj}\n")
                f.write(f"color {color}, {obj}\n")
                f.write(f"align {obj}, {dname}_AF2\n")

            f.write(f"\n# Group {dname}\n")
            group_members = f"{dname}_AF2 "
            group_members += " ".join(
                f'{dname}_{t.replace("uniform_", "u")}' for t in TESTS)
            f.write(f"group {dname}_group, {group_members}\n")

        f.write("\n# Display settings\n")
        f.write("show cartoon, all\n")
        f.write("set cartoon_trace_atoms, 1\n")
        f.write("set cartoon_tube_radius, 0.2\n")
        f.write("bg_color white\n")
        f.write("set ray_opaque_background, 1\n")

    print(f"  PyMOL script: {script_path}")


if __name__ == '__main__':
    main()
