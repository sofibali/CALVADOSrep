#!/usr/bin/env python3
"""
Option 3: Generate initial structure with domain3 regions repositioned.

This script:
1. Creates chain PDB files for CALVADOS to read domain restraints
2. Creates an initial positions PDB where:
   - Domain3 of each chain is translated and rotated to target positions
   - Linker residues between domain2 and domain3 are interpolated smoothly
   - Residues before domain2 end and after domain3 (C-terminal tail) are kept
3. Uses standard CALVADOS run.py (no custom forces needed)

The approach:
- FOXP4's domain3 orientation is used as the reference
- FOXP1's domain3 is rotated (Kabsch algorithm) to match FOXP4's orientation
- Both domain3 COMs are placed at target positions (separated along a direction)
- Linker residues (between domain2 and domain3) are interpolated along a
  parabolic arc to maintain physically viable bond lengths (~0.38 nm)
- C-terminal residues after domain3 move rigidly with domain3
"""

import os
import sys
import numpy as np
import subprocess

# Add CALVADOS to path
calvados_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, calvados_root)

from calvados.cfg import Config, Components

# ============================================================================
# CONFIGURATION
# ============================================================================

cwd = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(cwd)
input_dir = f'{parent_dir}/input'

# System name
sysname = 'option3_reposition'

# Chain mapping: CALVADOS name -> PDB chain ID
chain_mapping = {
    'FOXP4': 'A',
    'FOXP1': 'B',
}

# Domain definitions (1-based, inclusive)
domains = {
    'FOXP4': [(118, 199), (299, 375), (451, 551)],
    'FOXP1': [(118, 199), (299, 375), (451, 551)],
}

# Domain3 positioning parameters
DOMAIN3_SEPARATION = 5.0       # nm between domain3 COMs
SEPARATION_DIRECTION = np.array([1.0, 0.0, 0.0])  # Along x-axis
BOND_LENGTH = 3.8              # Angstrom (0.38 nm) - typical CA-CA distance

# Restraint parameters
INTRA_DOMAIN_K = 2000.0   # kJ/mol/nm^2

# Simulation parameters
BOX_SIZE = 50  # nm
N_SAVE = 500
N_FRAMES = 1000

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def parse_pdb_ca_atoms(pdb_file):
    """Parse CA atom coordinates from PDB file.
    Returns: {(chain_id, resnum): np.array([x, y, z])} in Angstrom.
    """
    ca_atoms = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line:
                chain = line[21]
                resnum = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ca_atoms[(chain, resnum)] = np.array([x, y, z])
    return ca_atoms


def parse_pdb_resnames(pdb_file):
    """Parse residue names for CA atoms.
    Returns: {(chain_id, resnum): resname_3letter}
    """
    resnames = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line:
                chain = line[21]
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                resnames[(chain, resnum)] = resname
    return resnames


def extract_chain_pdb(complex_pdb, chain_id, output_pdb):
    """Extract a single chain from complex PDB (original coordinates)."""
    with open(complex_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM') and line[21] == chain_id:
                f_out.write(line)
        f_out.write("TER\n")
        f_out.write("END\n")
    print(f"  Extracted chain {chain_id} -> {output_pdb}")


def kabsch_rotation(P, Q):
    """
    Compute rotation matrix R that best aligns Q onto P.
    Both P and Q must be centered at origin (COM subtracted).
    Returns R such that Q @ R.T ~ P (minimizes RMSD).
    """
    H = Q.T @ P
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ sign_matrix @ U.T
    return R


def interpolate_linker(start_pos, end_pos, n_residues, bond_length=3.8):
    """
    Create physically viable linker positions between two endpoints.
    Uses a parabolic arc with arc-length parameterization to ensure
    uniform bond lengths (~bond_length between consecutive residues).

    All coordinates in Angstrom.

    Args:
        start_pos: position of last domain2 residue (Angstrom)
        end_pos: position of first domain3 residue (Angstrom)
        n_residues: number of linker residues to place
        bond_length: target CA-CA bond length (Angstrom)

    Returns:
        Array of shape (n_residues, 3) with linker positions (Angstrom)
    """
    if n_residues == 0:
        return np.zeros((0, 3))

    direction = end_pos - start_pos
    distance = np.linalg.norm(direction)

    if distance < 1e-6:
        direction = np.array([1.0, 0.0, 0.0])
        distance = 0.01

    direction_unit = direction / distance

    # n_residues linker residues => n_residues+1 segments
    n_segments = n_residues + 1
    total_desired_length = n_segments * bond_length

    if distance >= total_desired_length:
        print(f"    WARNING: Linker stretched! Distance {distance/10:.2f} nm "
              f"> max extension {total_desired_length/10:.2f} nm")
        positions = np.zeros((n_residues, 3))
        for i in range(n_residues):
            t = (i + 1) / n_segments
            positions[i] = start_pos + t * direction
        return positions

    # Find perpendicular direction for the arc
    if abs(np.dot(direction_unit, [1, 0, 0])) < 0.9:
        perp = np.cross(direction_unit, [1, 0, 0])
    else:
        perp = np.cross(direction_unit, [0, 1, 0])
    perp = perp / np.linalg.norm(perp)

    # Helper: compute a point on the parabolic arc at parameter t in [0, 1]
    def arc_point(t, h):
        sag = 4 * h * t * (1 - t)
        return start_pos + t * direction + sag * perp

    # Helper: compute total arc length using a fine discretization
    n_fine = 2000  # fine resolution for arc-length computation
    def compute_arc_data(h):
        """Compute fine arc points, cumulative arc lengths."""
        ts = np.linspace(0, 1, n_fine + 1)
        points = np.array([arc_point(t, h) for t in ts])
        diffs = np.diff(points, axis=0)
        seg_lengths = np.sqrt(np.sum(diffs**2, axis=1))
        cum_lengths = np.concatenate([[0], np.cumsum(seg_lengths)])
        return ts, points, cum_lengths

    # Binary search for arc height h that gives the desired total path length
    h_low, h_high = 0.0, total_desired_length * 2.0
    for _ in range(100):
        h = (h_low + h_high) / 2
        _, _, cum = compute_arc_data(h)
        total_arc = cum[-1]
        if total_arc < total_desired_length:
            h_low = h
        else:
            h_high = h
        if abs(total_arc - total_desired_length) < 0.01:
            break

    # Sample positions at equal arc-length intervals
    ts_fine, points_fine, cum_lengths = compute_arc_data(h)
    total_arc = cum_lengths[-1]

    # Target arc-length positions: equally spaced from 0 to total_arc
    # Position 0 = start_pos, position total_arc = end_pos
    # Linker residues at positions: bond_length, 2*bond_length, ..., n_residues*bond_length
    positions = np.zeros((n_residues, 3))
    for i in range(n_residues):
        target_s = (i + 1) * bond_length
        # Find the fine index where cumulative arc length >= target_s
        idx = np.searchsorted(cum_lengths, target_s)
        if idx >= len(cum_lengths):
            idx = len(cum_lengths) - 1
        # Linear interpolation between fine points
        if idx == 0:
            positions[i] = points_fine[0]
        else:
            s0 = cum_lengths[idx - 1]
            s1 = cum_lengths[idx]
            frac = (target_s - s0) / (s1 - s0) if s1 > s0 else 0.0
            positions[i] = points_fine[idx - 1] + frac * (points_fine[idx] - points_fine[idx - 1])

    # Print diagnostics
    all_points = [start_pos] + [positions[i] for i in range(n_residues)] + [end_pos]
    bls = [np.linalg.norm(all_points[i+1] - all_points[i])
           for i in range(len(all_points)-1)]
    print(f"    Linker: {n_residues} residues, endpoint distance {distance/10:.2f} nm")
    print(f"    Arc height: {h/10:.2f} nm")
    print(f"    Bond lengths: min={min(bls)/10:.3f} nm, "
          f"max={max(bls)/10:.3f} nm, "
          f"mean={np.mean(bls)/10:.3f} nm")

    return positions


def create_initial_positions_pdb_repositioned(complex_pdb, chain_mapping, domains,
                                               box_size, output_pdb,
                                               domain3_separation, separation_direction):
    """
    Create a CG PDB with domain3 regions repositioned.

    For each chain:
    1. Residues before linker: keep at centered-in-box positions
    2. Domain3 residues: translate (and rotate for FOXP1) to target position
    3. C-terminal residues after domain3: move rigidly with domain3
    4. Linker residues (domain2_end+1 to domain3_start-1): interpolated arc

    Args:
        complex_pdb: path to the complex PDB file
        chain_mapping: dict {name: chain_id}
        domains: dict {name: [(start, end), ...]} with domain ranges
        box_size: box size in nm
        output_pdb: output PDB path
        domain3_separation: separation between domain3 COMs in nm
        separation_direction: unit vector for separation direction
    """
    print("\n--- Computing domain3 repositioning ---")

    # Parse all CA atoms and residue names
    ca_atoms = parse_pdb_ca_atoms(complex_pdb)
    resnames = parse_pdb_resnames(complex_pdb)

    # Calculate complex COM (protein chains only)
    protein_chains = list(chain_mapping.values())
    all_coords = [c for (ch, _), c in ca_atoms.items() if ch in protein_chains]
    complex_com = np.mean(all_coords, axis=0)

    # Box center in Angstrom
    box_center = np.array([box_size * 10 / 2.0] * 3)

    # Separation in Angstrom
    sep_dir = separation_direction / np.linalg.norm(separation_direction)
    sep_angstrom = domain3_separation * 10.0

    # Target COM positions for domain3
    target_coms = {
        'FOXP4': box_center + sep_dir * sep_angstrom / 2.0,
        'FOXP1': box_center - sep_dir * sep_angstrom / 2.0,
    }

    # --- Step 1: Get domain3 centered coordinates for both chains ---
    domain3_centered = {}  # {name: centered coords (COM subtracted)}
    domain3_coms_centered = {}  # {name: COM in centered-in-box frame}

    for name in chain_mapping:
        chain_id = chain_mapping[name]
        d3_start, d3_end = domains[name][2]  # domain3 range

        # Get domain3 coords, center in box
        d3_coords = []
        for rnum in range(d3_start, d3_end + 1):
            key = (chain_id, rnum)
            if key in ca_atoms:
                d3_coords.append(ca_atoms[key] - complex_com + box_center)
        d3_coords = np.array(d3_coords)
        d3_com = np.mean(d3_coords, axis=0)

        domain3_centered[name] = d3_coords - d3_com  # centered at origin
        domain3_coms_centered[name] = d3_com

        print(f"  {name} domain3 ({d3_start}-{d3_end}): {len(d3_coords)} residues, "
              f"COM = [{d3_com[0]:.1f}, {d3_com[1]:.1f}, {d3_com[2]:.1f}] A")

    # --- Step 2: Compute rotation to align FOXP1 domain3 with FOXP4 ---
    P = domain3_centered['FOXP4']  # reference (target shape)
    Q = domain3_centered['FOXP1']  # source (to be rotated)

    if len(P) == len(Q):
        R = kabsch_rotation(P, Q)
        # Verify: RMSD after alignment
        Q_aligned = Q @ R.T
        rmsd = np.sqrt(np.mean(np.sum((P - Q_aligned)**2, axis=1)))
        print(f"  Kabsch alignment RMSD: {rmsd/10:.3f} nm")
    else:
        print(f"  WARNING: Domain3 sizes differ ({len(P)} vs {len(Q)}), using principal axis alignment")
        R = np.eye(3)

    print(f"  Target FOXP4 domain3 COM: {target_coms['FOXP4']/10} nm")
    print(f"  Target FOXP1 domain3 COM: {target_coms['FOXP1']/10} nm")

    # Chain letters for PDB output (one per component)
    chain_letters = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')

    # --- Step 3: Process each chain and build new coordinates ---
    with open(output_pdb, 'w') as f:
        atom_num = 1
        chain_idx = 0

        for name in chain_mapping:
            chain_id = chain_mapping[name]
            pdb_chain_letter = chain_letters[chain_idx]
            res_num = 1  # reset residue numbering per chain
            d2_end = domains[name][1][1]     # last residue of domain2
            d3_start = domains[name][2][0]   # first residue of domain3
            d3_end = domains[name][2][1]     # last residue of domain3

            # Get all CA atoms for this chain, sorted
            chain_cas = [(rnum, ca_atoms[(chain_id, rnum)])
                         for (c, rnum) in ca_atoms if c == chain_id]
            chain_cas.sort(key=lambda x: x[0])
            resnums = [rnum for rnum, _ in chain_cas]

            # Center all coords in box
            coords = np.array([coord - complex_com + box_center for _, coord in chain_cas])

            # Find indices for domain boundaries
            d2_end_idx = None
            d3_start_idx = None
            d3_end_idx = None

            for idx, rnum in enumerate(resnums):
                if rnum == d2_end:
                    d2_end_idx = idx
                if rnum == d3_start:
                    d3_start_idx = idx
                if rnum == d3_end:
                    d3_end_idx = idx

            if d2_end_idx is None or d3_start_idx is None or d3_end_idx is None:
                print(f"  WARNING: Could not find domain boundaries for {name}, using original positions")
                # Write unchanged
                for idx, (rnum, _) in enumerate(chain_cas):
                    rn = resnames.get((chain_id, rnum), 'ALA')
                    c = coords[idx]
                    f.write(f"ATOM  {atom_num:5d}  CA  {rn:3s} {pdb_chain_letter}{res_num:4d}    "
                            f"{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}"
                            f"  1.00  0.00           C\n")
                    atom_num += 1
                    res_num += 1
                f.write("TER\n")
                chain_idx += 1
                continue

            # Indices for each region
            pre_linker_indices = list(range(0, d2_end_idx + 1))
            linker_indices = list(range(d2_end_idx + 1, d3_start_idx))
            domain3_indices = list(range(d3_start_idx, d3_end_idx + 1))
            post_domain3_indices = list(range(d3_end_idx + 1, len(resnums)))

            print(f"\n  {name}: Processing regions:")
            print(f"    Pre-linker: residues {resnums[0]}-{resnums[d2_end_idx]} "
                  f"({len(pre_linker_indices)} residues) - unchanged")
            print(f"    Linker: residues {resnums[d2_end_idx+1] if linker_indices else 'N/A'}-"
                  f"{resnums[d3_start_idx-1] if linker_indices else 'N/A'} "
                  f"({len(linker_indices)} residues) - interpolated")
            print(f"    Domain3: residues {d3_start}-{d3_end} "
                  f"({len(domain3_indices)} residues) - repositioned")
            print(f"    C-terminal: {len(post_domain3_indices)} residues - moved with domain3")

            # --- Transform domain3 ---
            d3_coords_original = coords[domain3_indices]
            d3_com_original = np.mean(d3_coords_original, axis=0)
            target_com = target_coms[name]

            if name == 'FOXP4':
                # Reference chain: translate only
                new_d3_coords = d3_coords_original - d3_com_original + target_com
            else:
                # Rotate to match FOXP4 orientation, then translate
                centered = d3_coords_original - d3_com_original
                rotated = centered @ R.T
                new_d3_coords = rotated + target_com

            # --- Transform C-terminal residues (move rigidly with domain3) ---
            if post_domain3_indices:
                post_coords_original = coords[post_domain3_indices]
                if name == 'FOXP4':
                    new_post_coords = post_coords_original - d3_com_original + target_com
                else:
                    centered = post_coords_original - d3_com_original
                    rotated = centered @ R.T
                    new_post_coords = rotated + target_com

            # --- Build new coordinate array ---
            new_coords = coords.copy()

            # Domain3
            for i, idx in enumerate(domain3_indices):
                new_coords[idx] = new_d3_coords[i]

            # C-terminal
            if post_domain3_indices:
                for i, idx in enumerate(post_domain3_indices):
                    new_coords[idx] = new_post_coords[i]

            # --- Interpolate linker ---
            if linker_indices:
                linker_start_pos = new_coords[d2_end_idx]      # end of domain2
                linker_end_pos = new_coords[d3_start_idx]       # start of domain3
                print(f"    Interpolating linker for {name}:")
                linker_positions = interpolate_linker(
                    linker_start_pos, linker_end_pos,
                    len(linker_indices), bond_length=BOND_LENGTH
                )
                for i, idx in enumerate(linker_indices):
                    new_coords[idx] = linker_positions[i]

            # --- Write all residues for this chain ---
            for idx, (rnum, _) in enumerate(chain_cas):
                rn = resnames.get((chain_id, rnum), 'ALA')
                c = new_coords[idx]
                f.write(f"ATOM  {atom_num:5d}  CA  {rn:3s} {pdb_chain_letter}{res_num:4d}    "
                        f"{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}"
                        f"  1.00  0.00           C\n")
                atom_num += 1
                res_num += 1

            f.write("TER\n")
            chain_idx += 1

        f.write("END\n")

    print(f"\nCreated repositioned initial positions PDB with {atom_num-1} atoms -> {output_pdb}")


# ============================================================================
# MAIN SCRIPT
# ============================================================================

if __name__ == '__main__':

    # Paths
    complex_pdb = f'{input_dir}/COMPLEX_WITHDNA.pdb'
    output_dir = cwd
    sim_path = f'{output_dir}/{sysname}'
    initial_pdb = f'{sim_path}/restart.pdb'

    # Create directories
    subprocess.run(f'mkdir -p {sim_path}', shell=True)
    subprocess.run(f'mkdir -p {output_dir}/input', shell=True)
    subprocess.run(f'mkdir -p {output_dir}/data', shell=True)

    # Copy required input files
    subprocess.run(f'cp {input_dir}/residues_CALVADOS3.csv {output_dir}/input/', shell=True)
    subprocess.run(f'cp {input_dir}/proteins.fasta {output_dir}/input/', shell=True)

    # Write domains.yaml with FOXP1 entry (parent input/ has FOX, not FOXP1)
    domains_yaml = f'{output_dir}/input/domains.yaml'
    with open(domains_yaml, 'w') as f:
        f.write("FOXP4:\n")
        f.write("- [118,199]\n")
        f.write("- [299,375]\n")
        f.write("- [451,551]\n")
        f.write("FOXP1:\n")
        f.write("- [118,199]\n")
        f.write("- [299,375]\n")
        f.write("- [451,551]\n")

    print("=" * 60)
    print("OPTION 3: Repositioned domain3 with interpolated linkers")
    print("=" * 60)
    print(f"  Domain3 separation: {DOMAIN3_SEPARATION} nm")
    print(f"  Separation direction: {SEPARATION_DIRECTION}")
    print(f"  Box size: {BOX_SIZE} nm")

    print("\n" + "=" * 60)
    print("Step 1: Extract individual chain PDBs (for restraint calculation)")
    print("=" * 60)

    for name, pdb_chain in chain_mapping.items():
        output_pdb = f'{output_dir}/input/{name}.pdb'
        extract_chain_pdb(complex_pdb, pdb_chain, output_pdb)

    print("\n" + "=" * 60)
    print("Step 2: Create repositioned initial positions PDB")
    print("=" * 60)

    create_initial_positions_pdb_repositioned(
        complex_pdb, chain_mapping, domains, BOX_SIZE, initial_pdb,
        DOMAIN3_SEPARATION, SEPARATION_DIRECTION
    )

    print("\n" + "=" * 60)
    print("Step 3: Creating CALVADOS configuration files")
    print("=" * 60)

    residues_file = f'{output_dir}/input/residues_CALVADOS3.csv'

    config = Config(
        sysname=sysname,
        box=[BOX_SIZE, BOX_SIZE, BOX_SIZE],
        temp=293,
        ionic=0.15,
        pH=7.0,
        topol='random',
        wfreq=N_SAVE,
        steps=N_FRAMES * N_SAVE,
        platform='CPU',
        restart='pdb',
        frestart='restart.pdb',
        verbose=True,
    )

    analyses = f'''
from calvados.analysis import calc_com_traj, calc_contact_map

chainid_dict = dict(FOXP4=0, FOXP1=1)
calc_com_traj(path="{sim_path}",sysname="{sysname}",output_path="{output_dir}/data",
              residues_file="{residues_file}",chainid_dict=chainid_dict,start=100)
calc_contact_map(path="{sim_path}",sysname="{sysname}",output_path="{output_dir}/data",
                 chainid_dict=chainid_dict)
'''

    config.write(sim_path, name='config.yaml', analyses=analyses)

    components = Components(
        molecule_type='protein',
        nmol=1,
        restraint=True,
        charge_termini='both',
        fresidues=residues_file,
        ffasta=f'{output_dir}/input/proteins.fasta',
        fdomains=f'{output_dir}/input/domains.yaml',
        pdb_folder=f'{output_dir}/input',
        restraint_type='harmonic',
        k_harmonic=INTRA_DOMAIN_K,
        cutoff_restr=0.9,
        use_com=True,
        colabfold=1,
        ext_restraint=True,
    )

    components.add(name='FOXP4', restraint=True)
    components.add(name='FOXP1', restraint=True)

    components.write(sim_path, name='components.yaml')

    print(f"\nConfiguration files written to: {sim_path}/")
    print(f"  - config.yaml")
    print(f"  - components.yaml")
    print(f"  - restart.pdb (repositioned domain3)")

    print("\n" + "=" * 60)
    print("SETUP COMPLETE")
    print("=" * 60)
    print(f"\nTo run the simulation:")
    print(f"  cd {sim_path}")
    print(f"  python run.py")
    print(f"\nTo visualize the initial structure:")
    print(f"  pymol {initial_pdb}")
    print(f"\nKey settings:")
    print(f"  - Domain3 separation: {DOMAIN3_SEPARATION} nm along {SEPARATION_DIRECTION}")
    print(f"  - Linker residues interpolated with parabolic arc")
    print(f"  - Intra-domain restraints: k = {INTRA_DOMAIN_K} kJ/mol/nm^2")
    print(f"  - NO inter-chain restraints")
