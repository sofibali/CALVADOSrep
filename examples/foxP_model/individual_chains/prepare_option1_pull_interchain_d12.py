#!/usr/bin/env python3
"""
Option 1 + interchain restraints for domain1 & domain2.

This script combines the domain3 pulling force (from option1_pull_force)
with interchain harmonic restraints between FOXP4 and FOXP1 for domain1
and domain2 only. Domain3 is NOT restrained inter-chain since it is
being pulled to target positions.

Steps:
1. Extract chain PDBs for CALVADOS domain restraints
2. Create initial positions PDB (complex centered in box)
3. Compute domain3 target positions (translated + Kabsch-rotated)
4. Compute interchain restraints for domain1 and domain2 from complex PDB
5. Generate config/components YAML and custom run.py with pulling force
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
sysname = 'option1_pull_interchain_d12'

# Chain mapping: CALVADOS name -> PDB chain ID in COMPLEX_WITHDNA.pdb
chain_mapping = {
    'FOXP4': 'A',
    'FOXP1': 'B',
}

# Domain definitions (1-based, inclusive)
domains = {
    'FOXP4': [(118, 199), (299, 375), (451, 551)],
    'FOXP1': [(118, 199), (299, 375), (451, 551)],
}

# Domain3 pulling parameters
DOMAIN3_SEPARATION = 5.0       # nm between domain3 COMs
SEPARATION_DIRECTION = np.array([1.0, 0.0, 0.0])  # Along x-axis
PULL_K = 20.0                  # kJ/mol/nm^2 - pulling force constant

# Restraint parameters
INTRA_DOMAIN_K = 2000.0        # kJ/mol/nm^2 - intra-domain restraints
INTER_CHAIN_K = 700.0          # kJ/mol/nm^2 - inter-chain restraints for d1 & d2
CUTOFF_RESTR = 0.9             # nm - cutoff for interchain restraint pairs

# Simulation parameters
BOX_SIZE = 50  # nm
N_SAVE = 500
N_FRAMES = 1000

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def parse_pdb_ca_atoms(pdb_file):
    """Parse CA atom coordinates from PDB file."""
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
    """Parse residue names for CA atoms."""
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
    """Extract a single chain from complex PDB."""
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
    Both P and Q must be centered at origin.
    Returns R such that Q @ R.T ~ P.
    """
    H = Q.T @ P
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ sign_matrix @ U.T
    return R


def create_initial_positions_pdb(complex_pdb, chain_mapping, box_size, output_pdb):
    """Create a CG PDB with chains positioned in the box (standard centering)."""
    ca_atoms = parse_pdb_ca_atoms(complex_pdb)
    resnames = parse_pdb_resnames(complex_pdb)

    protein_chains = list(chain_mapping.values())
    all_coords = [c for (ch, _), c in ca_atoms.items() if ch in protein_chains]
    complex_com = np.mean(all_coords, axis=0)
    box_center = np.array([box_size * 10 / 2.0] * 3)

    chain_letters = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')

    with open(output_pdb, 'w') as f:
        atom_num = 1
        chain_idx = 0

        for name in chain_mapping:
            chain_id = chain_mapping[name]
            pdb_chain_letter = chain_letters[chain_idx]
            res_num = 1
            chain_cas = [(rnum, ca_atoms[(chain_id, rnum)])
                         for (c, rnum) in ca_atoms if c == chain_id]
            chain_cas.sort(key=lambda x: x[0])

            for rnum, coord in chain_cas:
                new_coord = coord - complex_com + box_center
                rn = resnames.get((chain_id, rnum), 'ALA')
                f.write(f"ATOM  {atom_num:5d}  CA  {rn:3s} {pdb_chain_letter}{res_num:4d}    "
                        f"{new_coord[0]:8.3f}{new_coord[1]:8.3f}{new_coord[2]:8.3f}"
                        f"  1.00  0.00           C\n")
                atom_num += 1
                res_num += 1

            f.write("TER\n")
            chain_idx += 1

        f.write("END\n")

    print(f"  Created initial positions PDB with {atom_num-1} atoms -> {output_pdb}")
    return complex_com, box_center


def compute_domain3_targets(complex_pdb, chain_mapping, domains, box_size,
                            domain3_separation, separation_direction):
    """Compute target positions for domain3 residues of both chains."""
    ca_atoms = parse_pdb_ca_atoms(complex_pdb)

    protein_chains = list(chain_mapping.values())
    all_coords = [c for (ch, _), c in ca_atoms.items() if ch in protein_chains]
    complex_com = np.mean(all_coords, axis=0)
    box_center_ang = np.array([box_size * 10 / 2.0] * 3)

    sep_dir = separation_direction / np.linalg.norm(separation_direction)
    sep_ang = domain3_separation * 10.0

    target_coms = {
        'FOXP4': box_center_ang + sep_dir * sep_ang / 2.0,
        'FOXP1': box_center_ang - sep_dir * sep_ang / 2.0,
    }

    domain3_centered = {}
    for name in chain_mapping:
        chain_id = chain_mapping[name]
        d3_start, d3_end = domains[name][2]
        d3_coords = []
        for rnum in range(d3_start, d3_end + 1):
            key = (chain_id, rnum)
            if key in ca_atoms:
                d3_coords.append(ca_atoms[key] - complex_com + box_center_ang)
        d3_coords = np.array(d3_coords)
        d3_com = np.mean(d3_coords, axis=0)
        domain3_centered[name] = d3_coords - d3_com

    # Kabsch rotation: align FOXP1 domain3 to FOXP4 domain3
    P = domain3_centered['FOXP4']
    Q = domain3_centered['FOXP1']

    if len(P) == len(Q):
        R = kabsch_rotation(P, Q)
        Q_aligned = Q @ R.T
        rmsd = np.sqrt(np.mean(np.sum((P - Q_aligned)**2, axis=1)))
        print(f"  Kabsch alignment RMSD: {rmsd/10:.3f} nm")
    else:
        R = np.eye(3)
        print(f"  WARNING: Domain3 sizes differ, no rotation applied")

    # Compute absolute bead offsets for each chain
    chain_sizes = {}
    bead_offset = 0
    for name in chain_mapping:
        chain_id = chain_mapping[name]
        n_beads = sum(1 for (c, _) in ca_atoms if c == chain_id)
        chain_sizes[name] = (bead_offset, n_beads)
        bead_offset += n_beads

    # Build target list
    targets = []
    for name in chain_mapping:
        chain_id = chain_mapping[name]
        d3_start, d3_end = domains[name][2]
        chain_offset, chain_nbeads = chain_sizes[name]
        chain_resnums = sorted([rnum for (c, rnum) in ca_atoms if c == chain_id])
        target_com = target_coms[name]

        d3_idx = 0
        for local_idx, rnum in enumerate(chain_resnums):
            if d3_start <= rnum <= d3_end:
                bead_index = chain_offset + local_idx
                if name == 'FOXP4':
                    target_pos = domain3_centered['FOXP4'][d3_idx] + target_com
                else:
                    target_pos = (domain3_centered['FOXP1'][d3_idx] @ R.T) + target_com
                target_nm = target_pos / 10.0
                targets.append((bead_index, target_nm[0], target_nm[1], target_nm[2]))
                d3_idx += 1

    print(f"  Computed {len(targets)} target positions")
    print(f"  FOXP4 domain3: beads {targets[0][0]}-{targets[len(P)-1][0]}")
    print(f"  FOXP1 domain3: beads {targets[len(P)][0]}-{targets[-1][0]}")

    return targets


def compute_interchain_restraints(complex_pdb, chain_mapping, domains,
                                  cutoff_nm, k_force, domain_indices=(0, 1)):
    """
    Compute interchain restraint pairs between FOXP4 and FOXP1 for specified domains.

    Only considers residues within domain1 and domain2 (by default).
    Pairs within cutoff distance are added as restraints.

    Args:
        domain_indices: tuple of domain indices to include (0=domain1, 1=domain2)

    Returns:
        list of (name1, resid1, name2, resid2, dist_nm, k)
    """
    ca_atoms = parse_pdb_ca_atoms(complex_pdb)
    names = list(chain_mapping.keys())
    name1, name2 = names[0], names[1]
    chain1, chain2 = chain_mapping[name1], chain_mapping[name2]

    # Build sets of residues in the selected domains
    domain_residues_1 = set()
    domain_residues_2 = set()
    for di in domain_indices:
        d_start, d_end = domains[name1][di]
        domain_residues_1.update(range(d_start, d_end + 1))
        d_start, d_end = domains[name2][di]
        domain_residues_2.update(range(d_start, d_end + 1))

    # Get CA coords for domain residues
    coords1 = {rnum: ca_atoms[(chain1, rnum)]
                for rnum in domain_residues_1 if (chain1, rnum) in ca_atoms}
    coords2 = {rnum: ca_atoms[(chain2, rnum)]
                for rnum in domain_residues_2 if (chain2, rnum) in ca_atoms}

    cutoff_ang = cutoff_nm * 10.0
    restraints = []

    for r1, c1 in coords1.items():
        for r2, c2 in coords2.items():
            dist_ang = np.linalg.norm(c1 - c2)
            if dist_ang <= cutoff_ang:
                dist_nm = dist_ang / 10.0
                restraints.append((name1, r1, name2, r2, dist_nm, k_force))

    domain_names = [f"domain{di+1}" for di in domain_indices]
    print(f"  Domains included: {', '.join(domain_names)}")
    print(f"  Cutoff: {cutoff_nm} nm")
    print(f"  {name1} domain residues: {len(coords1)}")
    print(f"  {name2} domain residues: {len(coords2)}")
    print(f"  Interchain restraint pairs found: {len(restraints)}")

    return restraints


def write_interchain_restraints_file(restraints, output_file):
    """Write interchain restraints in CALVADOS custom restraints format."""
    with open(output_file, 'w') as f:
        for name1, r1, name2, r2, dist_nm, k in restraints:
            f.write(f"{name1} 1 {r1} | {name2} 1 {r2} | {dist_nm:.4f} {k}\n")
    print(f"  Saved {len(restraints)} interchain restraints -> {output_file}")


def write_targets_file(targets, output_file):
    """Write target positions to a text file."""
    with open(output_file, 'w') as f:
        f.write("# bead_index target_x_nm target_y_nm target_z_nm\n")
        for bead_idx, x, y, z in targets:
            f.write(f"{bead_idx} {x:.6f} {y:.6f} {z:.6f}\n")
    print(f"  Saved {len(targets)} target positions -> {output_file}")


def write_custom_run_py(output_path, pull_k, targets_file):
    """Write a custom run.py that adds the domain3 pulling force."""
    run_code = f'''#!/usr/bin/env python3
"""
Custom run script with domain3 pulling force + interchain restraints for d1 & d2.

This script:
1. Builds the CALVADOS system normally (bonds, intra-domain restraints,
   non-bonded forces, AND interchain restraints via custom_restraints)
2. Adds a CustomExternalForce that pulls domain3 residues toward target positions
3. Runs the simulation
"""

import os
import numpy as np
import openmm
from openmm import unit
from yaml import safe_load
from calvados import sim

# Pulling force parameters
PULL_K = {pull_k}  # kJ/mol/nm^2
TARGETS_FILE = '{targets_file}'


def add_domain3_pulling_force(mysim, pull_k, targets_file):
    """
    Add a CustomExternalForce that pulls domain3 residues toward target positions.
    """
    targets = []
    with open(f'{{mysim.path}}/{{targets_file}}', 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            bead_idx = int(parts[0])
            x0, y0, z0 = float(parts[1]), float(parts[2]), float(parts[3])
            targets.append((bead_idx, x0, y0, z0))

    pull_expr = 'k_pull * periodicdistance(x, y, z, x0, y0, z0)^2'
    pull_force = openmm.CustomExternalForce(pull_expr)
    pull_force.addGlobalParameter('k_pull',
                                   pull_k * unit.kilojoules_per_mole / unit.nanometer**2)
    pull_force.addPerParticleParameter('x0')
    pull_force.addPerParticleParameter('y0')
    pull_force.addPerParticleParameter('z0')

    for bead_idx, x0, y0, z0 in targets:
        pull_force.addParticle(bead_idx, [
            x0 * unit.nanometer,
            y0 * unit.nanometer,
            z0 * unit.nanometer
        ])

    mysim.system.addForce(pull_force)

    print(f"\\nAdded domain3 pulling force:")
    print(f"  k_pull = {{pull_k}} kJ/mol/nm^2")
    print(f"  {{len(targets)}} particles with target positions")
    print(f"  Force group: {{pull_force.getForceGroup()}}")


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--path', nargs='?', default='.', const='.', type=str)
    parser.add_argument('--config', nargs='?', default='config.yaml', const='config.yaml', type=str)
    parser.add_argument('--components', nargs='?', default='components.yaml', const='components.yaml', type=str)

    args = parser.parse_args()

    path = args.path
    fconfig = args.config
    fcomponents = args.components

    with open(f'{{path}}/{{fconfig}}', 'r') as stream:
        config = safe_load(stream)

    with open(f'{{path}}/{{fcomponents}}', 'r') as stream:
        components = safe_load(stream)

    # Build the system normally (includes interchain restraints via custom_restraints)
    mysim = sim.Sim(path, config, components)
    mysim.build_system()

    # Add the domain3 pulling force AFTER system build, BEFORE simulation
    add_domain3_pulling_force(mysim, PULL_K, TARGETS_FILE)

    # Run the simulation
    mysim.simulate()
'''
    run_file = f'{output_path}/run.py'
    with open(run_file, 'w') as f:
        f.write(run_code)
    print(f"  Wrote custom run.py -> {run_file}")


# ============================================================================
# MAIN SCRIPT
# ============================================================================

if __name__ == '__main__':

    # Paths
    complex_pdb = f'{input_dir}/COMPLEX_WITHDNA.pdb'
    output_dir = cwd
    sim_path = f'{output_dir}/{sysname}'
    initial_pdb = f'{sim_path}/restart.pdb'
    targets_file = 'domain3_targets.txt'
    restraints_file = f'{output_dir}/input/inter_chain_restraints_d12.txt'

    # Create directories
    subprocess.run(f'mkdir -p {sim_path}', shell=True)
    subprocess.run(f'mkdir -p {output_dir}/input', shell=True)
    subprocess.run(f'mkdir -p {output_dir}/data', shell=True)

    # Copy required input files
    subprocess.run(f'cp {input_dir}/residues_CALVADOS3.csv {output_dir}/input/', shell=True)
    subprocess.run(f'cp {input_dir}/proteins.fasta {output_dir}/input/', shell=True)

    # Write domains.yaml with FOXP1 entry
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
    print("OPTION 1 + INTERCHAIN RESTRAINTS (domain1 & domain2)")
    print("=" * 60)
    print(f"  Domain3 separation: {DOMAIN3_SEPARATION} nm")
    print(f"  Separation direction: {SEPARATION_DIRECTION}")
    print(f"  Pull force constant: {PULL_K} kJ/mol/nm^2")
    print(f"  Interchain restraint k: {INTER_CHAIN_K} kJ/mol/nm^2")
    print(f"  Interchain cutoff: {CUTOFF_RESTR} nm")
    print(f"  Box size: {BOX_SIZE} nm")

    print("\n" + "=" * 60)
    print("Step 1: Extract individual chain PDBs")
    print("=" * 60)

    for name, pdb_chain in chain_mapping.items():
        output_pdb = f'{output_dir}/input/{name}.pdb'
        extract_chain_pdb(complex_pdb, pdb_chain, output_pdb)

    print("\n" + "=" * 60)
    print("Step 2: Create initial positions PDB (standard centering)")
    print("=" * 60)

    create_initial_positions_pdb(complex_pdb, chain_mapping, BOX_SIZE, initial_pdb)

    print("\n" + "=" * 60)
    print("Step 3: Compute domain3 target positions")
    print("=" * 60)

    targets = compute_domain3_targets(
        complex_pdb, chain_mapping, domains, BOX_SIZE,
        DOMAIN3_SEPARATION, SEPARATION_DIRECTION
    )
    write_targets_file(targets, f'{sim_path}/{targets_file}')

    print("\n" + "=" * 60)
    print("Step 4: Compute interchain restraints (domain1 & domain2 only)")
    print("=" * 60)

    interchain_restraints = compute_interchain_restraints(
        complex_pdb, chain_mapping, domains,
        cutoff_nm=CUTOFF_RESTR, k_force=INTER_CHAIN_K,
        domain_indices=(0, 1)  # domain1 and domain2 only, NOT domain3
    )
    write_interchain_restraints_file(interchain_restraints, restraints_file)

    print("\n" + "=" * 60)
    print("Step 5: Creating CALVADOS configuration files")
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
        # Interchain restraints
        custom_restraints=True,
        custom_restraint_type='harmonic',
        fcustom_restraints=restraints_file,
    )

    config.write(sim_path, name='config.yaml', analyses='')

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
        cutoff_restr=CUTOFF_RESTR,
        use_com=True,
        colabfold=1,
        ext_restraint=True,
    )

    components.add(name='FOXP4', restraint=True)
    components.add(name='FOXP1', restraint=True)

    components.write(sim_path, name='components.yaml')

    print("\n" + "=" * 60)
    print("Step 6: Writing custom run.py with pulling force")
    print("=" * 60)

    write_custom_run_py(sim_path, PULL_K, targets_file)

    print("\n" + "=" * 60)
    print("SETUP COMPLETE")
    print("=" * 60)
    print(f"\nFiles written to: {sim_path}/")
    print(f"  - config.yaml")
    print(f"  - components.yaml")
    print(f"  - restart.pdb (standard centered positions)")
    print(f"  - {targets_file} (domain3 target positions in nm)")
    print(f"  - run.py (custom: pulling force + interchain restraints)")
    print(f"\nInterchain restraints file:")
    print(f"  - {restraints_file}")
    print(f"  - {len(interchain_restraints)} pairs (domain1 & domain2 only)")
    print(f"  - k = {INTER_CHAIN_K} kJ/mol/nm^2")
    print(f"\nTo run the simulation:")
    print(f"  cd {sim_path}")
    print(f"  python run.py")
