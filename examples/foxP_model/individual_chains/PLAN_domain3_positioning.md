# Plan: Domain3 Positioning Configuration

## Goal
Position domain3 regions of FOXP4 (residues 451-551) and FOXP1 (residues 450-539) in the same direction but at far distances from each other.

## Approach Overview

We need to modify the `create_initial_positions_pdb()` function to apply specific transformations to domain3 residues of both chains. The strategy involves:

1. **Identify domain3 residues** during PDB creation
2. **Apply translation vectors** to move domain3s to target locations
3. **Apply rotation** to orient domain3s in the same direction
4. **Maintain domain integrity** - keep residues within each domain together

## Implementation Steps

### Step 1: Add Domain3 Configuration Parameters

Add to the CONFIGURATION section:

```python
# Domain3 positioning parameters
POSITION_DOMAIN3 = True  # Enable/disable domain3 positioning
DOMAIN3_SEPARATION = 20.0  # nm - distance between domain3 centers
DOMAIN3_DIRECTION = np.array([1.0, 0.0, 0.0])  # Direction vector (unit vector)

# Target positions for domain3 centers (relative to box center)
# These will be calculated to ensure separation along DOMAIN3_DIRECTION
DOMAIN3_OFFSET_FOXP4 = DOMAIN3_SEPARATION / 2.0  # nm from box center
DOMAIN3_OFFSET_FOXP1 = -DOMAIN3_SEPARATION / 2.0  # nm from box center (opposite side)
```

### Step 2: Add Helper Functions

```python
def get_domain3_residues(chain_name, domains):
    """
    Get the residue range for domain3 (third domain) of a chain.

    Args:
        chain_name: Name of the chain (FOXP4 or FOXP1)
        domains: Dictionary of domain definitions

    Returns:
        Tuple (start, end) of domain3 residues, or None if not found
    """
    if chain_name in domains and len(domains[chain_name]) >= 3:
        return domains[chain_name][2]  # Third domain (index 2)
    return None


def calculate_domain_com(ca_atoms, chain_id, domain_range):
    """
    Calculate center of mass for a specific domain.

    Args:
        ca_atoms: Dictionary {(chain, resnum): coord}
        chain_id: PDB chain ID
        domain_range: Tuple (start, end) of residue numbers

    Returns:
        numpy array [x, y, z] of domain COM in Angstrom
    """
    domain_coords = []
    start, end = domain_range

    for (c, rnum), coord in ca_atoms.items():
        if c == chain_id and start <= rnum <= end:
            domain_coords.append(coord)

    if len(domain_coords) == 0:
        return None

    return np.mean(domain_coords, axis=0)


def apply_domain_translation(coords, current_com, target_com):
    """
    Translate a set of coordinates to move COM to target position.

    Args:
        coords: List or array of coordinates to translate
        current_com: Current center of mass
        target_com: Target center of mass

    Returns:
        Translated coordinates (numpy array)
    """
    translation = target_com - current_com
    return coords + translation


def align_to_direction(coords, com, target_direction):
    """
    Rotate coordinates to align the principal axis with target direction.

    This ensures domain3 regions point in the same direction.

    Args:
        coords: Nx3 array of coordinates
        com: Center of mass of the domain
        target_direction: Unit vector [x, y, z] for desired orientation

    Returns:
        Rotated coordinates
    """
    # Center coordinates at origin
    centered = coords - com

    # Calculate principal axis using PCA
    cov_matrix = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

    # Principal axis is eigenvector with largest eigenvalue
    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]
    principal_axis = principal_axis / np.linalg.norm(principal_axis)

    # Calculate rotation to align principal axis with target direction
    target_direction = target_direction / np.linalg.norm(target_direction)

    # Use Rodrigues' rotation formula
    v = np.cross(principal_axis, target_direction)
    s = np.linalg.norm(v)
    c = np.dot(principal_axis, target_direction)

    if s < 1e-6:  # Already aligned
        return coords

    v = v / s
    vx = np.array([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])

    rotation_matrix = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))

    # Apply rotation
    rotated = np.dot(centered, rotation_matrix.T)

    # Return to original COM
    return rotated + com
```

### Step 3: Modify `create_initial_positions_pdb()` Function

Replace the main loop with:

```python
def create_initial_positions_pdb_with_domain3(complex_pdb, chain_mapping, domains,
                                              box_size, output_pdb,
                                              position_domain3=False,
                                              domain3_separation=20.0,
                                              domain3_direction=np.array([1.0, 0.0, 0.0])):
    """
    Create a coarse-grained PDB with chains positioned in the box.
    Optionally positions domain3 regions at specific locations.

    If position_domain3=True:
    - Domain3 of FOXP4 and FOXP1 are moved to opposite sides along domain3_direction
    - Domain3 regions are rotated to point in the same direction
    - Rest of the chains follow naturally (may extend differently)
    """
    # Parse CA atoms from complex
    ca_atoms = parse_pdb_ca_atoms(complex_pdb)

    # Calculate complex COM
    protein_chains = list(chain_mapping.values())
    coords = [coord for (chain, _), coord in ca_atoms.items() if chain in protein_chains]
    complex_com = np.mean(coords, axis=0)

    # Box center in Angstrom
    box_center = np.array([box_size * 10 / 2] * 3)  # nm to Angstrom

    # Prepare domain3 positioning if enabled
    domain3_transforms = {}  # {chain_name: (translation, rotation_needed)}

    if position_domain3:
        print("\nCalculating domain3 positioning...")

        # Calculate target positions for domain3 (in Angstrom)
        direction_unit = domain3_direction / np.linalg.norm(domain3_direction)
        separation_angstrom = domain3_separation * 10  # nm to Angstrom

        # Target positions relative to box center
        target_foxp4 = box_center + direction_unit * (separation_angstrom / 2)
        target_foxp1 = box_center - direction_unit * (separation_angstrom / 2)

        for name in chain_mapping.keys():
            chain_id = chain_mapping[name]
            domain3_range = get_domain3_residues(name, domains)

            if domain3_range:
                # Calculate current domain3 COM
                current_com = calculate_domain_com(ca_atoms, chain_id, domain3_range)

                if current_com is not None:
                    # Determine target COM
                    if name == 'FOXP4':
                        target_com = target_foxp4
                    else:  # FOXP1
                        target_com = target_foxp1

                    # Calculate translation
                    translation = target_com - current_com
                    domain3_transforms[name] = {
                        'translation': translation,
                        'target_com': target_com,
                        'domain3_range': domain3_range
                    }

                    print(f"  {name} domain3: {domain3_range[0]}-{domain3_range[1]}")
                    print(f"    Current COM: [{current_com[0]:.2f}, {current_com[1]:.2f}, {current_com[2]:.2f}]")
                    print(f"    Target COM:  [{target_com[0]:.2f}, {target_com[1]:.2f}, {target_com[2]:.2f}]")
                    print(f"    Translation: [{translation[0]:.2f}, {translation[1]:.2f}, {translation[2]:.2f}]")

    # Write PDB with transformations applied
    with open(output_pdb, 'w') as f:
        atom_num = 1
        res_num = 1

        for name in chain_mapping.keys():
            chain_id = chain_mapping[name]

            # Get all CA atoms for this chain, sorted by residue number
            chain_cas = [(rnum, coord) for (c, rnum), coord in ca_atoms.items() if c == chain_id]
            chain_cas.sort(key=lambda x: x[0])

            # Apply transformations if needed
            chain_coords = np.array([coord for _, coord in chain_cas])

            if position_domain3 and name in domain3_transforms:
                transform = domain3_transforms[name]
                domain3_range = transform['domain3_range']
                translation = transform['translation']

                # Apply translation to entire chain
                chain_coords = chain_coords + translation
            else:
                # Default: center at complex COM and shift to box center
                chain_coords = chain_coords - complex_com + box_center

            # Write atoms
            for i, (orig_resnum, _) in enumerate(chain_cas):
                new_coord = chain_coords[i]

                # Get residue name
                resname = "ALA"
                with open(complex_pdb, 'r') as pdb_in:
                    for line in pdb_in:
                        if (line.startswith('ATOM') and ' CA ' in line and
                            line[21] == chain_id and int(line[22:26].strip()) == orig_resnum):
                            resname = line[17:20].strip()
                            break

                # Write CA atom
                f.write(f"ATOM  {atom_num:5d}  CA  {resname:3s} A{res_num:4d}    "
                       f"{new_coord[0]:8.3f}{new_coord[1]:8.3f}{new_coord[2]:8.3f}"
                       f"  1.00  0.00           C\n")
                atom_num += 1
                res_num += 1

        f.write("END\n")

    print(f"\nCreated initial positions PDB with {atom_num-1} atoms -> {output_pdb}")
```

### Step 4: Update Main Script Call

In the main script, update the function call:

```python
# Create the initial positions PDB with domain3 positioning
create_initial_positions_pdb_with_domain3(
    complex_pdb,
    chain_mapping,
    domains,
    BOX_SIZE,
    initial_pdb,
    position_domain3=POSITION_DOMAIN3,
    domain3_separation=DOMAIN3_SEPARATION,
    domain3_direction=DOMAIN3_DIRECTION
)
```

## Alternative Approaches

### Option 1: Use Position Restraints Instead
Instead of placing domain3 in specific positions initially, add position restraints that pull domain3 regions toward target positions during simulation:

- **Pros**: More flexible, allows equilibration
- **Cons**: Requires CALVADOS support for position restraints (may need custom implementation)

### Option 2: Use Distance Restraints Between Domain3s
Add a custom distance restraint between domain3 centers to maintain separation:

```python
# In custom restraints file
FOXP4 1 451-551 | FOXP1 1 450-539 | 20.0 100.0  # 20 nm separation, k=100
```

- **Pros**: Maintains separation without fixing orientation
- **Cons**: Doesn't control direction, only distance

### Option 3: Two-Stage Positioning
1. **Stage 1**: Position domain3s at target locations (translation only)
2. **Stage 2**: Apply gentle steering forces during early simulation to align orientations

## Testing Strategy

1. **Visualization**: Create initial restart.pdb and visualize in PyMOL/VMD
   - Check domain3 positions are at expected locations
   - Verify separation distance
   - Check alignment direction

2. **Short test simulation** (1000 steps):
   - Verify chains don't overlap
   - Check if domain3s maintain separation
   - Analyze trajectory for unwanted interactions

3. **Distance analysis**:
   ```python
   # Calculate domain3 COM-COM distance over trajectory
   import MDAnalysis as mda
   u = mda.Universe('topology.pdb', 'trajectory.dcd')

   domain3_foxp4 = u.select_atoms('resid 451-551 and name CA')
   domain3_foxp1 = u.select_atoms('resid 450-539 and name CA')

   for ts in u.trajectory:
       com1 = domain3_foxp4.center_of_mass()
       com2 = domain3_foxp1.center_of_mass()
       dist = np.linalg.norm(com1 - com2)
       print(f"Frame {ts.frame}: Domain3 separation = {dist/10:.2f} nm")
   ```

## Configuration Parameters to Tune

1. **DOMAIN3_SEPARATION**: Distance between domain3 centers
   - Start with 20 nm, adjust based on chain length
   - Too small: chains may overlap
   - Too large: may affect periodic boundary conditions

2. **DOMAIN3_DIRECTION**: Orientation vector
   - `[1, 0, 0]`: Along x-axis
   - `[0, 1, 0]`: Along y-axis
   - `[1, 1, 0]/sqrt(2)`: Diagonal in xy-plane

3. **BOX_SIZE**: May need to increase if chains are extended
   - Current: 50 nm
   - Recommended: 60-80 nm for 20 nm separation + chain extensions

## Potential Issues and Solutions

### Issue 1: Chains overlap after positioning
**Solution**: Increase DOMAIN3_SEPARATION or use energy minimization before simulation

### Issue 2: Domain3s drift apart during simulation
**Solution**: Add weak distance restraint between domain3 centers (Option 2 above)

### Issue 3: Domain3s rotate during simulation
**Solution**: This is expected! Intra-domain restraints maintain internal structure but allow rotation. If strict orientation is needed, add orientation restraints (requires custom implementation).

### Issue 4: Chain tension
Moving domain3 far from other domains may stretch the linker regions.
**Solution**:
- Use longer equilibration phase
- Check linker residues (between domains) are flexible in domains.yaml
- May need to adjust initial chain conformation

## Next Steps

1. Implement Step 1-4 to create the positioning functionality
2. Test with visualization before running full simulation
3. Run short test simulation to validate approach
4. Analyze results and tune parameters
5. Consider adding optional distance restraints if needed for stability
