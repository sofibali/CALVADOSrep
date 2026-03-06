#!/usr/bin/env python3
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
PULL_K = 20.0  # kJ/mol/nm^2
TARGETS_FILE = 'domain3_targets.txt'


def add_domain3_pulling_force(mysim, pull_k, targets_file):
    """
    Add a CustomExternalForce that pulls domain3 residues toward target positions.
    """
    targets = []
    with open(f'{mysim.path}/{targets_file}', 'r') as f:
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

    print(f"\nAdded domain3 pulling force:")
    print(f"  k_pull = {pull_k} kJ/mol/nm^2")
    print(f"  {len(targets)} particles with target positions")
    print(f"  Force group: {pull_force.getForceGroup()}")


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

    with open(f'{path}/{fconfig}', 'r') as stream:
        config = safe_load(stream)

    with open(f'{path}/{fcomponents}', 'r') as stream:
        components = safe_load(stream)

    # Build the system normally (includes interchain restraints via custom_restraints)
    mysim = sim.Sim(path, config, components)
    mysim.build_system()

    # Add the domain3 pulling force AFTER system build, BEFORE simulation
    add_domain3_pulling_force(mysim, PULL_K, TARGETS_FILE)

    # Run the simulation
    mysim.simulate()
