#!/usr/bin/env python3
"""
Run script for proteins-only simulation (FOXP4, FOX, chain_C)
Quick simulation with highly restrained domains
"""

import os
from calvados import sim

# Paths - sim.run joins path with fconfig/fcomponents, so use relative names
cwd = os.path.dirname(os.path.abspath(__file__))
path = cwd
fconfig = 'config.yaml'
fcomponents = 'components.yaml'

# Run simulation
sim.run(path=path, fconfig=fconfig, fcomponents=fcomponents)
