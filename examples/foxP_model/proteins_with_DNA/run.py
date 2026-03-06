#!/usr/bin/env python3
"""
Run script for proteins + DNA simulation
DNA modeled as coarse-grained RNA
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
