#!/usr/bin/env python3
"""HyRes run with CALVADOS-equivalent domain restraints.

argv: start.pdb  complex.psf  restraints.txt  gpu_id  seed

Adapted from HyRes_GPU/run_OpenMM/run.nvt.py with three changes:
  * intra-domain harmonic CA-CA restraints are added (the repo's
    restraints/folded_restriants.py is a TODO stub with a syntax error and is
    wired into nothing), matching CALVADOS: same domains.yaml, cutoff 0.9 nm,
    k = 700 kJ/mol/nm^2;
  * minimise BEFORE setting velocities -- the repo order NaNs immediately here;
  * temperature 293 K and a per-replicate thermostat seed.
"""
from sys import argv
from openmm.app import *
from openmm import *
from openmm import unit
import numpy as np

pdb_file, psf_file, restr_file, gpu_id, seed = argv[1], argv[2], argv[3], argv[4], int(argv[5])
top_inp, param_inp = 'top_hyres_GPU.inp', 'param_hyres_GPU.inp'

dt = 0.004*unit.picoseconds
total_step = 125000000          # 500 ns, matching CALVADOS
temperature = 293*unit.kelvin
log_freq, dcd_freq = 250000, 12500     # 50 ps/frame, matching CALVADOS
c_ion, er = 0.15, 20.0
eps_hb, sigma_hb = 2.0*unit.kilocalorie_per_mole, 0.29*unit.nanometer
r_cut = 1.8*unit.nanometer
friction, freq = 0.1/unit.picosecond, 25
lx = 40.0*unit.nanometer
a = Vec3(lx, 0.0, 0.0); b = Vec3(0.0, lx, 0.0); c = Vec3(0.0, 0.0, lx)

exec(open('_ff_body.py').read())        # force-field construction, verbatim from the repo

# ---- CALVADOS-equivalent intra-domain restraints -------------------------
restr = HarmonicBondForce()
restr.setUsesPeriodicBoundaryConditions(True)
n = 0
for line in open(restr_file):
    if line.startswith('#') or not line.strip():
        continue
    i, j, r0, k = line.split()
    restr.addBond(int(i), int(j), float(r0)*unit.nanometer,
                  float(k)*unit.kilojoule_per_mole/unit.nanometer**2)
    n += 1
system.addForce(restr)
print(f'added {n} intra-domain harmonic restraints (CALVADOS-equivalent)')

integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setRandomNumberSeed(seed)
plat = Platform.getPlatformByName('CUDA')
simulation = Simulation(top, system, integrator, plat,
                        {'Precision': 'mixed', 'DeviceIndex': gpu_id})
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy(maxIterations=10000)
print('E after minimisation:', simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.context.setVelocitiesToTemperature(temperature, seed)
simulation.reporters.append(DCDReporter('system.dcd', dcd_freq))
simulation.reporters.append(StateDataReporter('system.log', log_freq, step=True, time=True,
                            progress=True, totalSteps=total_step, temperature=True))
print('# NVT simulation')
simulation.step(total_step)
