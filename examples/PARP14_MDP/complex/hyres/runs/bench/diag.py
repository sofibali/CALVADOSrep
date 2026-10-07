import numpy as np
from openmm.app import *
from openmm import *
from openmm import unit
from scipy.spatial import cKDTree
pdb = PDBFile('complex_start.pdb')
X = np.array(pdb.positions.value_in_unit(unit.nanometer))
d,_ = cKDTree(X).query(X, k=2)
print(f'{len(X)} particles; closest pair {d[:,1].min()*10:.3f} A; '
      f'{(d[:,1]<0.05).sum()} pairs closer than 0.5 A')
psf = CharmmPsfFile('complex.psf'); psf.setBox(40*unit.nanometer,40*unit.nanometer,40*unit.nanometer)
params = CharmmParameterSet('top_hyres_GPU.inp','param_hyres_GPU.inp')
system = psf.createSystem(params, nonbondedMethod=CutoffPeriodic, constraints=HBonds)
integ = LangevinMiddleIntegrator(293*unit.kelvin, 0.1/unit.picosecond, 0.004*unit.picoseconds)
ctx = Context(system, integ, Platform.getPlatformByName('CUDA'), {'Precision':'mixed','DeviceIndex':'0'})
ctx.setPositions(pdb.positions)
e0 = ctx.getState(getEnergy=True).getPotentialEnergy()
print('E before minimisation:', e0)
LocalEnergyMinimizer.minimize(ctx, 10*unit.kilojoule_per_mole/unit.nanometer, 5000)
e1 = ctx.getState(getEnergy=True).getPotentialEnergy()
print('E after  minimisation:', e1)
