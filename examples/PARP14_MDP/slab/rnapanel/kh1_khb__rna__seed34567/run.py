import os
from calvados import sim

path = os.path.dirname(os.path.abspath(__file__))
sim.run(path=path, fconfig='config.yaml', fcomponents='components.yaml')
