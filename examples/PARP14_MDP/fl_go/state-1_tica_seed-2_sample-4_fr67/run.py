from calvados import sim
from argparse import ArgumentParser

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument('--path', nargs='?', default='.', const='.', type=str)
    p.add_argument('--config', nargs='?', default='config.yaml', const='config.yaml', type=str)
    p.add_argument('--components', nargs='?', default='components.yaml', const='components.yaml', type=str)
    a = p.parse_args()
    sim.run(path=a.path, fconfig=a.config, fcomponents=a.components)
