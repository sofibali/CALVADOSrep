from calvados import sim
from argparse import ArgumentParser

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--path', default='.')
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--components', default='components.yaml')
    args = parser.parse_args()
    sim.run(path=args.path, fconfig=args.config, fcomponents=args.components)
