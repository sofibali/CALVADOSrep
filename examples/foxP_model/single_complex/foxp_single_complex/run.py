from calvados import sim
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--path',nargs='?', default='.', const='.', type=str)
    parser.add_argument('--config',nargs='?', default='config.yaml', const='config.yaml', type=str)
    parser.add_argument('--components',nargs='?', default='components.yaml', const='components.yaml', type=str)

    args = parser.parse_args()

    path = args.path
    fconfig = args.config
    fcomponents = args.components

    sim.run(path=path,fconfig=fconfig,fcomponents=fcomponents)

from calvados.analysis import save_conf_prop

save_conf_prop(path="/home/sbali/CALVADOS/examples/foxP_model/single_complex/foxp_single_complex", name="foxp_single_complex",
               residues_file="/home/sbali/CALVADOS/examples/foxP_model/single_complex/input/residues_CALVADOS3.csv",
               output_path="/home/sbali/CALVADOS/examples/foxP_model/single_complex/data",
               start=100, is_idr=False, select='all')
