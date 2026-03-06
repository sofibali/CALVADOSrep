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

from calvados.analysis import calc_com_traj, calc_contact_map

chainid_dict = dict(FOXP4 = 0, FOX = 1)
calc_com_traj(path="/home/sbali/CALVADOS/examples/foxP_model/FOXP4_FOX_chain_C",sysname="FOXP4_FOX_chain_C",output_path="data",residues_file="/home/sbali/CALVADOS/examples/foxP_model/input/residues_CALVADOS3.csv",chainid_dict=chainid_dict,start=100)
calc_contact_map(path="/home/sbali/CALVADOS/examples/foxP_model/FOXP4_FOX_chain_C",sysname="FOXP4_FOX_chain_C",output_path="data",chainid_dict=chainid_dict)
