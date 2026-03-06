#!/usr/bin/env python
"""
Prepare all-atom trajectory movies for FOXP model simulations.

Steps for each simulation:
1. Subsample the CG DCD trajectory (stride=20 → ~50 frames)
2. Run cg2all --full on the subsampled trajectory
3. Output: viz/movie/allatom.dcd + allatom_top.pdb

Usage:
    conda run -n CALVADOS python prepare_movie_trajectories.py

    # Custom stride:
    conda run -n CALVADOS python prepare_movie_trajectories.py --stride 10

    # Single simulation only:
    conda run -n CALVADOS python prepare_movie_trajectories.py --sim option1_pull_interchain_d12
"""

import os
import sys
import argparse
import numpy as np

# --- Compatibility shim for mdtraj < 1.11 (missing bfactors attribute) ---
import mdtraj
if not hasattr(mdtraj.Trajectory, 'bfactors'):
    mdtraj.Trajectory.bfactors = property(
        lambda self: getattr(self, '_bfactors', np.zeros((self.n_frames, self.n_atoms))),
        lambda self, val: setattr(self, '_bfactors', val),
    )
import mdtraj as md

BASE = os.path.dirname(os.path.abspath(__file__))
INDIVIDUAL = os.path.join(BASE, 'individual_chains')

SIMULATIONS = {
    'option1_pull_interchain_d12': {
        'sysname': 'option1_pull_interchain_d12',
        'top': 'restart.pdb',
    },
    'two_chains_coils_interchain': {
        'sysname': 'two_chains_coils_interchain',
        'top': 'restart.pdb',
    },
}


def subsample_trajectory(sim_dir, sysname, top_file, stride, output_dir):
    """Subsample CG trajectory and save as DCD + topology PDB."""
    top_path = os.path.join(sim_dir, top_file)
    dcd_path = os.path.join(sim_dir, f'{sysname}.dcd')

    print(f"\nLoading CG trajectory: {dcd_path}")
    traj = md.load(dcd_path, top=top_path)
    print(f"  Total frames: {traj.n_frames}, Atoms: {traj.n_atoms}")

    # Subsample
    sub = traj[::stride]
    print(f"  Subsampled to {sub.n_frames} frames (stride={stride})")

    os.makedirs(output_dir, exist_ok=True)
    sub_dcd = os.path.join(output_dir, 'subsampled_cg.dcd')
    sub_top = os.path.join(output_dir, 'cg_top.pdb')

    sub.save_dcd(sub_dcd)
    sub[0].save_pdb(sub_top)
    print(f"  Saved: {sub_dcd}")
    print(f"  Saved: {sub_top}")

    return sub_top, sub_dcd, sub.n_frames


def reconstruct_allatom(cg_top, cg_dcd, output_dir, device='cpu'):
    """Run cg2all full trajectory reconstruction."""
    # Import here so the script can still show help without cg2all
    import torch
    sys.path.insert(0, BASE)
    from recunstruct_all_atom import _load_cg2all_model, reconstruct_full_trajectory

    cg_model_name = 'ResidueBasedModel'  # use_com=true in components.yaml
    print(f"\nReconstructing all-atom trajectory...")
    print(f"  CG model: {cg_model_name}")
    print(f"  Device: {device}")

    device = torch.device(device)
    model, config, cg_model_class = _load_cg2all_model(cg_model_name, device)

    out_dcd = os.path.join(output_dir, 'allatom.dcd')
    out_pdb = os.path.join(output_dir, 'allatom_top.pdb')

    reconstruct_full_trajectory(cg_top, cg_dcd, out_dcd, out_pdb, model, config, cg_model_class, device)

    print(f"  All-atom DCD:      {out_dcd}")
    print(f"  All-atom topology: {out_pdb}")
    return out_pdb, out_dcd


def main():
    parser = argparse.ArgumentParser(description='Prepare all-atom trajectory movies')
    parser.add_argument('--stride', type=int, default=20,
                        help='Subsample stride (default: 20, gives ~50 frames from 1000)')
    parser.add_argument('--sim', type=str, default=None,
                        help='Process single simulation (default: all)')
    parser.add_argument('--device', type=str, default='cpu',
                        help='Device for cg2all (cpu or cuda)')
    parser.add_argument('--skip-reconstruct', action='store_true',
                        help='Only subsample, skip cg2all reconstruction')
    args = parser.parse_args()

    sims_to_process = {args.sim: SIMULATIONS[args.sim]} if args.sim else SIMULATIONS

    for name, info in sims_to_process.items():
        sim_dir = os.path.join(INDIVIDUAL, name)
        output_dir = os.path.join(sim_dir, 'viz', 'movie')

        print(f"\n{'='*60}")
        print(f"Processing: {name}")
        print(f"{'='*60}")

        cg_top, cg_dcd, n_frames = subsample_trajectory(
            sim_dir, info['sysname'], info['top'], args.stride, output_dir
        )

        if not args.skip_reconstruct:
            reconstruct_allatom(cg_top, cg_dcd, output_dir, args.device)
        else:
            print("  Skipping reconstruction (--skip-reconstruct)")

    print(f"\n{'='*60}")
    print("Done! Next steps:")
    print("  1. Open PyMOL scripts to view movies:")
    print(f"     pymol {os.path.join(BASE, 'pymol_movie_option1_pull.pml')}")
    print(f"     pymol {os.path.join(BASE, 'pymol_movie_coils_interchain.pml')}")
    print("  2. Or render headless:")
    print(f"     pymol -cq {os.path.join(BASE, 'pymol_movie_option1_pull.pml')}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
