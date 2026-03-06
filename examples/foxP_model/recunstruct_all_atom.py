#!/usr/bin/env python
"""
Reconstruct all-atom structures from CALVADOS CG trajectory using cg2all.

Uses the ResidueBasedModel (for COM-based CG) or CalphaBasedModel (for CA-based CG)
depending on the use_com setting in components.yaml.

Usage:
    # Specific frames (0-indexed):
    python recunstruct_all_atom.py --frames 0 22 44

    # First, middle, last frames (default):
    python recunstruct_all_atom.py

    # Every 10th frame:
    python recunstruct_all_atom.py --stride 10

    # Use -1 for the last frame:
    python recunstruct_all_atom.py --frames 0 100 -1

    # Use a different trajectory:
    python recunstruct_all_atom.py --top top.pdb --traj viz/subsampled.dcd --frames 0 10 20

    # Full trajectory reconstruction (outputs DCD + PDB topology):
    python recunstruct_all_atom.py --full --traj viz/subsampled.dcd
"""

import os
import sys
import argparse
import tempfile
import numpy as np

# --- Compatibility shim for mdtraj < 1.11 (missing bfactors attribute) ---
# cg2all expects mdtraj.Trajectory.bfactors which only exists in dev builds.
# This adds it at runtime without modifying any installed packages.
import mdtraj
if not hasattr(mdtraj.Trajectory, 'bfactors'):
    mdtraj.Trajectory.bfactors = property(
        lambda self: getattr(self, '_bfactors', np.zeros((self.n_frames, self.n_atoms))),
        lambda self, val: setattr(self, '_bfactors', val),
    )
# --- End compatibility shim ---

import mdtraj as md


def get_cg_model(components_yaml='components.yaml'):
    """Determine CG model from components.yaml (use_com -> ResidueBasedModel)."""
    try:
        import yaml
        with open(components_yaml) as f:
            config = yaml.safe_load(f)
        use_com = config.get('defaults', {}).get('use_com', False)
    except (FileNotFoundError, KeyError):
        use_com = False
    return 'ResidueBasedModel' if use_com else 'CalphaBasedModel'


def _load_cg2all_model(cg_model_name, device):
    """Load cg2all model and return (model, config, cg_model_class)."""
    import torch
    import cg2all.lib.libcg
    import cg2all.lib.libmodel
    from cg2all.lib.libconfig import MODEL_HOME

    model_type_map = {
        'CalphaBasedModel': 'CalphaBasedModel',
        'CA': 'CalphaBasedModel',
        'ResidueBasedModel': 'ResidueBasedModel',
        'RES': 'ResidueBasedModel',
    }
    model_type = model_type_map.get(cg_model_name, cg_model_name)

    ckpt_fn = MODEL_HOME / f"{model_type}.ckpt"
    if not ckpt_fn.exists():
        cg2all.lib.libmodel.download_ckpt_file(model_type, ckpt_fn, fix_atom=False)

    ckpt = torch.load(ckpt_fn, map_location=device)
    config = ckpt["hyper_parameters"]

    cg_model_class = getattr(cg2all.lib.libcg, config["cg_model"])

    config = cg2all.lib.libmodel.set_model_config(config, cg_model_class, flattened=False)
    model = cg2all.lib.libmodel.Model(config, cg_model_class, compute_loss=False)

    state_dict = ckpt["state_dict"]
    for key in list(state_dict):
        state_dict[".".join(key.split(".")[1:])] = state_dict.pop(key)
    model.load_state_dict(state_dict)
    model = model.to(device)
    model.set_constant_tensors(device)
    model.eval()
    return model, config, cg_model_class


def reconstruct_frame(input_pdb, output_pdb, model, config, cg_model_class, device):
    """Run cg2all on a single CG PDB frame using the Python API."""
    import torch
    import dgl
    from cg2all.lib.libdata import (
        PredictionData,
        create_trajectory_from_batch,
        standardize_atom_name,
    )
    from cg2all.lib.libter import patch_termini

    input_s = PredictionData(
        input_pdb,
        cg_model_class,
        topology_map=None,
        dcd_fn=None,
        radius=config.globals.radius,
        chain_break_cutoff=1.0,
        is_all=False,
        fix_atom=False,
        batch_size=1,
    )
    loader = dgl.dataloading.GraphDataLoader(input_s, batch_size=1, num_workers=0, shuffle=False)
    batch = next(iter(loader)).to(device)

    with torch.no_grad():
        R = model.forward(batch)[0]["R"]

    traj_s, ssbond_s = create_trajectory_from_batch(batch, R)
    output = patch_termini(traj_s[0])
    output.save(output_pdb)
    return True


def reconstruct_full_trajectory(top_pdb, traj_dcd, output_dcd, output_pdb, model, config, cg_model_class, device):
    """Run cg2all on the full trajectory using the Python API."""
    import torch
    import tqdm
    import dgl
    from cg2all.lib.libdata import (
        PredictionData,
        create_topology_from_data,
        standardize_atom_name,
    )
    from cg2all.lib.libter import patch_termini

    input_s = PredictionData(
        top_pdb,
        cg_model_class,
        topology_map=None,
        dcd_fn=traj_dcd,
        radius=config.globals.radius,
        chain_break_cutoff=1.0,
        is_all=False,
        fix_atom=False,
        batch_size=1,
    )
    n_frame0 = input_s.n_frame0
    unitcell_lengths = input_s.cg.unitcell_lengths
    unitcell_angles = input_s.cg.unitcell_angles

    loader = dgl.dataloading.GraphDataLoader(input_s, batch_size=1, num_workers=0, shuffle=False)

    xyz = []
    for batch in tqdm.tqdm(loader, total=len(loader)):
        batch = batch.to(device)
        with torch.no_grad():
            R = model.forward(batch)[0]["R"].cpu().detach().numpy()
            mask = batch.ndata["output_atom_mask"].cpu().detach().numpy()
            xyz.append(R[mask > 0.0])

    xyz = np.array(xyz)
    top, atom_index = create_topology_from_data(batch)
    xyz = xyz[:, atom_index]
    traj = mdtraj.Trajectory(
        xyz=xyz, topology=top,
        unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles,
    )
    output = patch_termini(traj)
    output.save(output_dcd)
    output[-1].save(output_pdb)
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Reconstruct all-atom structures from CALVADOS CG trajectory using cg2all'
    )
    parser.add_argument('--top', default='top.pdb', help='CG topology file (default: top.pdb)')
    parser.add_argument('--traj', default='traj.dcd', help='CG trajectory file (default: traj.dcd)')
    parser.add_argument('--frames', type=int, nargs='+', default=None,
                        help='Frame indices to reconstruct (use -1 for last). Default: first, middle, last')
    parser.add_argument('--stride', type=int, default=None,
                        help='Reconstruct every Nth frame')
    parser.add_argument('--full', action='store_true',
                        help='Reconstruct full trajectory as DCD (instead of individual PDBs)')
    parser.add_argument('--cg', default=None,
                        help='CG model (default: auto-detect from components.yaml)')
    parser.add_argument('--device', default='cpu', help='Device for cg2all (cpu or cuda)')
    parser.add_argument('--output', default='viz/allatom_frames', help='Output directory')
    args = parser.parse_args()

    import torch
    device = torch.device(args.device)

    # Determine CG model
    cg_model_name = args.cg or get_cg_model()
    print(f"CG model: {cg_model_name}")

    # Load model once (reused for all frames)
    print("Loading cg2all model...")
    model, config, cg_model_class = _load_cg2all_model(cg_model_name, device)
    print("  Model loaded.")

    # Full trajectory mode
    if args.full:
        os.makedirs(args.output, exist_ok=True)
        out_dcd = os.path.join(args.output, 'allatom.dcd')
        out_pdb = os.path.join(args.output, 'allatom_top.pdb')
        print(f"Reconstructing full trajectory: {args.traj}")
        success = reconstruct_full_trajectory(
            args.top, args.traj, out_dcd, out_pdb, model, config, cg_model_class, device
        )
        if success:
            print(f"All-atom trajectory: {out_dcd}")
            print(f"All-atom topology:   {out_pdb}")
        return

    # Load trajectory
    print(f"Loading trajectory: {args.traj}")
    traj = md.load(args.traj, top=args.top)
    n_frames = traj.n_frames
    print(f"  Total frames: {n_frames}")
    print(f"  Chains: {traj.n_chains}")
    for c in traj.top.chains:
        print(f"    Chain {c.index}: {c.n_residues} residues")

    # Determine which frames to reconstruct
    if args.stride:
        frame_indices = list(range(0, n_frames, args.stride))
    elif args.frames:
        frame_indices = [f if f >= 0 else n_frames + f for f in args.frames]
    else:
        frame_indices = [0, n_frames // 2, n_frames - 1]

    print(f"\nReconstructing {len(frame_indices)} frames: {frame_indices}")
    os.makedirs(args.output, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp_dir:
        for frame_idx in frame_indices:
            if frame_idx < 0 or frame_idx >= n_frames:
                print(f"  Frame {frame_idx} out of range (0-{n_frames-1}), skipping")
                continue

            # Extract CG frame as PDB
            cg_pdb = os.path.join(tmp_dir, f'frame_{frame_idx:06d}_cg.pdb')
            traj[frame_idx].save_pdb(cg_pdb)

            # Reconstruct all-atom
            out_pdb = os.path.join(args.output, f'frame_{frame_idx:06d}_AA.pdb')
            print(f"  Frame {frame_idx}...", end='', flush=True)
            try:
                success = reconstruct_frame(cg_pdb, out_pdb, model, config, cg_model_class, device)
                print(f" -> {out_pdb}")
            except Exception as e:
                print(f" FAILED: {e}")

    print(f"\nDone! All-atom PDBs in: {args.output}/")


if __name__ == '__main__':
    main()
