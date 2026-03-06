#!/usr/bin/env python
"""
Extract frames from CALVADOS trajectory as PDB files.

This uses mdtraj (same as CALVADOS) to extract frames from the CG trajectory.
Output is coarse-grained (CA atoms only) - CALVADOS does not include backmapping.

Usage:
    python extract_frames.py --top top.pdb --traj trajectory.dcd --frames 0 -1 100
    python extract_frames.py --top top.pdb --traj trajectory.dcd --all --stride 100
"""

import os
import argparse
import mdtraj as md


def main():
    parser = argparse.ArgumentParser(description='Extract frames from CALVADOS trajectory')
    parser.add_argument('--top', type=str, required=True, help='Topology file (top.pdb)')
    parser.add_argument('--traj', type=str, required=True, help='Trajectory file (DCD)')
    parser.add_argument('--frames', type=int, nargs='+', default=None,
                        help='Frame indices to extract (use -1 for last)')
    parser.add_argument('--all', action='store_true', help='Extract all frames')
    parser.add_argument('--stride', type=int, default=1, help='Stride for --all mode')
    parser.add_argument('--output', type=str, default='frames', help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    print(f"Loading trajectory: {args.traj}")
    print(f"Topology: {args.top}")
    traj = md.load(args.traj, top=args.top)

    n_frames = traj.n_frames
    print(f"Total frames: {n_frames}")
    print(f"Atoms per frame: {traj.n_atoms}")

    if args.all:
        # Extract all frames with stride
        frames = list(range(0, n_frames, args.stride))
        print(f"\nExtracting {len(frames)} frames (stride={args.stride})...")
    elif args.frames:
        # Extract specific frames
        frames = [f if f >= 0 else n_frames + f for f in args.frames]
        print(f"\nExtracting frames: {frames}")
    else:
        # Default: first, middle, last
        frames = [0, n_frames // 2, n_frames - 1]
        print(f"\nExtracting default frames (first, middle, last): {frames}")

    for frame_idx in frames:
        if 0 <= frame_idx < n_frames:
            output_file = os.path.join(args.output, f"frame_{frame_idx:06d}.pdb")
            traj[frame_idx].save_pdb(output_file)
            print(f"  Frame {frame_idx} -> {output_file}")
        else:
            print(f"  WARNING: Frame {frame_idx} out of range, skipping")

    print(f"\nDone! PDB files saved to: {args.output}/")


if __name__ == '__main__':
    main()
