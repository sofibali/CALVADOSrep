#!/usr/bin/env python
"""
Drop re-simulated frames from a slab trajectory after a preemption.

THE PROBLEM
-----------
`sim.py` writes a DCD frame every `wfreq` steps but only checkpoints every
`steps/10`. When `run_slab_opportunistic.sh` yields a GPU it SIGTERMs the run
mid-interval, so the DCD already contains frames for steps the checkpoint never
captured. On resume, `DCDReporter(..., append=True)` does not rewind the file:
the leg restarts from the checkpoint, re-simulates that window with a different
random stream, and appends a second trajectory through it.

The result is an abandoned branch sitting in the middle of the file. It does
not bias c_sat -- those frames are still valid equilibrium samples -- but frame
index no longer maps to simulation time and the repeated window is counted
twice.

THE FIX
-------
Before resuming, truncate the DCD back to the last checkpoint. Then the resumed
leg appends onto a clean boundary and nothing is ever duplicated.

Which frame is the boundary? `slab_steps.py prepare` sets
`logfreq == steps/10`, i.e. exactly the checkpoint interval, so the last line in
`{sysname}.log` is written at the same step the checkpoint was saved at. Valid
frames are therefore `last_logged_step // wfreq`. Anything beyond that belongs
to the abandoned branch.

This only rewrites the DCD header's frame count and truncates the file, so it
is O(1) regardless of trajectory size -- no re-reading 200 MB.

USAGE
-----
    python trim_dcd.py <run_dir>            # report only
    python trim_dcd.py <run_dir> --apply    # actually truncate

Run it only between legs, never against a live run.
"""
import argparse
import os
import struct
import sys

import yaml


def dcd_layout(path, natoms):
    """(offset of first frame, bytes per frame, NSET field, frames present)."""
    with open(path, 'rb') as fh:
        blk = struct.unpack('<i', fh.read(4))[0]
        magic = fh.read(4)
        if magic != b'CORD':
            raise ValueError(f'{path}: not a DCD (magic {magic!r})')
        nset = struct.unpack('<i', fh.read(4))[0]
        fh.seek(0)
        fh.read(4 + blk + 4)                    # header block
        t = struct.unpack('<i', fh.read(4))[0]  # title block
        fh.read(t + 4)
        n = struct.unpack('<i', fh.read(4))[0]  # natoms block
        fh.read(n + 4)
        start = fh.tell()
    # unit-cell record (4+48+4) + one Fortran record per coordinate axis
    frame = 56 + 3 * (8 + 4 * natoms)
    size = os.path.getsize(path)
    if (size - start) % frame:
        raise ValueError(
            f'{path}: {(size - start) % frame} trailing bytes -- frame size '
            f'{frame} does not divide the data. Refusing to touch it.')
    return start, frame, nset, (size - start) // frame


def last_logged_step(path):
    last = 0
    if not os.path.isfile(path):
        return 0
    with open(path, errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                last = max(last, int(float(line.split('\t')[0])))
            except (ValueError, IndexError):
                continue
    return last


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir')
    ap.add_argument('--apply', action='store_true',
                    help='truncate; without this it only reports')
    args = ap.parse_args()

    d = os.path.abspath(args.run_dir)
    meta = yaml.safe_load(open(os.path.join(d, 'slab_meta.yaml')))
    cfg = yaml.safe_load(open(os.path.join(d, 'config.yaml')))
    sysname, natoms = meta['sysname'], int(meta['beads'])
    wfreq = int(float(cfg['wfreq']))

    dcd = os.path.join(d, f'{sysname}.dcd')
    log = os.path.join(d, f'{sysname}.log')
    if not os.path.isfile(dcd):
        print(f'  [trim] {sysname}: no production DCD yet, nothing to do')
        return 0

    start, frame, nset, present = dcd_layout(dcd, natoms)
    step = last_logged_step(log)
    keep = step // wfreq
    excess = present - keep

    if excess <= 0:
        print(f'  [trim] {sysname}: {present} frames, {keep} expected -- clean')
        return 0

    pct = 100.0 * excess / present
    print(f'  [trim] {sysname}: {present} frames on disk, {keep} backed by the '
          f'checkpoint (step {step:,} / wfreq {wfreq:,})')
    print(f'  [trim] {excess} re-simulated frame(s) to drop ({pct:.2f}% of the file)')

    if not args.apply:
        print('  [trim] report only -- pass --apply to truncate')
        return 0

    with open(dcd, 'r+b') as fh:
        fh.seek(8)
        fh.write(struct.pack('<i', keep))       # NSET
        fh.truncate(start + keep * frame)
    # re-read to confirm
    _, _, nset2, present2 = dcd_layout(dcd, natoms)
    ok = (present2 == keep == nset2)
    print(f'  [trim] truncated to {present2} frames, NSET={nset2} '
          f'{"OK" if ok else "-- MISMATCH"}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
