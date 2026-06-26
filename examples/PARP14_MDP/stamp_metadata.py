#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write a metadata.json into a simulation-set folder so the analysis scripts can
analyze it via `--sim-folder PATH` with no other arguments.

A metadata.json records which full-length PARP14 domain UNITS the construct
contains (and the dcd basename). The analysis scripts read it to remap active
sites / domain boundaries into the construct's residue numbering.

Usage
-----
  # Stamp every existing named set (fl, md, core, mka, norrm, noart, md3art,
  # fl_optimized) that is present on disk — units come from sim_registry:
  python stamp_metadata.py --all-named

  # Stamp one new construct folder, giving its FL domain units explicitly:
  python stamp_metadata.py mynewconstruct --units md1l1 md2 md3

  # Stamp a folder whose basename already matches a known named set:
  python stamp_metadata.py md            # units inferred from the 'md' set

Valid units: rrm1 rrm2 rrm3 kh1-kh6 kh7a md1l1 md2 md3 khb-kh8 wwe art
"""

import os
import sys
from argparse import ArgumentParser

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim_registry as reg


def stamp(folder, units=None):
    """Write metadata.json into `folder`; return the written dict."""
    folder = os.path.abspath(folder)
    if not os.path.isdir(folder):
        raise FileNotFoundError(folder)
    resolved = reg._resolve_units(folder, units)
    sysname = reg.detect_sysname(folder)
    meta = reg.write_metadata(folder, resolved, sysname=sysname)
    print(f"  {folder}\n    units={meta['units']} sysname={meta['sysname']} "
          f"sites={meta['sites']}")
    return meta


def main():
    ap = ArgumentParser(description='Write metadata.json into simulation folders')
    ap.add_argument('folders', nargs='*', help='Folders to stamp')
    ap.add_argument('--units', nargs='+', default=None,
                    help='FL domain units for the folder(s) (e.g. md1l1 md2 md3). '
                         'Omit to infer from a matching named set.')
    ap.add_argument('--all-named', action='store_true',
                    help='Stamp every named set folder present on disk.')
    args = ap.parse_args()

    print("Writing metadata.json ...")
    n = 0
    if args.all_named:
        for set_key, info in reg.SETS.items():
            folder = os.path.join(reg.CWD, set_key)
            if os.path.isdir(folder):
                stamp(folder, units=info['units'])
                n += 1
    for folder in args.folders:
        stamp(folder, units=args.units)
        n += 1

    if n == 0:
        ap.error('nothing to do: pass folder paths and/or --all-named')
    print(f"Done ({n} folder(s)).")


if __name__ == '__main__':
    main()
