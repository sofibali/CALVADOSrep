#!/usr/bin/env python
"""
Two 2-domain contiguous constructs, to test whether the extreme compaction
seen in md_full (MD1L1-MD2-MD3, 599 res -- roughly HALF the inter-domain
distance of the same 3 domains embedded in mka_full's longer chain) is a
general "short isolated fragment" effect or specific to the 3-domain case:

  md1_md2   FL 790-1193  (MD1L1-MD2, 404 res)
  md2_md3   FL 1004-1388 (MD2-MD3, 385 res)

Same convention as prepare_extra_full_constructs.py: sliced from the shared
FL AF2 model, 5 replicates x 1 us, harmonic restraints only (no KH domains
involved here, so no Go-model restraint needed).

    python prepare_md_pair_constructs.py
    bash md1_md2/run_all.sh parallel; bash md2_md3/run_all.sh parallel
"""
from prepare_extra_full_constructs import prepare_one

CONSTRUCTS = [
    {'key': 'md1_md2', 'fl_start': 790, 'fl_end': 1193, 'box': 70,
     'units': ('md1l1', 'md2'), 'needs_go': False},
    {'key': 'md2_md3', 'fl_start': 1004, 'fl_end': 1388, 'box': 70,
     'units': ('md2', 'md3'), 'needs_go': False},
]


def main():
    for spec in CONSTRUCTS:
        prepare_one(spec)


if __name__ == "__main__":
    main()
