#!/usr/bin/env python3
"""Centre a two-chain frame and remove the periodic-boundary artefact.

THE PROBLEM
-----------
CALVADOS and HyRes write each chain WHOLE, so no chain is split across the box
edge -- but nothing keeps the two chains in the same periodic image. A pair that
is 0.4 nm apart under the minimum-image convention can be drawn 39 nm apart,
with one chain outside the cell entirely. Distances computed with `box=` are
correct either way; the picture is not.

THE FIX
-------
1. Translate so PARP9's centre of mass sits at the box centre.
2. Shift DTX3L AS A RIGID WHOLE by whichever periodic image vector minimises the
   inter-chain distance.

Shifting the partner as a whole -- rather than wrapping atom by atom -- is what
keeps the chain intact. Per-atom wrapping is the usual way this goes wrong: it
puts individual residues on opposite faces of the box and shreds the cartoon.

After this the drawn inter-chain distance equals the minimum-image distance, so
what you see is what was measured.
"""
import numpy as np
from MDAnalysis.analysis.distances import distance_array


def center_pair(pos, n_first, box_nm):
    """pos: (N,3) Angstrom, chain 1 = pos[:n_first]. Returns centred copy."""
    L = box_nm * 10.0
    p = pos.copy().astype(np.float64)
    a, b = p[:n_first], p[n_first:]
    p -= a.mean(0) - np.array([L/2]*3)          # PARP9 COM to the box centre
    a, b = p[:n_first], p[n_first:]
    # whole-chain image shift on the partner
    best, shift = None, np.zeros(3)
    boxv = np.array([L, L, L, 90., 90., 90.], np.float32)
    ref = distance_array(a.astype(np.float32), b.astype(np.float32), box=boxv).min()
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                s = np.array([dx, dy, dz]) * L
                d = np.linalg.norm(
                    a[:, None, :] - (b + s)[None, ::7, :], axis=2).min()
                if best is None or d < best:
                    best, shift = d, s
    p[n_first:] += shift
    return p, float(ref), float(best)


def annotate(pos, n_first, box_nm):
    """True if the pair needed an image shift to be drawn adjacent."""
    _, ref, got = center_pair(pos, n_first, box_nm)
    return abs(got/10.0 - ref/10.0) > 0.01


def wrap_into_box(pos, n_first, box_nm):
    """Shift each chain WHOLE so its centroid lies inside the box.

    For display only. CALVADOS's grid placement puts molecules at grid points
    that can sit on a box corner or face -- correct under periodic boundaries,
    but it draws as a chain hanging outside the cell (and PyMOL then clips it).
    Shifting each chain by a whole box vector preserves every internal distance
    and the inter-chain minimum-image distance; it only chooses which image to
    draw. Per-atom wrapping would split the chain and is never used here.
    """
    import numpy as np
    L = box_nm * 10.0
    p = pos.copy().astype(np.float64)
    for sl in (slice(0, n_first), slice(n_first, len(p))):
        c = p[sl].mean(0)
        p[sl] -= np.floor(c / L) * L
    return p
