#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Render the PARP14 conformational-state surface-property gallery for one
all-atom representative-frame state PDB:

  - domain-colored cartoon (reference view)
  - electrostatic-potential surface (front + back): real Poisson-Boltzmann
    electrostatics via pdb2pqr + APBS (SBgrid), not PyMOL's built-in
    util.protein_vacuum_esp -- that has a reproducible bug that makes it
    unusable in this scripted batch context (see prep_apbs.py / git history)
  - Kyte-Doolittle hydrophobic surface (front + back)
  - per-domain face classification surface: active-site face vs back
    face (dot-product projection, same method as
    sim_analysis/figure_sasa_faces.py::classify_face) for MD1L1/MD2/MD3/ART;
    generic 2-means geometric split (scipy) for domains with no defined
    active site (RRM1-3, KH1-KH6, KH7a, KHb-KH8, WWE)
  - active-site pocket close-ups (one dir per site present in this
    construct): electrostatics, hydrophobicity, discrete residue-chemistry
    classification, H-bond donor/acceptor atom mapping

Every image is saved separately (PNG). Optionally also saves a per-state
.pse with everything this script built, for later combination into a
multi-state session by make_surface_gallery_session.py.

Run inside the pymol-render conda env (has pymol + numpy + scipy); SBgrid
(pdb2pqr + apbs) is invoked internally as a subprocess, no conda env needed
for it:
    conda run -n pymol-render python make_surface_gallery.py \\
        --set-key fl_optimized \\
        --pdb representative_frames/2026-06-25/fl_optimized_ca25_tica_allatom/state_1_allatom.pdb \\
        --label state_1 \\
        --out-dir figures/06_surface_gallery/fl_optimized_ca25_tica/state_1 \\
        --pse-out figures/06_surface_gallery/fl_optimized_ca25_tica/state_1/state_1.pse
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.cluster.vq import kmeans2

import pymol
from pymol import cmd
pymol.finish_launching(['pymol', '-cq'])

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surface_gallery_lib as sgl   # noqa: E402
import sim_registry as sr           # noqa: E402

IMG_SIZE = (900, 900)
FACE_COLORS_ACTIVE = {'active': '#27ae60', 'back': '#7f8c8d'}
SITE_TO_DOMAIN = {'MD1': 'MD1L1', 'MD2': 'MD2', 'MD3': 'MD3', 'ART': 'ART'}
DOMAIN_TO_SITE = {v: k for k, v in SITE_TO_DOMAIN.items()}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _hex_to_rgb(hexstr):
    hexstr = hexstr.lstrip('#')
    return tuple(int(hexstr[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


_COLOR_CACHE = {}


def hexcol(hexstr):
    """Register a hex color with PyMOL under a sanitized name (cmd.color()
    does not accept raw hex strings, only built-ins or set_color names)
    and return that name."""
    if hexstr not in _COLOR_CACHE:
        name = 'c_' + hexstr.lstrip('#')
        cmd.set_color(name, list(_hex_to_rgb(hexstr)))
        _COLOR_CACHE[hexstr] = name
    return _COLOR_CACHE[hexstr]


def setup_common():
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('antialias', 2)
    cmd.set('ray_trace_mode', 0)
    # Quality 0 looks fine at the whole-molecule scale these renders are
    # used at, and keeps ray-trace time reasonable across ~30 surfaces/state.
    cmd.set('surface_quality', 0)
    cmd.set('specular', 0.15)


def get_residue_table(obj):
    """[(resi:int, resn:str)] for polymer CA atoms, N->C order."""
    space = {'out': []}
    cmd.iterate(f'{obj} and polymer and name CA',
                'out.append((int(resi), resn))', space=space)
    return space['out']


def show_only(obj):
    """Hide every other object's representations so only `obj`'s geometry
    appears in the next render -- multiple leftover copies (esp/hydro/face/
    etc.) would otherwise all render on top of each other.

    Uses cmd.hide() (representation-level) on every OTHER object, never
    cmd.disable()/cmd.enable() (object-level) and never on `obj` itself.
    util.protein_vacuum_esp's derived "_e_pot" object is permanently broken
    the moment cmd.disable() touches it -- reproduced repeatedly, and a
    later cmd.enable() does not undo it. Representation-level hide/show is
    safe for it (and for everything else); this keeps the whole script on
    one safe pattern instead of a special case just for ESP. Caller is
    responsible for cmd.show()-ing whatever representation `obj` needs."""
    for name in cmd.get_names('objects'):
        if name != obj:
            cmd.hide('everything', name)


def snap(view_sel, path, buffer=3):
    cmd.orient(view_sel)
    cmd.zoom(view_sel, buffer=buffer)
    cmd.ray(*IMG_SIZE)
    cmd.png(str(path), dpi=200)
    log(f"  saved {path.name}")


def color_domains(obj, set_key):
    cmd.color('gray80', obj)
    for dname, (color_name, rgb, (cs, ce)) in sgl.get_domain_ranges_for_set(set_key).items():
        cmd.set_color(color_name, list(rgb))
        cmd.color(color_name, f'{obj} and resi {cs}-{ce}')


def show_catalytic_spheres(obj, set_key):
    all_cat = []
    for site in sr.SITE_NAMES:
        cat, _ = sgl.get_pocket_for_set(set_key, site)
        all_cat += cat
    if all_cat:
        sel = f'{obj} and resi {"+".join(map(str, all_cat))} and name CA'
        cmd.show('spheres', sel)
        cmd.color('yellow', sel)
        cmd.set('sphere_scale', 0.7, sel)
    return all_cat


def render_domain_cartoon(obj, set_key, out_dir):
    log("Rendering domain cartoon...")
    show_only(obj)
    cmd.hide('everything', obj)
    cmd.show('cartoon', obj)
    color_domains(obj, set_key)
    show_catalytic_spheres(obj, set_key)
    snap(obj, out_dir / 'domain_cartoon.png', buffer=5)


def ensure_esp_surface(obj, dx_path):
    """Create (once) a surface object colored by a real APBS electrostatic
    potential map (kT/e), returning its name ('{obj}_esp').

    Uses pdb2pqr + APBS (prep_apbs.py, run beforehand as a subprocess) --
    NOT PyMOL's util.protein_vacuum_esp, which has a reproducible internal
    bug making it unusable in this scripted batch context (its derived
    surface object silently fails to render, with no exception and no
    visible warning, depending on unrelated prior cmd.set_color()/
    cmd.disable() calls elsewhere in the same session -- not worth chasing
    further when a standard, independent tool does the same job reliably)."""
    esp_obj = f'{obj}_esp'
    if esp_obj not in cmd.get_names('objects'):
        cmd.create(esp_obj, obj)
        map_obj = f'{esp_obj}_map'
        cmd.load(dx_path, map_obj)
        cmd.hide('everything', esp_obj)
        cmd.show('surface', esp_obj)
        cmd.ramp_new(f'{esp_obj}_ramp', map_obj, [-5, 0, 5], ['red', 'white', 'blue'])
        cmd.set('surface_color', f'{esp_obj}_ramp', esp_obj)
    return esp_obj


def render_electrostatic(obj, dx_path, out_dir):
    log("Rendering APBS electrostatic-potential surface...")
    esp_obj = ensure_esp_surface(obj, dx_path)
    show_only(esp_obj)
    snap(esp_obj, out_dir / 'electrostatic_surface_front.png')
    cmd.turn('y', 180)
    cmd.ray(*IMG_SIZE)
    cmd.png(str(out_dir / 'electrostatic_surface_back.png'), dpi=200)
    log("  saved electrostatic_surface_back.png")
    return esp_obj


def render_hydrophobic(obj, out_dir):
    log("Rendering Kyte-Doolittle hydrophobic surface...")
    tmp = f'{obj}_hydro'
    if tmp not in cmd.get_names('objects'):
        cmd.create(tmp, obj)
        cmd.hide('everything', tmp)
        cmd.show('surface', tmp)
        for resi, resn in get_residue_table(tmp):
            kd = sgl.KYTE_DOOLITTLE.get(resn, 0.0)
            cmd.alter(f'{tmp} and resi {resi}', f'b={kd}')
        cmd.spectrum('b', 'cyan_white_orange', tmp, -4.5, 4.5)
    show_only(tmp)
    snap(tmp, out_dir / 'hydrophobic_surface_front.png')
    cmd.turn('y', 180)
    cmd.ray(*IMG_SIZE)
    cmd.png(str(out_dir / 'hydrophobic_surface_back.png'), dpi=200)
    log("  saved hydrophobic_surface_back.png")
    return tmp


def classify_active_face(obj, domain_range, catalytic_resids):
    """Dot-product projection: (CA - domain_COM) . unit(active_COM - domain_COM).
    Same concept as sim_analysis/figure_sasa_faces.py::classify_face."""
    d_start, d_end = domain_range
    coords = {}
    for r in range(d_start, d_end + 1):
        sel = f'{obj} and resi {r} and name CA'
        if cmd.count_atoms(sel) == 1:
            coords[r] = np.array(cmd.get_atom_coords(sel))
    if not coords:
        return {}
    d_com = np.mean(list(coords.values()), axis=0)
    act_coords = [coords[r] for r in catalytic_resids if r in coords]
    if not act_coords:
        return {}
    act_com = np.mean(act_coords, axis=0)
    dir_vec = act_com - d_com
    norm = np.linalg.norm(dir_vec)
    if norm < 1e-6:
        return {}
    dir_unit = dir_vec / norm
    return {r: ('active' if np.dot(c - d_com, dir_unit) > 0 else 'back')
            for r, c in coords.items()}


def classify_generic_face(obj, domain_range, seed):
    """k-means (k=2) split of a domain's CA coordinates into two
    geometric surface faces, for domains with no defined active site."""
    d_start, d_end = domain_range
    resids, coords = [], []
    for r in range(d_start, d_end + 1):
        sel = f'{obj} and resi {r} and name CA'
        if cmd.count_atoms(sel) == 1:
            resids.append(r)
            coords.append(cmd.get_atom_coords(sel))
    if len(resids) < 4:
        return {}
    coords = np.array(coords)
    rng = np.random.RandomState(seed)
    _, labels = kmeans2(coords, 2, seed=seed, minit='++')
    return {r: ('faceA' if lbl == 0 else 'faceB')
            for r, lbl in zip(resids, labels)}


def render_face_classification(obj, set_key, out_dir):
    log("Classifying + rendering domain faces...")
    tmp = f'{obj}_face'
    cmd.create(tmp, obj)
    cmd.hide('everything', tmp)
    cmd.show('surface', tmp)
    cmd.color('gray90', tmp)

    domains = sgl.get_domain_ranges_for_set(set_key)
    for i, (dname, (color_name, base_rgb, (cs, ce))) in enumerate(domains.items()):
        site = DOMAIN_TO_SITE.get(dname)
        faces = {}
        if site:
            cat, _ = sgl.get_pocket_for_set(set_key, site)
            if cat:
                faces = classify_active_face(tmp, (cs, ce), cat)
                for r, f in faces.items():
                    cmd.color(hexcol(FACE_COLORS_ACTIVE[f]), f'{tmp} and resi {r}')
        if not faces:
            faces = classify_generic_face(tmp, (cs, ce), seed=42 + i)
            # color_name is already a sanitized identifier (e.g. 'kh16_col');
            # dname itself may contain hyphens ('KH1-KH6') which PyMOL color
            # names can't hold.
            colA, colB = f'{color_name}_faceA', f'{color_name}_faceB'
            cmd.set_color(colA, list(base_rgb))
            muted = tuple(min(1.0, 0.55 + 0.45 * c) for c in base_rgb)
            cmd.set_color(colB, list(muted))
            for r, f in faces.items():
                cmd.color(colA if f == 'faceA' else colB, f'{tmp} and resi {r}')

    show_only(tmp)
    snap(tmp, out_dir / 'face_classification_front.png', buffer=5)
    cmd.turn('y', 180)
    cmd.ray(*IMG_SIZE)
    cmd.png(str(out_dir / 'face_classification_back.png'), dpi=200)
    log("  saved face_classification_back.png")
    return tmp


def render_pocket_views(obj, set_key, site, out_dir, dx_path):
    cat, pocket = sgl.get_pocket_for_set(set_key, site)
    if not pocket:
        return
    log(f"Rendering {site} pocket close-ups ({len(pocket)} residues)...")
    pocket_sel_str = "+".join(map(str, pocket))
    site_dir = out_dir / f'pocket_{site}'
    site_dir.mkdir(parents=True, exist_ok=True)

    # 1. Electrostatics: a fresh surface object per site, not a reuse of the
    #    whole-molecule one -- reusing it here rendered blank once enough
    #    other cmd.create()/cmd.color() activity (hydrophobic, face
    #    classification) had happened first. Skip entirely if APBS prep
    #    failed for this state.
    if dx_path:
        esp_obj_site = f'{obj}_esp_{site}'
        cmd.create(esp_obj_site, obj)
        cmd.load(dx_path, f'{esp_obj_site}_map')
        cmd.hide('everything', esp_obj_site)
        cmd.show('surface', esp_obj_site)
        cmd.ramp_new(f'{esp_obj_site}_ramp', f'{esp_obj_site}_map', [-5, 0, 5], ['red', 'white', 'blue'])
        cmd.set('surface_color', f'{esp_obj_site}_ramp', esp_obj_site)
        show_only(esp_obj_site)
        snap(f'{esp_obj_site} and resi {pocket_sel_str}', site_dir / 'electrostatic.png', buffer=4)
    else:
        log(f"  electrostatic SKIPPED for pocket_{site} (no APBS map)")

    # 2. Hydrophobicity: reuse the whole-molecule hydrophobic object if present.
    tmp_h = f'{obj}_hydro'
    if tmp_h not in cmd.get_names('objects'):
        cmd.create(tmp_h, obj)
        cmd.hide('everything', tmp_h)
        cmd.show('surface', tmp_h)
        for resi, resn in get_residue_table(tmp_h):
            cmd.alter(f'{tmp_h} and resi {resi}', f'b={sgl.KYTE_DOOLITTLE.get(resn, 0.0)}')
        cmd.spectrum('b', 'cyan_white_orange', tmp_h, -4.5, 4.5)
    show_only(tmp_h)
    snap(f'{tmp_h} and resi {pocket_sel_str}', site_dir / 'hydrophobicity.png', buffer=4)

    # 3. Discrete residue-chemistry classification (cartoon+sticks; a surface
    #    would hide the per-residue class detail at this zoom).
    tmp_c = f'{obj}_class_{site}'
    cmd.create(tmp_c, obj)
    cmd.hide('everything', tmp_c)
    cmd.show('cartoon', tmp_c)
    cmd.show('sticks', f'{tmp_c} and resi {pocket_sel_str} and not name C+N+O')
    cmd.color('gray70', tmp_c)
    for resi, resn in get_residue_table(tmp_c):
        if resi in pocket:
            _, color = sgl.RESIDUE_CLASS.get(resn, ('nonpolar', '#dcdcdc'))
            cmd.color(hexcol(color), f'{tmp_c} and resi {resi}')
    if cat:
        cat_sel = f'{tmp_c} and resi {"+".join(map(str, cat))} and name CA'
        cmd.show('spheres', cat_sel)
        cmd.set('sphere_scale', 0.5, cat_sel)
        cmd.color('black', cat_sel)
    show_only(tmp_c)
    snap(f'{tmp_c} and resi {pocket_sel_str}', site_dir / 'residue_class.png', buffer=4)

    # 4. H-bond donor/acceptor atom mapping.
    tmp_hb = f'{obj}_hbond_{site}'
    cmd.create(tmp_hb, obj)
    cmd.hide('everything', tmp_hb)
    cmd.show('cartoon', tmp_hb)
    cmd.color('gray80', tmp_hb)
    cmd.show('sticks', f'{tmp_hb} and resi {pocket_sel_str} and not name C+N+O')
    cmd.color('wheat', f'{tmp_hb} and resi {pocket_sel_str}')
    for resi, resn in get_residue_table(tmp_hb):
        if resi not in pocket:
            continue
        if resn != 'PRO':
            sel = f'{tmp_hb} and resi {resi} and name N'
            cmd.color(hexcol(sgl.DONOR_COLOR), sel)
            cmd.show('spheres', sel)
        sel = f'{tmp_hb} and resi {resi} and name O'
        cmd.color(hexcol(sgl.ACCEPTOR_COLOR), sel)
        cmd.show('spheres', sel)
        for atom_name, role in sgl.SIDECHAIN_HBOND_ATOMS.get(resn, []):
            sel = f'{tmp_hb} and resi {resi} and name {atom_name}'
            if cmd.count_atoms(sel) == 1:
                cmd.color(hexcol(sgl.HBOND_ROLE_COLOR[role]), sel)
                cmd.show('spheres', sel)
    cmd.set('sphere_scale', 0.35, f'{tmp_hb} and resi {pocket_sel_str}')
    show_only(tmp_hb)
    snap(f'{tmp_hb} and resi {pocket_sel_str}', site_dir / 'hbond_donor_acceptor.png', buffer=4)


def prep_apbs_maps(pdb_path, work_dir):
    """Run prep_apbs.py (pdb2pqr + APBS, via SBgrid) as a subprocess.
    Returns (centered_pdb_path, dx_path), or (None, None) if it failed --
    callers should skip electrostatics and continue with everything else
    rather than retry (SBgrid tools are already the stable path)."""
    script = Path(__file__).resolve().parent / 'prep_apbs.py'
    result = subprocess.run(
        [sys.executable, str(script), '--pdb', str(pdb_path), '--out-dir', str(work_dir)],
        capture_output=True, text=True)
    if result.returncode != 0:
        log(f"  APBS prep FAILED: {result.stderr.strip()[-500:]}")
        return None, None
    centered_pdb = dx = None
    for line in result.stdout.splitlines():
        if line.startswith('centered_pdb='):
            centered_pdb = line.split('=', 1)[1]
        elif line.startswith('dx='):
            dx = line.split('=', 1)[1]
    return centered_pdb, dx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set-key', required=True,
                     help="Registry set key for FL<->construct remapping "
                          "(e.g. fl_optimized, md, core, mka, norrm, noart)")
    ap.add_argument('--pdb', required=True)
    ap.add_argument('--label', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--pse-out', default=None)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    log(f"=== {args.label} ({args.set_key}) -> {out_dir} ===")

    log("Preparing APBS electrostatic potential map (pdb2pqr + APBS)...")
    centered_pdb, dx_path = prep_apbs_maps(args.pdb, out_dir / '_apbs_work')
    load_pdb = centered_pdb or args.pdb
    if not dx_path:
        log("  electrostatics will be SKIPPED for this state (APBS prep failed)")

    setup_common()
    cmd.load(load_pdb, 'mol')
    cmd.remove('solvent')
    log(f"Loaded {cmd.count_atoms('mol')} atoms")

    try:
        render_domain_cartoon('mol', args.set_key, out_dir)
    except Exception as e:
        log(f"  domain_cartoon FAILED: {e}")

    if dx_path:
        try:
            render_electrostatic('mol', dx_path, out_dir)
        except Exception as e:
            log(f"  electrostatic FAILED: {e}")

    try:
        render_hydrophobic('mol', out_dir)
    except Exception as e:
        log(f"  hydrophobic FAILED: {e}")

    try:
        render_face_classification('mol', args.set_key, out_dir)
    except Exception as e:
        log(f"  face_classification FAILED: {e}")

    for site in sr.SITE_NAMES:
        try:
            render_pocket_views('mol', args.set_key, site, out_dir, dx_path)
        except Exception as e:
            log(f"  pocket {site} FAILED: {e}")

    if args.pse_out:
        Path(args.pse_out).parent.mkdir(parents=True, exist_ok=True)
        cmd.save(args.pse_out)
        log(f"Saved session: {args.pse_out}")

    log(f"=== DONE {args.label} in {time.time() - t0:.0f}s ===")


if __name__ == '__main__':
    main()
