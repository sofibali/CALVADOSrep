#!/usr/bin/env python
"""
Prepare BindCraft targets + hotspots for PARP14 macrodomain (MD1-MD3) binder design.

Design goals (4):
  md1_block, md2_block, md3_block  -> block a single MD's ADP-ribose pocket
  md2md3_lock                      -> clamp the (transient) MD2-MD3 closed interface

Structural sources (2 per goal -> 8 BindCraft runs):
  *_af3 : coordinates from the clash-free AF3 md1l1_md2_md3 model (sim-validated
          arrangement; MD1-MD2 are permanently associated in the trajectories).
  *_sim : the CALVADOS-observed conformer for that goal, built by grafting AF3
          all-atom domains onto the backmapped CA frame (per-domain Kabsch on the
          structured core), then RIGID-BODY DECLASHED -- interpenetrating domains
          (a one-bead-per-residue artifact) are separated along their COM axis to
          van-der-Waals contact while preserving each domain's orientation.

Each target is truncated to the design-relevant domains to keep the AF2 multimer
prediction in BindCraft tractable:
  md1_block  : MD1+MD2  (local 1-404)   -- MD1 pocket; MD2 is its constant partner
  md2_block  : MD1+MD2+MD3 (1-586)      -- MD2 is central, flanked both sides
  md3_block  : MD2+MD3  (216-586)
  md2md3_lock: MD2+MD3  (216-586)

Construct local numbering 1-586 = MD1L1(FL 790-1004)+MD2(FL 1005-1193)+MD3(FL 1207-1388).
Residue numbers are PRESERVED (not renumbered) so hotspot numbers are consistent
across truncations and sources.

Run in the CALVADOS env:
    /home/sbali/miniconda3/envs/CALVADOS/bin/python prepare_targets.py
"""
import os
import csv
import json
import copy
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PDBIO
from Bio.PDB import Structure, Model, Chain
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.SVDSuperimposer import SVDSuperimposer

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT     = "/home/sbali/CALVADOS/examples/PARP14_MDP"
HERE     = os.path.join(ROOT, "bindcraft_md")
BACKMAP  = os.path.join(ROOT, "backmapped_states")
AF3_CIF  = "/home/sbali/CALVADOS/parp14/alphafold_outputs/md1l1_md2_md3/md1_md2_md3_model.cif"
SASA_CSV = os.path.join(ROOT, "data", "sasa_face_per_residue.csv")
BC_REPO  = "/home/sbali/Projects/BindCraft"

TARGET_DIR   = os.path.join(HERE, "targets")
SETTINGS_DIR = os.path.join(HERE, "settings")
DESIGN_ROOT  = os.path.join(HERE, "designs")
for d in (TARGET_DIR, SETTINGS_DIR):
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# Construct geometry
# ---------------------------------------------------------------------------
BLOCKS = {"MD1": (1, 215), "MD2": (216, 404), "MD3": (405, 586)}   # atom-assignment
CORES  = {"MD1": (2, 189), "MD2": (214, 401), "MD3": (414, 585)}   # rigid-fit core

def fl_to_local(fl):
    if 790 <= fl <= 1004:  return fl - 789
    if 1005 <= fl <= 1193: return fl - 789
    if 1207 <= fl <= 1388: return fl - 802
    return None

CATALYTIC_FL = {"MD1": [831, 923, 962],
                "MD2": [1035, 1046, 1134, 1171],
                "MD3": [1248, 1259, 1330, 1371]}
SASA_DOMAIN  = {"MD1": "MD1L1", "MD2": "MD2", "MD3": "MD3"}

def block_of(resid):
    for k, (lo, hi) in BLOCKS.items():
        if lo <= resid <= hi:
            return k
    return None

# ---------------------------------------------------------------------------
# Structure IO
# ---------------------------------------------------------------------------
def af3_residues():
    """Fresh list of AF3 all-atom residues (in order, 586), renumbered 1..586."""
    st = MMCIFParser(QUIET=True).get_structure("af3", AF3_CIF)
    m = next(st.get_models())
    res = [r for ch in m for r in ch if r.id[0] == " "][:586]
    res = copy.deepcopy(res)
    for i, r in enumerate(res, start=1):
        r.id = (" ", i, " ")
    ca = np.array([r["CA"].get_coord() for r in res])
    return ca, res

def conformer_ca(pdb_path):
    st = PDBParser(QUIET=True).get_structure("c", pdb_path)
    m = next(st.get_models())
    res = [r for ch in m for r in ch if r.id[0] == " "][:586]
    return np.array([r["CA"].get_coord() for r in res])

def assemble(residues, keep_lo, keep_hi, out_pdb):
    """Write residues with local number in [keep_lo, keep_hi] as chain A."""
    s = Structure.Structure("t"); m = Model.Model(0); c = Chain.Chain("A")
    for r in residues:
        if keep_lo <= r.id[1] <= keep_hi:
            c.add(r)
    m.add(c); s.add(m)
    io = PDBIO(); io.set_structure(s); io.save(out_pdb)

# ---------------------------------------------------------------------------
# Sim target: graft AF3 domains onto backmapped CA frame
# ---------------------------------------------------------------------------
def graft(af3_ca, af3_res, conf_ca):
    sup = SVDSuperimposer()
    rms = {}
    for name, (blo, bhi) in BLOCKS.items():
        clo, chi = CORES[name]
        sup.set(conf_ca[clo - 1:chi], af3_ca[clo - 1:chi])   # ref=conf, mov=af3 core
        sup.run()
        rot, tran = sup.get_rotran(); rms[name] = sup.get_rms()
        for r in af3_res[blo - 1:bhi]:
            for a in r:
                a.set_coord(np.dot(a.get_coord(), rot) + tran)
    return rms

# ---------------------------------------------------------------------------
# Rigid-body declash: separate interpenetrating domains to vdW contact
# ---------------------------------------------------------------------------
def declash(residues, present_domains, fixed="MD2", clash_cut=2.6, step=0.3, max_iter=200):
    """Translate each non-fixed domain outward along its COM->fixed-COM axis until
    no heavy atom is within clash_cut of another domain. Preserves orientation."""
    def dom_atoms(dom):
        lo, hi = BLOCKS[dom]
        return [a for r in residues if lo <= r.id[1] <= hi for a in r if a.element != "H"]
    def dom_com(dom):
        xs = np.array([a.get_coord() for a in dom_atoms(dom)])
        return xs.mean(0)
    if fixed not in present_domains:
        fixed = present_domains[0]
    movers = [d for d in present_domains if d != fixed]
    fixed_com = dom_com(fixed)
    for dom in movers:
        axis = dom_com(dom) - fixed_com
        n = np.linalg.norm(axis)
        axis = axis / n if n > 1e-6 else np.array([1.0, 0, 0])
        others = [d for d in present_domains if d != dom]
        other_atoms = [a for d in others for a in dom_atoms(d)]
        for _ in range(max_iter):
            ns = NeighborSearch(dom_atoms(dom) + other_atoms)
            clash = False
            for a, b in ns.search_all(clash_cut, "A"):
                da, db = block_of(a.get_parent().id[1]), block_of(b.get_parent().id[1])
                if da != db and da in present_domains and db in present_domains:
                    if dom in (da, db):
                        clash = True; break
            if not clash:
                break
            for a in dom_atoms(dom):
                a.set_coord(a.get_coord() + axis * step)
    # final min inter-domain heavy-atom distance
    allat = [a for d in present_domains for a in dom_atoms(d)]
    ns = NeighborSearch(allat); mind = 99.0
    for a, b in ns.search_all(8.0, "A"):
        da, db = block_of(a.get_parent().id[1]), block_of(b.get_parent().id[1])
        if da != db:
            mind = min(mind, np.linalg.norm(a.get_coord() - b.get_coord()))
    return mind

# ---------------------------------------------------------------------------
# Hotspots
# ---------------------------------------------------------------------------
def load_sasa():
    with open(SASA_CSV) as fh:
        return list(csv.DictReader(fh))

def block_hotspots(md, sasa, rsa_cut=0.15):
    local = set()
    for fl in CATALYTIC_FL[md]:
        loc = fl_to_local(fl)
        if loc: local.add(loc)
    for r in sasa:
        if r["domain"] != SASA_DOMAIN[md] or not r["is_pocket_in"]:
            continue
        try: rsa = float(r["rsa"])
        except ValueError: continue
        if rsa >= rsa_cut and r.get("face", "") != "back":
            loc = fl_to_local(int(r["resid"]))
            if loc: local.add(loc)
    return sorted(local)

def clamp_hotspots(residues, mdA, mdB, per_side=7, cut=10.0):
    """Contact-patch hotspots for a clamp: on the grafted+declashed state pose, the
    closest `per_side` residues of each domain core to the other domain (the surfaces
    facing each other across the cleft). Returns sorted local residue numbers."""
    def core_ca(md):
        lo, hi = CORES[md]
        return [(r.id[1], r["CA"].get_coord()) for r in residues if lo <= r.id[1] <= hi and "CA" in r]
    A, B = core_ca(mdA), core_ca(mdB)
    Aid = np.array([x[0] for x in A]); Axyz = np.array([x[1] for x in A])
    Bid = np.array([x[0] for x in B]); Bxyz = np.array([x[1] for x in B])
    D = np.linalg.norm(Axyz[:, None] - Bxyz[None], axis=2) / 1.0  # Angstrom
    aMin, bMin = D.min(1), D.min(0)
    aSel = [int(Aid[i]) for i in np.argsort(aMin)[:per_side] if aMin[i] < cut * 10]
    bSel = [int(Bid[i]) for i in np.argsort(bMin)[:per_side] if bMin[i] < cut * 10]
    return sorted(set(aSel + bSel))

def interface_hotspots(af3_ca, mdA, mdB, cut=8.0, max_each=8):
    # use structured CORES only: excludes the artificial MD2/MD3 sequence junction
    # (FL linker 1194-1206 deleted in the construct) and flexible termini, so we
    # capture genuine domain-body packing rather than chain-adjacency artifacts.
    aLo, aHi = CORES[mdA]; bLo, bHi = CORES[mdB]
    A = np.arange(aLo - 1, aHi); B = np.arange(bLo - 1, bHi)
    d = np.linalg.norm(af3_ca[A][:, None] - af3_ca[B][None], axis=2)
    contacts = d < cut
    ac, bc = contacts.sum(1), contacts.sum(0)
    aSel = [int(A[i]) + 1 for i in np.argsort(-ac)[:max_each] if ac[i] > 0]
    bSel = [int(B[i]) + 1 for i in np.argsort(-bc)[:max_each] if bc[i] > 0]
    return sorted(set(aSel + bSel))

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
def write_settings(name, pdb, hotspots, lengths, n_final=8):
    cfg = {"design_path": os.path.join(DESIGN_ROOT, name) + "/",
           "binder_name": name, "starting_pdb": pdb, "chains": "A",
           "target_hotspot_residues": ",".join(map(str, hotspots)),
           "lengths": lengths, "number_of_final_designs": n_final}
    out = os.path.join(SETTINGS_DIR, name + ".json")
    with open(out, "w") as fh:
        json.dump(cfg, fh, indent=4)
    return out

# ---------------------------------------------------------------------------
# Goal table
# ---------------------------------------------------------------------------
GOALS = {
    "md1_block":   dict(kind="block", md="MD1", keep=(1, 404),   conf="md_MD1_exposed_allatom.pdb", present=["MD1","MD2"],       lengths=[70,130]),
    "md2_block":   dict(kind="block", md="MD2", keep=(1, 586),   conf="md_MD2_exposed_allatom.pdb", present=["MD1","MD2","MD3"], lengths=[70,130]),
    "md3_block":   dict(kind="block", md="MD3", keep=(216, 586), conf="md_MD3_exposed_allatom.pdb", present=["MD2","MD3"],       lengths=[70,130]),
}
# Interface "locking" is handled by the CLAMP targets below (enforce a chosen TICA-state
# juxtaposition), which supersede the old contact-frequency-based lock approach.

def main():
    sasa = load_sasa()
    af3_ca_ref, _ = af3_residues()

    rows = []
    for name, g in GOALS.items():
        lo, hi = g["keep"]

        if g["kind"] == "block":
            # well-defined pocket -> build BOTH sources with the same hotspots
            hot = block_hotspots(g["md"], sasa); desc = f"BLOCK {g['md']} pocket"

            _, af3_res = af3_residues()
            af3_pdb = os.path.join(TARGET_DIR, f"{name}_af3.pdb")
            assemble(af3_res, lo, hi, af3_pdb)
            write_settings(f"{name}_af3", af3_pdb, hot, g["lengths"])

            af3_ca2, af3_res2 = af3_residues()
            rms = graft(af3_ca2, af3_res2, conformer_ca(os.path.join(BACKMAP, g["conf"])))
            mind = declash(af3_res2, g["present"])
            sim_pdb = os.path.join(TARGET_DIR, f"{name}_sim.pdb")
            assemble(af3_res2, lo, hi, sim_pdb)
            write_settings(f"{name}_sim", sim_pdb, hot, g["lengths"])

            print(f"[{name}] {desc} | keep {lo}-{hi} | {len(hot)} hotspots")
            print(f"    af3: {af3_pdb}")
            print(f"    sim: {sim_pdb}  (graft {','.join(f'{k}:{v:.2f}' for k,v in rms.items())}; declash min {mind:.2f} A)")
            print(f"    hotspots(local): {hot}\n")
            for src, pdb in (("af3", af3_pdb), ("sim", sim_pdb)):
                rows.append(dict(target=f"{name}_{src}", goal=desc, source=src,
                                 keep=f"{lo}-{hi}", n_hotspots=len(hot),
                                 hotspots=" ".join(map(str, hot)),
                                 lengths=f"{g['lengths'][0]}-{g['lengths'][1]}", pdb=pdb))
        else:
            # LOCK: the interface exists ONLY in the sim compact pose (AF3 has the
            # MD2/MD3 bodies apart). Build sim-only; hotspots from the actual sim
            # interface (core-core contacts on the declashed pose).
            mdA, mdB = g["pair"]; desc = f"LOCK {mdA}-{mdB} (transient; sim-only)"
            af3_ca2, af3_res2 = af3_residues()
            rms = graft(af3_ca2, af3_res2, conformer_ca(os.path.join(BACKMAP, g["conf"])))
            mind = declash(af3_res2, g["present"])
            sim_ca = np.array([r["CA"].get_coord() for r in af3_res2])
            hot = interface_hotspots(sim_ca, mdA, mdB, cut=10.0)
            sim_pdb = os.path.join(TARGET_DIR, f"{name}_sim.pdb")
            assemble(af3_res2, lo, hi, sim_pdb)
            write_settings(f"{name}_sim", sim_pdb, hot, g["lengths"])
            print(f"[{name}] {desc} | keep {lo}-{hi} | {len(hot)} hotspots")
            print(f"    sim: {sim_pdb}  (graft {','.join(f'{k}:{v:.2f}' for k,v in rms.items())}; declash min {mind:.2f} A)")
            print(f"    hotspots(local): {hot}\n")
            rows.append(dict(target=f"{name}_sim", goal=desc, source="sim",
                             keep=f"{lo}-{hi}", n_hotspots=len(hot),
                             hotspots=" ".join(map(str, hot)),
                             lengths=f"{g['lengths'][0]}-{g['lengths'][1]}", pdb=sim_pdb))

    # -----------------------------------------------------------------------
    # CLAMP targets: enforce a NEW interface from TICA states (md_ca25_tica 05-20).
    # Sim-pose-only (the pose IS the state). Hotspots = the state's contact patch on
    # BOTH domains -> binder bridges the cleft to clamp that juxtaposition.
    # state1 MD2-MD3 is open (no facing surface) -> skipped.
    # -----------------------------------------------------------------------
    STATE_DIR = os.path.join(ROOT, "representative_frames", "2026-05-20", "md_ca25_tica")
    CLAMP_KEEP = {("MD1", "MD2"): (1, 404), ("MD2", "MD3"): (216, 586)}
    CLAMP_PRESENT = {("MD1", "MD2"): ["MD1", "MD2"], ("MD2", "MD3"): ["MD2", "MD3"]}
    CLAMPS = [
        ("state3", "state_3.pdb", [("MD1", "MD2"), ("MD2", "MD3")]),
        ("state5", "state_5.pdb", [("MD1", "MD2"), ("MD2", "MD3")]),
        ("state1", "state_1.pdb", [("MD1", "MD2")]),   # MD2-MD3 open -> skip
    ]
    for sname, sfile, prs in CLAMPS:
        spath = os.path.join(STATE_DIR, sfile)
        for (mdA, mdB) in prs:
            af3_ca2, af3_res2 = af3_residues()
            graft(af3_ca2, af3_res2, conformer_ca(spath))
            present = CLAMP_PRESENT[(mdA, mdB)]
            mind = declash(af3_res2, present)
            lo, hi = CLAMP_KEEP[(mdA, mdB)]
            hot = [h for h in clamp_hotspots(af3_res2, mdA, mdB) if lo <= h <= hi]
            name = f"clamp_{mdA}{mdB}_{sname}".replace("MD", "md").lower()
            pdb = os.path.join(TARGET_DIR, name + ".pdb")
            assemble(af3_res2, lo, hi, pdb)
            write_settings(name, pdb, hot, [90, 150])
            desc = f"CLAMP {mdA}-{mdB} @ {sname} (enforce new interface)"
            print(f"[{name}] {desc} | keep {lo}-{hi} | declash min {mind:.2f} A")
            print(f"    hotspots(local): {hot}\n")
            rows.append(dict(target=name, goal=desc, source="sim-state",
                             keep=f"{lo}-{hi}", n_hotspots=len(hot),
                             hotspots=" ".join(map(str, hot)), lengths="90-150", pdb=pdb))

    with open(os.path.join(HERE, "targets_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"Wrote {len(rows)} targets -> targets_summary.csv")


if __name__ == "__main__":
    main()
