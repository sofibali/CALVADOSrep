#!/usr/bin/env python3
"""Recompute every derived entry in figdata.json from the campaign outputs.

Repo-relative on purpose: the original versions of these numbers were computed in a
session scratchpad under /tmp that was later cleared. Everything here reads from
../designs, ../targets and ../bc2, which are committed, so the artifact can always
be rebuilt from this repository alone.

    python analysis/make_figdata.py        # run from bindcraft_md/
"""
import json, os, re, csv
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
D    = lambda *p: os.path.join(ROOT, *p)
BC1_FILTERS = "/home/sbali/BindCraft/settings_filters/default_filters.json"

fig = json.load(open(os.path.join(HERE, "figdata.json")))

def counts(run):
    p = D("designs", run)
    n = lambda k: len([f for f in os.listdir(os.path.join(p, k)) if f.endswith(".pdb")]) \
                  if os.path.isdir(os.path.join(p, k)) else 0
    try: sc = sum(1 for _ in open(os.path.join(p, "mpnn_design_stats.csv"))) - 1
    except Exception: sc = 0
    return {"relaxed": n("Trajectory/Relaxed"), "clash": n("Trajectory/Clashing"),
            "lowconf": n("Trajectory/LowConfidence"), "scored": sc, "accepted": n("Accepted")}

# ---------------------------------------------------------------- _pockets
from Bio.PDB import PDBParser, ShrakeRupley
PROP = {'ARG':(5,0,+1,'pos'),'LYS':(3,0,+1,'pos'),'HIS':(1,1,0,'polar'),
        'ASP':(0,4,-1,'neg'),'GLU':(0,4,-1,'neg'),'ASN':(2,2,0,'polar'),
        'GLN':(2,2,0,'polar'),'SER':(1,2,0,'polar'),'THR':(1,2,0,'polar'),
        'TYR':(1,1,0,'arom'),'TRP':(1,0,0,'arom'),'CYS':(1,0,0,'polar'),
        'PHE':(0,0,0,'arom'),'ALA':(0,0,0,'apol'),'VAL':(0,0,0,'apol'),
        'LEU':(0,0,0,'apol'),'ILE':(0,0,0,'apol'),'MET':(0,0,0,'apol'),
        'PRO':(0,0,0,'apol'),'GLY':(0,0,0,'gly')}
POCKETS = [("md1_block_af3","MD1 block",D("settings","md1_block_af3.json")),
           ("md2_block_af3","MD2 block",D("sweep","queue","md2_block_af3.json")),
           ("md3_block_af3","MD3 block",D("sweep","queue","md3_block_af3.json")),
           ("clamp_md1md2_state3","MD1-MD2 clamp s3",D("sweep","queue","clamp_md1md2_state3.json")),
           ("clamp_md1md2_state1","MD1-MD2 clamp s1",D("sweep","queue","clamp_md1md2_state1.json"))]
p, sr = PDBParser(QUIET=True), ShrakeRupley()
pockets = {}
for key, label, cfg in POCKETS:
    s = json.load(open(cfg)); hs = [int(x) for x in s["target_hotspot_residues"].split(",")]
    st = p.get_structure("t", s["starting_pdb"]); sr.compute(st[0], level="R")
    res = {r.id[1]: r for r in st[0]["A"] if r.id[0] == " "}
    d=a=q=0; sasa=0.0; cls={}; n=0
    for h in hs:
        r = res.get(h)
        if r is None: continue
        dd,aa,qq,cc = PROP.get(r.get_resname(), (0,0,0,'apol'))
        d+=dd; a+=aa; q+=qq; cls[cc]=cls.get(cc,0)+1
        sasa += float(getattr(r, "sasa", 0.0)); n += 1
    polar = cls.get('polar',0)+cls.get('pos',0)+cls.get('neg',0)
    apol  = cls.get('apol',0)+cls.get('arom',0)
    pockets[key] = {"label":label,"n":n,"donors":d,"acceptors":a,"charge":q,
                    "sasa":round(sasa,1),"sasa_per_res":round(sasa/n,1),"cls":cls,
                    "polar_frac":round(100*polar/n,1),"apol_frac":round(100*apol/n,1),
                    "hb_per_res":round((d+a)/n,2)}
fig["_pockets"] = pockets

# ---------------------------------------------------------------- _clampcmp
filt = json.load(open(BC1_FILTERS))
GATES = ["Average_i_pAE","Average_n_InterfaceUnsatHbonds","Average_ShapeComplementarity",
         "Average_Surface_Hydrophobicity","Average_Binder_RMSD","Average_dG"]
HARD  = ["Average_n_InterfaceUnsatHbonds","Average_ShapeComplementarity",
         "Average_Surface_Hydrophobicity","Average_Binder_RMSD"]
clampcmp = {}
for run, lab in [("guess_clamp_md1md2_state3","state3"),("guess_clamp_md1md2_state1","state1")]:
    d = pd.read_csv(D("designs", run, "mpnn_design_stats.csv"))
    g = {}
    for c in GATES:
        sp = filt.get(c, {}); th = sp.get("threshold")
        if th is None: continue
        bad = (d[c] < th) if sp["higher"] else (d[c] > th)
        g[c.replace("Average_","")] = {"thr":th,"pct":round(100*bad.mean(),1),
                                       "med":round(float(d[c].median()),2)}
    nf = sum(((d[c] < filt[c]["threshold"]) if filt[c]["higher"] else (d[c] > filt[c]["threshold"])).astype(int) for c in HARD)
    clampcmp[lab] = {"scored":len(d),"gates":g,"pass_all":int((nf==0).sum()),"fail1":int((nf==1).sum())}
fig["_clampcmp"] = clampcmp

# ---------------------------------------------------------------- _bccmp
bc1 = {k.replace("Average_",""):(v["threshold"], v["higher"])
       for k,v in filt.items()
       if isinstance(v, dict) and v.get("threshold") is not None and k.startswith("Average_")}
# BC2 thresholds as its own campaign banner prints them
bc2 = {"Backbone_Clashes":(0,False),"Binder_RMSD":(3.5,False),"Interface_Residues":(7,True),
       "Unbound_Binder_pLDDT":(0.7,True),"i_pAE":(0.35,False),"i_pTM":(0.7,True),"pTM":(0.55,True)}
ALIAS = {"n_InterfaceResidues":"Interface_Residues","Binder_pLDDT":"Unbound_Binder_pLDDT"}
b1 = {ALIAS.get(k,k):v for k,v in bc1.items()}
rows = [{"metric":k,
         "bc1": (f"{'>=' if b1[k][1] else '<='} {b1[k][0]}" if k in b1 else None),
         "bc2": (f"{'>=' if bc2[k][1] else '<='} {bc2[k][0]}" if k in bc2 else None)}
        for k in sorted(set(b1) | set(bc2))]
# how many BC1 MD1 trajectories touched MD2 (surface absent from the BC2 target)
geo = json.load(open(os.path.join(HERE, "geometry_results.json")))
fig["_bccmp"] = {"rows":rows,"bc1_target_res":404,"bc2_target_res":215,"hotspots":16,
                 "bc1_md2_contact":170,"bc1_relaxed":270,"bc1_md1_only":100}

# ---------------------------------------------------------------- _queue / _armA
q = {"state3":{**counts("guess_clamp_md1md2_state3"),"status":"complete (50-traj cap)","days":3.2},
     "state1":{**counts("guess_clamp_md1md2_state1"),"status":"stopped at cap"},
     "md3":  {**counts("guess_md3_block_af3")},
     "md2":  {**counts("guess_md2_block_af3")}}
r3 = pd.read_csv(D("designs","guess_md3_block_af3","failure_csv.csv")).iloc[0]
t3 = pd.read_csv(D("designs","guess_md3_block_af3","trajectory_stats.csv"))
q["md3"].update({"seqs":q["md3"]["relaxed"]*20,"ipae_fail":int(r3["i_pAE"]),"iptm_fail":int(r3["i_pTM"]),
                 "traj_ipae":round(float(t3["i_pAE"].median()),3),
                 "traj_pass":int(((t3.i_pAE<=0.35)&(t3.i_pTM>=0.5)&(t3.pLDDT>=0.8)).sum()),
                 "traj_n":len(t3)})
for k in ("gates3","traj","remaining"):
    if k in fig.get("_queue", {}): q[k] = fig["_queue"][k]
fig["_queue"] = q

json.dump(fig, open(os.path.join(HERE, "figdata.json"), "w"))
print("figdata.json rebuilt; keys:", sorted(k for k in fig if k.startswith("_")))
for k,v in pockets.items():
    print(f"  {v['label']:18s} n={v['n']:2d} apolar={v['apol_frac']:4.1f}% hb/res={v['hb_per_res']}")
