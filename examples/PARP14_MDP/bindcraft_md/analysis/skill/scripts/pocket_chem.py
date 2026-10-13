#!/usr/bin/env python3
"""Pre-screen a hotspot patch BEFORE spending GPU-days on it.

    python pocket_chem.py <settings.json> [<settings.json> ...]
    python pocket_chem.py --pdb target.pdb --hotspots 112,115,119 [--chain A]

Reads target_hotspot_residues + starting_pdb from a BindCraft settings file.
Apolar fraction of the patch is the single strongest predictor of designability
observed so far: see the calibration table printed at the end.
"""
import argparse, json, os

# (donors, acceptors, formal charge, class)
PROP = {'ARG':(5,0,+1,'pos'),'LYS':(3,0,+1,'pos'),'HIS':(1,1,0,'polar'),
        'ASP':(0,4,-1,'neg'),'GLU':(0,4,-1,'neg'),'ASN':(2,2,0,'polar'),
        'GLN':(2,2,0,'polar'),'SER':(1,2,0,'polar'),'THR':(1,2,0,'polar'),
        'TYR':(1,1,0,'arom'),'TRP':(1,0,0,'arom'),'CYS':(1,0,0,'polar'),
        'PHE':(0,0,0,'arom'),'ALA':(0,0,0,'apol'),'VAL':(0,0,0,'apol'),
        'LEU':(0,0,0,'apol'),'ILE':(0,0,0,'apol'),'MET':(0,0,0,'apol'),
        'PRO':(0,0,0,'apol'),'GLY':(0,0,0,'gly')}

CALIBRATION = """
  calibration (PARP14 macrodomain campaign, BindCraft 1 + predict_initial_guess)
    apolar %   patch                 relaxed  scored  accepted
      43.8     MD1 block                 12      51       9
      42.9     MD1-MD2 clamp state3      50     197       5
      38.1     MD2 block                 50      16       0
      28.6     MD1-MD2 clamp state1      28      49       0
      28.6     MD3 block                 50       0       0
  Behaves like a THRESHOLD between 38.1% and 42.9%, not a gradient."""


def analyse(pdb, hotspots, chain, label):
    from Bio.PDB import PDBParser, ShrakeRupley
    st = PDBParser(QUIET=True).get_structure("t", pdb)
    ShrakeRupley().compute(st[0], level="R")
    res = {r.id[1]: r for r in st[0][chain] if r.id[0] == " "}

    d = a = q = 0; sasa = 0.0; cls = {}; n = 0; missing = []
    for h in hotspots:
        r = res.get(h)
        if r is None:
            missing.append(h); continue
        dd, aa, qq, cc = PROP.get(r.get_resname(), (0, 0, 0, 'apol'))
        d += dd; a += aa; q += qq
        cls[cc] = cls.get(cc, 0) + 1
        sasa += float(getattr(r, "sasa", 0.0)); n += 1

    if not n:
        print(f"\n{label}: NO hotspot matched the structure -- "
              "numbering mismatch (UniProt vs construct?)")
        print(f"  requested {hotspots}\n  available {min(res)}-{max(res)} on chain {chain}")
        return

    apol = cls.get('apol', 0) + cls.get('arom', 0)
    polar = cls.get('polar', 0) + cls.get('pos', 0) + cls.get('neg', 0)
    af = 100 * apol / n

    print(f"\n{'='*58}\n{label}\n{'='*58}")
    if missing:
        print(f"  !! {len(missing)} hotspot(s) absent from structure: {missing}")
    print(f"  hotspots matched     {n}")
    print(f"  apolar fraction      {af:.1f}%   (apol+arom {apol} / {n})")
    print(f"  polar fraction       {100*polar/n:.1f}%")
    print(f"  composition          " + ", ".join(f"{k}:{v}" for k, v in sorted(cls.items())))
    print(f"  net charge           {q:+d}")
    print(f"  H-bond donors/acc    {d}/{a}   ({(d+a)/n:.2f} per residue)")
    print(f"  patch SASA           {sasa:.0f} A^2   ({sasa/n:.0f} per residue)")

    if af >= 42:   v = "GOOD -- at or above every patch that yielded binders (>=42.9)"
    elif af > 38:  v = "MARGINAL -- in the untested gap between 38.1 and 42.9"
    else:          v = "POOR -- at or below every patch that produced zero binders"
    print(f"\n  VERDICT  {v}")
    if (d + a) / n > 2.5:
        print(f"  WARNING  {(d+a)/n:.2f} H-bond partners per residue: high burden. "
              "Expect\n           n_InterfaceUnsatHbonds (threshold 4) to be the limiter.")
    if abs(q) >= 4:
        print(f"  WARNING  net charge {q:+d}: binders will be driven to the opposite charge,")
        print( "           which tends to raise Surface_Hydrophobicity failures.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("settings", nargs="*")
    ap.add_argument("--pdb"); ap.add_argument("--hotspots"); ap.add_argument("--chain", default="A")
    a = ap.parse_args()
    if a.pdb:
        analyse(a.pdb, [int(x) for x in a.hotspots.split(",")], a.chain,
                f"{os.path.basename(a.pdb)}  {a.hotspots}")
    for s in a.settings:
        cfg = json.load(open(s))
        ch = (cfg.get("chains") or "A").split(",")[0].strip()
        analyse(cfg["starting_pdb"],
                [int(x) for x in cfg["target_hotspot_residues"].split(",")],
                ch, os.path.basename(s))
    print(CALIBRATION)


if __name__ == "__main__":
    main()
