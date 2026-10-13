#!/usr/bin/env python3
"""Where did a BindCraft campaign's designs die?

    python triage.py <run_dir> [<run_dir> ...] [--filters default_filters.json]

<run_dir> is a BindCraft 1 design_path: the folder holding Trajectory/, Accepted/,
mpnn_design_stats.csv, failure_csv.csv. Prints the true funnel (including aborts
that never produce a file), then which filter actually rejected each scored design.
"""
import argparse, json, os, re

DEF_FILTERS = "/home/sbali/BindCraft/settings_filters/default_filters.json"
# failure_csv.csv is ONE flat row; these are its trajectory-stage columns.
EARLY = ["Trajectory_logits_pLDDT", "Trajectory_softmax_pLDDT", "Trajectory_one-hot_pLDDT"]
AF2_GATES = ["pLDDT", "pTM", "i_pTM", "pAE", "i_pAE", "Clashes", "Binder_Energy_Score",
             "Surface_Hydrophobicity", "ShapeComplementarity", "dG", "dSASA",
             "n_InterfaceResidues", "n_InterfaceHbonds", "n_InterfaceUnsatHbonds",
             "Interface_Helix%", "Interface_BetaSheet%", "Interface_Loop%",
             "Binder_Helix%", "Binder_BetaSheet%", "Binder_Loop%",
             "InterfaceAAs", "Hotspot_RMSD", "Binder_pLDDT", "Binder_pTM",
             "Binder_pAE", "Binder_RMSD"]


def npdb(d):
    return len([f for f in os.listdir(d) if f.endswith(".pdb")]) if os.path.isdir(d) else 0


def base_metric(col):
    """Average_i_pAE -> i_pAE ; 3_ShapeComplementarity -> ShapeComplementarity"""
    return re.sub(r"^(Average|[1-5])_", "", col)


def triage(run, filters_path, verbose=False):
    import pandas as pd
    print(f"\n{'='*70}\n{os.path.basename(os.path.normpath(run))}\n{'='*70}")

    relaxed = npdb(os.path.join(run, "Trajectory", "Relaxed"))
    clash   = npdb(os.path.join(run, "Trajectory", "Clashing"))
    lowconf = npdb(os.path.join(run, "Trajectory", "LowConfidence"))
    acc     = npdb(os.path.join(run, "Accepted"))

    fc = os.path.join(run, "failure_csv.csv")
    f = pd.read_csv(fc).iloc[0].to_dict() if os.path.exists(fc) else {}
    geti = lambda k: int(f.get(k, 0) or 0)

    early = sum(geti(k) for k in EARLY)
    attempts = early + clash + lowconf + relaxed

    print("TRAJECTORY STAGE  (hallucination)")
    print(f"  attempts                        {attempts:>6}")
    print(f"    aborted logits/softmax/onehot {early:>6}   in failure_csv ONLY -- no folder, no file")
    print(f"    clashing                      {clash:>6}")
    print(f"    low confidence                {lowconf:>6}   = low final pLDDT ({geti('Trajectory_final_pLDDT')})"
          f" + too-few-contacts ({geti('Trajectory_Contacts')})")
    print(f"    relaxed -> MPNN               {relaxed:>6}")
    if attempts:
        print(f"  trajectory yield {100*relaxed/attempts:.0f}%"
              f"   (folder-only count would say {100*relaxed/max(1,clash+lowconf+relaxed):.0f}% -- overstated)")

    seqs = relaxed * 20  # num_seqs default; counters below are per-sequence events
    print("\nMPNN REDESIGN -> AF2 RE-PREDICTION  (per-SEQUENCE rejections, not per-design)")
    rej = [(g, geti(g)) for g in AF2_GATES if geti(g)]
    if rej:
        for g, v in sorted(rej, key=lambda x: -x[1]):
            pct = f"{100*v/seqs:5.1f}% of ~{seqs}" if seqs else "n/a"
            print(f"  {g:<24} {v:>6}   ({pct} seqs)")
    else:
        print("  (none recorded)")

    mp = os.path.join(run, "mpnn_design_stats.csv")
    d = pd.read_csv(mp) if os.path.exists(mp) and os.path.getsize(mp) > 0 else None
    scored = 0 if d is None else len(d)
    print(f"\n  reached Rosetta scoring   {scored}")
    print(f"  ACCEPTED                  {acc}")

    if not scored:
        if relaxed:
            print("\n  VERDICT: nothing survived AF2 re-prediction. The wall is confidence")
            print("           (i_pAE/i_pTM), not interface chemistry. Try predict_initial_guess.")
        return

    filt = json.load(open(filters_path))
    cols = {c: s for c, s in filt.items()
            if isinstance(s, dict) and s.get("threshold") is not None and c in d.columns
            and pd.api.types.is_numeric_dtype(d[c])}

    # per-design failure mask, aggregated to the base metric
    bymetric, perdesign = {}, pd.Series(0, index=d.index)
    for c, s in cols.items():
        bad = (d[c] < s["threshold"]) if s["higher"] else (d[c] > s["threshold"])
        bad = bad.fillna(False)
        m = base_metric(c)
        bymetric[m] = bymetric.get(m, pd.Series(False, index=d.index)) | bad
    for m, bad in bymetric.items():
        perdesign += bad.astype(int)

    print(f"\nFILTERS  ({len(cols)} numeric gates over Average_ and per-model 1_-5_)")
    print(f"  {'metric':<26} {'thr':>7}  {'fail':>10}  median")
    for m, bad in sorted(bymetric.items(), key=lambda x: -x[1].sum()):
        thr = filt.get(f"Average_{m}", {}).get("threshold")
        col = f"Average_{m}"
        med = d[col].median() if col in d.columns and pd.api.types.is_numeric_dtype(d[col]) else float("nan")
        flag = "  <<<" if bad.sum() and bad.sum() == max(b.sum() for b in bymetric.values()) else ""
        print(f"  {m:<26} {str(thr):>7}  {int(bad.sum()):>4}/{len(d)} {100*bad.mean():>4.0f}%  {med:>8.2f}{flag}")

    clean = int((perdesign == 0).sum())
    print(f"\n  failing 0 gates {clean}   exactly 1 {int((perdesign==1).sum())}"
          f"   2+ {int((perdesign>=2).sum())}")
    if clean != acc:
        print(f"  note: {clean} pass all numeric gates but {acc} are in Accepted/ --"
              " the rest lost to\n        non-numeric checks (InterfaceAAs composition, duplicate sequence).")
    worst = max(bymetric, key=lambda m: bymetric[m].sum())
    if bymetric[worst].sum():
        near = [m for m, b in bymetric.items() if b.sum() and m != worst]
        print(f"\n  VERDICT: limited by {worst}"
              f"{'  (then ' + ', '.join(sorted(near, key=lambda m: -bymetric[m].sum())[:2]) + ')' if near else ''}")
        print(f"  Metrics failing 0% are NOT the problem -- do not tune them.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--filters", default=DEF_FILTERS)
    a = ap.parse_args()
    for r in a.runs:
        triage(r, a.filters)


if __name__ == "__main__":
    main()
