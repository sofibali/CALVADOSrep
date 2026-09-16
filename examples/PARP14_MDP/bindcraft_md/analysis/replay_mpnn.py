"""
Replay the MPNN + AF2-validation stage on ALREADY-GENERATED trajectories.

Hallucination (the expensive step, ~45 min/trajectory) is reused from disk, so this
answers the i_pAE question in hours instead of days -- and on a far bigger sample.

For each sampled trajectory it generates MPNN sequences ONCE, then predicts those SAME
sequences under both validation conditions:
    guess=OFF -> predict_initial_guess False  (reproduces the runs that failed)
    guess=ON  -> predict_initial_guess True   (the hardtarget fix)
That pairing is what makes the comparison strong: identical backbone, identical sequence,
only the validation protocol differs.

No filtering -- every metric is recorded, including the failures the real runs discard.

Usage:
  python replay_mpnn.py <target> [n_trajectories] [num_seqs]
"""
import sys, os, glob, re, csv, random

sys.path.insert(0, "/home/sbali/BindCraft")
os.chdir("/home/sbali/BindCraft")

import pandas as pd
from functions.generic_utils import load_json_settings, load_af2_models, perform_advanced_settings_check
from functions.colabdesign_utils import mpnn_gen_sequence
from colabdesign import mk_afdesign_model, clear_mem
from colabdesign.shared.utils import copy_dict

TARGET   = sys.argv[1] if len(sys.argv) > 1 else "md1_block_af3"
N_TRAJ   = int(sys.argv[2]) if len(sys.argv) > 2 else 40
NUM_SEQS = int(sys.argv[3]) if len(sys.argv) > 3 else 5

BASE = "/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md"
OUT  = f"/tmp/claude-64170/-home-sbali-BindCraft/adecb958-709c-481a-8e03-db4ccc10a3b9/scratchpad/replay_{TARGET}.csv"
SETTINGS = f"{BASE}/settings/{TARGET}.json"
FILTERS  = "/home/sbali/BindCraft/settings_filters/default_filters.json"
ADVANCED = "/home/sbali/BindCraft/settings_advanced/default_4stage_multimer.json"

target_settings, advanced_settings, filters = load_json_settings(SETTINGS, FILTERS, ADVANCED)
advanced_settings = perform_advanced_settings_check(advanced_settings, "/home/sbali/BindCraft")
advanced_settings["num_seqs"] = NUM_SEQS
design_models, prediction_models, multimer_validation = load_af2_models(advanced_settings["use_multimer_design"])

stats = pd.read_csv(f"{BASE}/designs/{TARGET}/trajectory_stats.csv")
iface = dict(zip(stats["Design"], stats["InterfaceResidues"]))

pdbs = sorted(glob.glob(f"{BASE}/designs/{TARGET}/Trajectory/Relaxed/*.pdb"))
random.seed(0)
random.shuffle(pdbs)
picked = []
for p in pdbs:
    name = os.path.basename(p)[:-4]
    if name in iface and isinstance(iface[name], str) and iface[name].strip():
        picked.append((name, p))
    if len(picked) >= N_TRAJ:
        break
print(f"[{TARGET}] replaying {len(picked)} trajectories x {NUM_SEQS} seqs, both conditions", flush=True)

# ---- PHASE 1: MPNN sequences (mpnn_gen_sequence calls clear_mem, so do these first) ----
jobs = []
for i, (name, pdb) in enumerate(picked, 1):
    length = int(re.search(r"_l(\d+)_s", name).group(1))
    try:
        seqs = mpnn_gen_sequence(pdb, "B", iface[name], advanced_settings)
        uniq = []
        for s in seqs["seq"]:
            s = s[-length:]
            if s not in uniq:
                uniq.append(s)
        jobs.append(dict(name=name, pdb=pdb, length=length, seqs=uniq))
        print(f"  [{i}/{len(picked)}] {name} len={length} -> {len(uniq)} seqs", flush=True)
    except Exception as e:
        print(f"  [{i}/{len(picked)}] {name} MPNN FAILED: {e}", flush=True)
print(f"MPNN done: {sum(len(j['seqs']) for j in jobs)} sequences total", flush=True)

# ---- PHASE 2: AF2 validation under both conditions ----
fh = open(OUT, "w", newline="")
w = csv.writer(fh)
w.writerow(["target","trajectory","length","seq_idx","guess","model","pLDDT","pTM","i_pTM","pAE","i_pAE"])

for guess in (False, True):
    jobs.sort(key=lambda j: j["length"])       # group by length -> one compile per length
    cur_len = None
    model = None
    for j in jobs:
        if j["length"] != cur_len:
            clear_mem()
            model = mk_afdesign_model(protocol="binder",
                                      num_recycles=advanced_settings["num_recycles_validation"],
                                      data_dir=advanced_settings["af_params_dir"],
                                      use_multimer=multimer_validation,
                                      use_initial_guess=guess, use_initial_atom_pos=False)
            cur_len = j["length"]
        if guess:
            model.prep_inputs(pdb_filename=j["pdb"], chain="A", binder_chain="B",
                              binder_len=j["length"], use_binder_template=True,
                              rm_target_seq=advanced_settings["rm_template_seq_predict"],
                              rm_target_sc=advanced_settings["rm_template_sc_predict"],
                              rm_template_ic=True)
        else:
            model.prep_inputs(pdb_filename=target_settings["starting_pdb"],
                              chain=target_settings["chains"], binder_len=j["length"],
                              rm_target_seq=advanced_settings["rm_template_seq_predict"],
                              rm_target_sc=advanced_settings["rm_template_sc_predict"])
        for si, seq in enumerate(j["seqs"]):
            for m in prediction_models:
                try:
                    model.predict(seq=seq, models=[m],
                                  num_recycles=advanced_settings["num_recycles_validation"],
                                  verbose=False)
                    L = copy_dict(model.aux["log"])
                    w.writerow([TARGET, j["name"], j["length"], si, int(guess), m+1,
                                round(L["plddt"],3), round(L["ptm"],3), round(L["i_ptm"],3),
                                round(L["pae"],3), round(L["i_pae"],3)])
                except Exception as e:
                    w.writerow([TARGET, j["name"], j["length"], si, int(guess), m+1,
                                "ERR","ERR","ERR","ERR",str(e)[:60]])
            fh.flush()
        print(f"  guess={int(guess)} {j['name']} done", flush=True)
fh.close()
print("wrote", OUT, flush=True)
