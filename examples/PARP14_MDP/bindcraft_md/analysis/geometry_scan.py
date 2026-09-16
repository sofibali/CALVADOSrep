import glob, os, json, re
from Bio.PDB import PDBParser, Selection
from scipy.spatial import cKDTree
import numpy as np

parser = PDBParser(QUIET=True)

TARGETS = {
    "md1_block_af3": {
        "base": "/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md/designs/md1_block_af3",
        "hotspot_groups": {"MD1": [33,37,39,41,42,43,44,46,47,134,135,136,138,172,173]},
    },
    "clamp_md1md2_state3": {
        "base": "/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md/designs/clamp_md1md2_state3",
        "hotspot_groups": {
            "MD1": [71,102,103,106,107,109,110],
            "MD2": [351,354,357,395,396,399,400],
        },
    },
}

OUT = "/tmp/claude-64170/-home-sbali-BindCraft/adecb958-709c-481a-8e03-db4ccc10a3b9/scratchpad/geometry_results.json"

def min_dist(binder_coords, group_coords):
    if len(group_coords) == 0 or len(binder_coords) == 0:
        return None
    tree = cKDTree(group_coords)
    d, _ = tree.query(binder_coords, k=1)
    return float(d.min())

results = {}
for tname, cfg in TARGETS.items():
    base = cfg["base"]
    rows = []
    for bucket in ["Relaxed", "Clashing", "LowConfidence"]:
        folder = os.path.join(base, "Trajectory", bucket)
        pdbs = sorted(glob.glob(os.path.join(folder, "*.pdb")))
        for pdb_path in pdbs:
            name = os.path.basename(pdb_path).replace(".pdb", "")
            m = re.search(r"_l(\d+)_s", name)
            length = int(m.group(1)) if m else None
            try:
                structure = parser.get_structure(name, pdb_path)
                model = structure[0]
                if 'A' not in model or 'B' not in model:
                    continue
                binder_atoms = Selection.unfold_entities(model['B'], 'A')
                binder_coords = np.array([a.coord for a in binder_atoms])

                target_res_by_num = {r.id[1]: r for r in model['A'] if r.id[0] == ' '}

                dists = {}
                for gname, resnums in cfg["hotspot_groups"].items():
                    group_coords = []
                    for rn in resnums:
                        res = target_res_by_num.get(rn)
                        if res is None:
                            continue
                        for atom in res:
                            group_coords.append(atom.coord)
                    group_coords = np.array(group_coords) if group_coords else np.empty((0,3))
                    d = min_dist(binder_coords, group_coords)
                    dists[gname] = d

                rows.append({
                    "design": name,
                    "bucket": bucket,
                    "length": length,
                    "dists": dists,
                })
            except Exception as e:
                rows.append({"design": name, "bucket": bucket, "length": length, "error": str(e)})
    results[tname] = rows
    print(f"{tname}: processed {len(rows)} trajectories")

with open(OUT, "w") as f:
    json.dump(results, f)
print("wrote", OUT)
