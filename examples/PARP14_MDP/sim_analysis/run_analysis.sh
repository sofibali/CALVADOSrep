#!/bin/bash
# Run the full PARP14 trajectory-analysis suite on a simulation set OR any folder.
# Everything routes through sim_registry.py, so it works on ANY construct without
# editing code -- pass a named --set or a --sim-folder PATH (+ --units if the folder
# has no metadata.json).
#
#   bash run_analysis.sh --set fl_optimized
#   bash run_analysis.sh --sim-folder md_full --units md1l1 md2 md3
#   WORKERS=24 bash run_analysis.sh --sim-folder some/new/construct
#
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
PY="${PY:-python}"
WORKERS="${WORKERS:-16}"
if [ "$#" -eq 0 ]; then
  echo "usage: run_analysis.sh (--set NAME | --sim-folder PATH [--units u1 u2 ...])"
  exit 1
fi

echo "==== analyze_all: conf-prop, dmap, cmap, fnc, energy, wcn, active-sites, accessibility ===="
$PY analyze_all.py "$@" --workers "$WORKERS"

echo "==== lysine contacts ===="
$PY analyze_lys_contacts.py "$@"

echo "==== DONE. Per-residue/array results in data/ ===="
echo "Figures (read data/, run as needed):"
echo "  python figure_accessibility.py | figure_sasa_faces.py | figure_md_distances.py"
echo "  figure_domain_comparisons.py | figure_lys_exposed_persistence.py | figure_kh_domains.py"
