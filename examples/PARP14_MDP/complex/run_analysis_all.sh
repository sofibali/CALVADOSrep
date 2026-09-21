#!/bin/bash
# Full analysis chain over every completed set. Used to recover after the
# per-set watcher failed to fire (see analysis/watch.log, which stayed empty).
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
PY="conda run -n calvados python"
echo "=== START $(date -Is) ==="
rm -f analysis/data/raw.npy analysis/domain_face_contacts.csv
echo "--- stage 1: per-replicate observables (all sets) ---"
$PY analyze_binding.py 2>&1 | grep -v LIBCIFPP
echo "--- stage 2: dissociation table + figures ---"
$PY report2_binding.py 2>&1 | grep -v LIBCIFPP
echo "--- stage 3: contact maps + difference maps ---"
$PY figmaps_binding.py 2>&1 | grep -v LIBCIFPP
echo "--- stage 4: episode lifetimes + extend triage ---"
$PY episode_stats.py 2>&1 | grep -v LIBCIFPP
echo "--- stage 5: domain/face separation (all sets) ---"
$PY figure_face_separation.py 2>&1 | grep -v LIBCIFPP
echo "=== DONE $(date -Is) ==="
