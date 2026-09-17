#!/usr/bin/env bash
# Regenerate any *_lys_sasd.csv whose header and rows disagree in width.
#
# Sets processed before the domain-column fix wrote an 11-column header over
# 8-column rows, so every downstream read saw shifted columns and reported zero
# inter-domain crosslinks -- indistinguishable from a real negative result.
# Waits for any in-flight sweep to finish first so the two do not collide.
set -u
PY=/home/sbali/miniconda3/envs/CALVADOS/bin/python
TOP=${TOP:-500}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
cd "$(dirname "$0")"

# Wait for any in-flight sweep. Match only the python worker, not this
# script: pgrep -f matches full command lines and this file's own text
# contains the pattern, so a broader match makes the script wait on itself.
while pgrep -f "[a]nalyze_lys_contacts.py --set" >/dev/null; do sleep 20; done

BAD=()
for f in ../data/*_lys_sasd.csv; do
  [ -e "$f" ] || continue
  h=$(head -1 "$f" | awk -F, '{print NF}')
  r=$(sed -n '2p' "$f" | awk -F, '{print NF}')
  [ "$h" = "$r" ] || BAD+=("$(basename "$f" _lys_sasd.csv)")
done

if [ ${#BAD[@]} -eq 0 ]; then echo "all *_lys_sasd.csv well-formed"; exit 0; fi
echo "regenerating ${#BAD[@]} set(s): ${BAD[*]}"
for s in "${BAD[@]}"; do
  echo ""; echo "######## $s ########"
  $PY -u analyze_lys_contacts.py --set "$s" --sasd --sasd-top "$TOP" 2>&1 \
    | grep -E "SASD|replicates completed|Traceback|ERROR"
done
echo ""; echo "=== revalidating ==="
for f in ../data/*_lys_sasd.csv; do
  h=$(head -1 "$f" | awk -F, '{print NF}'); r=$(sed -n '2p' "$f" | awk -F, '{print NF}')
  s=$(basename "$f" _lys_sasd.csv)
  [ "$h" = "$r" ] && echo "  OK      $s" || echo "  BROKEN  $s"
done
