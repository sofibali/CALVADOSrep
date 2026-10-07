#!/bin/bash
# all-atom PDB -> HyRes PDB, the full corrected pipeline.
#   $1 input all-atom pdb   $2 output tag   $3 4-char segid
set -eu
P=/home/sbali/miniconda3/envs/hyres/bin/python
D=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
A=$D/../../HyRes_GPU/at2hyres
S=/tmp/claude-64170/-home-sbali-CALVADOS-examples-PARP14-MDP-complex/d20d78d0-35be-4d3b-9d0a-a0d0ae134789/scratchpad/charmmify.py
in=$1; tag=$2; seg=$3
$P $D/addH.py "$in" "$D/${tag}_H.pdb" >/dev/null
# CHARMM names the backbone amide hydrogen HN; OpenMM/PDB writes H
awk '{ if (substr($0,1,4)=="ATOM" && substr($0,13,4)==" H  ") print substr($0,1,12) " HN " substr($0,17); else print }' \
    "$D/${tag}_H.pdb" > "$D/${tag}_Hn.pdb"
$P $S "$D/${tag}_Hn.pdb" "$D/${tag}_cg.pdb" "$seg" >/dev/null
$P $A/pdbfix_res.py "$D/${tag}_cg.pdb" "$D/${tag}_fix.pdb" 0 >/dev/null 2>&1
$P $A/at2hyres.py "$D/${tag}_fix.pdb" "$D/${tag}_hyres.pdb" >/dev/null 2>&1
echo "  ${tag}_hyres.pdb: $(grep -c '^ATOM' "$D/${tag}_hyres.pdb") particles"
