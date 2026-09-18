#!/bin/bash
# Analyse each binding set as soon as run_campaign.sh finishes it.
#
# run_campaign.sh touches logs/{set}.complete when a set's 20 replicates are done.
# This watcher picks those up, runs stage 1 (per-replicate observables, merged
# into analysis/data/raw.npy), then regenerates the tables and figures over every
# set analysed so far. Exits once the campaign is complete and nothing is pending.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
PY="conda run -n calvados python"
mkdir -p logs analysis

while true; do
    did_work=0
    for marker in logs/*.complete; do
        [ -e "$marker" ] || continue
        s=$(basename "$marker" .complete)
        [ -e "logs/${s}.analyzed" ] && continue

        echo "[$(date -Is)] analysing $s"
        if $PY analyze_binding.py "$s" 2>&1 | grep -v LIBCIFPP; then
            $PY report2_binding.py  2>&1 | grep -v LIBCIFPP
            $PY figmaps_binding.py  2>&1 | grep -v LIBCIFPP
            $PY episode_stats.py    2>&1 | grep -v LIBCIFPP
            $PY figure_face_separation.py "$s" 2>&1 | grep -v LIBCIFPP
            touch "logs/${s}.analyzed"
            echo "[$(date -Is)] $s analysed"
        else
            echo "[$(date -Is)] $s analysis FAILED - will retry next pass"
        fi
        did_work=1
    done

    if grep -q 'CAMPAIGN COMPLETE' campaign5.log 2>/dev/null; then
        pending=0
        for marker in logs/*.complete; do
            [ -e "$marker" ] || continue
            s=$(basename "$marker" .complete)
            [ -e "logs/${s}.analyzed" ] || pending=1
        done
        [ "$pending" -eq 0 ] && { echo "[$(date -Is)] all sets analysed, watcher exiting"; break; }
    fi

    [ "$did_work" -eq 0 ] && sleep 120
done
