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
        [ -e "logs/${s}.failed" ] && continue

        echo "[$(date -Is)] analysing $s"
        # Run each stage and test the STAGE's exit code, not a pipeline's.
        # `cmd | grep -v LIBCIFPP` returns grep's status, so a crashed stage that
        # printed a traceback looked like success (grep matched the traceback
        # lines -> rc 0) and the set was marked analysed with no data; a silent
        # stage looked like failure and, with did_work=1 skipping the sleep,
        # span in a tight retry loop. Both actually happened.
        ok=1
        for stage in "analyze_binding.py $s" "report2_binding.py" \
                     "figmaps_binding.py" "episode_stats.py" \
                     "figure_face_separation.py $s"; do
            log=$(mktemp)
            # shellcheck disable=SC2086
            $PY $stage > "$log" 2>&1; rc=$?
            grep -v LIBCIFPP "$log"; rm -f "$log"
            if [ "$rc" -ne 0 ]; then
                echo "[$(date -Is)] $s: stage '${stage%% *}' exited $rc - not marking analysed"
                ok=0; break
            fi
        done
        if [ "$ok" -eq 1 ]; then
            touch "logs/${s}.analyzed"
            echo "[$(date -Is)] $s analysed"
        else
            touch "logs/${s}.failed"
            echo "[$(date -Is)] $s analysis FAILED - see logs/${s}.failed"
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
