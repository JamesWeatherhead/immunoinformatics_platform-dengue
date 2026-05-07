#!/usr/bin/env bash
# Periodically commit pipeline outputs from cluster back to GitHub fork.
#
# Runs in a tmux session on hsv-server. Every 30 min, stages outputs/,
# logs/, and any new docs, commits if anything changed, and pushes.
#
# Uses the gh auth set up in /home/pathinformatics/.config/gh/ (PAT
# was set up via gh auth login --with-token earlier).
#
# Usage on cluster:
#   tmux new -s periodic_commit -d "bash scripts/dengue/cluster_periodic_commit.sh"
set -uo pipefail

REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/periodic_commit.log
mkdir -p "$(dirname "$LOG")"

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

cd "$REPO"
log "periodic_commit starting"

iter=0
while true; do
    iter=$((iter + 1))
    log "iteration $iter"
    cd "$REPO"

    # Stage anything new under outputs/ or docs/dengue/ or logs/
    git add outputs/ docs/dengue/ scripts/dengue/ 2>>"$LOG" || true

    # Check if there's anything to commit
    if git diff --cached --quiet; then
        log "  no changes to commit"
    else
        n_files=$(git diff --cached --name-only | wc -l)
        msg="auto(cluster): pipeline outputs iteration $iter ($n_files files)"
        if git commit -m "$msg" 2>>"$LOG"; then
            log "  committed: $msg"
            if git push origin master 2>>"$LOG"; then
                log "  pushed to origin/master"
            else
                log "  PUSH FAILED; will retry next iteration"
            fi
        else
            log "  COMMIT FAILED"
        fi
    fi

    # Sleep 30 min
    sleep 1800
done
