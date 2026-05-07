#!/usr/bin/env bash
# Single command for James to run when he wakes up: shows complete cluster state.
#
# Local usage (from your Mac):
#   source ~/.claude/secrets.env && \
#   sshpass -p "$SSH_CLUSTER_PASS" ssh hsv-server bash /data/james/ica-dengue/immunoinformatics_platform-dengue/scripts/dengue/morning_status.sh
#
# Or copy this file local and run with sshpass on it.

REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
SIF_PATH=/data/james/ica-dengue/sif_images
LOGS=/data/james/ica-dengue/logs

cyan() { printf '\033[36m%s\033[0m\n' "$*"; }
green() { printf '\033[32m%s\033[0m\n' "$*"; }
red() { printf '\033[31m%s\033[0m\n' "$*"; }
yellow() { printf '\033[33m%s\033[0m\n' "$*"; }

cyan "=================================================================="
cyan " ICA-Dengue cluster morning status"
cyan " $(date)"
cyan "=================================================================="
echo

cyan "[1] tmux sessions running:"
if tmux ls 2>/dev/null; then :; else red "  no tmux sessions"; fi
echo

cyan "[2] SIF images built ($(ls "$SIF_PATH" 2>/dev/null | wc -l) files):"
ls -lah "$SIF_PATH" 2>/dev/null | tail -n +2 | awk '{printf "  %-12s %s\n", $5, $NF}' || red "  no SIF cache"
echo

cyan "[3] Latest 10 commits in fork:"
cd "$REPO" 2>/dev/null && git log --oneline -10 || red "  cant cd to repo"
echo

cyan "[4] Pipeline outputs (any new TSVs?):"
if [[ -d /data/james/ica-dengue/outputs ]]; then
    find /data/james/ica-dengue/outputs -type f -name "*.tsv" 2>/dev/null | head -20 | while read f; do
        printf "  %s (%s)\n" "$f" "$(stat -c %s "$f" 2>/dev/null || echo ?)b"
    done
else
    yellow "  outputs dir does not exist yet"
fi
echo

cyan "[5] Auto-dispatch log (last 30 lines):"
[[ -f "$LOGS/auto_dispatch.log" ]] && tail -30 "$LOGS/auto_dispatch.log" || yellow "  no auto_dispatch.log yet"
echo

cyan "[6] Periodic-commit log (last 10 lines):"
[[ -f "$LOGS/periodic_commit.log" ]] && tail -10 "$LOGS/periodic_commit.log" || yellow "  no periodic_commit.log yet"
echo

cyan "[7] Disk usage:"
df -h /data | tail -1
echo

cyan "[8] What to do based on state:"
if grep -q "phase 3: nextflow complete" "$LOGS/auto_dispatch.log" 2>/dev/null; then
    green "  -> B-cell pipeline finished. Tell Claude: 'morning - bcell done, dispatch full pipeline'"
elif grep -q "phase 3: dispatching" "$LOGS/auto_dispatch.log" 2>/dev/null; then
    yellow "  -> B-cell pipeline still running. Check back in a few hours."
elif grep -q "phase 2 FAILED" "$LOGS/auto_dispatch.log" 2>/dev/null; then
    red "  -> SIF build failure overnight. Tell Claude: 'morning - SIF failure, debug'"
elif tmux has-session -t sif_builds 2>/dev/null; then
    yellow "  -> SIF builds still running. Check sif_builds tmux."
else
    yellow "  -> Indeterminate state. Tell Claude: 'morning - unclear, investigate'"
fi
echo

cyan "=================================================================="
