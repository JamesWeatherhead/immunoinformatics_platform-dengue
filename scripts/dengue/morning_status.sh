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

cyan "[6] Master-pipeline log (last 30 lines):"
[[ -f "$LOGS/master_pipeline.log" ]] && tail -30 "$LOGS/master_pipeline.log" || yellow "  no master_pipeline.log yet"
echo

cyan "[7] AlphaFold dispatch log (last 15 lines):"
[[ -f "$LOGS/alphafold_dispatch.log" ]] && tail -15 "$LOGS/alphafold_dispatch.log" || yellow "  no alphafold_dispatch.log yet"
echo

cyan "[8] NetMHCpan build log (last 15 lines):"
[[ -f "$LOGS/netmhcpan_build.log" ]] && tail -15 "$LOGS/netmhcpan_build.log" || yellow "  no netmhcpan_build.log yet"
echo

cyan "[9] Periodic-commit log (last 10 lines):"
[[ -f "$LOGS/periodic_commit.log" ]] && tail -10 "$LOGS/periodic_commit.log" || yellow "  no periodic_commit.log yet"
echo

cyan "[10] Disk usage:"
df -h /data | tail -1
echo

cyan "[11] Manuscript draft status:"
if [[ -f /data/james/ica-dengue/outputs/manuscript/dengue_npj_vaccines_draft.md ]]; then
    green "  -> draft EXISTS at outputs/manuscript/dengue_npj_vaccines_draft.md"
    n_tbd=$(grep -c "TBD" /data/james/ica-dengue/outputs/manuscript/dengue_npj_vaccines_draft.md 2>/dev/null || echo 0)
    if (( n_tbd > 0 )); then
        yellow "  -> $n_tbd TBD placeholders remain (need either pipeline rerun or manual fill)"
    else
        green "  -> all placeholders filled, ready for paperclip enrichment + author review"
    fi
else
    yellow "  -> draft not yet generated"
fi
echo

cyan "[12] Figures:"
fig_dir=/data/james/ica-dengue/outputs/figures
if [[ -d "$fig_dir" ]]; then
    n_png=$(ls "$fig_dir"/*.png 2>/dev/null | wc -l)
    green "  -> $n_png figures rendered in $fig_dir"
    ls "$fig_dir"/*.png 2>/dev/null | xargs -I {} basename {} | sed 's/^/    /'
else
    yellow "  -> no figures yet"
fi
echo

cyan "[13] What to do based on state:"
if [[ -f "$REPO/.git/refs/tags/v1.0-dengue-results" ]] || (cd "$REPO" && git tag -l 2>/dev/null | grep -q "v1.0-dengue-results"); then
    green "  -> v1.0-dengue-results TAGGED. Pipeline + figures + manuscript draft all done."
    green "  -> Pull from origin to your local: git pull && git checkout v1.0-dengue-results"
    green "  -> Then run paperclip locally to enrich citations:"
    green "       cd ~/Desktop/dengue-fork && bash scripts/dengue/run_paperclip_enrichment.sh"
elif grep -q "master_pipeline COMPLETE" "$LOGS/master_pipeline.log" 2>/dev/null; then
    green "  -> Master pipeline complete; tag should appear in 'git tag -l'"
elif grep -q "phase 4: launching master_pipeline" "$LOGS/auto_dispatch.log" 2>/dev/null; then
    yellow "  -> Master pipeline running. Check master_pipeline.log for current phase."
elif grep -q "phase 3: nextflow complete" "$LOGS/auto_dispatch.log" 2>/dev/null; then
    yellow "  -> B-cell complete; master_pipeline should be starting any minute."
elif grep -q "phase 3: dispatching" "$LOGS/auto_dispatch.log" 2>/dev/null; then
    yellow "  -> B-cell still running. Master pipeline will auto-launch when done."
elif grep -q "phase 2 FAILED" "$LOGS/auto_dispatch.log" 2>/dev/null; then
    red "  -> SIF build failure. Tell Claude: 'morning - SIF failure, debug'"
elif tmux has-session -t sif_builds 2>/dev/null; then
    yellow "  -> SIF builds still running."
else
    yellow "  -> Indeterminate state. Tell Claude: 'morning - unclear, investigate'"
fi
echo

cyan "=================================================================="
