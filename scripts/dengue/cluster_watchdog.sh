#!/usr/bin/env bash
# Cluster-side watchdog: detects known failure patterns in master_pipeline
# logs and applies fixes autonomously over multi-day operation.
#
# Designed to run forever in its own tmux. Polls every 15 min. Each known
# failure pattern has a remediation. Unknown failures are logged loudly
# but the watchdog continues so partial progress is preserved.
#
# Usage:
#   tmux new -s watchdog -d "bash scripts/dengue/cluster_watchdog.sh; sleep 604800"

set -uo pipefail

REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/watchdog.log
SIF_PATH=/data/james/ica-dengue/sif_images
OUTPUTS=/data/james/ica-dengue/outputs

mkdir -p "$(dirname "$LOG")"
log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

cd "$REPO"
log "watchdog starting (PID $$)"

# Heartbeat detection: if a tmux session is supposed to be alive but its
# log hasn't grown in N minutes, consider it stuck.
last_logsize_master=0
master_stuck_count=0

iter=0
while true; do
    iter=$((iter + 1))
    log "iter $iter: checking system state"

    # ----- Auto_dispatch state -----
    if tmux has-session -t auto_dispatch 2>/dev/null; then
        if grep -q "phase 4: launching master_pipeline" /data/james/ica-dengue/logs/auto_dispatch.log 2>/dev/null; then
            log "  auto_dispatch: chained to master_pipeline OK"
        elif grep -q "phase 3: nextflow complete" /data/james/ica-dengue/logs/auto_dispatch.log 2>/dev/null; then
            # B-cell finished but master_pipeline not launched yet (newer auto_dispatch needed)
            if ! tmux has-session -t master_pipeline 2>/dev/null; then
                log "  REMEDIATION: B-cell done but master_pipeline not running; launching directly"
                cd "$REPO"
                tmux new -s master_pipeline -d "bash $REPO/scripts/dengue/cluster_master_pipeline.sh; sleep 604800"
            fi
        fi
    fi

    # ----- Master_pipeline heartbeat -----
    if tmux has-session -t master_pipeline 2>/dev/null; then
        cur_logsize=$(stat -c %s /data/james/ica-dengue/logs/master_pipeline.log 2>/dev/null || echo 0)
        if (( cur_logsize == last_logsize_master )); then
            master_stuck_count=$((master_stuck_count + 1))
            if (( master_stuck_count > 8 )); then  # 8 * 15 min = 2 hours
                log "  WARN: master_pipeline log unchanged for 2h; possible stuck phase"
                log "  REMEDIATION: capturing pane state"
                tmux capture-pane -t master_pipeline -p | tail -30 >> "$LOG"
            fi
        else
            master_stuck_count=0
        fi
        last_logsize_master=$cur_logsize
    fi

    # ----- Known failure pattern: missing SIF -----
    # Detect lines like: "FATAL:   could not open image /path/to/foo.sif: ..."
    # grep --no-filename avoids the bug where multiple files prefix matches with "<path>:"
    missing_sif=$(grep --no-filename -hoE "could not open image /data/james/ica-dengue/sif_images/[^[:space:]:]+\.sif" \
                  /data/james/ica-dengue/logs/master_pipeline.log \
                  /data/james/ica-dengue/logs/auto_dispatch.log \
                  /data/james/ica-dengue/logs/alphafold_dispatch.log 2>/dev/null \
                  | sed -E 's/^could not open image //' | sort -u | head -3)
    if [[ -n "$missing_sif" ]]; then
        log "  detected missing SIFs: $missing_sif"
        for sif in $missing_sif; do
            # Strict: only treat exact /data/james/ica-dengue/sif_images/X.sif paths
            if [[ ! "$sif" =~ ^/data/james/ica-dengue/sif_images/[A-Za-z0-9._-]+\.sif$ ]]; then
                log "    skipping malformed match: $sif"
                continue
            fi
            bn=$(basename "$sif")
            if [[ "$bn" == "blast_latest.sif" ]] || [[ "$bn" == "clustalomega_latest.sif" ]]; then
                log "    $bn: workflow gating fixed in 0a1f93e; not building"
            elif [[ "$bn" == "alphafold.sif" ]] && [[ ! -f "$sif" ]]; then
                log "    $bn: critical, attempting docker-daemon rebuild"
                if ! tmux has-session -t af_sif_rebuild 2>/dev/null; then
                    tmux new -s af_sif_rebuild -d "apptainer build $sif docker-daemon://alphafold:2.3.2 2>&1 | tee /data/james/ica-dengue/logs/alphafold_sif_rebuild.log; sleep 7200"
                fi
            fi
        done
    fi

    # ----- Known failure pattern: nextflow workdir permission -----
    if grep -qE "Permission denied" /data/james/ica-dengue/logs/*.log 2>/dev/null; then
        log "  WARN: detected permission errors; check ownership of workdirs"
    fi

    # ----- Known failure pattern: disk full -----
    df_pct=$(df -P /data | awk 'NR==2 {print int($5)}' | tr -d '%')
    if (( df_pct > 90 )); then
        log "  CRITICAL: /data disk is ${df_pct}% full; pruning old work directories"
        # Remove oldest nextflow work dirs older than 12h
        find "$REPO" -maxdepth 1 -type d -name "work_*" -mtime +0.5 -exec du -sh {} \; 2>/dev/null | head -5 | tee -a "$LOG"
        # NOT actually deleting until manual review; just reporting
    fi

    # ----- Pipeline progress markers -----
    log "  state markers:"
    for marker in \
        "$OUTPUTS/dengue_bcell/.master_pipeline_bcell_complete" \
        "$OUTPUTS/alphafold_e_dimers/.master_pipeline_af_complete" \
        "$OUTPUTS/dengue_tcell/.master_pipeline_tcell_complete"; do
        if [[ -f "$marker" ]]; then
            log "    [DONE] $(basename "$(dirname "$marker")")"
        fi
    done

    # ----- Tag detection -----
    cd "$REPO"
    if git tag -l 2>/dev/null | grep -q "^v1.0-dengue-results$"; then
        log "  v1.0-dengue-results TAGGED. Pipeline considered complete."
        log "  watchdog will continue polling for safety but no remediations needed."
    fi

    # ----- AlphaFold dispatch progress -----
    af_pdbs=$(find "$OUTPUTS/alphafold_e_dimers" -name "ranked_0.pdb" 2>/dev/null | wc -l)
    log "  alphafold ranked_0.pdb count: $af_pdbs / 4 expected"

    # ----- Sleep 15 min -----
    sleep 900
done
