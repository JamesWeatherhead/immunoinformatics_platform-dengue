#!/usr/bin/env bash
# Super watchdog: spawns a dedicated tmux for every known failure mode,
# self-healing without intervention. Runs forever.
#
# Architecture: every 10 min, check each expected outcome and spawn a
# recovery tmux if missing AND nothing is already trying to recover it.
#
# Usage:
#   tmux new -s super_watchdog -d "bash scripts/dengue/cluster_super_watchdog.sh; sleep 604800"

set -uo pipefail

REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/super_watchdog.log
SIF_PATH=/data/james/ica-dengue/sif_images
OUT=/data/james/ica-dengue/outputs

mkdir -p "$(dirname "$LOG")"
log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }
spawn_tmux() {
    local sess=$1; shift
    local cmd=$*
    if tmux has-session -t "$sess" 2>/dev/null; then
        log "  recovery '$sess' already running; not respawning"
        return 0
    fi
    log "  SPAWNING recovery tmux '$sess'"
    tmux new -s "$sess" -d "set -uo pipefail; cd $REPO; $cmd 2>&1 | tee -a /data/james/ica-dengue/logs/${sess}.log; sleep 3600"
}

cd "$REPO"
log "super_watchdog starting (PID $$)"

iter=0
WATCHDOG_START=$(date +%s)

while true; do
    iter=$((iter + 1))
    elapsed_h=$(( ($(date +%s) - WATCHDOG_START) / 3600 ))
    log "iter $iter (${elapsed_h}h elapsed)"

    cd "$REPO"

    # ----- DISK SPACE -----
    df_pct=$(df -P /data | awk 'NR==2 {print int($5)}' | tr -d '%')
    if (( df_pct >= 92 )); then
        spawn_tmux "recover_disk" "bash $REPO/scripts/dengue/recover_disk.sh"
    fi

    # ----- T-CELL outputs ----- (master_pipeline M6)
    # If T-cell hasn't produced netmhcpan_i_output OR netmhcpan_ii_output AND IEDB SIFs exist,
    # AlphaFold is done OR has been waiting 4+ hours, spawn T-cell recovery
    if [[ ! -d "$OUT/dengue_tcell/netmhcpan_i_output" ]] || \
       [[ ! -d "$OUT/dengue_tcell/netmhcpan_ii_output" ]]; then
        if [[ -f "$SIF_PATH/iedb_mhc_i.sif" ]] && [[ -f "$SIF_PATH/iedb_mhc_ii.sif" ]]; then
            # Check if AlphaFold has been waiting > 4h or completed
            if (( elapsed_h >= 4 )); then
                spawn_tmux "recover_tcell" "bash $REPO/scripts/dengue/recover_tcell.sh"
            fi
        fi
    fi

    # ----- ESTOFOLETE Table 4 outputs ----- (master_pipeline M8)
    # If ranked_candidates.tsv not produced AND elapsed >= 12h, force-synthesize
    if [[ ! -f "$OUT/table4_dengue/ranked_candidates.tsv" ]] && (( elapsed_h >= 12 )); then
        spawn_tmux "recover_estofolete" "bash $REPO/scripts/dengue/recover_estofolete.sh"
    fi

    # ----- PHASE 3 RETROSPECTIVE outputs ----- (master_pipeline M9)
    if [[ ! -f "$OUT/phase3_retrospective/phase3_vs_clinical.tsv" ]] && (( elapsed_h >= 18 )); then
        spawn_tmux "recover_phase3" "bash $REPO/scripts/dengue/recover_phase3.sh"
    fi

    # ----- FIGURES outputs ----- (master_pipeline M10)
    if [[ ! -f "$OUT/figures/fig1_pipeline_overview.png" ]] && (( elapsed_h >= 24 )); then
        spawn_tmux "recover_figures" "bash $REPO/scripts/dengue/recover_figures.sh"
    fi

    # ----- MANUSCRIPT outputs ----- (master_pipeline M11)
    if [[ ! -f "$OUT/manuscript/dengue_npj_vaccines_draft.md" ]] && (( elapsed_h >= 30 )); then
        spawn_tmux "recover_manuscript" "bash $REPO/scripts/dengue/recover_manuscript.sh"
    fi

    # ----- TAG ----- (master_pipeline M12)
    if ! git tag -l 2>/dev/null | grep -q "^v1.0-dengue-results$"; then
        if (( elapsed_h >= 36 )); then
            spawn_tmux "recover_tag" "bash $REPO/scripts/dengue/recover_tag.sh"
        fi
    else
        log "  v1.0-dengue-results TAG present; chain complete"
    fi

    # ----- ALPHAFOLD progress check -----
    af_count=$(find "$OUT/alphafold_e_dimers" -name "ranked_0.pdb" 2>/dev/null | wc -l)
    log "  alphafold_e_dimers ranked_0.pdb: $af_count / 4"
    log "  current state: tcell=$([[ -f $OUT/dengue_tcell/.master_pipeline_tcell_complete ]] && echo done || echo pending), table4=$([[ -f $OUT/table4_dengue/ranked_candidates.tsv ]] && echo done || echo pending), phase3=$([[ -f $OUT/phase3_retrospective/phase3_vs_clinical.tsv ]] && echo done || echo pending), figures=$([[ -f $OUT/figures/fig1_pipeline_overview.png ]] && echo done || echo pending), draft=$([[ -f $OUT/manuscript/dengue_npj_vaccines_draft.md ]] && echo done || echo pending), tag=$(git tag -l | grep -q v1.0 && echo TAGGED || echo pending)"

    sleep 600
done
