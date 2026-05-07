#!/usr/bin/env bash
# Overnight auto-dispatch for the dengue extension pipeline.
#
# Runs in a tmux session on hsv-server. Monitors the `sif_builds` tmux
# session. When SIF builds complete, fires the B-cell pipeline run
# (CDHIT + EpiDope only — BepiPred is skipped tonight because the IEDB
# B-cell container needs IEDB-licensed tarballs that are not on the
# cluster) on the 4 DENV reference polyproteins.
#
# Logs to logs/auto_dispatch.log so morning_status.sh shows the timeline.
#
# Usage on cluster:
#   tmux new -s auto_dispatch -d "bash scripts/dengue/cluster_auto_dispatch.sh; sleep 172800"
set -uo pipefail

REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/auto_dispatch.log
SIF_PATH=/data/james/ica-dengue/sif_images
mkdir -p "$(dirname "$LOG")"

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

cd "$REPO"
log "auto_dispatch starting"
log "repo HEAD: $(git rev-parse --short HEAD)"
log "SIF cache before: $(ls "$SIF_PATH" 2>/dev/null | wc -l) files"

# Phase 1: wait for SIF builds tmux to finish
log "phase 1: waiting for sif_builds tmux to finish"
poll_count=0
while tmux has-session -t sif_builds 2>/dev/null; do
    poll_count=$((poll_count + 1))
    pane_lines=$(tmux capture-pane -t sif_builds -p 2>/dev/null | tail -3 | tr '\n' ' ' | head -c 200)
    log "  poll $poll_count: $pane_lines"
    if echo "$pane_lines" | grep -q '=== all builds done ==='; then
        log "  detected all-builds-done banner"
        break
    fi
    sleep 60
done
log "phase 1 complete; SIF cache after: $(ls "$SIF_PATH" 2>/dev/null | wc -l) files"
ls -la "$SIF_PATH" >> "$LOG" 2>&1

# Phase 2: figure out what we can run
# CDHIT only is always possible. EpiDope is the only B-cell scorer
# currently buildable (BepiPred container needs IEDB tarballs that don't
# exist on cluster).
need=("cd-hit-4-8-1.sif" "cdhit-to-tsv.sif")
missing=0
for sif in "${need[@]}"; do
    if [[ ! -f "$SIF_PATH/$sif" ]]; then
        log "  MISSING required SIF: $sif"
        missing=1
    fi
done
if (( missing )); then
    log "phase 2 FAILED: minimum SIFs not built; aborting auto-dispatch"
    log "morning task: investigate /tmp/build_*.log for which container failed"
    log "auto_dispatch exiting"
    exit 2
fi

epidope_arg="--epidope no"
if [[ -f "$SIF_PATH/epidope_v0.3.sif" ]] || [[ -f "$SIF_PATH/cuda-epidope.sif" ]]; then
    epidope_arg="--epidope yes"
    log "  epidope SIF present, enabling"
else
    log "  epidope SIF NOT present, will run CDHIT-only"
fi
log "phase 2 complete: required SIFs present, $epidope_arg"

# Phase 3: dispatch B-cell pipeline run on dengue references
log "phase 3: dispatching B-cell pipeline (CDHIT + EpiDope; BepiPred skipped due to IEDB container build blocker)"

bcell_outdir=/data/james/ica-dengue/outputs/dengue_bcell
rm -rf "$bcell_outdir" "$REPO/work_bcell"
mkdir -p "$bcell_outdir"

export SINGULARITY_SIF_PATH="$SIF_PATH"
log "dispatch command: nextflow run ... --bepipred no $epidope_arg --cdhit_input_proteins yes (T-cell stages off)"

cd "$REPO"
~/nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    -profile hyperplane_singularity \
    -c host/nextflow_config/nextflow.config \
    --protein_file ref_fastas/dengue/dengue_polyproteins_multiseq.fasta \
    --epitope_output_folder "$bcell_outdir" \
    --cdhit_input_proteins yes \
    --bepipred no \
    $epidope_arg \
    --dc_bcell no \
    --netmhcpani no --netmhcpanii no \
    --consolidate_epitopes no --tcrpmhc no --jessev no \
    -work-dir "$REPO/work_bcell" \
    >> "$LOG" 2>&1 &
NF_PID=$!
log "nextflow PID: $NF_PID; will tail status every 60s"

# Wait up to 4 hours and report status periodically
for i in $(seq 1 240); do
    sleep 60
    if kill -0 "$NF_PID" 2>/dev/null; then
        last3=$(tail -3 "$LOG" 2>/dev/null | tr '\n' ' ' | head -c 200)
        if (( i % 5 == 0 )); then
            log "  ($i min) nextflow still running"
        fi
    else
        log "  nextflow exited (after $i min)"
        break
    fi
done

if kill -0 "$NF_PID" 2>/dev/null; then
    log "phase 3: nextflow still running after 4 hours; auto_dispatch script done watching"
else
    log "phase 3: nextflow complete"
    log "outputs:"
    find "$bcell_outdir" -type f >> "$LOG" 2>&1
fi

# Phase 4: chain into master_pipeline for everything downstream
# (AlphaFold, DiscoTope, T-cell, JessEV, Table 4 mapping, Phase 3 retro,
#  figures, manuscript, commit, tag). master_pipeline.sh is idempotent and
#  resumes from wherever it last stopped.
log "phase 4: launching master_pipeline.sh in tmux session 'master_pipeline'"
if [[ -x "$REPO/scripts/dengue/cluster_master_pipeline.sh" ]]; then
    if tmux has-session -t master_pipeline 2>/dev/null; then
        log "  master_pipeline tmux session already exists; not relaunching"
    else
        tmux new -s master_pipeline -d "bash $REPO/scripts/dengue/cluster_master_pipeline.sh; sleep 604800"
        log "  launched master_pipeline tmux"
    fi
else
    log "  cluster_master_pipeline.sh missing or not executable; downstream phases will not run"
fi

log "auto_dispatch script done; master_pipeline now driving multi-day chain"
