#!/usr/bin/env bash
# Overnight auto-dispatch for the dengue extension pipeline.
#
# Runs in a tmux session on hsv-server. Monitors the `sif_builds` tmux
# session. When SIF builds complete, fires the B-cell-only pipeline run
# (CDHIT + BepiPred + EpiDope) on the 4 DENV reference polyproteins.
#
# Skips dc_bcell tonight because DiscoTope requires pre-computed AlphaFold
# PDBs that we run separately (AlphaFold dispatch is a morning task; takes
# its own GPU schedule).
#
# Logs everything to logs/auto_dispatch.log so morning_status.sh can show
# the timeline at a glance.
#
# Usage on cluster:
#   tmux new -s auto_dispatch -d "bash scripts/dengue/cluster_auto_dispatch.sh"
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
while tmux has-session -t sif_builds 2>/dev/null; do
    pane_lines=$(tmux capture-pane -t sif_builds -p 2>/dev/null | tail -3)
    log "  sif_builds tail: $(echo "$pane_lines" | tr '\n' ' ' | head -c 200)"
    if echo "$pane_lines" | grep -q '=== all builds done ==='; then
        log "  detected all builds done banner"
        break
    fi
    sleep 60
done
log "phase 1 complete; SIF cache after: $(ls "$SIF_PATH" 2>/dev/null | wc -l) files"
ls -la "$SIF_PATH" >> "$LOG" 2>&1

# Phase 2: verify minimum SIFs present for B-cell pipeline
need=("cd-hit-4-8-1.sif" "cdhit-to-tsv.sif" "iedb-b-cell-31.sif" "bepipred-to-tsv.sif")
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
    exit 2
fi
log "phase 2 complete: minimum required SIFs present"

# Phase 3: dispatch B-cell pipeline run on dengue references
log "phase 3: dispatching B-cell pipeline (CDHIT + BepiPred [+ EpiDope if SIF built])"

# Detect what we can actually run based on available SIFs
epidope_arg="--epidope no"
if [[ -f "$SIF_PATH/epidope_v0.3.sif" ]] || [[ -f "$SIF_PATH/cuda-epidope.sif" ]]; then
    epidope_arg="--epidope yes"
    log "  epidope SIF present, enabling"
fi

bcell_outdir=/data/james/ica-dengue/outputs/dengue_bcell
rm -rf "$bcell_outdir" "$REPO/work_bcell"
mkdir -p "$bcell_outdir"

export SINGULARITY_SIF_PATH="$SIF_PATH"
log "dispatch command:"
log "  nextflow run ... --bepipred yes $epidope_arg --cdhit_input_proteins yes (T-cell stages off)"

cd "$REPO"
~/nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    -profile hyperplane_singularity \
    -c host/nextflow_config/nextflow.config \
    --protein_file ref_fastas/dengue/dengue_polyproteins_multiseq.fasta \
    --epitope_output_folder "$bcell_outdir" \
    --cdhit_input_proteins yes \
    --bepipred yes \
    $epidope_arg \
    --dc_bcell no \
    --netmhcpani no --netmhcpanii no \
    --consolidate_epitopes no --tcrpmhc no --jessev no \
    -work-dir "$REPO/work_bcell" \
    >> "$LOG" 2>&1 &
NF_PID=$!
log "nextflow PID: $NF_PID; tailing for 30 min then leaving in background"

# Wait up to 30 min and report status periodically
for i in $(seq 1 30); do
    sleep 60
    if kill -0 "$NF_PID" 2>/dev/null; then
        log "  ($i/30) nextflow still running"
    else
        log "  nextflow exited (after $i min). exit code captured later."
        break
    fi
done

if kill -0 "$NF_PID" 2>/dev/null; then
    log "phase 3: nextflow continuing in background; exit auto_dispatch"
else
    log "phase 3: nextflow complete"
    log "outputs:"
    find "$bcell_outdir" -type f >> "$LOG" 2>&1
fi

log "auto_dispatch script done"
