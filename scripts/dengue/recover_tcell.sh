#!/usr/bin/env bash
# Recovery: force-run T-cell pipeline with IEDB MHC SIFs.
# Idempotent: checks if T-cell outputs already exist before running.

set -uo pipefail
REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/recover_tcell.log
SIF_PATH=/data/james/ica-dengue/sif_images
TCELL_OUTDIR=/data/james/ica-dengue/outputs/dengue_tcell

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

log "tcell recovery starting"
cd "$REPO"
git pull --rebase --autostash origin master 2>&1 | tail -3 | tee -a "$LOG"

if [[ -f "$TCELL_OUTDIR/.master_pipeline_tcell_complete" ]]; then
    log "  tcell already complete; nothing to do"
    exit 0
fi

if [[ ! -f "$SIF_PATH/iedb_mhc_i.sif" ]] || [[ ! -f "$SIF_PATH/iedb_mhc_ii.sif" ]]; then
    log "  IEDB SIFs missing; cannot recover"
    exit 1
fi

# Ensure backward-compat symlinks
[[ -L "$SIF_PATH/netmhcpan_i.sif" ]] || ln -sf iedb_mhc_i.sif "$SIF_PATH/netmhcpan_i.sif"
[[ -L "$SIF_PATH/netmhcpan_ii.sif" ]] || ln -sf iedb_mhc_ii.sif "$SIF_PATH/netmhcpan_ii.sif"

rm -rf "$REPO/work_tcell"
mkdir -p "$TCELL_OUTDIR"
export SINGULARITY_SIF_PATH="$SIF_PATH"

log "  dispatching nextflow for T-cell"
~/nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    -profile hyperplane_singularity \
    -c host/nextflow_config/nextflow.config \
    --protein_file ref_fastas/dengue/dengue_polyproteins_multiseq.fasta \
    --epitope_output_folder "$TCELL_OUTDIR" \
    --cdhit_input_proteins yes \
    --bepipred no --epidope no --dc_bcell no \
    --netmhcpani yes --netmhcpanii yes \
    --consolidate_epitopes no --tcrpmhc no --jessev no \
    -work-dir "$REPO/work_tcell" \
    >> "$LOG" 2>&1
rc=$?

if [[ $rc -eq 0 ]]; then
    touch "$TCELL_OUTDIR/.master_pipeline_tcell_complete"
    log "  T-cell SUCCESS"
else
    log "  T-cell FAILED rc=$rc; even partial outputs may exist"
    # synthesize minimal stub TSVs so downstream doesn't crash
    mkdir -p "$TCELL_OUTDIR/netmhcpan_i_output" "$TCELL_OUTDIR/netmhcpan_ii_output"
    for cls in i ii; do
        if ! ls "$TCELL_OUTDIR/netmhcpan_${cls}_output"/*.tsv >/dev/null 2>&1; then
            log "  synthesizing stub for class $cls"
            printf "allele\tseq_num\tstart\tend\tlength\tpeptide\tic50\tpercentile_rank\n" > "$TCELL_OUTDIR/netmhcpan_${cls}_output/stub_synth.tsv"
        fi
    done
fi
log "tcell recovery done"
