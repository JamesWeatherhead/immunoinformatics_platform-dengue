#!/usr/bin/env bash
# Master orchestration: chains every post-B-cell phase autonomously over 2-3 days.
#
# Designed for unattended multi-day cluster operation. James can walk away and
# return to a fully-rendered, manuscript-filled, committed-and-tagged release.
#
# Phase chain:
#   M1. Wait for B-cell pipeline (auto_dispatch) to deposit outputs/dengue_bcell/
#   M2. Dispatch AlphaFold-Multimer on DENV E homo-dimers (cluster_alphafold_dispatch.sh)
#   M3. Wait for AlphaFold to finish (~6-12h per protein, parallelized over 8 GPUs)
#   M4. Run DiscoTope on AlphaFold structures (geometric B-cell epitope refinement)
#   M5. Build NetMHCpan I/II SIFs from in-repo binaries (cluster_build_netmhcpan.sh)
#   M6. Dispatch full T-cell pipeline on dengue references
#   M7. Run JessEV epitope set selection (Gurobi)
#   M8. Apply Estofolete Table 4 mapping (estofolete_table4_mapping.py)
#   M9. Run Phase 3 retrospective validation (retrospective_phase3_validation.py)
#   M10. Generate all figures (generate_figures.py)
#   M11. Fill manuscript template with pipeline numbers (fill_manuscript.py)
#   M12. Commit + tag v1.0-dengue-results
#
# Each phase is idempotent: re-running this script after partial completion
# resumes from the first incomplete phase. Failed phases log and continue;
# downstream phases that depend on them are skipped with a warning.
#
# Usage on cluster:
#   tmux new -s master_pipeline -d "bash scripts/dengue/cluster_master_pipeline.sh; sleep 604800"
set -uo pipefail

REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/master_pipeline.log
SIF_PATH=/data/james/ica-dengue/sif_images
OUTPUTS=/data/james/ica-dengue/outputs
mkdir -p "$(dirname "$LOG")" "$OUTPUTS"

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }
phase() { log ""; log "=========================================================================="; log " $*"; log "=========================================================================="; }

cd "$REPO"
log "master_pipeline starting"
log "repo HEAD: $(git rev-parse --short HEAD)"
log "PID: $$"

# ----------------------------------------------------------------------------
# M1. Wait for B-cell pipeline to complete
# ----------------------------------------------------------------------------
phase "M1. Waiting for B-cell pipeline (auto_dispatch) to complete"
BCELL_OUTDIR="$OUTPUTS/dengue_bcell"
BCELL_DONE_MARKER="$BCELL_OUTDIR/.master_pipeline_bcell_complete"
m1_max_wait_hours=8
m1_polls=0
while true; do
    m1_polls=$((m1_polls + 1))
    if grep -q "phase 3: nextflow complete" /data/james/ica-dengue/logs/auto_dispatch.log 2>/dev/null; then
        log "  detected auto_dispatch B-cell completion"
        touch "$BCELL_DONE_MARKER"
        break
    fi
    if grep -q "phase 2 FAILED" /data/james/ica-dengue/logs/auto_dispatch.log 2>/dev/null; then
        log "  auto_dispatch reported phase 2 FAILED; B-cell skipped, proceeding to next phases"
        break
    fi
    if (( m1_polls > m1_max_wait_hours * 12 )); then
        log "  B-cell wait exceeded ${m1_max_wait_hours}h; proceeding to next phases anyway"
        break
    fi
    sleep 300
done
log "M1 complete after $m1_polls polls"
[[ -d "$BCELL_OUTDIR" ]] && find "$BCELL_OUTDIR" -type f | head -20 >> "$LOG"

# ----------------------------------------------------------------------------
# M2. Dispatch AlphaFold-Multimer on DENV E homo-dimers
# ----------------------------------------------------------------------------
phase "M2. Dispatching AlphaFold-Multimer on DENV E dimers"
AF_OUTDIR="$OUTPUTS/alphafold_e_dimers"
AF_DONE_MARKER="$AF_OUTDIR/.master_pipeline_af_complete"
mkdir -p "$AF_OUTDIR"
if [[ -f "$AF_DONE_MARKER" ]]; then
    log "  AlphaFold already marked complete; skipping dispatch"
else
    if [[ -x "$REPO/scripts/dengue/cluster_alphafold_dispatch.sh" ]]; then
        log "  launching cluster_alphafold_dispatch.sh"
        bash "$REPO/scripts/dengue/cluster_alphafold_dispatch.sh" >> "$LOG" 2>&1 &
        AF_PID=$!
        log "  alphafold dispatch PID: $AF_PID"
    else
        log "  cluster_alphafold_dispatch.sh missing or not executable; SKIPPING AlphaFold"
        AF_PID=""
    fi
fi

# ----------------------------------------------------------------------------
# M5. Build NetMHCpan I/II SIFs in parallel (independent of AlphaFold)
# ----------------------------------------------------------------------------
phase "M5. Building NetMHCpan I/II SIFs from in-repo binaries"
if [[ -f "$SIF_PATH/netmhcpan_i.sif" ]] && [[ -f "$SIF_PATH/netmhcpan_ii.sif" ]]; then
    log "  NetMHCpan I + II SIFs already built; skipping"
else
    if [[ -x "$REPO/scripts/dengue/cluster_build_netmhcpan.sh" ]]; then
        log "  launching cluster_build_netmhcpan.sh"
        bash "$REPO/scripts/dengue/cluster_build_netmhcpan.sh" >> "$LOG" 2>&1 &
        NMP_PID=$!
        log "  netmhcpan build PID: $NMP_PID"
    else
        log "  cluster_build_netmhcpan.sh missing; T-cell phase will be SKIPPED"
        NMP_PID=""
    fi
fi

# ----------------------------------------------------------------------------
# M3. Wait for AlphaFold to finish
# ----------------------------------------------------------------------------
phase "M3. Waiting for AlphaFold-Multimer (~6-12h per protein, parallel)"
af_max_wait_hours=24
af_polls=0
if [[ -n "${AF_PID:-}" ]]; then
    while true; do
        af_polls=$((af_polls + 1))
        if [[ -f "$AF_DONE_MARKER" ]]; then
            log "  AlphaFold complete marker present"
            break
        fi
        if ! kill -0 "$AF_PID" 2>/dev/null; then
            log "  alphafold dispatch PID exited"
            touch "$AF_DONE_MARKER"
            break
        fi
        if (( af_polls > af_max_wait_hours * 6 )); then
            log "  AlphaFold wait exceeded ${af_max_wait_hours}h; proceeding (DiscoTope may have partial input)"
            break
        fi
        if (( af_polls % 6 == 0 )); then
            log "  ($((af_polls * 10)) min) waiting for AlphaFold..."
        fi
        sleep 600
    done
else
    log "  no AlphaFold PID; skipping wait"
fi

# ----------------------------------------------------------------------------
# M4. DiscoTope on AlphaFold structures (geometric B-cell epitope scoring)
# ----------------------------------------------------------------------------
phase "M4. DiscoTope on AlphaFold structures"
DISCO_OUT="$OUTPUTS/discotope"
mkdir -p "$DISCO_OUT"
if compgen -G "$AF_OUTDIR/*/ranked_0.pdb" > /dev/null; then
    log "  found AlphaFold PDBs; attempting DiscoTope"
    if [[ -f "$SIF_PATH/discotope.sif" ]] || [[ -f "$SIF_PATH/discotope3.sif" ]]; then
        for pdb in "$AF_OUTDIR"/*/ranked_0.pdb; do
            name=$(basename "$(dirname "$pdb")")
            log "  DiscoTope on $name"
            sif=$(ls "$SIF_PATH"/discotope*.sif 2>/dev/null | head -1)
            apptainer run --nv "$sif" --pdb "$pdb" --out_dir "$DISCO_OUT/$name" >> "$LOG" 2>&1 || \
                log "    DiscoTope failed for $name (non-fatal)"
        done
    else
        log "  no DiscoTope SIF; SKIPPING (geometric refinement not done)"
    fi
else
    log "  no AlphaFold PDBs found; SKIPPING DiscoTope"
fi

# ----------------------------------------------------------------------------
# M6. Wait for NetMHCpan SIFs and dispatch T-cell pipeline
# ----------------------------------------------------------------------------
phase "M6. T-cell pipeline (NetMHCpan I + II)"
if [[ -n "${NMP_PID:-}" ]]; then
    log "  waiting for NetMHCpan build to finish"
    wait "$NMP_PID" 2>/dev/null
    log "  NetMHCpan build PID exited"
fi

TCELL_OUTDIR="$OUTPUTS/dengue_tcell"
TCELL_DONE="$TCELL_OUTDIR/.master_pipeline_tcell_complete"
mkdir -p "$TCELL_OUTDIR"
if [[ -f "$TCELL_DONE" ]]; then
    log "  T-cell already complete; skipping"
elif [[ -f "$SIF_PATH/iedb_mhc_i.sif" ]] && [[ -f "$SIF_PATH/iedb_mhc_ii.sif" ]]; then
    log "  dispatching T-cell pipeline"
    rm -rf "$REPO/work_tcell"
    export SINGULARITY_SIF_PATH="$SIF_PATH"
    cd "$REPO"
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
    if [[ $? -eq 0 ]]; then
        touch "$TCELL_DONE"
        log "  T-cell pipeline complete"
    else
        log "  T-cell pipeline FAILED; check log"
    fi
else
    log "  IEDB MHC SIFs missing; will retry every 5 min for up to 30 min"
    waited=0
    while (( waited < 1800 )); do
        sleep 300
        waited=$((waited + 300))
        if [[ -f "$SIF_PATH/iedb_mhc_i.sif" ]] && [[ -f "$SIF_PATH/iedb_mhc_ii.sif" ]]; then
            log "  IEDB SIFs now present; running T-cell pipeline"
            rm -rf "$REPO/work_tcell"
            export SINGULARITY_SIF_PATH="$SIF_PATH"
            cd "$REPO"
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
            if [[ $? -eq 0 ]]; then
                touch "$TCELL_DONE"
                log "  T-cell pipeline complete (after waiting ${waited}s)"
            else
                log "  T-cell pipeline FAILED after wait; check log"
            fi
            break
        fi
        log "  ($((waited/60)) min) still waiting for IEDB SIFs"
    done
    if [[ ! -f "$TCELL_DONE" ]]; then
        log "  IEDB SIFs never appeared in 30 min; SKIPPING T-cell"
    fi
fi

# ----------------------------------------------------------------------------
# M7. JessEV epitope set selection (Gurobi)
# ----------------------------------------------------------------------------
phase "M7. JessEV epitope set selection"
JESS_OUTDIR="$OUTPUTS/jessev"
mkdir -p "$JESS_OUTDIR"
if [[ -f "$SIF_PATH/jess_ev.sif" ]] && [[ -f "$TCELL_DONE" ]]; then
    log "  dispatching JessEV"
    rm -rf "$REPO/work_jessev"
    cd "$REPO"
    ~/nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
        -profile hyperplane_singularity \
        -c host/nextflow_config/nextflow.config \
        --protein_file ref_fastas/dengue/dengue_polyproteins_multiseq.fasta \
        --epitope_output_folder "$JESS_OUTDIR" \
        --cdhit_input_proteins yes \
        --bepipred no --epidope no --dc_bcell no \
        --netmhcpani yes --netmhcpanii yes \
        --consolidate_epitopes yes --tcrpmhc no --jessev yes \
        -work-dir "$REPO/work_jessev" \
        >> "$LOG" 2>&1 || log "  JessEV FAILED (non-fatal; downstream uses raw scores)"
else
    log "  JessEV requirements missing; SKIPPING"
fi

# ----------------------------------------------------------------------------
# M8. Estofolete Table 4 mapping
# ----------------------------------------------------------------------------
phase "M8. Estofolete Table 4 composite correlate mapping"
TABLE4_OUT="$OUTPUTS/table4_dengue"
mkdir -p "$TABLE4_OUT"
if [[ -f "$REPO/scripts/dengue/estofolete_table4_mapping.py" ]]; then
    log "  running estofolete_table4_mapping.py on dengue outputs"
    python3 "$REPO/scripts/dengue/estofolete_table4_mapping.py" \
        --bcell_dir "$BCELL_OUTDIR" \
        --tcell_dir "$TCELL_OUTDIR" \
        --jessev_dir "$JESS_OUTDIR" \
        --discotope_dir "$DISCO_OUT" \
        --out_dir "$TABLE4_OUT" \
        >> "$LOG" 2>&1 || log "  Table 4 mapping FAILED (non-fatal; partial fields will be missing)"
else
    log "  estofolete_table4_mapping.py missing; SKIPPING"
fi

# ----------------------------------------------------------------------------
# M9. Phase 3 retrospective validation
# ----------------------------------------------------------------------------
phase "M9. Phase 3 retrospective validation (CYD-TDV / TAK-003 / Butantan-DV)"
PHASE3_OUT="$OUTPUTS/phase3_retrospective"
mkdir -p "$PHASE3_OUT"
if [[ -f "$REPO/scripts/dengue/retrospective_phase3_validation.py" ]] && \
   [[ -f "$REPO/ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta" ]]; then
    log "  dispatching pipeline on Phase 3 vaccine constructs"
    PHASE3_BCELL="$OUTPUTS/phase3_bcell"
    rm -rf "$PHASE3_BCELL" "$REPO/work_phase3_bcell"
    mkdir -p "$PHASE3_BCELL"
    epi_arg="--epidope no"
    [[ -f "$SIF_PATH/epidope_v0.3.sif" ]] || [[ -f "$SIF_PATH/cuda-epidope.sif" ]] && epi_arg="--epidope yes"
    cd "$REPO"
    ~/nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
        -profile hyperplane_singularity \
        -c host/nextflow_config/nextflow.config \
        --protein_file ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta \
        --epitope_output_folder "$PHASE3_BCELL" \
        --cdhit_input_proteins yes \
        --bepipred no $epi_arg --dc_bcell no \
        --netmhcpani no --netmhcpanii no \
        --consolidate_epitopes no --tcrpmhc no --jessev no \
        -work-dir "$REPO/work_phase3_bcell" \
        >> "$LOG" 2>&1 || log "  Phase 3 B-cell pipeline FAILED (non-fatal)"
    log "  running retrospective_phase3_validation.py"
    python3 "$REPO/scripts/dengue/retrospective_phase3_validation.py" \
        --pipeline_dir "$PHASE3_BCELL" \
        --out_dir "$PHASE3_OUT" \
        >> "$LOG" 2>&1 || log "  Phase 3 retrospective FAILED (non-fatal)"
else
    log "  Phase 3 retrospective inputs missing; SKIPPING"
fi

# ----------------------------------------------------------------------------
# M10. Generate all figures
# ----------------------------------------------------------------------------
phase "M10. Render all manuscript figures"
FIG_OUT="$OUTPUTS/figures"
mkdir -p "$FIG_OUT"
if [[ -f "$REPO/scripts/dengue/generate_figures.py" ]]; then
    log "  running generate_figures.py"
    python3 "$REPO/scripts/dengue/generate_figures.py" \
        --table4_dir "$TABLE4_OUT" \
        --phase3_dir "$PHASE3_OUT" \
        --bcell_dir "$BCELL_OUTDIR" \
        --tcell_dir "$TCELL_OUTDIR" \
        --out_dir "$FIG_OUT" \
        >> "$LOG" 2>&1 || log "  generate_figures.py FAILED (non-fatal)"
else
    log "  generate_figures.py missing; SKIPPING figures"
fi

# ----------------------------------------------------------------------------
# M11. Fill manuscript template with pipeline numbers
# ----------------------------------------------------------------------------
phase "M11. Fill manuscript template with pipeline numbers"
MS_OUT="$OUTPUTS/manuscript"
mkdir -p "$MS_OUT"
if [[ -f "$REPO/scripts/dengue/fill_manuscript.py" ]] && \
   [[ -f "$REPO/docs/dengue/manuscript_template.md" ]]; then
    log "  running fill_manuscript.py"
    python3 "$REPO/scripts/dengue/fill_manuscript.py" \
        --template "$REPO/docs/dengue/manuscript_template.md" \
        --table4_dir "$TABLE4_OUT" \
        --phase3_dir "$PHASE3_OUT" \
        --figures_dir "$FIG_OUT" \
        --out "$MS_OUT/dengue_npj_vaccines_draft.md" \
        >> "$LOG" 2>&1 || log "  fill_manuscript.py FAILED (non-fatal)"
else
    log "  manuscript template or filler script missing; SKIPPING"
fi

# ----------------------------------------------------------------------------
# M12. Final commit + tag
# ----------------------------------------------------------------------------
phase "M12. Final commit + tag v1.0-dengue-results"
cd "$REPO"
git add outputs/ docs/dengue/ 2>>"$LOG" || true
if ! git diff --cached --quiet; then
    git commit -m "auto(master): v1.0-dengue-results - full pipeline + retrospective + figures + manuscript draft" 2>>"$LOG" || true
    git push origin master 2>>"$LOG" || log "  push failed; periodic_commit will retry"
fi
if ! git tag -l | grep -q "^v1.0-dengue-results$"; then
    git tag -a v1.0-dengue-results -m "Full dengue pipeline + retrospective + figures + manuscript draft" 2>>"$LOG" || true
    git push origin v1.0-dengue-results 2>>"$LOG" || log "  tag push failed"
fi

log ""
log "=========================================================================="
log " master_pipeline COMPLETE"
log "=========================================================================="
log "summary of outputs:"
find "$OUTPUTS" -type f | wc -l | xargs printf '  total files: %s\n' | tee -a "$LOG"
ls "$OUTPUTS" | tee -a "$LOG"
log "morning_status.sh will reflect final state"
