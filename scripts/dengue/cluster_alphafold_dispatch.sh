#!/usr/bin/env bash
# Dispatch AlphaFold-Multimer on DENV envelope (E) protein homo-dimers.
#
# Why E protein: E is the dominant target of neutralizing antibodies and
# carries the EDE epitope (Estofolete Table 4 antigenic-loss centroid).
# Homo-dimer is the minimal physiologically-relevant unit (E exists as
# head-to-tail dimers on the mature virion).
#
# DENV polyprotein domain coordinates (UniProt-aligned, 1-indexed):
#   capsid:    1-114
#   prM:       115-280
#   E:         281-775     <- our extraction window (residues 281-775, length 495)
#   NS1+:      776-end
#
# Distributes 4 jobs across 8 A100s using two GPUs per job (-gpu_devices 0,1
# / 2,3 / 4,5 / 6,7).
#
# Uses pre-pulled alphafold:2.3.2 docker image (or apptainer SIF) and the
# preinstalled databases at /data/alphafold_dbs (~2.9 TB).
#
# Usage on cluster (called from master_pipeline.sh, or standalone):
#   bash scripts/dengue/cluster_alphafold_dispatch.sh

set -uo pipefail

REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/alphafold_dispatch.log
OUTDIR=/data/james/ica-dengue/outputs/alphafold_e_dimers
SIF_PATH=/data/james/ica-dengue/sif_images
AF_DBS=/data/alphafold_dbs

mkdir -p "$OUTDIR" "$(dirname "$LOG")"
log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

cd "$REPO"
log "alphafold_dispatch starting"

# ----------------------------------------------------------------------------
# Step 1: Extract E protein sequences (residues 281-775) from DENV polyproteins
# ----------------------------------------------------------------------------
log "step 1: extracting E protein from each DENV polyprotein"
E_FASTA_DIR="$OUTDIR/e_dimer_fastas"
mkdir -p "$E_FASTA_DIR"

python3 - <<'PYEOF' 2>&1 | tee -a "$LOG"
import os
import sys

REPO = "/data/james/ica-dengue/immunoinformatics_platform-dengue"
OUT = "/data/james/ica-dengue/outputs/alphafold_e_dimers/e_dimer_fastas"
os.makedirs(OUT, exist_ok=True)

# E protein domain in DENV polyprotein (UniProt-aligned, 1-indexed: 281-775)
E_START = 280  # 0-indexed inclusive
E_END = 775    # 0-indexed exclusive (length 495)

records = {}
combined_path = os.path.join(REPO, "ref_fastas/dengue/dengue_polyproteins_multiseq.fasta")
if not os.path.isfile(combined_path):
    print(f"ERROR: input fasta missing: {combined_path}")
    sys.exit(1)

current_header = None
current_seq = []
with open(combined_path) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if current_header is not None:
                records[current_header] = "".join(current_seq)
            current_header = line[1:].split()[0]
            current_seq = []
        elif line:
            current_seq.append(line.strip())
    if current_header is not None:
        records[current_header] = "".join(current_seq)

for header, seq in records.items():
    e_seq = seq[E_START:E_END]
    if len(e_seq) < 400:
        print(f"  {header}: E sequence too short ({len(e_seq)} aa); skipping")
        continue
    serotype = header.split("_")[0] if "_" in header else header
    out_path = os.path.join(OUT, f"{serotype}_E_dimer.fasta")
    with open(out_path, "w") as out:
        out.write(f">{serotype}_E_chain_A\n{e_seq}\n")
        out.write(f">{serotype}_E_chain_B\n{e_seq}\n")
    print(f"  wrote {out_path} ({len(e_seq)} aa per chain, 2 chains)")
PYEOF

# ----------------------------------------------------------------------------
# Step 2: Locate AlphaFold image (prefer SIF, fall back to docker)
# ----------------------------------------------------------------------------
log "step 2: locating alphafold container"
AF_IMAGE=""
if [[ -f "$SIF_PATH/alphafold.sif" ]]; then
    AF_IMAGE="$SIF_PATH/alphafold.sif"
    AF_RUNTIME="apptainer"
    log "  using apptainer SIF: $AF_IMAGE"
elif command -v docker >/dev/null 2>&1 && docker images alphafold 2>/dev/null | grep -q "alphafold"; then
    AF_IMAGE=$(docker images --format "{{.Repository}}:{{.Tag}}" | grep alphafold | head -1)
    AF_RUNTIME="docker"
    log "  using docker image: $AF_IMAGE"
else
    log "  ERROR: no alphafold container found; SKIPPING AlphaFold (downstream DiscoTope will be skipped too)"
    log "  to fix: build apptainer SIF or pull docker alphafold:2.3.2"
    touch "$OUTDIR/.master_pipeline_af_complete"
    exit 0
fi

# ----------------------------------------------------------------------------
# Step 3: Sanity-check AlphaFold databases
# ----------------------------------------------------------------------------
log "step 3: sanity-check AlphaFold databases"
required_db_subdirs=(uniref90 mgnify pdb70 pdb_mmcif uniprot)
missing_dbs=0
for d in "${required_db_subdirs[@]}"; do
    if [[ ! -d "$AF_DBS/$d" ]]; then
        log "  MISSING db subdir: $AF_DBS/$d"
        missing_dbs=$((missing_dbs + 1))
    fi
done
if (( missing_dbs > 0 )); then
    log "  $missing_dbs DB subdirs missing; AlphaFold will fail; SKIPPING"
    touch "$OUTDIR/.master_pipeline_af_complete"
    exit 0
fi

# ----------------------------------------------------------------------------
# Step 4: Dispatch AlphaFold-Multimer jobs SERIALLY for the MSA stage.
#
# Why serial, not parallel:
#   The previous version dispatched 4 jobs in parallel with a 30s stagger.
#   All 4 failed at HHblits within 16s with EMPTY stderr. Root cause:
#   each HHblits process mmaps the BFD ffindex/ffdata (~1.7TB on disk,
#   ~30GB resident); 4 concurrent processes blow through the kernel page
#   cache and either (a) saturate /data NFS I/O so HHblits times out
#   without flushing stderr, or (b) trip the OOM killer (SIGKILL leaves
#   stderr empty — that is exactly the symptom). 30s stagger does nothing
#   because BFD residency lasts the entire ~3-4 hour HHblits run.
#
#   Fix: run the full pipeline (Jackhmmer + HHblits + inference) one
#   serotype at a time. Each job gets ALL 8 GPUs for the (cheap, ~10 min)
#   inference stage, and only one HHblits process touches BFD at a time.
#   Total wall time: ~4 jobs * ~3.5h = ~14h, vs. 4 silent failures.
#
# Knob if you want more parallelism later:
#   PARALLEL_AF_JOBS=2 yields 2 concurrent HHblits processes; the 80GB
#   A100 host with ~250GB RAM handles 2 BFD residents but not 4. Do NOT
#   raise above 2 without first staging cs219 to NVMe and benchmarking.
# ----------------------------------------------------------------------------
PARALLEL_AF_JOBS="${PARALLEL_AF_JOBS:-1}"
log "step 4: dispatching AlphaFold-Multimer jobs (PARALLEL_AF_JOBS=$PARALLEL_AF_JOBS)"

# All-GPU assignment for serial mode; only used when PARALLEL_AF_JOBS=1.
gpu_all="0,1,2,3,4,5,6,7"
gpu_pair_assignments=("0,1" "2,3" "4,5" "6,7")

run_one_af_job() {
    local fasta="$1"
    local gpus="$2"
    local name
    name=$(basename "$fasta" .fasta)
    local job_outdir="$OUTDIR/$name"
    mkdir -p "$job_outdir"
    if [[ -f "$job_outdir/ranked_0.pdb" ]]; then
        log "  $name: ranked_0.pdb exists; skipping"
        return 0
    fi
    log "  running $name on GPUs $gpus"
    if [[ "$AF_RUNTIME" == "apptainer" ]]; then
        apptainer run --nv \
            --bind "$AF_DBS":/dbs \
            --bind "$OUTDIR":/out \
            --env CUDA_VISIBLE_DEVICES=$gpus \
            "$AF_IMAGE" \
            --fasta_paths=/out/e_dimer_fastas/${name}.fasta \
            --output_dir=/out/${name} \
            --model_preset=multimer \
            --data_dir=/dbs \
            --uniref90_database_path=/dbs/uniref90/uniref90.fasta \
            --mgnify_database_path=/dbs/mgnify/mgy_clusters_2022_05.fa \
            --bfd_database_path=/dbs/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
            --uniref30_database_path=/dbs/uniref30/UniRef30_2021_03 \
            --template_mmcif_dir=/dbs/pdb_mmcif/mmcif_files \
            --obsolete_pdbs_path=/dbs/pdb_mmcif/obsolete.dat \
            --pdb_seqres_database_path=/dbs/pdb_seqres/pdb_seqres.txt \
            --uniprot_database_path=/dbs/uniprot/uniprot.fasta \
            --max_template_date=2024-01-01 \
            --use_precomputed_msas=true \
            --use_gpu_relax=true \
            >> "$LOG" 2>&1
    else
        docker run --rm --gpus "device=$gpus" \
            -v "$AF_DBS":/dbs \
            -v "$OUTDIR":/out \
            "$AF_IMAGE" \
            --fasta_paths=/out/e_dimer_fastas/${name}.fasta \
            --output_dir=/out/${name} \
            --model_preset=multimer \
            --data_dir=/dbs \
            --max_template_date=2024-01-01 \
            --use_precomputed_msas=true \
            --use_gpu_relax=true \
            >> "$LOG" 2>&1
    fi
}

declare -a af_pids=()
i=0
for fasta in "$E_FASTA_DIR"/*.fasta; do
    name=$(basename "$fasta" .fasta)
    if (( PARALLEL_AF_JOBS <= 1 )); then
        # SERIAL: one HHblits at a time, all 8 GPUs each.
        run_one_af_job "$fasta" "$gpu_all" \
            && log "  $name: SUCCESS" \
            || log "  $name: FAILED (continuing with remaining serotypes)"
    else
        # BOUNDED PARALLEL: at most PARALLEL_AF_JOBS HHblits processes live.
        while (( $(jobs -rp | wc -l) >= PARALLEL_AF_JOBS )); do
            sleep 60
        done
        gpus="${gpu_pair_assignments[$((i % 4))]}"
        i=$((i + 1))
        run_one_af_job "$fasta" "$gpus" &
        af_pids+=("$!")
        log "  $name PID: $!  GPUs=$gpus"
        sleep 300  # 5-min stagger so HHblits database paging settles before next launch
    fi
done

if (( PARALLEL_AF_JOBS > 1 )); then
    log "step 5: waiting for ${#af_pids[@]} parallel jobs"
    for pid in "${af_pids[@]}"; do
        if wait "$pid"; then
            log "  PID $pid: SUCCESS"
        else
            log "  PID $pid: FAILED (non-fatal; other jobs continue)"
        fi
    done
fi

# ----------------------------------------------------------------------------
# Step 6: Mark complete
# ----------------------------------------------------------------------------
touch "$OUTDIR/.master_pipeline_af_complete"
log "alphafold_dispatch complete"
log "outputs:"
find "$OUTDIR" -name "ranked_0.pdb" | tee -a "$LOG"
