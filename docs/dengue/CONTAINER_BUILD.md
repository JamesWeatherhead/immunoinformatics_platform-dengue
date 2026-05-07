# Building the SIF (Singularity/Apptainer) images for the pipeline

The Versiani et al. pipeline runs each step in its own container, mapped via
`-profile hyperplane_singularity` in `host/nextflow_config/nextflow.config`.
The profile expects pre-built `.sif` files at `$SINGULARITY_SIF_PATH`.

## What's already built on hsv-server (Docker images, need conversion to SIF)

The following Pete-namespace Docker images are present on hsv-server and
can be converted to .sif with `apptainer build`:

| Docker image | SIF target name | Maps to process |
|---|---|---|
| `pmccaffrey6/jess_ev:latest` | `jess_ev_latest.sif` | (JessEV; gated by --jessev) |
| `pmccaffrey6/cuda-epidope:latest` | `epidope_v0.3.sif` | EPIDOPE |
| `pmccaffrey6/polisher:latest` | (not in profile; legacy) | - |
| `flomock/epidope:v0.3` | (alternate EPIDOPE) | EPIDOPE |
| `alphafold:2.3.2` | `alphafold_2.3.2.sif` | (AlphaFold predict; not gated by profile mappings, see DISCOTOPE) |

## What needs to be built from this repo's containers/ directory

Each subdirectory under `containers/` contains a Dockerfile that builds one
of the pipeline tools. To convert to .sif, build the Docker image first
then convert. Example for `cd-hit`:

```bash
export SINGULARITY_SIF_PATH=/data/james/ica-dengue/sif_images
mkdir -p "$SINGULARITY_SIF_PATH"

cd containers/cd-hit
docker build -t cd-hit-4-8-1 .
apptainer build "$SINGULARITY_SIF_PATH/cd-hit-4-8-1.sif" docker-daemon://cd-hit-4-8-1:latest
```

Apply the same pattern to every container that the profile references. The
full mapping (target SIF name -> containers/ source dir):

| SIF target | containers/ source |
|---|---|
| `cd-hit-4-8-1.sif` | `containers/cd-hit/` |
| `cdhit-to-tsv.sif` | `containers/cdhit_to_tsv/` |
| `epidope_v0.3.sif` | `containers/epidope_03/` (or pull from `flomock/epidope:v0.3`) |
| `iedb-b-cell_latest.sif` | `containers/iedb_b_cell_31/` |
| `bepipred-to-tsv_latest.sif` | `containers/bepipred_to_tsv/` |
| `discotope_1.1.sif` | `containers/discotope_11/` |
| `discotope_to_tsv_latest.sif` | `containers/discotope_to_tsv/` |
| `iedb-mhc-i_latest.sif` | `containers/iedb_mhc_i_312/` |
| `iedb-mhc-ii_latest.sif` | `containers/iedb_mhc_ii_316/` |
| `netmhcpan_i_4.1.sif` | `containers/netmhcpan_i_4.1/` (NEEDS NetMHCpan-4.1 license tarball) |
| `netmhcpan_ii_4.1.sif` | `containers/netmhcpan_ii_4.1/` (NEEDS NetMHCIIpan-4.3 license tarball) |
| `split_fastas_latest.sif` | `containers/split_fastas/` |
| `tcrpmhc_1.0.sif` | `containers/tcr-pmhc-models/` |
| `pisa_server_latest.sif` | (PDB ePISA web service; may not need a container) |
| `blast_latest.sif` | (NCBI BLAST+; biocontainers/blast or build from Dockerfile) |
| `clustalomega_latest.sif` | (Clustal Omega; biocontainers/clustalo) |

## Build all in one pass

```bash
export SINGULARITY_SIF_PATH=/data/james/ica-dengue/sif_images
mkdir -p "$SINGULARITY_SIF_PATH"

cd /data/james/ica-dengue/immunoinformatics_platform-dengue

# Free containers (no licenses needed)
for c in cd-hit cdhit_to_tsv epidope_03 iedb_b_cell_31 bepipred_to_tsv \
         discotope_11 discotope_to_tsv iedb_mhc_i_312 iedb_mhc_ii_316 \
         split_fastas tcr-pmhc-models; do
    [ -d "containers/$c" ] || continue
    img=$(echo "$c" | tr '_' '-')
    echo "=== building $c -> $img ==="
    docker build -t "$img:latest" "containers/$c" || echo "  build failed for $c"
done

# Convert built Docker images to SIF
for c in cd-hit cdhit-to-tsv epidope-03 iedb-b-cell-31 bepipred-to-tsv \
         discotope-11 discotope-to-tsv iedb-mhc-i-312 iedb-mhc-ii-316 \
         split-fastas tcr-pmhc-models; do
    target_sif="$SINGULARITY_SIF_PATH/${c}.sif"
    [ -f "$target_sif" ] && { echo "exists: $target_sif"; continue; }
    apptainer build "$target_sif" "docker-daemon://${c}:latest" || echo "  convert failed for $c"
done

# License-gated containers (build only after license registration)
# These will fail until the NetMHCpan tarballs are placed in the build context:
for c in netmhcpan_i_4.1 netmhcpan_ii_4.1; do
    img=$(echo "$c" | tr '_' '-')
    echo "=== building $c -> $img (REQUIRES LICENSE TARBALL) ==="
    docker build -t "$img:latest" "containers/$c" 2>&1 | tail -3
done

# Pull pre-built upstream images that match expected SIF names
docker pull pmccaffrey6/jess_ev:latest
apptainer build "$SINGULARITY_SIF_PATH/jess_ev_latest.sif" docker-daemon://pmccaffrey6/jess_ev:latest
```

## Then dispatch with the right profile

```bash
export SINGULARITY_SIF_PATH=/data/james/ica-dengue/sif_images
~/nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    -profile hyperplane_singularity \
    -c host/nextflow_config/nextflow.config \
    --protein_file ref_fastas/dengue/dengue_polyproteins_multiseq.fasta \
    --epitope_output_folder /data/james/ica-dengue/outputs/dengue_main \
    --allele_target_region 'Brazil,Thailand,Vietnam,Philippines' \
    --cdhit_input_proteins yes --bepipred yes --epidope yes \
    --netmhcpani yes --netmhcpanii yes --dc_bcell yes \
    --consolidate_epitopes yes --tcrpmhc yes --jessev yes
```

## Smoke test that doesn't need licenses

For a fast plumbing check that doesn't need NetMHCpan or Gurobi:

```bash
~/nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    -profile hyperplane_singularity \
    -c host/nextflow_config/nextflow.config \
    --protein_file ref_fastas/dengue/dengue_polyproteins_multiseq.fasta \
    --epitope_output_folder /data/james/ica-dengue/outputs/dengue_smoke \
    --cdhit_input_proteins yes \
    --bepipred no --epidope no --netmhcpani no --netmhcpanii no --dc_bcell no \
    --consolidate_epitopes no --tcrpmhc no --jessev no
```

After the cd-hit SIF is built and `SINGULARITY_SIF_PATH` is set, this run
should succeed in seconds and produce `outputs/dengue_smoke/processed_fastas/`
plus a CDHIT cluster table.

## Disk planning for SIF storage

Each SIF is roughly the same size as the source Docker image:
- AlphaFold: 6-10 GB
- EpiDope (CUDA): 16 GB
- JessEV: 7 GB
- Polisher: 21 GB
- Smaller tools: 100 MB - 2 GB each

Total SIF cache: ~60-80 GB. /data has 1.6 TB free; comfortable.
