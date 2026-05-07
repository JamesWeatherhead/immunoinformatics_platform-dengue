# hsv-server cluster runbook for the dengue extension

## State as of v1.0-dengue-prep tag

Workspace: `/data/james/ica-dengue/immunoinformatics_platform-dengue/`

### Verified ready (2026-05-06)
- Fork cloned and at tag `v1.0-dengue-prep`
- Dengue inputs fetched: 4 reference polyproteins + **12 of 12** Phase 3 vaccine
  parent strain sequences (CYD-TDV DENV-1 recovered via nucleotide-translate fallback)
- Nextflow 23.10.1 working at `~/nextflow`
- Apptainer 1.1.7 ready
- Pre-built Docker images on cluster:
  - `pmccaffrey6/jess_ev:latest` (6.93 GB)
  - `pmccaffrey6/cuda-epidope:latest` (16.5 GB)
  - `pmccaffrey6/polisher:latest` (21 GB)
  - `flomock/epidope:v0.3` (7.58 GB)
  - `alphafold:2.3.2` / `alphafold_23:latest` (6.33 / 10.3 GB)
- AlphaFold MSA databases at `/data/alphafold_dbs/` (2.9 TB; bfd, mgnify, pdb_mmcif, params multimer v3)
- Conda envs available: boltz, cellbender, mari-app, hsv_ad

### Outstanding blockers (require user action)

1. **NetMHCpan-4.1 academic license**
   - Register: https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCpan&version=4.1
   - 24-48h email turnaround
   - Required by `containers/netmhcpan_i_4.1/`

2. **NetMHCIIpan-4.3 academic license** (separate registration; same site)

3. **Gurobi WLS license** for JessEV epitope selection
   - First ask: McCaffrey (he authored Versiani et al. and his lab license is what was leaked in upstream)
   - Else register: https://www.gurobi.com/academia/academic-program-and-licenses/
   - Set `WLSACCESSID`, `WLSSECRET`, `LICENSEID` env vars at runtime; do NOT commit

4. **Upstream pipeline bug** at line 1146 of `nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf`
   - References `bepipred_tsv_ch` even when `--bepipred no`
   - Workaround: run with all stages enabled (the intended use case), OR patch
   - Not blocking for actual full pipeline runs; only matters for partial smoke tests

5. **GitHub credentials on cluster** — `git push` from `/data/james/ica-dengue/...`
   currently fails because no auth helper is configured. To enable result commits
   from cluster:
   ```bash
   gh auth login            # follow prompts, paste a personal access token
   git config --global credential.helper store
   ```

### Dispatch when licenses arrive

```bash
cd /data/james/ica-dengue/immunoinformatics_platform-dengue

# Copy license files into the relevant container build contexts:
cp /home/pathinformatics/licenses/netMHCpan-4.1a.Linux.tar.gz containers/netmhcpan_i_4.1/
cp /home/pathinformatics/licenses/netMHCIIpan-4.3.Linux.tar.gz containers/netmhcpan_ii_4.1/

# Set Gurobi env:
export WLSACCESSID=...
export WLSSECRET=...
export LICENSEID=...

mkdir -p /data/james/ica-dengue/{outputs,logs}

# Run main pipeline on dengue references (4 polyproteins):
tmux new -s dengue_main -d "
  ~/nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    --protein_file ref_fastas/dengue/dengue_polyproteins_multiseq.fasta \
    --epitope_output_folder /data/james/ica-dengue/outputs/dengue_main \
    --allele_target_region 'Brazil,Thailand,Vietnam,Philippines,Sri_Lanka' \
    --cdhit_input_proteins yes --bepipred yes --epidope yes \
    --netmhcpani yes --netmhcpanii yes --dc_bcell yes \
    --consolidate_epitopes yes --tcrpmhc yes --jessev yes \
    -profile apptainer 2>&1 | tee /data/james/ica-dengue/logs/dengue_main_\$(date +%Y%m%d_%H%M).log
"

# Run on Phase 3 vaccine parent strains (12 polyproteins):
tmux new -s dengue_phase3 -d "
  ~/nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    --protein_file ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta \
    --epitope_output_folder /data/james/ica-dengue/outputs/dengue_phase3 \
    --allele_target_region 'Brazil,Thailand,Philippines' \
    --cdhit_input_proteins yes --bepipred yes --epidope yes \
    --netmhcpani yes --netmhcpanii yes --dc_bcell yes \
    --consolidate_epitopes yes --tcrpmhc yes --jessev yes \
    -profile apptainer 2>&1 | tee /data/james/ica-dengue/logs/dengue_phase3_\$(date +%Y%m%d_%H%M).log
"
```

Expected runtime per Versiani et al.: ~48 hours per run on the listed
hardware (256 cores, 1000 GB RAM, 4x A100 80GB). Dengue inputs are
much smaller than the alphavirus run (4 + 12 vs 53 proteomes); expect
~12-24 hours per run.

### After runs complete

```bash
# Map main pipeline output to Estofolete Table 4:
python3 scripts/dengue/estofolete_table4_mapping.py \
    --epitope-output-folder /data/james/ica-dengue/outputs/dengue_main \
    --outdir /data/james/ica-dengue/outputs/dengue_table4

# Run retrospective Phase 3 validation:
python3 scripts/dengue/retrospective_phase3_validation.py \
    --epitope-output-folder /data/james/ica-dengue/outputs/dengue_phase3 \
    --outdir /data/james/ica-dengue/outputs/retrospective_phase3

# Commit results back to fork (after gh auth login on cluster):
git add outputs/dengue_table4/ outputs/retrospective_phase3/
git commit -m "v1.0-dengue results: pipeline runs on dengue references and 3 Phase 3 vaccines"
git push origin master
git tag -a v1.0-dengue-results -m "Pipeline outputs from dengue + Phase 3 retrospective runs"
git push origin v1.0-dengue-results
```

### Disk planning
/data has 1.6 TB free. Pipeline outputs per protein ~50-100 GB with
intermediates. For 4 + 12 = 16 sequences, plan for ~800 GB total.
If pressure: clean intermediates after the table4 mapping script extracts
the per-construct scoring TSVs (the rest can be regenerated).
