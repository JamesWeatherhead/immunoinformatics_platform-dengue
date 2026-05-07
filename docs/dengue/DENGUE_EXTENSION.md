# Dengue extension of immunoinformatics_platform

This fork extends the published pan-alphavirus pipeline of Versiani et al.
(*Sci. Adv.* 2026; 12:eaeb2066, [doi:10.1126/sciadv.aeb2066](https://doi.org/10.1126/sciadv.aeb2066))
to dengue, mapping the pipeline's epitope outputs to the four-row composite
correlate of dengue vaccine efficacy proposed by Estofolete et al.
(*npj Vaccines* 2026; 11:68, [doi:10.1038/s41541-026-01400-4](https://doi.org/10.1038/s41541-026-01400-4)).

## What is added on top of the upstream pipeline

| Addition | Path | Purpose |
|---|---|---|
| DENV 1-4 reference polyproteins | `ref_fastas/dengue/` | UniProt SwissProt P17763, P29990, P27915, P09866 + concatenated multi-FASTA for pipeline ingestion |
| Phase 3 vaccine parent strain CDS sequences | `ref_fastas/dengue_phase3_vaccines/` | GenBank accessions for CYD-TDV, TAK-003, Butantan-DV parent strains, used for retrospective scoring |
| Fetch script | `scripts/dengue/fetch_dengue_data.sh` | Pulls the above from UniProt + GenBank with provenance documentation |
| Estofolete Table 4 mapping layer | `scripts/dengue/estofolete_table4_mapping.py` | Joins upstream pipeline outputs (DiscoTope, ESMFold pLDDT, NetMHCpan I/II, AFND-weighted coverage) into per-candidate Tier A/B scores |
| Retrospective Phase 3 validation | `scripts/dengue/retrospective_phase3_validation.py` | Runs the mapping layer across the 12 Phase 3 vaccine parent strain sequences and produces a comparison table against published clinical efficacy |
| Documentation | `docs/dengue/DENGUE_EXTENSION.md` (this file) | Methodology and provenance |

## What is NOT added (and why)

Estofolete Table 4 has four rows. This pipeline computes two of them
pre-trial from public sequence + database data:

| Tier 4 row | Pipeline coverage | Why |
|---|---|---|
| Tier A neutralization breadth | YES (Geometry axis: EDE-residue B-cell epitope coverage on predicted dimer) | Computable from sequence alone |
| Tier A avidity | NO (deferred) | Requires post-Phase-1 vaccinee serum + chaotrope-displacement ELISA |
| Tier B memory B-cell breadth | NO (deferred) | Requires post-Phase-1 vaccinee BCR-seq |
| Tier B CD8 polyfunctionality | YES (Equity axis: NetMHCpan I + II x AFND geographic weighting) | Computable from sequence + AFND alone (with the explicit caveat that the binding-prediction proxy is not the same construct as the measured ICS readout in Estofolete Table 4; see manuscript Methods) |

The two deferred axes are described in the manuscript Methods as the
prospective minimum-viable-cohort (MVC) validation pathway. They are
**not** implemented as code in this fork because no public data exists
that could feed them for a never-administered vaccine candidate. The
Versiani et al. upstream pipeline implements the same two pre-trial axes
(Geometry, Equity) for alphaviruses; the dengue extension is the
input-swap plus Estofolete-Table-4 mapping.

## How to run

### Prerequisites

1. Clone this fork
2. Install the upstream pipeline dependencies (Nextflow, Apptainer/Singularity or Docker, Java 11+)
3. Obtain academic licenses for NetMHCpan-4.1 and NetMHCIIpan-4.3
   (https://services.healthtech.dtu.dk/) — the binaries cannot be redistributed
4. Obtain a Gurobi license for the JessEV epitope-selection step
   (https://www.gurobi.com/academia/academic-program-and-licenses/)
5. AlphaFold/ESMFold/DiscoTope databases as required by the upstream pipeline

### Fetch dengue inputs

```bash
bash scripts/dengue/fetch_dengue_data.sh
```

This populates `ref_fastas/dengue/` and `ref_fastas/dengue_phase3_vaccines/`.

### Run the pipeline on dengue inputs

```bash
nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    --protein_file ref_fastas/dengue/dengue_polyproteins_multiseq.fasta \
    --epitope_output_folder /path/to/epitope_outputs/dengue_main \
    --allele_target_region "Brazil,Thailand,Vietnam,Philippines,Indonesia,Malaysia,Sri_Lanka,India,Mexico,Colombia,Brazil"
```

(Adjust `--allele_target_region` to the comma-separated AFND population
search terms appropriate to the endemic regions of interest.)

### Apply the Estofolete Table 4 mapping

```bash
python scripts/dengue/estofolete_table4_mapping.py \
    --epitope-output-folder /path/to/epitope_outputs/dengue_main \
    --thresholds-yaml       config/dengue_thresholds.yaml \
    --outdir                outputs/dengue_table4/
```

Outputs:
- `outputs/dengue_table4/ranked_candidates.tsv` — per-construct Tier A/B scores
- `outputs/dengue_table4/ranked_candidates_summary.md` — Markdown summary
- `outputs/dengue_table4/tier_breakdown.json` — per-construct provenance

### Retrospective Phase 3 validation

```bash
# Run pipeline on the 12 Phase 3 vaccine parent strain sequences:
nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    --protein_file ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta \
    --epitope_output_folder /path/to/epitope_outputs/dengue_phase3 \
    --allele_target_region "Brazil,Thailand,Philippines"

# Then map and compare to published efficacy:
python scripts/dengue/retrospective_phase3_validation.py \
    --epitope-output-folder /path/to/epitope_outputs/dengue_phase3 \
    --thresholds-yaml       config/dengue_thresholds.yaml \
    --outdir                outputs/retrospective_phase3/
```

This produces a Markdown comparison table showing whether the pipeline
recovers the relative immunogenicity ordering documented in the Phase 3
record (CYD-TDV pooled ~56% efficacy, TAK-003 4.5y ~61%, Butantan-DV
single-dose ~80%).

## Attribution

The pipeline core (containers, Nextflow workflow, allele frequency
machinery, JessEV optimization) is from the published upstream:

> Versiani AF, McCaffrey P, Ribeiro-Filho HV, Silva NIO, Lopes-de-Oliveira PS, Carrera J-P, Nogueira ML, Marques RE, Rossi SL, Vasilakis N. Integrated reiterative pipeline for rapid epitope-based pan-alphavirus vaccines. *Sci. Adv.* 2026; 12(11):eaeb2066. [doi:10.1126/sciadv.aeb2066](https://doi.org/10.1126/sciadv.aeb2066).

Code: [pmccaffrey6/immunoinformatics_platform](https://github.com/pmccaffrey6/immunoinformatics_platform).

The dengue extension and Estofolete Table 4 mapping layer are added in this
fork. The conceptual framework is from the Saivish et al. dengue
Perspective (in revision); the empirical companion the Table 4 mapping
implements is:

> Estofolete CF, Saivish MV, Nogueira ML, Vasilakis N. From promise to pitfalls: immunological lessons from dengue vaccines and their implications. *npj Vaccines* 2026; 11:68. [doi:10.1038/s41541-026-01400-4](https://doi.org/10.1038/s41541-026-01400-4).

## Security note

The upstream repo committed Gurobi license credentials in
`containers/jess_ev/src/gurobi.lic` and the `Dockerfile`'s `ENV`
variables. This fork has scrubbed those credentials in HEAD; the upstream
git history still contains them. The original credentials should be
considered compromised and rotated by their owner. See
`containers/jess_ev/README.md` for how to supply your own license at
runtime.
