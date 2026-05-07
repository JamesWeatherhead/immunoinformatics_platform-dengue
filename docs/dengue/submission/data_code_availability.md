# Data and code availability

## Code

All code underlying this Perspective is publicly available at:

- **Repository:** https://github.com/JamesWeatherhead/immunoinformatics_platform-dengue
- **Tag:** `v1.0-dengue-results` (frozen at submission; all figures and numbers in the manuscript correspond to outputs produced at this tag).
- **Pipeline pedigree:** fork of `pmccaffrey6/immunoinformatics_platform` (Versiani, McCaffrey et al. 2026 *Sci Adv* 12:eaeb2066). Upstream history is preserved; the fork's diff is contained to the dengue extension files documented in `docs/dengue/DENGUE_EXTENSION.md`.
- **Implementation:** Nextflow DSL2 with all stages containerised (Apptainer/Singularity).
- **License:** the source code retains the upstream license; the dengue extension files are released under the same terms.

## Data

All output data underlying every figure and number in the manuscript are committed to the same repository under the `outputs/` subtree, at git tag `v1.0-dengue-results`:

- `outputs/table4_dengue/ranked_candidates.tsv` (Table 1, per-serotype Tier A and Tier B sub-scores).
- `outputs/phase3_retrospective/` (Figure 3, Phase 3 efficacy regression).
- `outputs/figures/` (Figures 2-5 source PNG and SVG).
- `outputs/manuscript/` (auto-filled manuscript draft and FINAL).

Output data are released under **CC-BY-4.0**.

## Pipeline DOI

A persistent DOI for the dengue fork at tag `v1.0-dengue-results` will be deposited at Zenodo (https://zenodo.org) prior to acceptance, via the Zenodo-GitHub integration. The DOI will be added to the final published manuscript and to this document. (Pre-acceptance: the GitHub tag `v1.0-dengue-results` serves as the canonical reproducibility anchor.)

## Container images

Pre-built Apptainer SIF images for every pipeline stage (EpiDope, BepiPred-3.1, AlphaFold-Multimer, DiscoTope, NetMHCpan-4.1, NetMHCIIpan-4.1, IEDB Population Coverage, JessEV) are stored at `/data/james/ica-dengue/sif_images/` on the UTMB Pathology Informatics GPU cluster. SIF images are available on request to the corresponding authors, subject to the upstream tools' academic licenses (Gurobi/JessEV solver and NetMHCpan binaries are used under Peter McCaffrey's existing academic licenses; downstream users must obtain their own licenses for these stages). Container build recipes are documented in `docs/dengue/CONTAINER_BUILD.md`.

## External data sources

- **AFND HLA frequencies.** Allele Frequency Net Database release **2024-12**, accessed via http://www.allelefrequencies.net (Gonzalez-Galarza et al. 2020 *Nucleic Acids Res* 48:D783-D788). The 16 primary HLA-typing papers underlying the regional coverage analysis are documented at `docs/dengue/paperclip_research/hla_primary_measurements.md`.
- **IEDB.** Immune Epitope Database release **3.1.x** (Vita et al. 2019 *Nucleic Acids Res* 47:D339-D343); IEDB Population Coverage tool used per Bui et al. 2006 *BMC Bioinformatics* 7:153.
- **UniProt reference polyproteins.** DENV-1 P17763, DENV-2 P29990, DENV-3 P27915, DENV-4 P09866 (UniProt release current at submission).
- **GenBank Phase 3 vaccine constructs.** Accessions documented in `ref_fastas/dengue_phase3_vaccines/PROVENANCE.md`.
- **PDB.** DENV-2 E dimer reference 1OAN (Modis et al. 2003).

## Reproducibility

The full pipeline can be re-run by any reader with access to a Nextflow-capable cluster and the upstream tool licenses noted above. The command to reproduce all manuscript outputs from a clean checkout is documented in `docs/dengue/CLUSTER_RUNBOOK.md`. Citation audits and the response to prior reviewers are committed at `docs/dengue/citation_audit.md` and `docs/dengue/response_to_prior_reviewers.md` respectively.
