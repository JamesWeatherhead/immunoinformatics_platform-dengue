# Zenodo deposition manifest

## Title

ICA-SP: Sequence-only proxies of the Estofolete composite correlate for dengue vaccine candidate prioritisation (data + code)

## Version

`v1.0-dengue-results`

## License

- Data and outputs (TSV, JSON, FASTA, figures, supplementary tables): **CC-BY-4.0**
- Source code, Nextflow modules, Python and shell scripts, configuration files: **MIT License**

## Description (~200 words)

This deposit accompanies the npj Vaccines submission "ICA-SP: Sequence-only proxies of the Estofolete composite immune correlate for dengue vaccine candidate prioritisation." It packages the complete reproducibility bundle for an Immunogenicity-Compartment Analysis (ICA) of the four-row composite correlate proposed by Estofolete, Nogueira and colleagues (Estofolete et al., 2026, *npj Vaccines*, DOI 10.1038/s41541-026-01400-4), restricted to features computable from primary amino acid sequence alone.

The bundle includes (i) the four DENV reference polyproteins, three Phase 3 parent strains (CYD-TDV, TAK-003, Butantan-DV), and AFND release 2024-12 regional HLA frequency tables used as inputs; (ii) ranked candidate tables, EDE-restricted antigenic-loss matrices, HLA population coverage by region, and CD-HIT cluster outputs for all four serotypes; (iii) a Phase 3 retrospective harness mapping pipeline scores to clinical efficacy; (iv) ablation, permutation, sensitivity, and IEDB-validation supplementary tables (S1-S5); (v) all five main figures as PNG and SVG; and (vi) the dengue fork of the McCaffrey `immunoinformatics_platform` Nextflow pipeline that produced these outputs. The deposit is intended to support independent re-execution, audit of the Estofolete Table 4 mapping, and re-use of the EDE-restricted scoring window in downstream flavivirus vaccine prioritisation work.

## Keywords

dengue, vaccine, composite correlate, immunoinformatics, EDE, EDE1, EDE2, nextflow, pipeline, Estofolete, Phase 3 retrospective, HLA population coverage, AFND, IEDB, CD8 T-cell epitopes, B-cell epitopes, ICA, sequence-only proxy

## Related identifiers

- **Manuscript GitHub repository (this fork):** https://github.com/JamesWeatherhead/immunoinformatics_platform-dengue (`isSupplementTo`)
- **Upstream pipeline (Versiani, McCaffrey et al. 2026):** https://github.com/pmccaffrey6/immunoinformatics_platform (`isDerivedFrom`)
- **Estofolete et al. 2026 npj Vaccines (composite correlate framework):** https://doi.org/10.1038/s41541-026-01400-4 (`isSupplementTo`)
- **AFND release 2024-12 (HLA frequencies):** http://www.allelefrequencies.net (`references`)
- **IEDB (T-cell and B-cell epitope queries):** https://www.iedb.org (`references`)

## File inventory

```
outputs/manuscript/dengue_npj_vaccines_FINAL.md
outputs/manuscript/dengue_npj_vaccines_draft.md
outputs/manuscript/methods_reproducibility.md
outputs/manuscript/.claude_round1_log.txt
outputs/manuscript/supplementary/S1_ablation_studies.md
outputs/manuscript/supplementary/S1_ablation_studies.tsv
outputs/manuscript/supplementary/S2_per_construct_breakdown.md
outputs/manuscript/supplementary/S2_per_construct_breakdown.tsv
outputs/manuscript/supplementary/S3_citation_audit.md
outputs/manuscript/supplementary/S3_citation_audit.tsv
outputs/manuscript/supplementary/S4_sensitivity_analysis.md
outputs/manuscript/supplementary/S4_sensitivity_analysis.tsv
outputs/manuscript/supplementary/S5_iedb_validation.md
outputs/manuscript/supplementary/S5_iedb_validation.tsv
outputs/manuscript/supplementary/population_coverage_v2.md
outputs/figures/fig1_pipeline_overview.png
outputs/figures/fig1_pipeline_overview.svg
outputs/figures/fig2_per_serotype_tier_scores.png
outputs/figures/fig2_per_serotype_tier_scores.svg
outputs/figures/fig3_phase3_retrospective.png
outputs/figures/fig3_phase3_retrospective.svg
outputs/figures/fig4_hla_coverage_heatmap.png
outputs/figures/fig4_hla_coverage_heatmap.svg
outputs/figures/fig5_ede_loss_matrix.png
outputs/figures/fig5_ede_loss_matrix.svg
outputs/figures/fig6_cdhit_clusters.png
outputs/figures/fig6_cdhit_clusters.svg
outputs/phase3_retrospective/statistical_analyses.json
outputs/phase3_retrospective/phase3_vs_clinical.tsv
outputs/phase3_retrospective/table4_mapping/ranked_candidates.tsv
outputs/dengue_smoke_evidence/cdhit_clusters_dengue.clstr
outputs/dengue_smoke_evidence/cdhit_clusters_dengue.tsv
outputs/dengue_smoke_evidence/dengue_processed.fasta
outputs/dengue_smoke_evidence/join_tables/input_fasta_table.tsv
outputs/dengue_smoke_evidence/join_tables/Brazil_regional_allele_frequencies.tsv
outputs/table4_dengue/hla_coverage_by_region.tsv
outputs/table4_dengue/hla_coverage_by_region_v2.tsv
outputs/table4_dengue/hla_coverage_by_region_v2.json
outputs/table4_dengue/hla_coverage_by_region_v2.log
outputs/table4_dengue/ede_antigenic_loss_matrix.tsv
outputs/table4_dengue/ranked_candidates.tsv
```

The full source tree (Nextflow workflows, modules, Python helpers, container recipes, `ref_fastas/`, `conf/`, and the Phase 3 retrospective harness) is included as the GitHub release tarball at git tag `v1.0-dengue-results`.

## Citation

Weatherhead J., McCaffrey P., Weiskopf D., Estofolete C.F., Nogueira M.L., Vasilakis N. (2026). *ICA-SP: Sequence-only proxies of the Estofolete composite correlate for dengue vaccine candidate prioritisation (data + code)*. Zenodo. https://doi.org/[Zenodo DOI to be assigned on publication].
