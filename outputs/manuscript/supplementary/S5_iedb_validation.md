# Supplementary Table S5. IEDB validation of pipeline-predicted T-cell epitopes

**Source files.** Pipeline-predicted strong-binder counts will be
drawn from `outputs/table4_dengue/ranked_candidates.tsv`
(`n_mhci_strong_binders`, `n_mhcii_strong_binders` columns) joined per
HLA allele to a curated IEDB cohort pull. The IEDB cohort pull
(experimentally validated dengue T-cell epitopes filtered by
host-organism = Homo sapiens, assay = T-cell-positive, restriction =
listed HLA allele) is **not yet pulled** in this artifact.

**STATUS: pending IEDB validation cohort pull.**

This table is committed as a stub with zero-fill rows so the
manuscript figure-counting and supplementary index lines up
correctly. The 8 listed alleles are the load-bearing HLAs in the
pipeline's population-coverage calculation: 5 MHC-I (HLA-A02:01,
A24:02, B07:02, B35:01, B57:01, the canonical Estofolete Section 2.3
panel) and 3 MHC-II (DRB1*07:01, DRB1*15:01, DPB1*04:01, the three
most-prevalent MHC-II loci in the AFND Brazilian regional allele
frequency table at
`outputs/dengue_smoke_evidence/join_tables/Brazil_regional_allele_frequencies.tsv`).

**Caption.** When the IEDB pull is performed, this table will report
per-allele counts of pipeline-predicted strong binders
(`n_predicted`), curated IEDB experimentally-validated dengue
epitopes for the same allele (`n_iedb_known`), the intersection
(`true_positives`), and per-allele precision and recall. Aggregate
precision and recall (micro-averaged over alleles) will be added to
Section 3 of the manuscript as the only external orthogonal
validation of the T-cell tier.

**Action item.** Pull IEDB curated dengue T-cell epitope set (filter:
organism = Dengue virus, host = Homo sapiens, qualitative assay =
positive, restriction = listed HLA) and rerun the join against
`ranked_candidates.tsv` strong-binder peptides before manuscript
submission. Estimated wall-clock time: 2-4 hours given existing IEDB
API tooling.
