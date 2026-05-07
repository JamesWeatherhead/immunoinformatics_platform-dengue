# Figure captions: dengue ICA Perspective (npj Vaccines submission)

All figures rendered at 300 DPI in PNG, PDF (vector), and SVG. Source pipeline: fork of Versiani 2026 immunoinformatics_platform; outputs at `outputs/figures/publication_quality/`.

## Figure 1. ICA-dengue empirical pipeline overview

Schematic of the five-axis Immunogenicity-Compartment Analysis (ICA) pipeline forked from the Versiani 2026 immunoinformatics platform. Inputs are DENV 1-4 polyproteins; pre-trial axes (Geometry: AlphaFold-Multimer + DiscoTope + EpiDope B-cell scoring; Equity: NetMHCpan I/II + IEDB population coverage across Brazil, Southeast Asia, Africa, India) are computable from sequence alone. Post-Phase-1 axes (Time: MBC Hill diversity via Alakazam; Intent: TCR-pMHC polyfunctionality; Effector: dendritic-cell B-cell maturation and opsonophagocytosis) require trial data and constitute future work. Output is the Estofolete Table 4 ranked Tier A and Tier B candidate list. Color coding: blue = inputs/outputs, yellow = Geometry, green = Equity, pink = post-trial axes.

## Figure 2. Per-DENV-serotype Estofolete Table 4 composite correlate scores

Bar plot of Tier A (neutralization breadth and antibody avidity) and Tier B (memory B-cell breadth and CD8 polyfunctionality) composite scores for each DENV serotype on a 0-1 scale. Higher composite values indicate closer match to the immune signature predictive of Phase 3 efficacy. Tier A is plotted in blue, Tier B in orange. Source: `outputs/table4_dengue/ranked_candidates.tsv`.

## Figure 3. Retrospective Phase 3 validation: predicted composite scores vs published clinical efficacy

Scatter plot of pipeline-predicted composite Tier A+B scores (x-axis) against published Phase 3 efficacy estimates (y-axis) for the three completed Phase 3 dengue vaccine programs: CYD-TDV (Sridhar 2018 NEJM), TAK-003 (Tricou 2024 Lancet Glob Health), and Butantan-DV (Kallas 2024 NEJM). Dashed gray line is the linear fit. Pearson r is annotated upper-left; permutation p (Tier A only, label-shuffled null) is annotated upper-right. Subtitle notes n=3 vaccine programs, which precludes formal inferential statistics; values are descriptive. Source: `outputs/phase3_retrospective/phase3_vs_clinical.tsv` and `statistical_analyses.json`.

## Figure 4. HLA Class I population coverage of DENV-1-4 polyproteins across 16 dengue-endemic regions

Single-column heatmap of mean IEDB population coverage of dengue polyprotein-derived MHC-I-restricted epitopes across 16 dengue-endemic regions (Mexico, Brazil overall plus five Brazilian sub-regions, India, Colombia, Thailand, Vietnam, Peru, Tanzania, Senegal, Nigeria). Cells are colored on an RdYlGn diverging scale from 0.70 to 1.00. Regions with coverage below the 0.78 equity threshold (Senegal 0.7769; Nigeria 0.7552) are outlined in red and flagged with an asterisk; these two West African populations fall short of the equity-axis inclusion criterion, motivating downstream allele-supplementation work. Source: `outputs/table4_dengue/hla_coverage_by_region_v2.tsv` (n_populations_aggregated and total_sample_size in source TSV).

## Figure 5. EDE epitope antigenic-loss matrix

4x4 heatmap of antigenic loss between vaccine source serotype (rows) and cross-reactive target serotype (columns) for the envelope dimer epitope (EDE) region. Values range from 0 (full cross-reactivity) to 1 (no cross-reactivity); high off-diagonal loss indicates ADE risk while low loss indicates cross-protective breadth potential. Diagonal cells (self-recognition) anchor the scale. Geometry-axis input to the Tier A composite. Source: `outputs/table4_dengue/ede_antigenic_loss_matrix.tsv`.

## Figure 6. CDHIT redundancy clustering across vaccine-construct sequences

Bar plot of CDHIT cluster sizes, sorted descending, applied to the input vaccine-construct B-cell sequences. Cluster centroids are propagated to the downstream B-cell (DiscoTope, EpiDope) and T-cell (NetMHCpan) scoring modules; this clustering step removes redundant near-duplicate sequences that would otherwise inflate composite scores. Source: `outputs/dengue_bcell/cdhit_clusters.tsv` (with fallback to `outputs/dengue_smoke_evidence/cdhit_clusters_dengue.tsv`).
