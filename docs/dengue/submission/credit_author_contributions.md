# Author contributions (CRediT)

Author contributions are reported using the CRediT taxonomy (Contributor Roles Taxonomy, 14 categories).

- **Conceptualization.** C.F.E., M.L.N., N.V. Cassia F. Estofolete and Maurício L. Nogueira originated the four-row composite immune correlate framework on which this Perspective builds; Nikos Vasilakis defined the sequence-only proxy question and the Phase 3 retrospective scope.

- **Methodology.** J.W., P.M. James Weatherhead designed the Estofolete Table 4 mapping (Section 2.3), the EDE-restricted scoring window, and the rank-normalised Tier A+B composite; Peter McCaffrey contributed the underlying immunoinformatics pipeline methodology (Versiani, McCaffrey et al. 2026 *Sci Adv*) on which the fork rests.

- **Software.** J.W., P.M. James Weatherhead led all software development for the dengue fork, including the Estofolete Table 4 mapping module, Phase 3 retrospective harness, EDE-restriction flag, stdlib-only patches resolving the Apptainer numpy overlay, and the 5-axis composite-score module. Peter McCaffrey is the original author of the upstream `immunoinformatics_platform` (Nextflow DSL2 with Apptainer/Singularity containerisation).

- **Validation.** J.W., P.M. The Phase 3 retrospective (Section 2.4) and the per-tier ablation plus permutation test (`real(stats)` commit) were designed and implemented by J.W. with pipeline-level validation review by P.M.

- **Formal analysis.** J.W. The Pearson and Spearman analyses, leave-one-out cross-validation, and EDE cross-serotype antigenic-loss matrix (Figure 5) were performed by J.W.

- **Investigation.** J.W., D.W. James Weatherhead performed all pipeline runs and produced the outputs at git tag `v1.0-dengue-results`. Daniela Weiskopf contributed T-cell epitope domain expertise and curated the IEDB-derived dengue CD8 epitope set used in Tier B.2 scoring.

- **Resources.** P.M., D.W., N.V. Peter McCaffrey provided pipeline pedigree (the published `immunoinformatics_platform` and its container builds) and academic-license access (Gurobi/JessEV, NetMHCpan binaries). Daniela Weiskopf provided IEDB-curated dengue CD8 epitope tables and structural references for the Tier B.2 axis. Nikos Vasilakis provided UTMB compute infrastructure and reference proteome curation.

- **Data curation.** J.W. Curation of the four DENV reference polyproteins, the three Phase 3 vaccine parent strains, and the AFND release 2024-12 frequency tables was performed by J.W. (`ref_fastas/` subtree).

- **Writing — original draft.** J.W. James Weatherhead drafted the manuscript including Abstract, Introduction, Methods, Results, Discussion, and the response-to-prior-reviewers document.

- **Writing — review and editing.** C.F.E., M.L.N., D.W., P.M., N.V. All authors reviewed and edited the manuscript; Estofolete and Nogueira reviewed the Table 4 mapping for fidelity to the original composite; Weiskopf reviewed the T-cell sections; McCaffrey reviewed the pipeline-pedigree text; Vasilakis reviewed the framing and the response to prior reviewers.

- **Visualization.** J.W. James Weatherhead produced all five figures (Figures 1-5) and Table 1.

- **Supervision.** N.V. Nikos Vasilakis supervised the work end-to-end at UTMB.

- **Project administration.** N.V. Nikos Vasilakis administered the project, the UTMB-FAMERP-LJI-Pathology-Informatics author coordination, and the response-to-prior-reviewers cycle.

- **Funding acquisition.** N.V. Nikos Vasilakis acquired the funding under which this work was performed at UTMB.

All authors have read and approved the final submitted manuscript. Co-corresponding authors (J.W. and N.V.) accept responsibility for the integrity of the work as a whole.
