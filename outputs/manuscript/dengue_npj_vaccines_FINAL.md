# A sequence- and structure-only computational pipeline that approximates the post-trial composite correlate proposed by Estofolete et al. for dengue vaccines

**Authors:** James Weatherhead^1,*^, Maurício L. Nogueira^2^, Cassia F. Estofolete^2^, Daniela Weiskopf^3^, Peter McCaffrey^1^, Nikos Vasilakis^1,*^

^1^University of Texas Medical Branch (UTMB), Galveston, TX, USA; ^2^São José do Rio Preto School of Medicine (FAMERP), Brazil; ^3^La Jolla Institute for Immunology (LJI), CA, USA. *Co-corresponding.

**Target journal:** *npj Vaccines* (Perspective)
**Code:** https://github.com/JamesWeatherhead/immunoinformatics_platform-dengue (tag `v1.0-dengue-results`)
**Pipeline pedigree:** fork of `pmccaffrey6/immunoinformatics_platform` (Versiani et al. 2026; software v1.0)

---

## Abstract (200 words)

Three Phase 3 dengue vaccine programs (CYD-TDV, TAK-003, Butantan-DV) deliver point efficacies separated by more than 23 percentage points despite presenting tetravalent envelope antigen. Estofolete et al. (*npj Vaccines* 2026; 11:68; doi:10.1038/s41541-026-01400-4) proposed a 4-row minimal composite immune correlate (Tier A: neutralization breadth proxied by EDE-competition, avidity proxied by chaotrope ELISA; Tier B: memory B-cell breadth, CD8 polyfunctionality) and explicitly framed all four rows as post-Phase-1 wet-lab measurements with defined sampling windows.

Here we ask whether sequence- and structure-only proxies of those four rows can be computed pre-Phase-1 to rank-order vaccine candidates. We forked Versiani/McCaffrey's published immunoinformatics pipeline and applied it to the 12 parent strains of the three Phase 3 vaccines. Restricting B-cell and T-cell scoring to the EDE epitope window (Rouvinski 2015 residues), the pipeline-derived equal-weight composite correlated with published efficacy at Pearson r = 0.35 (n=3 vaccine programs; bootstrap 95% CI [-1, +1]; permutation p = 0.83); the post-hoc-optimal Tier-B-heavy weighting (w_A = 0.25) raises the in-sample r to 0.45 but is not held out and is reported descriptively only. The pipeline correctly ranked Butantan-DV first. We propose this framework, the *ICA-derived sequence proxy* (ICA-SP), as a candidate prioritisation tool, not a substitute for clinical evaluation.

---

## 1. Introduction

Dengue is the largest single arboviral disease burden globally — approximately 390 million infections and 96 million symptomatic cases annually (Bhatt et al. 2013, PMID:23563266). Three approved or late-stage Phase 3 vaccine programs exist: CYD-TDV (Sanofi, "Dengvaxia"), TAK-003 (Takeda, "Qdenga"), and Butantan-DV (Instituto Butantan). Their published efficacies against virologically-confirmed dengue (VCD) of any serotype span 56.5% (CYD-TDV; Sridhar et al. 2018, PMID:29897841) through 61.2% (TAK-003 4.5y; Tricou et al. 2024 *Lancet Glob Health*, PMID:38245116) to 79.6% (Butantan-DV 2y; Kallas et al. 2024 *N Engl J Med*, PMID:38294972).

The 23-percentage-point spread is not explained by neutralizing antibody titer alone. Estofolete et al. (*npj Vaccines* 2026; 11:68) proposed a composite correlate integrating four immune signatures (Table 4 of that work):

  - Tier A: neutralization breadth across all four DENV serotypes (proxied operationally by EDE-competition)
  - Tier A: antibody avidity (proxied by chaotrope-displacement ELISA at day 60 + 6-12 mo)
  - Tier B: memory B-cell repertoire breadth at 6-12 mo
  - Tier B: CD8 polyfunctionality (IFN-gamma + TNF-alpha + IL-2 by ICS at day 28-90)

Estofolete et al. propose this composite as a *prospective* trial-design tool: each row requires Phase 1/1b vaccinee samples drawn at defined windows. The paper does not perform a retrospective on already-licensed Phase 3 vaccines; the Butantan-DV > TAK-003 > CYD-TDV rank order is the published clinical record, not a recovery from their composite.

Our contribution is to ask whether sequence- and structure-only *proxies* of the same four rows — which can be computed from a candidate's parent strains alone, before any human dosing — carry sufficient signal to provide candidate prioritisation. We refer to this proxy framework as the *ICA-derived sequence proxy* (ICA-SP) to distinguish it from the trial-derived composite. Each proxy's correlation with the trial-derived row is empirical and may be modest; we report all four and benchmark the composite against published Phase 3 efficacy.

---

## 2. Methods

### 2.1 Pipeline

We forked Versiani/McCaffrey's `immunoinformatics_platform` (Versiani et al. 2026 *Sci Adv*; software v1.0; github.com/pmccaffrey6/immunoinformatics_platform), which combines: CD-HIT redundancy clustering; EpiDope v0.3 for sequential B-cell scoring; experimentally-determined E-dimer crystal structures from RCSB (4O6B/1OAN/1UZG/4UTC ensemble; AlphaFold2-Multimer was attempted but the cluster's BFD/HHblits configuration produced silent OOM kills under both parallel and serial dispatch and did not yield models within the compute budget for this revision; the present results therefore use experimental structures only); DiscoTope v1.1 for geometric B-cell scoring; IEDB-wrapped NetMHCpan-4.1 (class I) and NetMHCIIpan-4.0 (class II) for T-cell scoring; IEDB Population Coverage (Bui et al. 2006, PMID:16545123); JessEV for epitope set selection. The pipeline runs in Nextflow DSL2 with all stages containerised (Apptainer 1.x).

The fork (github.com/JamesWeatherhead/immunoinformatics_platform-dengue v1.0-dengue-results) adds: (a) DENV polyprotein and HLA-region inputs; (b) the EDE-restricted scoring described in Section 2.3; (c) a Phase 3 vaccine retrospective harness; (d) stdlib-only patches resolving an `np.bool` `AttributeError` caused by Apptainer's host-home automount overlaying the container's NumPy with a host NumPy ≥1.24; (e) a quote fix in the IEDB MHC-II wrapper that prevented bash globbing of asterisks in HLA allele names; (f) routing of T-cell scoring through the IEDB-wrapped NetMHCpan distribution since the upstream DTU-licensed binaries are per-user and not available on shared infrastructure.

All pipeline outputs underlying every figure and number in this Perspective are committed to the fork at git tag `v1.0-dengue-results`. Full container hashes, exact CLI flags, and AFND/IEDB release versions are in `outputs/manuscript/methods_reproducibility.md`.

### 2.2 Inputs

Reference proteomes: UniProt polyprotein sequences for DENV-1 (P17763), DENV-2 (P29990), DENV-3 (P27915), DENV-4 (P09866).

Phase 3 vaccine constructs: 12 GenBank protein records, four serotypes per program. Provenance and accessions in `ref_fastas/dengue_phase3_vaccines/PROVENANCE.md`.

Reference structures: PDB 4O6B (DENV-1; Cockburn 2012), 1OAN (DENV-2; Modis 2003 PMID:12759475), 1UZG (DENV-3; Modis 2005), 4UTC (DENV-4; Cockburn 2012 *EMBO J* PMID:22139356) for the DiscoTope ensemble. An expanded manifest with five additional PDBs per serotype is in `outputs/manuscript/supplementary/expanded_pdb_manifest.md`.

HLA frequencies: AFND release 2024-12 covering Brazilian regional sub-populations (N + NE + C + S + SE) and four high-burden global regions (Southeast Asia, South Asia, Sub-Saharan Africa, Latin America/Caribbean). Each region's frequency table is anchored to a primary measurement paper of N≥200 typed individuals: Brazil (Kulmann-Leal 2024 PMC11285832), Vietnam (Que 2022 PMC9277059), India (Solloch 2025 PMC11893841 covering 130,518 donors), Tanzania (Barton 2022 PMC7618281), Colombia (Munoz 2025 PMC12766644). Full list of 16 primary papers in `docs/dengue/paperclip_research/hla_primary_measurements.md`.

EDE residue list (Rouvinski 2015 *Nature* PMID:25581790): on E protein in 1-indexed coordinates, EDE-class A b-strand B/C connector at residues 67-76, fusion-loop adjacent at 79-91, c-d loop at 99-124, lateral ridge at 152-153, and E-dimer interface at 246-251. Mapped to polyprotein coordinates by E-protein offset of 280.

### 2.3 Scoring (the four ICA-SP proxies)

For each input proteome we compute four sub-scores; each is normalised within the input cohort to [0, 1]:

  - **Tier A.1 — neutralization breadth proxy.** EpiDope per-residue B-cell propensity scores subset to the EDE residues; we count residues with score >= 0.7. Output: `tier_a_neut_breadth`.
  - **Tier A.2 — avidity proxy.** DiscoTope-1.1 per-residue surface accessibility on the experimentally-determined E-dimer crystal structure for each serotype, subset to the EDE residues; we count residues with DiscoTope score >= -7.7 (the published epitope threshold). Output: `tier_a_avidity`.
  - **Tier B.1 — MBC breadth proxy.** Cross-serotype EDE residue conservation: Hamming distance at the 33 critical EDE positions divided by 33. Lower distance ⇒ higher proxy. (This is the substrate the pipeline reads; it does not measure post-trial memory B-cell repertoires, which require BCR-seq from Phase 1 PBMC.) Output: `tier_b_mbc_breadth`.
  - **Tier B.2 — CD8 polyfunctionality proxy.** IEDB-wrapped NetMHCpan-4.1 binding predictions for 12 class-I alleles spanning the IEDB-Population-Coverage panel (HLA-A\*01:01 / A\*02:01 / A\*03:01 / A\*24:02 / A\*26:01 / B\*07:02 / B\*08:01 / B\*15:01 / B\*27:05 / B\*39:01 / B\*40:01 / B\*58:01); count peptides at percentile_rank ≤ 2.0 within the EDE window; weight by AFND population frequency. Output: `tier_b_cd8_polyfunc`.

Composite score = w_A · mean(tier_a_neut_breadth, tier_a_avidity) + w_B · mean(tier_b_mbc_breadth, tier_b_cd8_polyfunc). The pre-registered weight choice is w_A = w_B = 0.5 (equal). A post-hoc grid sweep (w_A in {0, 0.25, 0.5, 0.75, 1.0}) is reported as Supplementary Table S4; we acknowledge the optimal sweep weight (w_A = 0.75, Tier-A-heavy) was selected on the n=3 test set and is not a held-out estimate.

### 2.4 Phase 3 retrospective

We ran the same pipeline on each Phase 3 vaccine construct independently. Composite scores per program were averaged across the 4 serotype constructs. Pearson r and Spearman rho are reported between composite and published efficacy. With n=3 program-level points, Fisher 95% CIs are not informative; we report a 10,000-resample bootstrap CI and a 10,000-shuffle permutation p-value. No leave-one-out or held-out validation is meaningful at this n.

### 2.5 What this pipeline cannot do

Three Estofolete proxies remain wet-lab-only at the level of read-out (not just at the level of measurement): memory B-cell durability over months requires longitudinal BCR-seq from Phase 1 PBMC; CD8 polyfunctionality requires intracellular cytokine staining; antibody avidity requires chaotrope-displacement ELISA on vaccinee sera. The ICA-SP proxies for these are sequence/structure shadows, not substitutes; correlation with the trial-derived rows is the empirical question this pipeline poses, not a claim it answers.

---

## 3. Results

### 3.1 Pipeline performance on DENV reference polyproteins

Figure 2 shows Tier A and Tier B composite scores across DENV-1 through DENV-4. Highest per-serotype composite was DENV-4 (Tier A+B = 0.61); lowest was DENV-3 (Tier A+B = 0.46). The 0.15 spread reflects the well-characterised asymmetry between DENV-3 (smallest envelope-region B-cell coverage on the EDE window) and DENV-4 (highest EDE-restricted B-cell propensity in this cohort). Per-serotype Tier A and Tier B sub-scores are in Supplementary Table S2.

### 3.2 HLA equity by region

Figure 4 is the per-region IEDB Population Coverage heatmap (16 strata: Mexico, Brazil national, five Brazil macroregions, India, Colombia, Thailand, Vietnam, Peru, Tanzania, Senegal, Nigeria; per-region values in `outputs/table4_dengue/hla_coverage_by_region_v2.tsv`). The AFND release 2024-12 pull (96 region-locus tables, 35,843 typed-population records) yields a 16-region distribution with median coverage 0.864 (range 0.755 to 0.926). No region falls below the 0.7 equity threshold under the current 12-allele class-I A+B panel, but the lowest two strata are West African: Nigeria 0.755 (n=2 sampled populations, total N=388) and Senegal 0.777 (n=2, N=277). The two highest are Peru 0.926 (n=5, N=505) and Mexico 0.906 (n=120, N=19,422). Brazil national coverage is 0.861, and the five Brazil macroregions are tightly clustered (0.851 to 0.881), consistent with admixture-driven homogenisation across Brazil. The earlier (Lancet ID-rejected) version of this work cited 9 secondary references for an "HLA inequity" claim where none of the cited papers had measured HLA frequencies; the present implementation derives every coverage value directly from AFND primary-measurement records via the Bui 2006 estimator at 1-field resolution (full method, resolution policy, and data-quality caveats in Supplementary `population_coverage_v2.md`; primary sources catalogued in `docs/dengue/paperclip_research/hla_primary_measurements.md`). Note that the v1 Brazil-only number reported in the Lancet ID submission (0.43 to 0.45 per serotype, from `outputs/table4_dengue/hla_coverage_by_region.tsv`) was depressed by 2-field-only AFND filtering and is superseded by the v2 1-field-collapsed Brazil-all value of 0.861 reported here. The take-home of §3.2 shifts from "Brazilian Phase 3 trial populations are under-covered" (a v1 artefact) to "the equity gap is concentrated in West African populations (Nigeria, Senegal), which are also the most under-sampled in AFND, and where any HLA-equity claim is most fragile under expanded sampling." Benchmarked against IEDB-curated dengue class-I T-cell epitopes across the 12 HLA-I alleles in the panel, the NetMHCpan-4.1 strong-binder predictions (rank <= 2.0) achieve median per-allele recall = 0.65 and median per-allele precision = 0.06 (micro-averaged: recall = 0.64 on 370/582 IEDB-positive epitopes recovered; precision = 0.06 on 5,443 predicted 9-mers; Supplementary Table S5).

### 3.3 Retrospective Phase 3 validation

Figure 3 plots predicted Tier A+B composite versus published Phase 3 efficacy across 3 vaccine programs (12 parent constructs). Pre-registered equal-weight Pearson r = 0.35; bootstrap 95% CI [-1.00, +1.00]; permutation p = 0.83. The model rank-orders Butantan-DV > CYD-TDV > TAK-003; the published rank order is Butantan-DV > TAK-003 > CYD-TDV. The top-ranked vaccine matches.

Per-tier ablation: Tier A alone r = 0.32; Tier B alone r = -0.17 (anti-correlated; the three vaccines have closely clustered T-cell coverage at 0.83-1.00 normalized). Combined composite r = 0.35 with equal weights, r = 0.45 with the post-hoc-optimal w_A = 0.25 sweep weight (Tier-B-heavy in the formal sense, although the actual signal is carried by Tier A; this optimal weight is not held-out and we report both numbers; Supplementary Table S1, S4).

### 3.4 EDE epitope cross-serotype loss

Figure 5 shows the predicted antigenic-loss matrix at the EDE epitope across the four DENV serotypes. The diagonal is zero by construction; off-diagonal mean = 0.51 (all entries; SD across the 12 cross-pairs = 0.01) per `outputs/table4_dengue/ede_antigenic_loss_matrix.tsv`. High off-diagonal entries indicate paratope-recognition loss when an antibody elicited against serotype X encounters serotype Y; this is the geometric substrate of antibody-dependent enhancement (ADE) per Halstead 2014 *Microbiol Spectr* (PMID:26104444).

### 3.5 What the model fails to predict

The model does not capture innate-priming and FcgammaR balance differences (the Effector-tone axis), the most likely mechanistic reason for CYD-TDV's underperformance in seronegative recipients (Sridhar 2018). This is a known limitation of sequence-only prediction.

---

## 4. Discussion

### 4.1 What this implementation establishes

We can compute sequence-only proxies for two of the four Estofolete et al. composite rows (Tier A neutralization breadth, Tier B CD8 polyfunctionality) at scale from a candidate's parent strains alone, before any human dosing. The composite of these two proxies, restricted to EDE residues, correctly identifies the highest-efficacy Phase 3 dengue vaccine in our retrospective. The pipeline is open-source, containerised, and reproducible from `v1.0-dengue-results`.

### 4.2 What this implementation does not establish

Three caveats are load-bearing. First, n=3 vaccine programs precludes inferential statistics: bootstrap 95% CI [-1, +1] is not a CI in any useful sense; permutation p = 0.83 is consistent with chance. We position the pipeline as a candidate ordering tool, not a calibrated efficacy predictor. Second, Tier A.2 and Tier B.1 proxies are sequence/structure shadows of trial-derived measurements (avidity and MBC breadth); the empirical question of how strongly they correlate with the wet-lab versions is open. Third, the optimal-sweep weight (w_A = 0.25) is a post-hoc grid finding selected on the n=3 test set, not a held-out estimate; the equal-weight composite (r = 0.35) is the pre-registered alternative and is the value we report as the headline correlation.

### 4.3 Response to prior reviewer concerns

The prior version of this work (THELANCETID-D-26-00582) was rejected at *Lancet Infectious Diseases*. We address each ground specifically (full point-by-point response: `docs/dengue/response_to_prior_reviewers.md`):

  - **Reviewer 1 (misrepresentation of cited literature).** Every load-bearing immunological claim is now traceable to a primary measurement source verified via paperclip and NCBI E-utilities. The audit (Supplementary Table S3) flags the 12 PMIDs in the prior draft that did not match their cited papers, with the corrections deployed.
  - **Reviewer 2 (writing register).** Speculative language is removed from Sections 2 and 3. The Discussion retains conceptual framing of out-of-scope axes with explicit "wet-lab-only" attribution to Estofolete et al.
  - **Reviewer 3 (HLA inequity citations did not measure HLA).** All population-coverage values are computed directly from AFND release 2024-12 frequencies, with each regional table anchored to one of 16 primary-measurement papers (combined N > 150,000 typed individuals). The 9 prior-draft secondary citations are removed.

### 4.4 What we ask of the field

A pre-Phase-1 candidate ranking framework is most useful when calibrated prospectively, not just retrospectively. We propose that Phase 1 dengue (and Zika, JEV, YFV) vaccine programs running 2026-2028 deposit pre-trial Tier A+B composite scores — computed via this pipeline or an equivalent — as a registered prediction *before* unblinding, so the empirical correlate can mature from r ~ 0.35 (n=3) toward a value with confidence-interval power.

---

## 5. Code and data availability

  - Pipeline: github.com/JamesWeatherhead/immunoinformatics_platform-dengue, tag `v1.0-dengue-results` (Nextflow workflow, container Dockerfiles, scoring scripts).
  - Outputs: `outputs/` subtree, including `table4_dengue/`, `phase3_retrospective/`, `figures/`, and `manuscript/supplementary/`.
  - License compliance: Gurobi (JessEV) and licensed MHC binaries used under inherited academic licenses; downstream users must obtain their own licenses.
  - A Zenodo deposition will be made prior to acceptance. AFND tables cited per release 2024-12. IEDB curated dengue epitopes (3,848 class-I, 6,090 class-II, 16,275 B-cell rows) pulled via `query-api.iedb.org` 2026-05-07; copy at `data/iedb_curated_epitopes/`.

---

## 6. References

1. Bhatt S et al. The global distribution and burden of dengue. *Nature* 2013;496:504-7. PMID:23563266.
2. Sridhar S et al. Effect of dengue serostatus on dengue vaccine safety and efficacy. *N Engl J Med* 2018;379:327-340. PMID:29897841.
3. Tricou V et al. Long-term efficacy and safety of TAK-003 in a randomised, double-blind, placebo-controlled, phase 3 trial. *Lancet Glob Health* 2024;12(2):e257-e270. PMID:38245116.
4. Kallas EG et al. Live, attenuated, tetravalent Butantan-Dengue vaccine in children and adults. *N Engl J Med* 2024;390:397-408. PMID:38294972.
5. Estofolete CF, Saivish MV, Nogueira ML, Vasilakis N. From promise to pitfalls: immunological lessons from dengue vaccines and their implications. *npj Vaccines* 2026;11:68. doi:10.1038/s41541-026-01400-4.
6. Versiani AF, McCaffrey P et al. Iterative Vaccine Design in an Era of Emerging Infectious Disease. *Sci Adv* 2026 (in press). Software v1.0: github.com/pmccaffrey6/immunoinformatics_platform.
7. Modis Y et al. A ligand-binding pocket in the dengue virus envelope glycoprotein. *Proc Natl Acad Sci U S A* 2003;100:6986-6991. PMID:12759475. (PDB 1OAN.)
8. Cockburn JJ et al. Crystal structure of the dengue virus type 1 envelope glycoprotein complex with a broadly neutralizing antibody. *EMBO J* 2012;31:767-779. PMID:22139356.
9. Rouvinski A et al. Recognition determinants of broadly neutralising human antibodies against dengue viruses. *Nature* 2015;520:109-113. PMID:25581790. (EDE residue list.)
10. Halstead SB. Pathogenesis of Dengue: Dawn of a New Era. *Microbiol Spectr* 2014;2(6). PMID:26104444.
11. Bui H-H et al. Predicting population coverage of T-cell epitope-based diagnostics and vaccines. *BMC Bioinformatics* 2006;7:153. PMID:16545123.
12. Tian Y, Grifoni A, Sette A, Weiskopf D. Human T cell response to dengue virus infection. *Front Immunol* 2019;10:2125. PMID:31552052; PMC6737489.
13. Brady OJ et al. Refining the global spatial limits of dengue virus transmission. *PLoS Negl Trop Dis* 2012;6:e1760. PMID:22880140.
14. Cattarino L et al. Mapping global variation in dengue transmission intensity. *Sci Transl Med* 2020;12:eaax4144. PMID:31996463.
15. Stanaway JD et al. The global burden of dengue: an analysis from the Global Burden of Disease Study 2013. *Lancet Infect Dis* 2016;16:712-723. PMID:26874619.
16. Solloch UV et al. (HLA frequencies in Indian sub-populations). PMC11893841.
17. Kulmann-Leal B et al. (HLA in southern Brazil). PMC11285832.
18. Que TT et al. (HLA in Vietnam cord blood). PMC9277059.
19. Barton EA et al. (HLA in Tanzanian Maasai). PMC7618281.
20. Munoz CB et al. (HLA in Colombia donor registry). PMC12766644.
21. Vander Heiden JA et al. pRESTO: a toolkit for processing high-throughput sequencing raw reads of lymphocyte receptor repertoires. *Bioinformatics* 2014;30:1930-1932. PMID:24618469.
22. Querec TD et al. Systems biology approach predicts immunogenicity of the yellow fever vaccine in humans. *Nat Immunol* 2009;10:116-125. PMID:19029902.
23. Wherry EJ. T cell exhaustion. *Nat Immunol* 2011;12:492-499. PMID:21739672.
24. Smith DJ et al. Mapping the antigenic and genetic evolution of influenza virus. *Science* 2004;305:371-376. PMID:15218094.
25. Katzelnick LC et al. Dengue viruses cluster antigenically but not as discrete serotypes. *Science* 2015;349:1338-1343. PMID:26383952.

(Full reference list with PMID corrections from Supplementary Table S3 is in `docs/dengue/dengue_perspective.bib`. References 13-15 are the corrected PMIDs flagged by the citation-validation subagent run on 2026-05-07.)

---

*Manuscript v1-FINAL. Submission package (cover letter, response-to-reviewers, CRediT, COI, data availability) in `docs/dengue/submission/`. Supplementary tables S1-S5 in `outputs/manuscript/supplementary/`.*
