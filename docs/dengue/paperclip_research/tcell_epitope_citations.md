# T-cell Epitope Mapping Citations (paperclip-verified)

Generated 2026-05-06 via paperclip CLI. All entries below were retrieved from
PMC, bioRxiv, medRxiv, or OpenAlex abstracts and confirmed by ID. Citations
keyed to: HLA-restricted dengue CD4/CD8 epitopes; polyfunctionality
(IFN-gamma + TNF-alpha + IL-2); cross-serotype reactivity.

## 1. Tian, Grifoni, Sette, Weiskopf 2019 (Frontiers in Immunology)
- ID: PMC6737489 | DOI: 10.3389/fimmu.2019.02125
- Title: "Human T Cell Response to Dengue Virus Infection"
- Measured: review of DENV-specific T cell epitope identification, phenotype,
  and function across all four serotypes; covers IEDB-curated epitopes.
- HLA: HLA-B*35:01, HLA-B*07:02, HLA-A*11:01, HLA-A*02:01, HLA-DRB1 panels.
- Immunogenic DENV proteins: NS3 and NS5 dominate CD8; capsid, NS3, NS5 dominate
  CD4. Highlights HLA-B*35-restricted NS3/NS5 epitopes as immunodominant.

## 2. Weiskopf and Sette 2014 (Frontiers in Immunology)
- ID: PMC3945531 | DOI: 10.3389/fimmu.2014.00093
- Title: "T-Cell Immunity to Infection with Dengue Virus in Humans"
- Measured: synthesis of HLA-class-I-restricted CD8 epitopes mapped in Sri
  Lankan cohort; argues protective rather than pathogenic role for T cells.
- HLA: HLA-A and HLA-B supertypes, with HLA-B*35-restricted responses prominent.
- Immunogenic proteins: NS3, NS4B, NS5 carry the dominant CD8 load; structural
  proteins (E, capsid) more prominent for CD4.

## 3. Adikari, Di Giallonardo, Leung, Grifoni, Sette, Weiskopf 2020 (Sci Rep)
- ID: PMC7687909 | DOI: 10.1038/s41598-020-77565-2
- Title: "Conserved epitopes with high HLA-I population coverage are targets
  of CD8+ T cells associated with high IFN-gamma responses against all dengue
  virus serotypes"
- Measured: IFN-gamma ELISpot mapping of 95 conserved CD8 peptides; identifies
  cross-serotype reactive epitopes with broad population coverage.
- HLA: HLA-A*01, A*02, A*03, A*11, A*24, B*07, B*35, B*44, B*58 supertypes.
- Immunogenic proteins: NS3, NS5, NS4B (CD8 dominance); cross-reactive across
  DENV-1/2/3/4 for the conserved subset.

## 4. Grifoni et al. 2022 (Vaccines, transcriptomics follow-up)
- ID: PMC9029181 | DOI: 10.3390/vaccines10040612
- Title: "Transcriptomics of Acute DENV-Specific CD8+ T Cells Does Not Support
  Qualitative Differences as Drivers of Disease Severity"
- Measured: scRNA-seq of HLA-multimer-sorted DENV-specific CD8 T cells from
  acute Sri Lankan patients across DF and DHF.
- HLA: HLA-B*35:01 multimers used to sort NS3/NS5-specific CD8 cells.
- Immunogenic proteins: NS3 epitope (NS3 526-534) and NS5 epitopes targeted;
  no qualitative transcriptional difference between mild and severe disease.

## 5. Wijeratne et al. 2019 (Clinical and Experimental Immunology, Sri Lanka)
- ID: PMC6842812 | DOI: 10.1111/cei.13363
- Title: "Association of dengue virus-specific polyfunctional T-cell responses
  with clinical disease severity in acute dengue infection"
- Measured: ICS for IFN-gamma + TNF-alpha + IL-2 + CD107a co-staining; SPICE
  polyfunctionality index.
- HLA: not HLA-typed; uses pooled NS3/NS5 megapools.
- Immunogenic proteins: NS3 and NS5 megapools drive polyfunctional CD4 and CD8.
  Polyfunctional responses enriched in DF over DHF, supporting protective T-cell
  hypothesis.

## 6. Tiberi et al. 2025 (Pathogens, megapool validation)
- ID: PMC12845116 | DOI: 10.3390/pathogens15010005
- Title: "Peptide MegaPools Approach to Evaluate the Dengue-Specific CD4 and
  CD8 T-Cell Response"
- Measured: AIM and ICS using DENV1-4 CD4 and CD8 megapools developed by the
  Sette/Weiskopf group; benchmarks vaccinees vs. natural infection.
- HLA: covers Class I and Class II supertypes via megapool design.
- Immunogenic proteins: capsid, NS1, NS3, NS5 contribute; cross-serotype CD8
  reactivity strongest against NS3/NS5 conserved cores.

## Notes on coverage gaps

- The canonical Weiskopf 2013 PNAS (PMID 23690582) and Grifoni 2019 J Virol
  papers were not returned as full-text hits in the PMC corpus, but their
  results are summarized in Tian/Grifoni/Sette/Weiskopf 2019 (PMC6737489) and
  Adikari 2020 (PMC7687909), both of which cite them directly.
- Rivino's Singapore/Vietnam cohort papers did not surface in paperclip; the
  Hamis et al. 2025 Singapore skin-resident CD8 paper (bio_699dca1d70f1,
  doi: 10.1101/2025.07.11.664360) is the closest current Singapore-cohort
  T-cell study indexed and should be considered when updating the manuscript.
- Friberg et al. 2011 (PMC3216538) provides the cross-serotype CD8 expansion
  result and is a useful adjunct citation for the cross-reactivity claim.
