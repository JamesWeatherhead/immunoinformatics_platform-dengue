# Response to prior reviewers

**Prior submission:** THELANCETID-D-26-00582, *The Lancet Infectious Diseases* (rejected).
**Current submission:** *npj Vaccines*, Perspective format.
**Manuscript file:** `outputs/manuscript/dengue_npj_vaccines_draft.md` (this repository).
**Compiled:** 2026-05-06.

We thank the three reviewers who evaluated the prior version of this work at *The Lancet Infectious Diseases*. Their critiques were specific, technically grounded, and have meaningfully improved the resubmission. We address each concern in turn below, with concrete file and line references into the new manuscript and a list of the changes made.

---

## 1. Reviewer 1: misrepresentation of cited literature

**Reviewer 1's concern (paraphrased verbatim from the decision letter):** "Several load-bearing immunological claims in the manuscript cite literature that does not, on inspection, support the specific quantity or mechanism attributed to it. The reader cannot reconstruct the chain of evidence between a stated claim and the cited paper. Each numerical or mechanistic claim must be traceable to a primary measurement of that quantity, not to a review or commentary that gestures at it."

**Our response:** We accept this critique without reservation. In the resubmission, every load-bearing immunological or epidemiological claim is traceable to a single primary measurement, and each citation has been re-verified against the cited paper's Methods and Results sections. A complete citation audit is now committed at `docs/dengue/citation_audit.md` and is referenced from the manuscript at line 261. Three concrete examples:

1. **Phase 3 efficacy point estimates** (Introduction, line 53). The prior draft attributed CYD-TDV, TAK-003, and Butantan-DV efficacies to a secondary review. We now cite the three primary trial reports directly: Sridhar et al. 2018 (PMID 30021080) which directly measured CYD-TDV 5-year efficacy of 56.5% in serostatus-stratified analysis; Tricou et al. 2024 (PMID 38598795) which directly measured TAK-003 5-year efficacy of 61.2%; and Kallas et al. 2024 (PMID 38320283) which directly measured Butantan-DV 2-year efficacy of 79.6%. These are the primary trial publications, not derivative reviews.

2. **Global dengue burden estimate** (Introduction, line 47). The prior draft cited a WHO factsheet that itself cites the burden figure. We now cite Bhatt et al. 2013 *Nature* (PMID 23563266) which directly modelled and measured the 390-million-infection global burden using cartographic and seroprevalence inputs.

3. **EDE epitope topology** (Methods, line 36; Results, line 211). The prior draft cited a review of dengue antibody responses. We now cite Rouvinski et al. 2015 *Nature* (PMID 25569347) which directly resolved the envelope dimer epitope (EDE) by crystal structure and is the primary measurement of the residue set we score against.

Beyond these three, every Phase 3 retrospective number in Section 3.3, every structural reference (e.g. PDB 1OAN for the DENV-2 E dimer at line 143), and every immunoinformatics-tool benchmark (EpiDope, BepiPred-3.1, NetMHCpan-4.1, NetMHCIIpan-4.1, DiscoTope, AlphaFold-Multimer, JessEV) is now cited to its primary methods publication rather than to a downstream application paper. The full mapping lives in `docs/dengue/citation_audit.md` and `docs/dengue/dengue_perspective.bib`.

---

## 2. Reviewer 2: writing register inappropriate for a methods paper

**Reviewer 2's concern (paraphrased verbatim):** "The Methods and Results sections frequently use language more suitable for a position paper or editorial than for a computational methods study. Statements such as 'this approach should transform' or 'these results suggest the field must rethink' are speculative and not supported by the evidence presented. A methods paper must describe what was computed and what was observed, in declarative, falsifiable terms, and reserve interpretive language for a clearly demarcated Discussion."

**Our response:** We agree. The Methods (Section 2, lines 93-172) and Results (Section 3, lines 176-222) of the resubmission have been rewritten in declarative register, with all speculative or hortatory language moved into Section 4 (Discussion) and explicitly labelled. Concrete examples of the language changes:

1. **Methods (Section 2.5, lines 161-172).** The prior draft contained language such as "the pipeline could revolutionise pre-Phase-1 candidate triage." This has been replaced by a numbered list titled "What this pipeline cannot do" (line 161) which states explicitly that Axis 3 (Time), Axis 4 (Intent), and Axis 5 (Effector tone) are out of scope for sequence-only prediction and require longitudinal BCR sequencing, intracellular cytokine staining, and human DC co-culture respectively.

2. **Results (Section 3.3, lines 198-205).** The prior draft framed the retrospective Pearson r as evidence that the pipeline "predicts" efficacy. The resubmission states the actual measured value (r = 0.45, n = 3, with explicit "[CI: n<4]" annotation indicating that a confidence interval is not computable at this sample size) and notes that the model's rank order diverges from the published clinical rank order. No claim of predictive validity is made.

3. **Discussion (Section 4.2, lines 237-249).** A new subsection titled "What this implementation does not establish" has been added, listing the three principal limitations in declarative form: axes 3-5 require trial data; the retrospective is underpowered (n = 3); sequence cannot capture immune compartment dynamics. Speculative reframings of these limitations as future strengths have been removed.

Throughout the manuscript, hedged or aspirational verbs ("may", "could", "should", "promises to") have been replaced where the underlying claim is empirical, and retained only in Section 4.4 (line 270) where they correctly describe a proposed future calibration exercise rather than a present finding.

---

## 3. Reviewer 3: HLA inequity citations did not measure HLA frequencies

**Reviewer 3's concern (paraphrased verbatim):** "The manuscript cites nine references in support of the claim that population coverage of HLA-restricted epitopes is inequitable across dengue-endemic regions. On inspection, none of the nine cited papers actually measure HLA allele frequencies in the populations under discussion; several are commentaries on health equity, and others discuss HLA in unrelated disease contexts. A claim about population HLA inequity must rest on primary measurements of HLA allele frequencies in the relevant populations, run through a transparent population-coverage formula. The current manuscript does neither."

**Our response:** We accept this criticism in full. It was the most consequential of the three reviewer concerns and has driven the most substantial change in the resubmission. Two changes follow.

**Change 1: the nine inappropriate citations have been removed.** They are no longer present in `docs/dengue/dengue_perspective.bib` and do not appear in the current manuscript text. The HLA-equity claim has been reduced from a literature-based assertion to a quantity that the pipeline computes directly.

**Change 2: a new corpus of 16 primary-measurement HLA papers has been assembled** to serve as the citation base for the equity calculation. This corpus is documented at `docs/dengue/paperclip_research/hla_primary_measurements.md` (compiled 2026-05-06 via paperclip CLI search across PMC, bioRxiv, and medRxiv). Every paper in the corpus directly genotypes HLA-A, -B, -C, -DRB1 (and frequently -DQB1, -DPB1, -DRB3/4/5) in N >= 200 individuals from a defined geographic population, uses IPD-IMGT/HLA-compliant nomenclature, and is interoperable with the Allele Frequency Net Database (AFND). The 16 papers are:

1. Kulmann-Leal et al. 2024, Southern Brazil (REDOME bone-marrow donors), PMC11285832.
2. Munoz et al. 2025, Colombia (11,576 blood donors, NGS, high-resolution), PMC12766644.
3. Hernandez-Mejia et al. 2023, Colombia (1,763 stem-cell donors), PMC9869256.
4. Del Angel-Pablo et al. 2020, Mexico (three urban admixed populations), PMC7168288.
5. Que et al. 2022, Vietnam Kinh (3,750 cord-blood units), PMC9277059.
6. Do et al. 2020, Vietnam Kinh (NGS, 4-field resolution), PMC7204072.
7. Satapornpong et al. 2020, Thailand (470 across 4 regional subpopulations), PMC7057685.
8. Puangpetch et al. 2015, Thailand (986 PCR-SSOP HLA-B), PMC4302987.
9. Yuliwulandari et al. 2024, Indonesian Malay (HLA-B 4-field), PMC10909668.
10. Ng et al. 2022, Singapore (multi-ethnic donors), PMC9873421.
11. Solloch et al. 2025, India (130,518 stem-cell donors across 8 subpopulations), PMC11893841.
12. Yadav et al. 2024, North India, PMC12392222.
13. Barton et al. 2022, Tanzanian Maasai (336 individuals, 6 loci), PMC7618281.
14. Rosario et al. 2025, multi-country East and West Africa (>1,000 individuals, 11 HLA loci, Uganda, Kenya, Tanzania, Nigeria), bioRxiv 2025.09.04.673795.
15. Banjoko et al. 2025, Eastern and Southern Africa (high-resolution Class I), PMC12222908.
16. Tshabalala et al. 2022, South Africa (Black, White, Indian, Mixed-ancestry, high-resolution), PMC8931603.

**Regional coverage achieved.** The 16 papers span the five regions reported in Figure 4 of the manuscript: Brazil (rows 1); Latin America and Caribbean (rows 2, 3, 4 covering Colombia and Mexico); Southeast Asia (rows 5-10 covering Vietnam, Thailand, Indonesia, Singapore); South Asia (rows 11-12 covering India); Sub-Saharan Africa (rows 13-16 covering Tanzania, Uganda, Kenya, Nigeria, Eastern and Southern Africa, South Africa). Sri Lanka, Bangladesh, and Pakistan did not yield primary measurement papers meeting the N >= 200, post-2015, PMC criteria; the Indian dataset (Solloch 2025) is the strongest South Asian proxy and stratifies by region.

**The IEDB Population Coverage formula now computed directly.** Section 2.2 of the manuscript (line 129) and Figure 4 (Section 3.2, line 192) compute population coverage from AFND release 2024-12 allele frequencies, populated with the 16 primary measurement papers above where AFND incorporates them, using the IEDB Population Coverage standard formula:

`P_covered(region) = 1 - prod_a (1 - f_a(region))^2`

where `f_a(region)` is the measured allele frequency of MHC allele `a` in the region, and the product runs over all alleles for which the predicted-binders set contains at least one strong binder. The formula assumes Hardy-Weinberg equilibrium at each locus and independence across loci, both of which are the standard IEDB assumptions and are documented as such in the manuscript Methods. The output is `P_covered` per region, as plotted in Figure 4 (line 192). For the Brazilian sub-region in which `P_covered` falls below 0.7 (line 193), the underlying allele-frequency vector traces back through AFND to the Kulmann-Leal 2024 (PMC11285832) REDOME cohort.

This change converts the equity claim from a literature gesture to a reproducible computation grounded in primary measurements, anchored in 16 papers with a combined N exceeding 150,000 typed individuals.

---

We thank the three reviewers again for the rigour of their reading of the prior submission. We believe the resubmission is materially stronger as a result and look forward to further evaluation at *npj Vaccines*.
