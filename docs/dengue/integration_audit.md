# Integration / consistency audit of dengue_npj_vaccines_FINAL.md

**Audit date:** 2026-05-07
**Audited file:** `outputs/manuscript/dengue_npj_vaccines_FINAL.md`
**Cross-checked against:** `methods_reproducibility.md`, supplementary S1-S5, `phase3_retrospective/statistical_analyses.json`, `phase3_retrospective/phase3_vs_clinical.tsv`, `table4_dengue/ranked_candidates.tsv`, `docs/dengue/response_to_prior_reviewers.md`, `docs/dengue/submission/cover_letter.md`, `docs/dengue/paperclip_research/*.md`.

Findings are listed in numbered order; severity is one of {critical, major, minor}.

---

## CRITICAL (number/claim disagreements between FINAL.md and canonical sources)

### 1. Pre-registered weights in FINAL conflict with methods_reproducibility.md and supplementary S1
- **Where:** `dengue_npj_vaccines_FINAL.md` Abstract (line 17), §2.3 (line 69), §3.3 (line 95), §4.2 (line 115)
- **What:** FINAL repeatedly states the pre-registered weights are `w_A = 0.5, w_B = 0.5` (equal) and that the post-hoc-optimal sweep weight is `w_A = 0.25` (i.e. Tier B heavy).
- **Conflict A (methods_reproducibility):** §2.6.2 of `methods_reproducibility.md` says the headline `r = 0.45` reported in §3.3 of the main text uses the **pre-registered 0.5/0.5 weights**, but the canonical JSON `statistical_analyses.json` reports `pearson_r_composite_avg = 0.3516` for the equal-weight composite (which in the JSON corresponds to the headline `composite_avg`). The 0.45 value does not appear in any artefact under `outputs/`. The actual headline value reproducible from the canonical JSON is **0.35**, not 0.45.
- **Conflict B (supplementary S1 / S4):** S1 (`S1_ablation_studies.tsv` row 9, `EDE_only / tier_A_heavy_w0.75`) and S4 (`weight_w_A=0.75`, `pearson_r=0.352`) both identify the post-hoc-optimal weight pair as **`w_A = 0.75, w_B = 0.25`** (Tier-A-heavy), with composite `r = 0.352`. FINAL.md inverts this and reports the optimum as **`w_A = 0.25, w_B = 0.75`** (Tier-B-heavy). This is contradicted by S1 row 9 caption ("Tier-A-heavy scheme [w_A = 0.75, w_B = 0.25] maximises composite r at 0.352"), by S4 weight_w_A sweep (peak at 0.75), and by S2 caption ("Tier-A-heavy [w_A = 0.75, w_B = 0.25] weighted sum used for the headline result").
- **Conflict C (methods_reproducibility §2.6.2 self-contradiction):** §2.6.2 names the optimum as `w_A = 0.25, w_B = 0.75` (matching FINAL.md) but contradicts the Tier-A vs Tier-B partial r values reported in `statistical_analyses.json` (Tier A r=+0.32, Tier B r=-0.17). Because Tier B is anti-correlated, the optimal weight should put **more** weight on Tier A, not less. The 0.25/0.75 assignment is internally impossible against the JSON.
- **Severity:** CRITICAL.
- **Suggested fix:** Reconcile to **w_A = 0.75, w_B = 0.25** as the post-hoc-optimum (per S1, S4, JSON), and either update headline `r` to **0.35** (equal-weight, pre-registered) or to **0.352** (post-hoc-optimum, w_A=0.75). The 0.45 value should be dropped or sourced from a stored artefact. (Applied below: change FINAL.md to consistent with the canonical artefacts.)

### 2. Headline Pearson r = 0.45 has no canonical source
- **Where:** Abstract line 17; §3.3 line 93-95; §4.2 line 115; §4.4 line 127; cover_letter.md.
- **What:** FINAL.md claims `Pearson r = 0.45` for the EDE-restricted composite. The canonical `statistical_analyses.json` reports `pearson_r_composite_avg = 0.3516`. S1 row 9 reports composite `r = 0.352` for the highest-correlation EDE-only weighting. The 0.45 number does not appear in any TSV/JSON under `outputs/`.
- **Severity:** CRITICAL.
- **Suggested fix:** Replace 0.45 with the canonical `0.35` (equal-weight pre-registered) value throughout, or reference a held-out artefact that justifies 0.45. Per commit `454166d` ("optimal-weight composite (0.25 Tier A + 0.75 Tier B) yields Pearson r = 0.45"), the 0.45 value was claimed as the post-hoc-optimal sweep, but this is in conflict with S1/S2/S4 which name the optimum at w_A=0.75.

### 3. Tier A/Tier B partial-r values inverted relative to JSON
- **Where:** §3.3 line 95.
- **What:** FINAL says "Tier A alone r = 0.32; Tier B alone r = -0.17". The JSON gives `pearson_r_tier_A_only = 0.324` and `pearson_r_tier_B_only = -0.173`. These match (modulo rounding). However, FINAL then says "Tier B (the three vaccines have nearly identical T-cell coverage at 0.97-1.00 normalized)". The JSON shows `tier_B_normalized = [0.829, 1.000, 0.857]` (CYD-TDV, TAK-003, Butantan-DV) -- a spread of 0.829 to 1.0, not 0.97 to 1.00.
- **Severity:** CRITICAL (the "0.97-1.00" claim is fabricated relative to the canonical JSON).
- **Suggested fix:** Change "0.97-1.00 normalized" to "0.83-1.00 normalized" or simply drop the parenthetical.

### 4. Tier A normalization values conflict with canonical JSON
- **Where:** §3.1 line 85.
- **What:** FINAL says "Highest composite was DENV-1 (Tier A+B = 0.95); lowest was DENV-2 (Tier A+B = 0.87). The 0.08 spread...". The canonical JSON gives `composite_avg = [0.914, 0.500, 0.929]` for [CYD-TDV, TAK-003, Butantan-DV] (program-level, not per-serotype). Per-serotype values from `ranked_candidates.tsv` give serotype Tier A scores of 0.484, 0.438, 0.391, 0.560 for DENV-1 to DENV-4, and Tier B of 0.09, 0.07, 0.07, 0.05 -- so Tier A+B sums of 0.574, 0.508, 0.461, 0.610. Highest is DENV-4 (not DENV-1) and lowest is DENV-3 (not DENV-2), with spread 0.149 (not 0.08). Even normalised to [0,1], no calculation produces 0.95/0.87.
- **Severity:** CRITICAL.
- **Suggested fix:** Replace per-serotype claim with the actual TSV values: highest DENV-4 (~0.61), lowest DENV-3 (~0.46), spread ~0.15. Or alternatively rephrase as program-level (CYD-TDV 0.872, TAK-003 0.750, Butantan-DV 0.893) but then it is not per-serotype. The earlier draft (`dengue_npj_vaccines_draft.md` line 181) had this correct: "Highest composite was DENV-4 (Tier A+B = 0.31); lowest was DENV-3 (Tier A+B = 0.23)".

### 5. EDE off-diagonal mean inconsistent with TSV
- **Where:** §3.4 line 99.
- **What:** FINAL says "off-diagonal mean = 0.40 (all entries; SD across the 12 cross-pairs = 0.05)". The canonical `ede_antigenic_loss_matrix.tsv` reports off-diagonal values 0.508, 0.513, 0.496, 0.508, 0.517, 0.500, 0.513, 0.517, 0.505, 0.496, 0.500, 0.505 -- mean ≈ **0.506**, SD ≈ **0.0078**. The earlier draft correctly reported "0.51 ± 0.01" (line 211).
- **Severity:** CRITICAL.
- **Suggested fix:** Change "off-diagonal mean = 0.40 ... SD = 0.05" to "off-diagonal mean = 0.51 (SD = 0.01)".

### 6. HLA equity claim contradicts the underlying TSV
- **Where:** §3.2 line 89.
- **What:** FINAL says "Coverage falls below 0.7 in 0 of 5 regions when the AFND-derived frequencies are computed directly". The canonical `hla_coverage_by_region.tsv` shows Brazil coverage of 0.454, 0.442, 0.442, 0.430 across DENV-1 to DENV-4 -- ALL below 0.7. The earlier draft reported "Coverage falls below 0.7 in 1 of the 5 regions (Brazil)" (line 192-193) which is more consistent with the available TSV (although the TSV only has Brazil rows; other regions are missing).
- **Severity:** CRITICAL.
- **Suggested fix:** Either (a) note that the TSV currently contains only Brazil rows (and 4/4 are below 0.7) so the "0 of 5 regions" claim has no canonical basis, or (b) restate as "Brazil coverage falls below 0.7 across all four serotypes".

---

## MAJOR (cross-section drift, missing files, misalignments)

### 7. Cover letter title vs FINAL.md title: MATCH
- **Where:** `cover_letter.md` line 10; `dengue_npj_vaccines_FINAL.md` line 1.
- **What:** Both use the title "A sequence- and structure-only computational pipeline that approximates the post-trial composite correlate proposed by Estofolete et al. for dengue vaccines." ✓ consistent.
- **Severity:** none (informational).

### 8. Reference to expanded_pdb_manifest.md and ede_residue_provenance.md not yet committed
- **Where:** §2.2 line 54 ("expanded manifest with five additional PDBs per serotype is in `outputs/manuscript/supplementary/expanded_pdb_manifest.md`"); methods_reproducibility.md GAP at line 175 ("docs/dengue/ede_residue_provenance.md before journal submission").
- **What:** Neither file exists in `outputs/manuscript/supplementary/` (only S1-S5 exist).
- **Severity:** MAJOR.
- **Suggested fix:** Either remove the cross-reference (which currently broken-links from FINAL.md), or commit the file before submission.

### 9. Reference to `data/iedb_curated_epitopes/` (not in repo as cited)
- **Where:** §5 line 136.
- **What:** "IEDB curated dengue epitopes (3,848 class-I, 6,090 class-II, 16,275 B-cell rows) pulled via `query-api.iedb.org` 2026-05-07; copy at `data/iedb_curated_epitopes/`."
- **Severity:** MAJOR (need to confirm path exists; commit `a3836ae` says "26,213 curated dengue epitope assays for cross-validation"; total = 3848+6090+16275 = 26,213 which matches commit message).
- **Suggested fix:** Verify the actual path matches the citation; if the corpus is at a different relative path, update the manuscript.

### 10. Estofolete reframing: §3.4 still references "the geometric substrate of antibody-dependent enhancement"
- **Where:** §3.4 line 99 cites "Halstead 2014 *Microbiol Spectr* (PMID:26104444)"; this PMID corresponds to "Halstead 2017 Annu Rev Virol" in the `structural_ade_citations.md` corrections table. The 2014 Microbiol Spectr is the correct citation; the citation in FINAL.md uses the correct PMID 26104444 — but reference list entry #10 is also "Halstead SB. *Microbiol Spectr* 2014;2(6). PMID:26104444." Consistent.
- **Severity:** MINOR — citation OK.

### 11. Estofolete reframing: Abstract and §2 discuss "all four rows as post-Phase-1 wet-lab measurements"; consistent with paperclip_research/estofolete_2026_table4.md
- **Where:** Abstract line 15, §1 line 32-34, §2.5 line 77, §4.2 line 115.
- **What:** FINAL.md correctly does NOT use the "Estofolete proposed 5 axes" framing or the "Estofolete recovered the rank ordering" language flagged in the paperclip research as drifts. It also correctly attributes the Butantan-DV > TAK-003 > CYD-TDV ordering to "the published clinical record, not a recovery from their composite" (line 32). ✓ consistent.
- **Severity:** none.
- **Note:** The earlier `dengue_npj_vaccines_draft.md` §4.1 (line 230) DID contain the drifted phrasing ("These two axes (Geometry, Equity) recover the Butantan-DV > TAK-003 > CYD-TDV rank order"). FINAL.md correctly removes this. ✓

### 12. n=3 caveat: present in Abstract, §3.3, §4.2 — consistent
- **Where:** Abstract line 17 ("n=3 vaccine programs; bootstrap 95% CI [-1, +1]; permutation p = 0.83"); §3.3 line 93 ("Pearson r = 0.45; bootstrap 95% CI [-1.00, +1.00]; permutation p = 0.83"); §4.2 line 115 ("n=3 vaccine programs precludes inferential statistics: bootstrap 95% CI [-1, +1] is not a CI in any useful sense; permutation p = 0.83 is consistent with chance").
- **Severity:** none — n=3 caveat is consistently present.
- **Note:** §4.4 also references "r ~ 0.45 (n=3)" (line 127). All four locations consistent.

### 13. Permutation p = 0.83 vs JSON value (no permutation_p_composite reported)
- **Where:** Abstract, §3.3, §4.2 all report p = 0.83.
- **What:** The canonical JSON only reports `permutation_p_tier_A = 1.0`. There is no `permutation_p_composite` field in the JSON, so the 0.83 value cannot be cross-checked against the canonical artefact.
- **Severity:** MAJOR (un-traceable but plausible).
- **Suggested fix:** Compute and store `permutation_p_composite` in `statistical_analyses.json`, or annotate FINAL.md to reference the script that produced 0.83.

### 14. n_constructs = 12 in JSON vs 12 GenBank protein records in §2.2
- **Where:** §2.2 line 52 ("12 GenBank protein records, four serotypes per program") and §3.3 line 93 ("3 vaccine programs (12 parent constructs)").
- **What:** Consistent. Each vaccine has 4 serotype components × 3 programs = 12. ✓
- **Severity:** none.

### 15. Reference list orphans (PMIDs in references not cited in body)
- **Where:** §6 references list lines 142-167.
- **What:** Audit each reference for body-text mention:
  - #1 Bhatt 2013 PMID 23563266 → cited in §1 line 23 ✓
  - #2 Sridhar 2018 PMID 29897841 → cited in §1 line 23, §3.5 line 103 ✓
  - #3 Tricou 2024 PMID 38245116 → cited in §1 line 23 ✓
  - #4 Kallas 2024 PMID 38294972 → cited in §1 line 23 ✓
  - #5 Estofolete 2026 → cited throughout ✓
  - #6 Versiani McCaffrey 2026 → cited in §2.1 line 42 ✓
  - #7 Modis 2003 PMID 12759475 → cited in §2.2 line 54 ✓
  - #8 Cockburn 2012 PMID 22139356 → cited in §2.2 line 54 ✓
  - #9 Rouvinski 2015 PMID 25581790 → cited in §2.2 line 58, §2.3 line 64 ✓
  - #10 Halstead 2014 PMID 26104444 → cited in §3.4 line 99 ✓
  - #11 Bui 2006 PMID 16545123 → cited in §2.1 line 42 ✓
  - **#12 Tian/Grifoni/Sette/Weiskopf 2019 PMID 31572359 → ORPHAN (not cited in body text)**
  - **#13 Brady 2012 PMID 22880140 → ORPHAN (not cited in body)**
  - **#14 Cattarino 2020 PMID 31996463 → ORPHAN**
  - **#15 Stanaway 2016 PMID 26874619 → ORPHAN**
  - #16-#20 HLA papers → cited collectively in §2.2 line 56 (Kulmann-Leal, Que, Solloch, Barton, Munoz appear by name); ✓
  - **#21 Vander Heiden 2014 PMID 24618469 (pRESTO) → ORPHAN (not cited in body)**
  - **#22 Querec 2009 PMID 19029902 → ORPHAN**
  - **#23 Wherry 2011 PMID 21739672 (T cell exhaustion) → ORPHAN**
  - **#24 Smith 2004 PMID 15218094 (antigenic cartography) → ORPHAN**
  - **#25 Katzelnick 2015 PMID 26383952 (dengue antigenic cartography) → ORPHAN**
- **Severity:** MAJOR. Eight orphan citations (refs 12, 13, 14, 15, 21, 22, 23, 24, 25) appear in the references list but are not cited in the body text.
- **Suggested fix:** Either remove the orphans from §6, or add body-text citations referring to them (e.g. §1 could cite Brady/Cattarino/Stanaway as additional burden refs; §4.2 could cite Wherry/Querec for T-cell exhaustion mechanism).

### 16. response_to_prior_reviewers.md vs FINAL.md §4.3 cross-walk
- **Where:** `response_to_prior_reviewers.md` and FINAL.md §4.3 line 117-123.
- **What:**
  - R1 ("misrepresentation of cited literature"): both documents cover this. Response document at lines 12-26 lists three concrete examples (Phase 3 efficacy attributions, Bhatt 2013, Rouvinski 2015) and references S3. FINAL §4.3 R1 states "Every load-bearing immunological claim is now traceable to a primary measurement source verified via paperclip and NCBI E-utilities. The audit (Supplementary Table S3) flags the 12 PMIDs in the prior draft that did not match their cited papers". S3 actually flags **multiple PMIDs** (4 are flagged `n` because they cite a review, plus 3 are flagged `n` because of indexing limitations, plus 3 are `NA`); the precise count of "12 PMIDs in the prior draft that did not match" is NOT supported by S3, which lists 23 total claims with mixed verification states.
  - **Inconsistency:** the "12 PMIDs" claim in FINAL §4.3 R1 line 121 has no canonical basis in S3. S3 has 23 rows of which 4 are flagged `n` because they cite a review, and 3 because of indexing. The total of "PMIDs that did not match the cited papers" in S3 is closer to 7, not 12.
  - The response_to_prior_reviewers.md cites different PMIDs than the FINAL.md/S3 (e.g. response document cites Sridhar 2018 PMID 30021080; FINAL and S3 cite Sridhar 2018 PMID 29897841). The PMID 29897841 is correct per `paperclip_research/phase3_efficacy_citations.md`, and 30021080 is wrong.
- **Severity:** MAJOR (number "12" is not traceable; PMIDs in response document are wrong).
- **Suggested fix:** (a) Change "12 PMIDs" to "the PMIDs in the prior draft" or "the citation drifts" to avoid an unsupported count. (b) Update response_to_prior_reviewers.md to use the corrected PMIDs (29897841, 38245116, 38294972) consistent with FINAL.md, S3, and `phase3_efficacy_citations.md`. (c) The response doc also cites Bhatt 2013 PMID 23563266 (correct) and "Rouvinski 2015 PMID 25569347" (WRONG -- correct is 25581790 per `structural_ade_citations.md` and FINAL §6 ref #9).

### 17. cover_letter.md "n = 3, with the explicit caveat that a confidence interval is not computable" — consistent with FINAL.md §3.3
- **Where:** cover_letter line 12.
- **What:** Cover letter says "Pearson r = 0.45 (n = 3, with the explicit caveat that a confidence interval is not computable at this sample size)". FINAL.md §3.3 says "bootstrap 95% CI [-1.00, +1.00]; permutation p = 0.83". These are not contradictory but the framings differ; cover letter is more cautious, FINAL is more technical. ✓
- **Severity:** minor; framings consistent.

### 18. response_to_prior_reviewers.md references `outputs/manuscript/dengue_npj_vaccines_draft.md` not the FINAL
- **Where:** response_to_prior_reviewers.md line 5 ("Manuscript file: `outputs/manuscript/dengue_npj_vaccines_draft.md`")
- **What:** The response document was last touched 2026-05-07 09:57 but still points to the draft (.md) rather than FINAL. After FINAL was created, response_to_prior_reviewers.md needs to update its cross-reference.
- **Severity:** MAJOR.
- **Suggested fix:** Change line 5 to point to `dengue_npj_vaccines_FINAL.md`.

### 19. Pipeline figures reference Figure 1 but body text never cites Figure 1
- **Where:** No `Figure 1` reference in FINAL.md body (search for "Figure" produces only Figure 2, 3, 4, 5).
- **What:** `outputs/figures/fig1_pipeline_overview.png` and `fig1_pipeline_overview.svg` exist but are not referenced in the manuscript body. CRediT (`credit_author_contributions.md` line 25) credits J.W. with "all five figures (Figures 1-5)".
- **Severity:** MAJOR.
- **Suggested fix:** Either add a Figure 1 reference (e.g. "Figure 1 shows the pipeline overview" near §2.1) or remove fig1 from the figure file set.

### 20. Pipeline figures fig6_cdhit_clusters not referenced anywhere
- **Where:** `outputs/figures/fig6_cdhit_clusters.png` and `.svg` exist but no Figure 6 in the manuscript.
- **Severity:** MINOR (orphan figure file).
- **Suggested fix:** Either remove from outputs or note as Supplementary Figure.

---

## MINOR (formatting / nomenclature drift)

### 21. CRediT references "5-axis composite-score module"; FINAL.md uses "4-row composite" / "ICA-SP" framework
- **Where:** `credit_author_contributions.md` line 9 ("5-axis composite-score module"); FINAL.md uses "four rows" or "ICA-SP" throughout.
- **What:** CRediT carries forward the 5-axis (Geometry / Equity / Time / Intent / Effector tone) framing that paperclip_research flagged as a drift. FINAL.md correctly drops it; CRediT does not.
- **Severity:** MINOR.
- **Suggested fix:** Update CRediT to "Estofolete Table 4 mapping module" or "4-row composite-score module".

### 22. methods_reproducibility.md cites Versiani as "Sci. Adv. 12:eaeb2066" but FINAL.md does not include the eLocator
- **Where:** `methods_reproducibility.md` line 17 ("Versiani et al. 2026 *Sci Adv* 12:eaeb2066"); FINAL.md ref #6 ("Versiani AF, McCaffrey P et al. ... Sci Adv 2026 (in press). Software v1.0: github.com/pmccaffrey6/immunoinformatics_platform").
- **What:** FINAL.md ref #6 says "in press" but methods_reproducibility.md (and the cover letter) cite the eLocator `12:eaeb2066`.
- **Severity:** MINOR.
- **Suggested fix:** Either harmonise FINAL ref #6 to include "12:eaeb2066" OR change methods_reproducibility / cover letter to "in press".

### 23. methods_reproducibility.md §2.5 cites Rouvinski PMID 25849258; FINAL ref #9 cites PMID 25581790
- **Where:** `methods_reproducibility.md` line 149 ("Rouvinski et al. 2015 (Nature 520:109-13, PMID 25849258)"); FINAL.md ref #9 ("Rouvinski A et al. ... PMID:25581790").
- **What:** Per `structural_ade_citations.md` (the citation-validation source), the correct Rouvinski 2015 PMID is **25581790** (not 25849258). methods_reproducibility.md uses the WRONG PMID. The error originated in the prior brief and was corrected by the citation audit.
- **Severity:** MINOR (auxiliary doc has the wrong PMID).
- **Suggested fix:** Update methods_reproducibility.md §2.5 PMID to 25581790.

### 24. FINAL.md title and methods_reproducibility.md cite different paper titles for the upstream
- **Where:** FINAL.md ref #6 says Versiani 2026 "Iterative Vaccine Design in an Era of Emerging Infectious Disease" is implied but the FINAL ref does NOT name the title; methods_reproducibility.md does not either; only `versiani_2026_pipeline.md` provides the title.
- **Severity:** MINOR.
- **Suggested fix:** Add the upstream title to FINAL ref #6 for completeness.

### 25. EDE residue list in FINAL §2.2 (line 58) vs methods_reproducibility.md §2.5 (line 152-160)
- **Where:** FINAL.md §2.2 (line 58) lists EDE residues as "67-76, 79-91, 99-124, 152-153, 246-251" (range 79-91 is 13 residues; range 99-124 is 26 residues; total = 10 + 13 + 26 + 2 + 6 = 57 residues).
- **What:** methods_reproducibility.md §2.5 (line 152-160) lists the canonical residue set as `range(67, 77) + range(99, 125) + [79, 80, 84, 86, 88, 91, 152, 153] + [246, 248, 249, 251]` = 10 + 26 + 8 + 4 = 48 residues. The "33 critical EDE positions" referenced in FINAL §2.3 line 66 ("Hamming distance at the 33 critical EDE positions") matches neither count (48 nor 57).
- **Severity:** MINOR (residue specifications and counts disagree across docs).
- **Suggested fix:** Reconcile residue counts. Either FINAL §2.2 should match the verbatim list in methods §2.5 (48 residues), or the "33 critical EDE positions" claim in FINAL §2.3 should be re-counted.

### 26. methods_reproducibility.md Versiani Sci Adv eLocator
- **Where:** methods_reproducibility.md line 17 cites "Versiani et al. 2026 *Sci Adv* 12:eaeb2066"; `versiani_2026_pipeline.md` notes the eLocator could not be programmatically verified and the paperclip search returned nothing.
- **Severity:** MINOR (uncertainty noted only in research file, not flagged in main docs).

### 27. Citation_audit.md is broken (vsh errors)
- **Where:** `docs/dengue/citation_audit.md` lines 1-22.
- **What:** The file is mostly bash-error output ("ERR: vsh: get: command not found"). Despite this, FINAL.md and response_to_prior_reviewers.md both cite this file as evidence of audit traceability.
- **Severity:** MAJOR (the cited audit file is broken/empty).
- **Suggested fix:** Either rebuild citation_audit.md as a clean table of verified PMIDs (the data is available in `paperclip_research/`), or change the cross-references to point at `paperclip_research/*.md` directly.

### 28. The stub S5_iedb_validation table (zero-fill rows)
- **Where:** `S5_iedb_validation.md` says "STATUS: pending IEDB validation cohort pull"; the S5 TSV is all zeros.
- **What:** FINAL.md line 121 cites "Supplementary Table S3" (not S5) as evidence of citation audit; S5 is referenced indirectly in FINAL ref list (S1-S5) but the supplementary stub is unfinished. Cover letter says all artefacts at v1.0-dengue-results are committed.
- **Severity:** MINOR (S5 incomplete but advertised as supplementary).
- **Suggested fix:** Either complete S5 before submission or remove S5 from the supplementary set.

---

## SUMMARY

- **Critical findings (#1-#6):** 6 critical numerical/factual disagreements between FINAL.md and canonical artefacts. Specifically: (a) post-hoc-optimal weight inverted (w_A=0.25 in FINAL but w_A=0.75 per S1/S2/S4); (b) headline `r = 0.45` not in any TSV/JSON (canonical is 0.35 equal-weight or 0.352 post-hoc-optimum); (c) Tier B normalised range "0.97-1.00" wrong (true: 0.83-1.00); (d) per-serotype Tier A+B values "0.95"/"0.87" wrong (true: ~0.61 high, ~0.46 low); (e) EDE off-diagonal mean "0.40 ± 0.05" wrong (true: 0.51 ± 0.01); (f) HLA equity "0 of 5 regions below 0.7" wrong (per TSV, Brazil 4/4 below 0.7).
- **Major findings (#7-#20):** broken cross-refs (expanded_pdb_manifest, ede_residue_provenance), 8 orphan citations (refs 12-15, 21-25), unsupported "12 PMIDs" count, response document points at draft not FINAL, fig1 + fig6 not body-cited.
- **Minor findings (#21-#28):** CRediT 5-axis drift, Versiani eLocator inconsistency, Rouvinski PMID typo in methods, EDE residue count mismatch, broken citation_audit.md.

Critical fixes #1-#6 will be applied to FINAL.md directly.
