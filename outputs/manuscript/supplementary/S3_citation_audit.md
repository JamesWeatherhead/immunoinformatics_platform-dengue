# Supplementary Table S3. Citation audit for load-bearing claims

**Source files.** Manuscript text from
`outputs/manuscript/dengue_npj_vaccines_draft.md`. Each load-bearing
empirical claim was extracted by hand and back-traced to its primary
or review citation. PMIDs were verified via paperclip (`paperclip
fetch <pmid>` against a corpus pulled 2026-05-06); preprints and
non-PubMed-indexed records (Versiani 2026 *Sci Adv*, JessEV
*Bioinformatics* 2020, Estofolete-Nogueira 2026 *npj Vaccines*) carry
DOI rather than PMID and are flagged `n` in the verified column where
paperclip could not retrieve full text on this pass.

**Caption.** Each row is a single empirical, methodological, or
parameter-choice claim that the manuscript depends on. The
`paperclip_verified_y_n` column reports whether the citation full text
was successfully retrieved and the claim quotation matched. Three
rows (Versiani 2026 *Sci Adv*, JessEV, ADE review citation) are
flagged `n` because the corresponding records are either too recent
for full-text indexing or hosted outside the paperclip-indexed
corpora (PMC, bioRxiv, medRxiv, arXiv); these require manual
verification before submission. Three further rows are
`NA_internal_review` or `NA_database_release` and refer to non-citable
sources (a prior reviewer comment, a database release tag); they are
listed for traceability. Four rows are flagged `n` despite having a
PMID because the citation is currently to a *review* article rather
than the primary measurement source; per Reviewer 1's prior critique,
these should be replaced with the primary papers before resubmission.

**Action items for resubmission.**
1. Replace ADE-review citation (PMID 21737702) with primary Halstead
   PMID (this is the load-bearing claim Reviewer 1 specifically called out).
2. Verify Versiani-McCaffrey 2026 *Sci Adv* full text against the
   pipeline-architecture description in Section 2.1 once available in PMC.
3. Pull a primary measurement of DENV-4 circulating diversity to
   replace review citation 28492867.
4. The seronegative-CYD-TDV-underperformance claim (PMID 29282091) is
   a primary trial paper but the FcgammaR-balance *mechanism* is
   inferred from a separate review; either cite both or attribute the
   mechanism more cautiously in Section 3.5.
