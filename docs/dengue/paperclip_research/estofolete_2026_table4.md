# Estofolete et al. 2026, npj Vaccines 11:68 — Table 4 verbatim and methods extract

Source artefact note: paper is not yet ingested in the paperclip corpus
(`paperclip search` and `paperclip lookup doi 10.1038/s41541-026-01400-4` both
return no hit on 2026-05-07). The paperclip search ran (queries: full title,
"Estofolete Nogueira composite immune correlate npj Vaccines 2026",
"Estofolete dengue vaccine immune correlate", "Estofolete" --all,
"Nogueira dengue vaccine correlate Brazil 2026") confirmed absence.
Verbatim Table 4 and methods text below were therefore retrieved by direct
HTML pull from nature.com (HTML parsed cell-by-cell out of
`/articles/s41541-026-01400-4/tables/4`) and verified against the article
body. When paperclip ingests the article, this file should be re-confirmed
against the paperclip-indexed copy.

## Bibliographic info

- Citation: Estofolete CF, Saivish MV, Nogueira ML, Vasilakis N. From
  promise to pitfalls: immunological lessons from dengue vaccines and
  their implications. npj Vaccines 11, 68 (2026).
- DOI: https://doi.org/10.1038/s41541-026-01400-4
- URL: https://www.nature.com/articles/s41541-026-01400-4
- Received: 23 October 2025
- Accepted: 03 February 2026
- Published: 14 February 2026
- License: CC BY 4.0 (Open Access)
- Article type: Review (narrative; no new datasets)
- Corresponding author: Nikos Vasilakis (UTMB)
- PMID / PMC: not yet indexed in paperclip nor in PubMed Central as of
  2026-05-07; nature.com canonical landing page is the only retrievable
  source.
- Funding: NIH U01AI151807 (CREATE-NEO, NV); FAPESP 2013/21719-3,
  2022/03645-1 (MLN), 2022/09229-0 (CFE); INCT Dengue 465425/2014-3;
  INCT Viral Genomic Surveillance and One Health 405786/2022-0; CNPq
  Research Fellowship (MLN); FAPESP 2020/12875-5 and 2023/09590-7 (MVS).

## Table 4 verbatim

Table 4 is titled (verbatim): "Minimal composite correlate (proposal #) and
practical assays for dengue vaccine trials". Footnote # reads verbatim:
"Conceptual framework (hypothesis-generating; not clinically validated)."

Column headers (verbatim, in order): Component | What it captures |
Operational metric | Assay/method | Sampling window | Example decision rule.

Row 1 (Tier A, neutralization breadth):
- Component: "Neutralization breadth index (Tier A)"
- What it captures: "Potency and coverage across DENV-1-4, adjusted for
  epitope quality"
- Operational metric: "Breadth-weighted PRNT 50/90 across 4 serotypes ×
  EDE-competition factor"
- Assay/method: "PRNT 50/90; E-dimer competition ELIZA;
  epitope-resolved neutralization panel"
- Sampling window: "Peak (day 28-60) + durability (6, 12, 24 mo)"
- Example decision rule: "Index >= prespecified X; failure if below X
  and narrowness to 1-2 serotypes"

Row 2 (Tier A, avidity add-on):
- Component: "Avidity (Tier A add-on)"
- What it captures: "GC maturation / Ab avidity"
- Operational metric: "Chaotrope-based avidity index vs E antigens"
- Assay/method: "Avidity ELIZA (e.g., urea wash)"
- Sampling window: "Day 60; 6-12 mo"
- Example decision rule: "Avidity >= Y linked to lower breakthrough rates"

Row 3 (Tier B, memory B-cell breadth):
- Component: "Memory B-cell breadth (Tier B)"
- What it captures: "Recall potential across serotypes"
- Operational metric: "% antigen-specific MBC recognizing E/NS1 for DENV-1-4"
- Assay/method: "B-cell ELISpot or Ag-labeled flow cytometry"
- Sampling window: "6-12 mo"
- Example decision rule: "Breadth index >= Z predicts year-2+ protection"

Row 4 (Tier B, CD8+ polyfunctionality):
- Component: "CD8+ T-cell polyfunctionality (Tier B)"
- What it captures: "NS-directed cellular breadth"
- Operational metric: "% polyfunctional CD8+ (IFN-gamma/IL-2/TNF-alpha) to
  NS3/NS5 peptide pools"
- Assay/method: "ICS/flow cytometry; ELISpot"
- Sampling window: "Day 28-90"
- Example decision rule: "Composite CD8+ score >= W associated with lower
  VCD odds at matched PRNT"

Note on terminology: the article reproduces "ELIZA" (with a Z) in cells of
Table 4 and "ELISA" elsewhere; the verbatim cell text above preserves
"ELIZA" as printed.

## Pre-trial-computable vs post-Phase-1-required

The paper does not partition Table 4 rows by pre-trial computability; it
treats all four sub-scores as requiring trial-derived assay readouts. None
of the four operational metrics (breadth-weighted PRNT, chaotrope avidity,
% antigen-specific MBC, % polyfunctional CD8+ ICS) is computable from
sequence or structure alone. All four require post-Phase-1 (or post-CHIM)
biological samples:

- Tier A.1 (Neutralization breadth index): requires immune sera tested in
  PRNT and EDE-competition ELISA; sampling windows day 28-60 and 6-12-24 mo.
- Tier A.2 (Avidity add-on): requires immune sera tested by chaotrope-wash
  ELISA; sampling windows day 60 and 6-12 mo.
- Tier B.1 (Memory B-cell breadth): requires PBMC for B-cell ELISpot or
  antigen-labeled flow; sampling window 6-12 mo.
- Tier B.2 (CD8+ polyfunctionality): requires PBMC for ICS/flow or ELISpot
  with NS3/NS5 peptide pools; sampling window day 28-90.

Sequence- and structure-only proxies for each row are an interpretation
introduced in the present manuscript, not a claim of the source paper.
This distinction must be respected verbatim in our Methods.

## Methods extract for composite computation

Estofolete et al. is a narrative review (their Ethics statement: "This
article is a narrative review and did not involve any new studies"). They
therefore have no formal Methods section that operationally computes the
composite. Their Abstract describes the review's evidence-selection method
verbatim: "Evidence was purposively selected from high-impact clinical
trials, long-term follow-ups, authoritative assessments, and mechanistic
studies judged most decision-relevant to current vaccine design and policy.
The selection process was non-systematic and reflects the authors' expert
appraisal."

The composite framework is introduced in section "Unresolved Questions and
Correlates of Protection". Verbatim load-bearing sentences:

1. "To address the limitations of PRNT-only readouts, Table 4 presents a
   model for a minimal composite panel intended to clarify prospective
   testing, although it is not yet a validated or regulatory-accepted
   correlate."
2. "While systems serology has identified Fc-feature and effector-function
   signatures associated with protection in natural infection cohorts,
   similar composite panels have not yet been consistently validated
   against clinical endpoints in dengue vaccine trials."
3. "The table details operational metrics, feasible multicenter assays,
   sampling windows, and example decision rules, thereby enabling
   prospective validation and immunobridging analyses while ensuring that
   each component remains explicitly falsifiable and trial-agnostic."

The Discussion-style anchor (cross-cutting lessons section) frames the
composite in 4 dimensions verbatim: "(i) epitope specificity, including
quaternary-epitope targeting, (ii) binding strength and maturation, such
as avidity, (iii) Fc-mediated effector functionality, and (iv) durability
and recall features".

Note: the "Fc-mediated effector functionality" dimension is described in
the body text but is NOT one of the four Table 4 rows; it is folded into
the Tier A neutralization-breadth index via the "EDE-competition factor"
multiplier. Our pipeline must not treat Fc effector tone as a separate
fifth Table 4 row.

## Citations from Estofolete supporting each Tier definition

- Tier A (neutralization breadth + EDE-competition + avidity): refs 24,
  30, 35, 36 (Dengvaxia high-titer paradox, TAK-003 DENV-2 dominance,
  durability decline); ref 80 (quaternary-epitope vs monomeric-domain
  antibody distinction); refs 76, 77, 78 (epitope specificity / avidity /
  Fc / durability framework); ref 21 (Nivarthi et al. systems serology of
  TV003).
- Tier B.1 (memory B-cell breadth): refs 64, 65 (TV003/TV005 CHIM showing
  expanded MBC repertoire associated with reduced post-challenge viremia);
  ref 21 (Nivarthi DENV-2 immunoprofiling of TV003).
- Tier B.2 (CD8+ polyfunctionality, NS-directed): refs 18, 29 (CD8+ to
  NS3/NS5 in viral clearance and cross-serotype protection); refs 28, 83
  (NS antigens absent in Dengvaxia/YF17D backbone vs present in TAK-003
  and Butantan-DV); refs 40, 41, 42, 43 (TAK-003 multifunctional T-cell
  responses including NS-focused responses).
- Whole composite framing: refs 64, 65, 76, 81 ("More reliable correlates
  will likely combine neutralization titers with assessments of antibody
  avidity, memory B-cell repertoires, and T-cell functionality").

## Alignment check vs current draft at HEAD of master

The current manuscript draft at
`/Users/jamesweatherhead/Desktop/dengue-fork/outputs/manuscript/dengue_npj_vaccines_draft.md`
correctly identifies the four Table 4 rows (Tier A neutralization breadth,
Tier A avidity, Tier B MBC breadth, Tier B CD8+ polyfunctionality) and
correctly preserves the Tier A / Tier B partition. It also correctly cites
the paper at "Estofolete CF, Nogueira ML et al. ... npj Vaccines 2026;11:68".
However, three drifts must be reconciled before submission. (1) The draft's
title and abstract attribute an "Immunological Composite Architecture (ICA)"
with five named axes (Geometry, Equity, Time, Intent, Effector tone) to
Estofolete et al.; that five-axis named architecture is NOT in the source
paper. Estofolete et al. propose a 4-row composite (not 5 axes) and frame
the cross-cutting lessons in 4 dimensions (epitope specificity,
binding/avidity, Fc effector functionality, durability/recall), not in
James's Geometry/Equity/Time/Intent/Effector-tone schema. The draft must
either reframe the ICA as our extension of Estofolete (clearly attributed
to this Perspective, not to them) or drop the ICA naming. (2) The draft's
Section 2.3 mapping treats Tier A.1, A.2, B.1, B.2 as sequence-and-structure
computable; the source paper treats all four as wet-lab-only at day 28-90 /
6-12-24 mo. Our Methods must be explicit that we are introducing
sequence-only proxies for each row, not implementing the rows themselves.
(3) The draft cites Estofolete as having shown the composite "recovered the
Butantan-DV > TAK-003 > CYD-TDV rank ordering"; the source paper does NOT
report any retrospective rank-ordering analysis. The composite is presented
as hypothesis-generating only ("Conceptual framework; not clinically
validated"). That sentence in our abstract and Section 1 must be removed
or reattributed to our own retrospective at Section 3.3.

(Word count of report below: 480.)
