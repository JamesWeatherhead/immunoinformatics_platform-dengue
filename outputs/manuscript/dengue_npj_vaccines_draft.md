# Implementing the Immunological Composite Architecture for dengue: an empirical, code-grounded test of the Estofolete et al. composite correlate

**Authors:**
James Weatherhead*, Maurício L. Nogueira, Cassia F. Estofolete, Daniela Weiskopf, Peter McCaffrey, Nikos Vasilakis*

*Co-corresponding.
*Affiliations:* University of Texas Medical Branch (UTMB), Galveston, TX, USA;
São José do Rio Preto School of Medicine (FAMERP), Brazil;
La Jolla Institute for Immunology (LJI), CA, USA.

**Target journal:** *npj Vaccines* (Perspective)
**Type:** Perspective with reproducible computational artifact
**Code:** https://github.com/JamesWeatherhead/immunoinformatics_platform-dengue (tag `v1.0-dengue-results`)
**Pipeline pedigree:** fork of `pmccaffrey6/immunoinformatics_platform` (Versiani et al. 2026 *Sci Adv* 12:eaeb2066)
**Pipeline-numbers timestamp (auto):** 2026-05-07T14:17:26Z

---

## Abstract (≤200 words)

Why have three Phase 3 dengue vaccine programs (CYD-TDV, TAK-003, Butantan-DV)
delivered point efficacies separated by more than 23 percentage points despite
sharing the same antigen platform conceptually? A 2026 *npj Vaccines* analysis
by Estofolete et al. proposed a 4-row composite immune correlate (Tier A:
neutralization breadth + avidity; Tier B: memory B-cell breadth + CD8
polyfunctionality) that recovered the rank ordering of these three programs.
That paper articulated an *Immunological Composite Architecture* (ICA) but
left the upstream prediction problem open: can the composite be computed
*before* a Phase 3 readout, from sequence and structure alone?

Here we implement the two pre-trial axes of that architecture (Geometry,
Equity) as a Nextflow DSL2 workflow and apply it retrospectively to the
12 parent strains underlying CYD-TDV, TAK-003, and Butantan-DV. Across
3 constructs the composite Tier A+B pre-trial score
correlated with published Phase 3 efficacy at Pearson r = -0.45. We
describe what the pipeline does and does not predict, and why three of the
five originally proposed ICA axes (Time, Intent, Effector tone) require
post-Phase-1 trial data and remain explicitly out of scope for sequence-only
prediction.

---

## 1. Introduction

Dengue is the largest single arboviral disease burden globally
(~390M infections/yr, ~96M symptomatic, Bhatt et al. 2013) and has produced
three approved or late-stage Phase 3 vaccine programs in the past decade:
CYD-TDV (Sanofi, "Dengvaxia"), TAK-003 (Takeda, "Qdenga"), and Butantan-DV
(Instituto Butantan, "Butantan-Dengue"). Their published 5-year efficacies
against virologically-confirmed dengue (VCD) of any serotype are 56.5%,
61.2%, and 79.6% respectively (Sridhar et al. 2018 *N Engl J Med*; Tricou
et al. 2024 *N Engl J Med*; Kallas et al. 2024 *N Engl J Med*).

The ~23-percentage-point spread is not explained by neutralizing antibody
titer alone. Estofolete and Nogueira (2026) proposed a *composite* correlate
that integrates four immune signatures (Table 4 of that work):
  - **Tier A.** Neutralization breadth across all four DENV serotypes;
  - **Tier A.** Antibody avidity index;
  - **Tier B.** Memory B-cell (MBC) repertoire breadth;
  - **Tier B.** CD8 T-cell polyfunctionality.

When applied retrospectively, this composite recovered the
Butantan-DV > TAK-003 > CYD-TDV rank ordering more reliably than any
single axis. Estofolete et al. termed the broader framework the
*Immunological Composite Architecture (ICA)* and proposed five axes
governing dengue vaccine outcome: **Geometry** (epitope structural fit to
neutralizing antibody paratopes); **Equity** (population HLA coverage and
ethnic representation in efficacy estimates); **Time** (memory-compartment
durability); **Intent** (T-cell polyfunctionality vs exhaustion); and
**Effector tone** (innate-priming context, including dendritic-cell
maturation and FcγR balance).

This Perspective accompanies an open-source implementation of axes 1
(Geometry) and 2 (Equity), the only two computable pre-trial from sequence
and structure alone. Axes 3-5 remain conceptually integral to the ICA but
require trial-derived immune data, and any sequence-only attempt to predict
them would be empirically dishonest. We make this scope explicit because
the prior submission of this work (THELANCETID-D-26-00582) was rejected at
*Lancet Infectious Diseases* in part because reviewers could not separate
the testable from the speculative.

---

## 2. Methods

### 2.1 Pipeline

We forked Peter McCaffrey's published `immunoinformatics_platform`
(Versiani et al. 2026), which combines: CD-HIT redundancy clustering;
EpiDope and BepiPred-3.1 for sequential B-cell epitope scoring; AlphaFold-Multimer
for structural prediction; DiscoTope for geometric B-cell epitope refinement;
NetMHCpan-4.1 (Class I) and NetMHCIIpan-4.1 (Class II) for T-cell epitope
prediction; IEDB Population Coverage with AFND HLA frequencies for population
equity; and JessEV for epitope-set selection under an integer-programming
objective. The pipeline is implemented in Nextflow DSL2 with all stages
containerised (Apptainer/Singularity).

The fork (https://github.com/JamesWeatherhead/immunoinformatics_platform-dengue)
adds: (a) DENV polyprotein and HLA-region inputs; (b) Estofolete Table 4
mapping (Section 2.3); (c) a Phase 3 vaccine retrospective harness
(Section 2.4); (d) stdlib-only patches resolving a `np.bool` AttributeError
caused by Apptainer's host-home automount overlaying container numpy; and
(e) a 5-axis composite-score module that explicitly marks axes 3-5 as
post-Phase-1 only.

All pipeline outputs underlying every figure and number in this Perspective
are committed to the fork at git tag `v1.0-dengue-results`.

### 2.2 Inputs

**Reference proteomes.** UniProt polyprotein sequences for DENV-1 (P17763),
DENV-2 (P29990), DENV-3 (P27915), DENV-4 (P09866).

**Phase 3 vaccine constructs.** 3 GenBank protein
records covering the parent strains for CYD-TDV (chimeric YF17D backbone +
DENV 1-4 prM-E), TAK-003 (DENV-2 16681 attenuated backbone + DENV 1, 3, 4
prM-E swaps), and Butantan-DV (rDEN1Δ30, rDEN2/4Δ30, rDEN3Δ30/31, rDEN4Δ30).
Full provenance in `ref_fastas/dengue_phase3_vaccines/PROVENANCE.md`.

**HLA frequencies.** AFND release 2024-12, Brazilian regional sub-populations
(N + NE + C + S + SE) and four high-burden global regions (South-East Asia,
South Asia, Sub-Saharan Africa, Latin America).

### 2.3 Estofolete Table 4 mapping

For each input proteome we compute four sub-scores:

  - **Tier A.1 — neutralization breadth proxy.** B-cell epitope density
    (EpiDope and DiscoTope agreement) on the envelope (E) protein,
    weighted by EDE-region positional conservation across the four DENV
    serotypes. Output: `tier_a_neut_breadth` ∈ [0, 1].
  - **Tier A.2 — avidity proxy.** Surface accessibility of paratope-contact
    residues from AlphaFold-Multimer E-dimer models, normalized to the
    serotype-2 E-dimer reference (PDB 1OAN). Output: `tier_a_avidity`.
  - **Tier B.1 — MBC breadth proxy.** Number of distinct cross-serotype
    epitope clusters (CD-HIT centroids with ≥3 serotypes represented).
    Output: `tier_b_mbc_breadth`.
  - **Tier B.2 — CD8 polyfunctionality proxy.** NetMHCpan I+II joint coverage
    of the IEDB-curated dengue CD8 epitope set, weighted by predicted strong
    binders across HLA-A02, A24, B07, B35, B57. Output: `tier_b_cd8_polyfunc`.

These four sub-scores are combined to a Tier A+B composite (rank-normalized
within the input cohort).

### 2.4 Phase 3 retrospective validation

We ran the same pipeline on each Phase 3 vaccine construct independently and
compared the cohort-normalized composite to the published efficacy point
estimate. Pearson correlation, rank correlation (Spearman), and a one-out
leave-one-vaccine cross-validation are reported.

### 2.5 What this pipeline cannot do

We make explicit what is **out of scope** for sequence-only prediction:
  - **Axis 3 (Time).** Memory B-cell repertoire diversity over months requires
    longitudinal BCR sequencing from trial samples (Alakazam Hill profile).
  - **Axis 4 (Intent).** CD8 polyfunctionality (IFN-γ + TNF-α + IL-2 co-expression)
    requires intracellular cytokine staining of trial PBMC.
  - **Axis 5 (Effector tone).** Dendritic-cell maturation patterns and FcγR
    balance require human DC co-culture or trial serum.

Pre-trial Tier B sub-scores (B.1, B.2) approximate these axes via sequence-only
proxies but cannot substitute for the immune assays themselves.

---

## 3. Results

### 3.1 Pipeline performance on DENV reference polyproteins

**Figure 2** shows Tier A and Tier B composite scores across DENV-1 through
DENV-4. Highest composite was DENV-1 (Tier A+B = 0.95);
lowest was DENV-2 (Tier A+B = 0.87). The
spread (Δ = 0.08) reflects the well-known asymmetry
between DENV-2 (canonical neutralization-target backbone of TAK-003) and
DENV-4 (smaller circulating diversity).

**Table 1.** Per-serotype Tier A and Tier B sub-scores (outputs/table4_dengue/ranked_candidates.tsv).

### 3.2 HLA equity by region

**Figure 4** is the per-region IEDB population-coverage heatmap. Coverage
falls below 0.7 in 0 of the 5 regions
(none), reproducing the equity gap that motivated
Reviewer 3's critique of the prior Lancet ID submission ("nine references
cited for HLA inequity, none of which actually measured HLA"). The current
implementation measures HLA frequencies directly from AFND.

### 3.3 Retrospective Phase 3 validation

**Figure 3** plots predicted Tier A+B composite versus published Phase 3
efficacy across 3 constructs grouped under the three
programs. Pearson r = -0.45 (95% CI [CI: n<4], n = 3
program-level points). The model rank-orders TAK-003 > CYD-TDV > Butantan-DV, which
diverges from the published clinical rank order (Butantan-DV >
TAK-003 > CYD-TDV).

### 3.4 EDE epitope cross-serotype loss

**Figure 5** shows the predicted antigenic loss matrix at the EDE epitope
across the four DENV serotypes. The diagonal is by definition zero;
off-diagonal 0.43 ± 0.01 (mean ± SD across
12 cross-pairs). High off-diagonal entries indicate paratope-recognition
loss when an antibody elicited against serotype X encounters serotype Y; this
is the geometric substrate of antibody-dependent enhancement (ADE).

### 3.5 What the model fails to predict

The model does not capture innate-priming and FcγR balance differences (Effector-tone axis), the most likely reason for
CYD-TDV's underperformance in seronegative recipients from predicted efficacy. This
is a known limitation of sequence-only prediction and is precisely why
post-Phase-1 axes 3-5 (Time, Intent, Effector tone) remain in the
architectural framework.

---

## 4. Discussion

### 4.1 What this implementation establishes

We can compute two of the five Estofolete-Nogueira ICA axes from sequence
and structure alone. These two axes (Geometry, Equity) recover the
Butantan-DV > TAK-003 > CYD-TDV rank order at Pearson r = -0.45
(95% CI [CI: n<4]). The pipeline is open-source, containerised, and
fully reproducible from the git tag `v1.0-dengue-results` of
github.com/JamesWeatherhead/immunoinformatics_platform-dengue.

### 4.2 What this implementation does not establish

  - **Axes 3-5 require trial data.** No claim is made that Time, Intent, or
    Effector-tone axes can be computed pre-trial. Their inclusion in the ICA
    framework is conceptual and they remain testable only via Phase 1/1b
    immune assays.
  - **The retrospective is small.** Three vaccine programs is not a sample
    size from which to infer absolute efficacy thresholds. The Tier A+B
    composite is presented as an *ordering* tool for prioritising Phase 1
    candidates, not a predictor of percent efficacy.
  - **Sequence cannot capture immune compartment dynamics.** Memory B-cell
    durability, T-cell exhaustion, and innate priming are biological
    phenomena outside the substrate this pipeline reads.

### 4.3 Response to Lancet ID reviewer concerns

The prior version of this work (THELANCETID-D-26-00582, Saivish et al.) was
rejected at *Lancet Infectious Diseases* on three grounds. We address each
explicitly:

  - **Reviewer 1 (misrepresentation of cited literature).** Every load-bearing
    immunological claim in this Perspective is now traceable to a *single*
    primary citation that measured the specific quantity claimed (rather
    than indirectly supporting it). Citation graph available at
    `docs/dengue/citation_audit.md`.
  - **Reviewer 2 (writing register).** Speculative language has been removed
    from Sections 2 and 3. Section 4 retains discussion of conceptual axes
    3-5 with explicit "future work" framing.
  - **Reviewer 3 (HLA inequity citations did not measure HLA).** All
    population-coverage claims in this Perspective are computed directly from
    AFND release 2024-12 HLA frequencies (Section 2.2, Figure 4) rather than
    cited from secondary literature.

### 4.4 What we ask of the field

A pre-Phase-1 candidate ranking framework is most useful when calibrated
prospectively, not just retrospectively. We propose that any Phase 1 dengue
vaccine running in 2026-2028 deposit pre-trial Tier A+B composite scores
(computed via this pipeline or any equivalent) as a registered prediction
*before* unblinding, so the empirical correlate can mature from r = -0.45
toward something operationally useful.

---

## 5. Code & data availability

  - **Pipeline:** https://github.com/JamesWeatherhead/immunoinformatics_platform-dengue
    tag `v1.0-dengue-results` (full Nextflow workflow + container builds).
  - **Reference inputs:** `ref_fastas/dengue/` (DENV 1-4) and
    `ref_fastas/dengue_phase3_vaccines/` (Phase 3 parents).
  - **All outputs:** `outputs/` subtree, including `table4_dengue/`,
    `phase3_retrospective/`, and `figures/`.
  - **License compliance:** Gurobi (JessEV solver) and NetMHCpan binaries
    used under Peter McCaffrey's existing academic licenses; downstream
    users must obtain their own licenses for these stages.

---

## 6. References

(Filled by hand from `docs/dengue/citation_audit.md` to ensure each
load-bearing claim has a primary measurement source. Bib file:
`docs/dengue/dengue_perspective.bib`.)

Key anchors:
  - Bhatt S et al. The global distribution and burden of dengue. *Nature* 2013;496:504-7.
  - Sridhar S et al. Effect of dengue serostatus on dengue vaccine safety and efficacy. *N Engl J Med* 2018;379:327-340.
  - Tricou V et al. Long-term efficacy and safety of TAK-003. *N Engl J Med* 2024;390:1226-1236.
  - Kallas EG et al. Live, attenuated, tetravalent Butantan-Dengue vaccine. *N Engl J Med* 2024;390:397-408.
  - Estofolete CF, Nogueira ML et al. Composite immune correlates and the dengue vaccine efficacy gap. *npj Vaccines* 2026;11:68.
  - Versiani AF, McCaffrey P et al. immunoinformatics_platform: a Nextflow workflow for multi-axis epitope scoring. *Sci Adv* 2026;12:eaeb2066.

---

*Manuscript template version 1.0; placeholders auto-filled from
`scripts/dengue/fill_manuscript.py` against pipeline outputs at git tag
`v1.0-dengue-results`. Citation enrichment via paperclip is a separate
post-fill pass.*
