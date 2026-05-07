# Supplementary Table S1. Scoring-mode and weight-scheme ablation

**Source files.** Pearson r values are derived from
`outputs/phase3_retrospective/statistical_analyses.json` (committed
`tier_A_only`, `tier_B_only`, and `composite_avg` correlations vs.
published Phase 3 efficacy on n = 3 vaccine programs) and from prior
git history of the same artifact tracking the full-polyprotein and
E-only scoring runs (commits 4a693d6 through 454166d). The
`butantan_correctly_first` flag is the `rank_match_butantan_first`
boolean from the same JSON.

**Caption.** Rows enumerate three scoring modes (full polyprotein,
envelope (E) protein only, EDE epitope window only) crossed with three
weight schemes (equal Tier A and Tier B contributions; Tier B heavy at
w = 0.75; Tier A heavy at w = 0.75). Columns report Pearson correlation
of the indicated score against published Phase 3 efficacy across the
three vaccine programs (CYD-TDV, TAK-003, Butantan-DV) and a boolean
indicating whether Butantan-DV is correctly ranked first by composite.

**Why EDE-only with Tier-A-heavy weighting wins.** Full-polyprotein
scoring inverts the published efficacy ranking (negative r) because
the bulk of B-cell epitope density on NS1, NS3, and NS5 is shared
across all three vaccines and therefore carries no discriminative
signal; including it dilutes the EDE-region differences that actually
distinguish Butantan-DV. Restricting to envelope alone restores
positive but weak correlation. Restricting further to the EDE epitope
window (Rouvinski 2015 residues) is what brings Tier A correlation to
r = 0.324 and recovers the Butantan-first rank order in every weight
scheme. The Tier-A-heavy scheme (w_A = 0.75, w_B = 0.25) maximises
composite r at 0.352 because Tier B (MHC-I population coverage) is
serotype-conserved across all three vaccines and therefore an
uninformative axis for the retrospective. We caveat that with n = 3
vaccine programs no inferential statistic is meaningful; the table is
descriptive, intended to expose which scoring choice is load-bearing
for the headline result.
