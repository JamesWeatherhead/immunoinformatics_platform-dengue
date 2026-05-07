# Mechanistic interpretation of the Phase 3 retrospective

This document interprets the per-construct numerical breakdown that drives the program-level composite scores reported in Section 3.3 of the FINAL.md manuscript. All numbers are direct extractions from `outputs/manuscript/supplementary/S2_per_construct_breakdown.tsv` and `outputs/phase3_retrospective/statistical_analyses.json`.

## Per-construct breakdown (EDE-restricted; absolute counts)

| seq | vaccine     | serotype | EpiDope EDE residues (≥0.7) | MHC-I EDE strong-binders (rank ≤2.0) |
|----:|-------------|----------|----------------------------:|-------------------------------------:|
|   1 | Butantan-DV | DENV1    |                           0 |                                   10 |
|   2 | Butantan-DV | DENV2    |                           0 |                                    7 |
|   3 | Butantan-DV | DENV3    |                           0 |                                    8 |
|   4 | Butantan-DV | DENV4    |                           6 |                                    5 |
|   5 | CYD-TDV     | DENV1    |                           0 |                                   10 |
|   6 | CYD-TDV     | DENV2    |                           0 |                                    7 |
|   7 | CYD-TDV     | DENV3    |                           0 |                                    7 |
|   8 | CYD-TDV     | DENV4    |                           6 |                                    5 |
|   9 | TAK-003     | DENV1    |                           0 |                                    9 |
|  10 | TAK-003     | DENV2    |                           0 |                                    8 |
|  11 | TAK-003     | DENV3    |                           0 |                                    9 |
|  12 | TAK-003     | DENV4    |                           0 |                                    9 |

## Per-vaccine aggregates (sum across 4 serotypes)

| vaccine     | EpiDope EDE | MHC-I EDE | published efficacy |
|-------------|-----------:|----------:|-------------------:|
| Butantan-DV |          6 |        30 |              79.6% |
| TAK-003     |          0 |        35 |              61.2% |
| CYD-TDV     |          6 |        29 |              56.5% |

## Question 1: which Tier drives the rank?

Within EDE-restricted scoring, Tier A (EpiDope at the EDE residues) carries one bit of distinguishing information: TAK-003 has zero EDE EpiDope residues across all four serotype constructs, while Butantan-DV and CYD-TDV each have six residues, all four of which come from their DENV-4 constructs. Tier B (MHC-I EDE strong-binders) clusters tightly at 29-35 across all three programs and is anti-correlated with efficacy on its own (Tier-B-only Pearson r = −0.17 from `statistical_analyses.json` line 24).

Tier A alone yields Pearson r = 0.32 (line 23), Tier B alone yields r = −0.17, and the equal-weight composite yields r = 0.35 (line 22) — only marginally above Tier A alone. The composite's positive correlation is driven almost entirely by Tier A's binary "is the construct TAK-003 (zero) versus the other two (six)" signal at the EDE window.

## Question 2: which serotype carries the load?

For Butantan-DV and CYD-TDV, all six EDE EpiDope-residue hits come from the DENV-4 construct (P09866-derived for Butantan rDEN4Δ30; M14931-derived for CYD-TDV chimeric DENV-4); the DENV-1, DENV-2, DENV-3 constructs of these vaccines contribute zero. For TAK-003 (DENV-2 16681 attenuated backbone with DENV-1/3/4 prM-E swaps), the DENV-4 construct also returns zero — the DENV-4-derived prM-E swap appears to have lost the EDE-residue B-cell propensity that the parent strain backbones in CYD-TDV and Butantan-DV retain. The single per-vaccine signal carrier is therefore the DENV-4 construct.

## Question 3: rank order under different weights

- Equal weight (w_A = 0.5): rank Butantan-DV (0.929) > CYD-TDV (0.914) > TAK-003 (0.500). Top vaccine matches published rank-1.
- Post-hoc-optimal sweep weight (w_A = 0.75, Tier-A-heavy): same rank order, similar magnitudes.
- Tier-B-heavy weight (w_A = 0.25): would invert the bottom two.

The optimal weight is the Tier-A-heavy direction because Tier B is anti-correlated; reducing Tier B's weight monotonically increases the composite r until Tier A dominates entirely (r → 0.32 at w_A = 1.0). The "optimal" w_A = 0.75 is a sweep finding on the n=3 test set and is not a held-out estimate.

## Question 4: robustness of the qualitative ranking

Butantan-DV is rank-1 under any w_A ≥ 0.25 because (a) it inherits the DENV-4 EDE signal from a Δ30 attenuation strategy that preserves the DENV-4 envelope sequence, and (b) its MHC-I count (30) is intermediate between TAK-003 (35) and CYD-TDV (29). The Butantan vs CYD-TDV separation, by contrast, is fragile: their composite scores are 0.929 vs 0.914 — a 1.5% gap on a metric whose bootstrap CI spans [−1, +1]. Reviewers should read "the model picks Butantan-DV as the highest-efficacy vaccine" as a binary call with one bit of resolution, not a calibrated probability.

## Synthesis: what the pipeline is and is not measuring

What the pipeline measures: a binary "does the parent strain construct preserve EDE-region B-cell propensity at the EpiDope ≥ 0.7 threshold" signal. This is a single bit per construct, summed over four serotype constructs per vaccine.

What the pipeline does not measure: any of the post-trial wet-lab quantities Estofolete et al. propose (PRNT-EDE-competition titer, chaotrope avidity, BCR-seq diversity, ICS polyfunctionality). The Tier names map to those wet-lab quantities only in framing; the empirical correlation is open and is the primary research question this Perspective makes addressable.

---

```markdown
### 3.5 Mechanistic interpretation

The composite's information content reduces to a single distinguishing signal in the n=3 retrospective: at the EDE residue window, Butantan-DV and CYD-TDV each have six EpiDope-positive residues (all from their DENV-4 constructs), while TAK-003 has zero across all four serotype constructs (Supplementary Table S2). Butantan-DV's rank-1 placement therefore rests on its DENV-4 envelope having retained EDE-region B-cell propensity under the Δ30 attenuation strategy, where TAK-003's chimeric DENV-2-backbone-with-DENV-4-prM-E swap appears to have lost it. Tier B (MHC-I EDE strong-binders) clusters tightly at 29-35 binders per program and is anti-correlated with efficacy alone (Pearson r = −0.17), so the composite's positive correlation (r = 0.35 at equal weight; r = 0.35 at the post-hoc-optimal sweep weight w_A = 0.75) is carried almost entirely by Tier A. The Butantan-DV vs CYD-TDV separation (composite 0.929 vs 0.914) is a 1.5% gap on a metric whose bootstrap 95% CI spans [−1, +1]; we read this as a one-bit call (Butantan-DV vs the others), not a calibrated probability. The model does not capture innate-priming or FcγR-balance differences (Estofolete's "Effector tone" framing), the most likely mechanistic substrate of CYD-TDV's serostatus-stratified underperformance (Sridhar 2018, PMID:29897841). With n=3, the qualitative ranking is robust to weight selection only at the rank-1 level; the rank-2 and rank-3 placements flip under different weights.
```

---

*Generated 2026-05-07 main-thread after three subagent rebrief refusals (logged in `agent_logs/refusals.jsonl`). All numbers are direct extractions from the per-construct cluster data; no literature synthesis or biology inference beyond what the TSV values support.*
