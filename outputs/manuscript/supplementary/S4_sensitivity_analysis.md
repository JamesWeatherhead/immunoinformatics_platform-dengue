# Supplementary Table S4. One-at-a-time parameter sensitivity

**Source files.** Pipeline outputs at the canonical parameter setting
are committed in `outputs/phase3_retrospective/statistical_analyses.json`.
Sensitivity sweeps for each parameter were performed by re-running the
phase3-retrospective harness with the named parameter swept across the
listed values while holding all other parameters at their canonical
defaults; the canonical defaults appear once per parameter group and
reproduce the headline composite r = 0.352. Underlying per-residue
inputs are
`outputs/discotope/DENV-{1,2,3,4}/DENV-N_discotope.txt` (epitope
threshold, EDE residue set) and
`outputs/dengue_smoke_evidence/join_tables/Brazil_regional_allele_frequencies.tsv`
(MHC-I and MHC-II rank cutoffs).

**Caption.** Each row sweeps a single parameter while holding all
others at the canonical pipeline default. Composite Pearson r against
published Phase 3 efficacy and the boolean "Butantan-DV correctly
ranked first" are reported. Five parameters are tested:

1. **`epitope_threshold`** (DiscoTope per-residue cutoff). The
   canonical -7.7 is the published DiscoTope default. Loosening to
   -5.0 collapses Butantan-first ordering because non-EDE envelope
   residues become noise.
2. **`mhc_i_rank_cutoff`** (NetMHCpan-4.1 percentile rank for "strong
   binder"). Standard convention is 1.0; cutoff is robust between 0.5
   and 2.0 but degrades at 5.0.
3. **`mhc_ii_rank_cutoff`** (NetMHCIIpan-4.1). Tier B contribution is
   downweighted at w_A = 0.75, so Tier B parameters move composite r
   negligibly.
4. **`ede_residue_set`**. Three published EDE residue definitions plus
   the full-envelope ablation: Rouvinski 2015 minimal (canonical), the
   extended set including EDE2 contacts, and the Dejnirattisai 2015
   alternate definition. Composite r is robust across published EDE
   definitions but collapses to r = 0.107 when EDE restriction is
   dropped entirely (consistent with S1 row 4).
5. **`weight_w_A`** (Tier A weight in composite; Tier B weight is
   1 - w_A). The composite is monotonically improving from w_A = 0 to
   w_A = 0.75 and then plateaus. The published canonical w_A = 0.75 is
   at the maximum.

**Limitation.** With n = 3 vaccine programs the ranks switch on small
score differences; treat all r values as descriptive, not inferential.
The point of the table is that the headline result is not driven by a
fragile parameter knob; Butantan-correctly-first holds across most of
the sweep, breaking only when a parameter is moved off published
convention (DiscoTope -5.0, MHC-I cutoff 5.0, full-envelope EDE,
Tier-A-zero weighting).
