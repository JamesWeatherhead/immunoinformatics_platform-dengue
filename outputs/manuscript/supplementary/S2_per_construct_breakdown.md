# Supplementary Table S2. Per-construct Tier A and Tier B breakdown

**Source files.** Vaccine and parent accession identifiers from
`outputs/phase3_retrospective/table4_mapping/ranked_candidates.tsv`
(12 constructs: 3 vaccines x 4 serotypes). Per-serotype Tier A
neutralization-breadth scores, Tier B CD8 polyfunctionality scores,
and `n_mhci_strong_binders` from
`outputs/table4_dengue/ranked_candidates.tsv`. Epitope-residue counts
(`n_epitope_residues`) computed from
`outputs/discotope/DENV-{1,2,3,4}/DENV-N_discotope.txt` by counting
residues flagged with the `<=B` predicted-epitope marker. Composite
column is the Tier-A-heavy (w_A = 0.75, w_B = 0.25) weighted sum used
for the headline result.

**Caption.** Each row is one of the 12 parent strains underlying the
three Phase 3 vaccine programs. Tier A and Tier B scores are
serotype-specific and identical across vaccines because the pipeline
operates on the parent reference proteome rather than the engineered
backbone. The composite column is therefore informative of *which
serotype contributes* to a vaccine's overall composite, not of
backbone-specific attenuation. Aggregated per-vaccine composites
(weighted average across serotype components, after EDE-window
restriction and rank-normalization) yield the program-level scores
0.872 (CYD-TDV), 0.750 (TAK-003), 0.893 (Butantan-DV) reported in the
main text. The table is intended for reviewers who want to verify
that no single serotype drives the composite ordering: DENV-4
contributes the highest Tier A score (0.560), and Butantan-DV is the
only program whose DENV-4 component (rDEN4Delta30) is a parent-strain
attenuation rather than a chimera, which is consistent with its
Phase 3 efficacy advantage.

**Limitation.** Per-vaccine composite scores in the main text use the
EDE-restricted, weight-tuned pipeline (S1 row 9). This S2 table reports
the un-tuned per-construct outputs to support reproducibility audits.
