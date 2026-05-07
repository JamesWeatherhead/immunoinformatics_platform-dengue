# Phase 3 dengue vaccine parent strain provenance

This directory holds GenBank protein records for the parent strains of the three
licensed Phase 3 dengue vaccine constructs. They are used by the retrospective
validation script (`scripts/dengue/retrospective_phase3_validation.py`) to
score CYD-TDV, TAK-003, and Butantan-DV through the pipeline and compare to
published clinical efficacy.

## Important: parent strain != licensed construct

The licensed vaccines carry attenuation mutations (3' UTR deletions in the
Butantan-DV components, PDK passage adaptations in TAK-003, YF17D chimerism
in CYD-TDV). For the Geometry and Equity axes that depend on the E protein
primary sequence, these mutations are not expected to substantively change
the predicted scores. For a stricter analysis, replace these FASTAs with
construct-specific sequences from each program's FDA Biologics License
Application supplement.

## Construct summary

| Vaccine | Sponsor | Backbone strategy | Parent strains (DENV-1/2/3/4) | Pivotal trial |
|---|---|---|---|---|
| CYD-TDV (Dengvaxia) | Sanofi Pasteur | YF17D chimeric | PUO-359/82, PUO-218, PaH881/88, 1228 | Capeding 2014, Hadinegoro 2015, Sridhar 2018 |
| TAK-003 (Qdenga) | Takeda | DENV-2 PDK-53 | 16007 PDK-13, 16681 PDK-53, 16562 PGMK-30, 1036 PDK-48 | Biswal 2019, Tricou 2024 |
| Butantan-DV / TV003 | NIAID / Butantan | rDEN attenuated tetravalent | WestPac-74, NewGuineaC (chimeric DENV-2/4), Sleman-78, Dominica-814669 | Kirkpatrick 2016, Kallas 2024, Nogueira 2024 |

## Citations

- Capeding MR et al. Lancet 2014;384:1358 (CYD-TDV Phase 3 Asia)
- Hadinegoro SR et al. NEJM 2015;373:1195 (CYD-TDV long-term)
- Sridhar S et al. NEJM 2018;379:327 (CYD-TDV serostatus stratified)
- Huang CY-H et al. J Virol 2013;87:7036 (TAK-003 backbone)
- Biswal S et al. NEJM 2019;381:2009 (TAK-003 Phase 3 1y)
- Tricou V et al. Lancet Glob Health 2024;12:e257 (TAK-003 4.5y)
- Whitehead SS et al. PLoS NTD 2017;11:e0005584 (TV003 Phase 1)
- Kirkpatrick BD et al. Sci Transl Med 2016;8:330ra36 (TV003 challenge)
- Kallas EG et al. NEJM 2024;390:397 (Butantan-DV Phase 3 single-dose)
- Nogueira ML et al. Lancet Infect Dis 2024;24 (Butantan-DV extended follow-up)
- Estofolete CF, Saivish MV, Nogueira ML, Vasilakis N. npj Vaccines 2026;11:68
