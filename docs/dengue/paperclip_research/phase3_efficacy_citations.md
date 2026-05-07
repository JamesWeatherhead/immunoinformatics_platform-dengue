# Phase 3 Dengue Vaccine Efficacy Citations

Verified efficacy numbers from PubMed E-utilities (verbatim from official abstracts) for the three landmark Phase 3 dengue vaccine trials. Numbers are exact; preserve decimal points (Lancet uses middle dot, NEJM uses period).

Source: NCBI E-utilities efetch on PubMed records, 2026-05-06.

Important corrections to original brief:
- Sridhar 2018 PMID is **29897841** (not 30067296)
- Tricou 2024 4.5-year paper is in **Lancet Global Health 12(2):e257-e270**, NOT NEJM 390:1226-1236 (which does not match this paper). PMID **38245116** (not 38477990)
- Kallas 2024 PMID is **38294972** (not 38294966)

---

## 1. Sridhar S et al. 2018 — CYD-TDV (Sanofi Dengvaxia), serostatus stratification

- **Citation:** N Engl J Med. 2018 Jul 26;379(4):327-340
- **DOI:** 10.1056/NEJMoa1800820
- **PMID:** 29897841
- **Vaccine:** CYD-TDV (Dengvaxia, Sanofi Pasteur)
- **Design:** Case-cohort post hoc reanalysis of 3 efficacy trials (CYD14, CYD15, CYD-TDV-23) using anti-NS1 IgG ELISA at month 13 to infer baseline serostatus
- **Population:** Ages 2-16y (and 9-16y subset). Trials: NCT00842530, NCT01983553, NCT01373281, NCT01374516
- **Follow-up:** 5 years
- **Funding:** Sanofi Pasteur

### Hospitalization for VCD (cumulative 5-year incidence; hazard ratio = vaccine vs. control)

| Subgroup | Age | Vaccine incidence | Control incidence | HR | 95% CI |
|---|---|---|---|---|---|
| Seronegative | 2-16y | 3.06% | 1.87% | **1.75** | 1.14 to 2.70 |
| Seronegative | 9-16y | 1.57% | 1.09% | **1.41** | 0.74 to 2.68 |
| Seropositive | 2-16y | 0.75% | 2.47% | **0.32** | 0.23 to 0.45 |
| Seropositive | 9-16y | 0.38% | 1.88% | **0.21** | 0.14 to 0.31 |

**Note:** This paper reports HRs for hospitalization, not classical "vaccine efficacy %" against symptomatic VCD. The abstract does not report per-serotype efficacy. Vaccine efficacy can be derived as (1 - HR), e.g., seropositive 9-16y VE = 79% (95% CI 69-86); seropositive 2-16y VE = 68% (95% CI 55-77). Severe VCD trends in same direction. Verify in main text Tables before citing derived VEs.

---

## 2. Tricou V et al. 2024 — TAK-003 (Takeda Qdenga), 4.5-year long-term efficacy

- **Citation:** Lancet Glob Health. 2024 Feb;12(2):e257-e270 *(NOT NEJM)*
- **DOI:** 10.1016/S2214-109X(23)00522-3
- **PMID:** 38245116
- **Vaccine:** TAK-003 (Qdenga, Takeda); 2 subcutaneous doses 3 months apart
- **Design:** Phase 3, double-blind, placebo-controlled, randomised 2:1 (TIDES / DEN-301; NCT02747927)
- **Population:** Ages 4-16y; 8 dengue-endemic countries (Brazil, Colombia, Dominican Republic, Nicaragua, Panama, Philippines, Sri Lanka, Thailand); 26 centres
- **Enrollment:** 20,099 randomly assigned (TAK-003 n=13,401; placebo n=6,698); 20,071 received >=1 dose; 18,257 (91.0%) completed ~4.5y follow-up
- **Cases:** 1,007 virologically confirmed dengue (placebo 560; TAK-003 447); 188 hospitalised (placebo 142; TAK-003 46) of 27,684 febrile illnesses
- **Follow-up:** ~4.5 years after second vaccination (cumulative from first vaccination)
- **Funding:** Takeda Vaccines

### Cumulative vaccine efficacy (overall and seronegative)

| Outcome | Population | VE % | 95% CI |
|---|---|---|---|
| Virologically confirmed dengue | Overall (safety set) | **61.2** | 56.0-65.8 |
| Hospitalised VCD | Overall | **84.1** | 77.8-88.6 |
| Virologically confirmed dengue | Baseline seronegative | **53.5** | 41.6-62.9 |
| Hospitalised VCD | Baseline seronegative | **79.3** | 63.5-88.2 |

### Per-serotype efficacy (qualitative; numeric CIs not in abstract)

- **Seropositive participants:** efficacy demonstrated against all four serotypes (DENV-1, DENV-2, DENV-3, DENV-4)
- **Seronegative participants:** efficacy against **DENV-1 and DENV-2 only**; **NOT against DENV-3**; DENV-4 incidence too low to evaluate

Per-serotype point estimates with CIs are reported in the main paper Tables, not the abstract. Pull from full text or Table 2 if needed for figure.

---

## 3. Kallas EG et al. 2024 — Butantan-DV (Instituto Butantan), single-dose

- **Citation:** N Engl J Med. 2024 Feb 1;390(5):397-408
- **DOI:** 10.1056/NEJMoa2301790
- **PMID:** 38294972
- **Vaccine:** Butantan-Dengue Vaccine (Butantan-DV); single-dose, live, attenuated, tetravalent
- **Design:** Phase 3, double-blind, randomised, placebo-controlled, ongoing (DEN-03-IB; NCT02406729; WHO ICTRP U1111-1168-8679)
- **Population:** Brazil; stratified by age cohorts: 2-6y, 7-17y, 18-59y
- **Enrollment:** 16,235 participants over 3-year enrollment period (Butantan-DV n=10,259; placebo n=5,976)
- **Follow-up:** 2 years per participant (5 years planned); efficacy against symptomatic VCD >28 days post-vaccination
- **Funding:** Instituto Butantan and others

### Vaccine efficacy at 2 years

| Subgroup | VE % | 95% CI |
|---|---|---|
| **Overall** | **79.6** | 70.0-86.3 |
| Baseline seronegative (no previous dengue) | **73.6** | 57.6-83.7 |
| Baseline seropositive (history of exposure) | **89.2** | 77.6-95.6 |
| Age 2-6y | **80.1** | 66.0-88.4 |
| Age 7-17y | **77.8** | 55.6-89.6 |
| Age 18-59y | **90.0** | 68.2-97.5 |

### Per-serotype efficacy

| Serotype | VE % | 95% CI |
|---|---|---|
| **DENV-1** | **89.5** | 78.7-95.0 |
| **DENV-2** | **69.6** | 50.8-81.5 |
| DENV-3 | not detected during follow-up | - |
| DENV-4 | not detected during follow-up | - |

### Safety

Solicited systemic vaccine- or placebo-related adverse events within 21 days: Butantan-DV 58.3% vs placebo 45.6%.

---

## Summary table for Phase 3 retrospective figure

| Vaccine | Trial | Headline VE (overall) | 95% CI | Follow-up | Population | Notes |
|---|---|---|---|---|---|---|
| CYD-TDV (Dengvaxia) | Sridhar 2018 (3-trial reanalysis) | HR 0.32 seropositive 2-16y / HR 1.75 seronegative 2-16y | 0.23-0.45 / 1.14-2.70 | 5 y | 2-16 y | HRs for hospitalization for VCD; not classical VE % in abstract |
| TAK-003 (Qdenga) | Tricou 2024 (DEN-301, 4.5 y) | 61.2% (any VCD) / 84.1% (hospitalised) | 56.0-65.8 / 77.8-88.6 | ~4.5 y post-dose 2 | 4-16 y | Seronegative VE 53.5%; no DENV-3 protection in seronegatives |
| Butantan-DV | Kallas 2024 (DEN-03-IB) | 79.6% (any VCD) | 70.0-86.3 | 2 y | 2-59 y | Single dose; DENV-3/4 not detected during follow-up |

---

## Verification provenance

All numbers extracted verbatim from PubMed via NCBI E-utilities efetch (2026-05-06). Paperclip CLI confirmed PMIDs/DOIs for all three papers in abstract corpus (oa_2807989886, oa_4391074691, oa_4391378947). Full-text PMC versions not available (NEJM and Lancet GH gated).

Cross-check before figure publication: pull main text Tables 2 and 3 from each paper for per-serotype CIs (Tricou) and any updated 5-y data for Kallas/Butantan-DV.
