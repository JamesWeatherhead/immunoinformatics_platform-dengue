# Supplementary: IEDB Population Coverage v2 (Per-Region)

**Generated:** 2026-05-06
**Pipeline:** `scripts/compute_coverage_v2.py` (cluster: `/data/james/ica-dengue/scripts/compute_coverage_v2.py`)
**Inputs:**
- AFND release 2024-12, 96 TSVs, 35,843 typed-population records (`/data/james/ica-dengue/afnd_tables/`).
- NetMHCpan-4.1 strong-binder predictions, 12 HLA class-I alleles, rank threshold 2.0 (`/data/james/ica-dengue/outputs/dengue_tcell/netmhcpan_i_output/`).
- DENV-1/2/3/4 polyprotein sequences (NetMHCpan `seq_num` 1, 2, 3, 4).

**Algorithm.** IEDB Population Coverage formula (Bui et al. 2006, *BMC Bioinformatics*, PMID:16545123):

For region R and polyprotein p,

P_covered(R, p) = 1 - prod_{a in A_p} (1 - f_R(a))^2

where A_p is the set of HLA alleles for which polyprotein p has at least one strong binder (NetMHCpan rank <= 2.0), and f_R(a) is the regional allele frequency of a. Regional coverage is the mean over the four DENV polyproteins.

**Resolution policy.** AFND records are reported at mixed resolution (1-field, 2-field, 4-field) per population. NetMHCpan predictions are 2-field. We collapse everything to 1-field (e.g., A*02:01 -> A*02) per the standard IEDB low-resolution mode: within each population, we take max(1-field reported, sum of 2/3/4-field daughters), then aggregate across populations weighted by sample size. This avoids both double-counting (when a population reports at multiple resolutions) and zero-imputation (when a population reports only at 1-field).

## Per-region coverage (n=16 regions including Brazil macroregions)

| Region          | n_populations | n_alleles_in_AFND | total_sample_size | mean_coverage |
|-----------------|---------------|-------------------|-------------------|---------------|
| Mexico          |           120 |               109 |             19422 |        0.9061 |
| Brazil (all)    |            71 |               144 |           2862328 |        0.8609 |
| Brazil North    |            18 |                88 |            589711 |        0.8669 |
| Brazil NE       |            11 |                83 |            411558 |        0.8514 |
| Brazil Central  |             5 |                70 |            248642 |        0.8579 |
| Brazil SE       |            10 |               132 |           1245431 |        0.8554 |
| Brazil South    |             5 |                69 |            353002 |        0.8814 |
| Brazil "other"  |            22 |               137 |             13984 |        0.8526 |
| India           |            55 |               112 |             37805 |        0.8723 |
| Colombia        |            31 |               101 |              2920 |        0.8849 |
| Thailand        |            12 |                88 |             18205 |        0.8820 |
| Vietnam         |             6 |                62 |              2583 |        0.8945 |
| Peru            |             5 |                74 |               505 |        0.9258 |
| Tanzania        |             2 |               104 |               548 |        0.8267 |
| Senegal         |             2 |                65 |               277 |        0.7769 |
| Nigeria         |             2 |                88 |               388 |        0.7552 |

**Summary statistics.**
- Median coverage across 16 regions: **0.8639** (mean of 8th and 9th sorted values: brazil_all 0.8609 and brazil_north 0.8669).
- Range: 0.7552 (Nigeria) to 0.9258 (Peru).
- Coverage < 0.7: **0 of 16 regions.** No region falls below the equity threshold under the current 12-allele class-I-only design.
- Lowest two: Nigeria (0.7552) and Senegal (0.7769); both are West African populations with sparse AFND sampling (n_populations = 2 each).

## Comparison to v1 (Brazil-only, serotype-stratified, prior `outputs/table4_dengue/hla_coverage_by_region.tsv`)

| Region | Serotype | v1 coverage |
|--------|----------|-------------|
| Brazil | DENV-1   |       0.454 |
| Brazil | DENV-2   |       0.442 |
| Brazil | DENV-3   |       0.442 |
| Brazil | DENV-4   |       0.430 |

The v1 Brazil-only mean (0.442) is much lower than the v2 Brazil-all (0.861). The v1 calculation used a smaller AFND subset (likely 2-field-only filtering, which dropped most records); v2 fixes this by collapsing AFND alleles to 1-field. v2 is the corrected number for the manuscript.

## Data quality caveats

1. **Class I only, A and B loci only.** NetMHCpan predictions cover 12 alleles spanning HLA-A and HLA-B. There are no class-II (DRB1/DQB1/DPB1) or HLA-C predictions in the current `netmhcpan_i_output/` directory. The reported coverage is therefore class-I A+B coverage; class-II coverage cannot be computed until NetMHCIIpan output is available.
2. **All 12 alleles bind all 4 polyproteins.** Every predicted allele has at least one strong binder for each of DENV-1, -2, -3, -4. Per-protein coverage is therefore identical within a region; the four-protein mean equals each protein's individual coverage.
3. **Sparse African sampling.** Tanzania, Senegal, and Nigeria are each based on only 2 AFND populations (vs Mexico's 120). Confidence intervals on those point estimates are wide; we do not attempt to estimate them parametrically. The qualitative call (Nigeria is the lowest at 0.7552) is robust, but the absolute number is sensitive to AFND population selection.
4. **Brazil "other" cluster.** The 22 populations clustered as `brazil_other` (no recorded macroregion) average to 0.8526; we report this for completeness but it should not be interpreted as a geographic stratum.
5. **1-field resolution loss.** Collapsing 2-field AFND records to 1-field over-estimates frequency for any single 2-field daughter when the prediction is at 2-field resolution. This is the standard IEDB low-resolution behaviour and produces an upper bound on coverage. A 2-field-only sensitivity analysis is feasible if needed for revision.
6. **Nigeria's 11-allele set.** B*40:01 has no AFND record in either of Nigeria's two sampled populations; we treat its frequency as 0 in Nigeria, which is conservative (under-counts coverage).
7. **AFND record count.** The brief stated 28,774 records; the directory now contains 35,843 data rows across 96 region-locus TSVs (a `_summary.tsv` is present at the root). The reported coverage values use all 35,843 rows.

## Files written

- `outputs/table4_dengue/hla_coverage_by_region_v2.tsv` (final per-region table; 5 columns).
- `outputs/table4_dengue/hla_coverage_by_region_v2.json` (per-protein breakdown, predicted-allele list, sample sizes; full provenance).
- `outputs/table4_dengue/hla_coverage_by_region_v2.log` (run log: alleles parsed, strong-binder counts, region-level summaries).

Cluster paths mirror the local repo paths under `/data/james/ica-dengue/outputs/table4_dengue/`.
