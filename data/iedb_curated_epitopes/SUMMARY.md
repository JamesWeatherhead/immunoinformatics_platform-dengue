# IEDB curated dengue epitopes
_Downloaded: 2026-05-07T15:04:33.653443Z from https://query-api.iedb.org_

Filter: parent_source_antigen_source_org_name ILIKE *Dengue*
(captures DENV-1, DENV-2, DENV-3, DENV-4, and "Dengue virus" generic)

| File | Endpoint | Rows |
|---|---|---|
| tcell_class_i.tsv | tcell_search?mhc_class=eq.I | 3848 |
| tcell_class_ii.tsv | tcell_search?mhc_class=eq.II | 6090 |
| bcell.tsv | bcell_search | 16275 |

## T-cell class I (3848 assay rows)

**Top 10 HLA alleles** (30 unique 4-digit):
  - HLA-B*07:02: 618
  - HLA-A*11:01: 497
  - HLA-A*02:01: 488
  - HLA-A*01:01: 287
  - HLA-B*35:01: 269
  - HLA-A*24:02: 153
  - HLA-B*08:01: 109
  - HLA-B*40:01: 84
  - HLA-B*58:01: 67
  - HLA-B*15:01: 50

**Top 10 source proteins** (53 unique):
  - polyprotein: 2653
  - Genome polyprotein: 322
  - Nonstructural protein NS3: 131
  - polyprotein [dengue virus type 4]: 70
  - polyprotein precursor: 69
  - viral polyprotein: 64
  - polyprotein [dengue virus type 3]: 60
  - nonstructural protein 3: 55
  - envelope protein: 52
  - NS3 protein: 50

**Top 5 PMIDs** (72 unique papers):
  - PMID 25056881: 363 rows
  - PMID 24190657: 359 rows
  - PMID 40926371: 240 rows
  - PMID 20421879: 144 rows
  - PMID 25320311: 109 rows

## T-cell class II (6090 assay rows)

**Top 10 HLA alleles** (38 unique 4-digit):
  - HLA-DRB1*07:01: 311
  - HLA-DRB1*04:01: 308
  - HLA-DRB1*01:01: 252
  - HLA-DRB1*08:02: 240
  - HLA-DRB1*15:01: 230
  - HLA-DRB1*11:01: 223
  - HLA-DRB1*03:01: 212
  - HLA-DRB1*09:01: 207
  - HLA-DRB1*13:01: 199
  - HLA-DRB1*04:03: 179

**Top 10 source proteins** (82 unique):
  - polyprotein: 2114
  - polyprotein, partial: 970
  - Genome polyprotein: 433
  - polyprotein [Dengue virus 3]: 295
  - polyprotein [Dengue virus 1]: 292
  - polyprotein [Dengue virus 4]: 199
  - envelope protein, partial: 161
  - envelope protein: 158
  - envelope glycoprotein: 141
  - Nonstructural protein NS3: 122

**Top 5 PMIDs** (41 unique papers):
  - PMID 29081779: 1931 rows
  - PMID 31333679: 1442 rows
  - PMID 26195744: 936 rows
  - PMID 24130917: 445 rows
  - PMID 27974563: 280 rows

## B-cell (16275 assay rows)


**Top 10 source proteins** (69 unique):
  - Genome polyprotein: 4541
  - polyprotein [Dengue virus 3]: 1806
  - polyprotein [Dengue virus 4]: 1300
  - polyprotein, partial [Dengue virus 2]: 1296
  - polyprotein: 1292
  - polyprotein [Dengue virus 1]: 1035
  - envelope protein: 1014
  - polyprotein [Dengue virus 2]: 846
  - envelope glycoprotein: 472
  - Envelope protein: 289

**Top 5 PMIDs** (163 unique papers):
  - PMID 34290707: 10386 rows
  - PMID 34732126: 600 rows
  - PMID 40753572: 296 rows
  - PMID 36513793: 282 rows
  - PMID 36560438: 235 rows
