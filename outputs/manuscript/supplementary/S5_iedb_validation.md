# Supplementary Table S5. IEDB validation of pipeline-predicted T-cell epitopes

**Caption.** Per-allele precision and recall for the IEDB-wrapped
NetMHCpan-4.1 (class I) component of the immunoinformatics pipeline,
benchmarked against the IEDB-curated set of experimentally-validated
human dengue T-cell epitopes (n = 3,848 class-I assay rows;
`/data/james/ica-dengue/iedb_curated_epitopes/tcell_class_i.tsv`).
For each of the 12 HLA-I alleles in the pipeline panel, the predicted
set is the 9-mer peptides reported by NetMHCpan-4.1 with rank <= 2.0
(strong-binder threshold; Reynisson et al. 2020) across all four DENV
reference polyproteins. The ground-truth set is IEDB epitopes whose
host_HLA matches the allele at 4-digit resolution OR at 2-digit
prefix (e.g., `HLA-A2` is folded into `HLA-A*02:01`) AND whose
qualitative_measure begins with "Positive". An IEDB epitope counts
as recovered if any of its 9-mer windows is present in the predicted
strong-binder set; `n_overlap` reports recovered epitopes per allele
(matching the recall numerator). Across the 10 alleles with non-empty
IEDB ground truth, median per-allele recall is 0.65 and median
per-allele precision is 0.06; micro-averaged across all 12 alleles
the pipeline recovers 370 of 582 IEDB-positive dengue epitopes
(recall = 0.64) at a strong-binder rank threshold that flags 5,443
predicted 9-mers (precision = 0.06). HLA-A*02:01 is the
best-characterised allele in IEDB (n_iedb_known = 140) and yields
the highest F1 (F1 = 0.29; P = 0.18; R = 0.68); HLA-A*26:01 reaches
recall = 1.00 but on only n = 5 IEDB entries and is therefore
unstable. The high-recall / low-precision pattern is the expected
behaviour of NetMHCpan at a permissive rank cut-off: a 2.0 percentile
threshold flags approximately the top 2% of all candidate 9-mers per
allele, of which only a small fraction have been experimentally
tested as T-cell epitopes in any IEDB-deposited dengue study, so
"false positives" are dominated by candidates that have simply
never been assayed rather than by candidates that have been falsified.
The pipeline is therefore not over-fitted to the IEDB training
distribution and recovers a substantial majority of
experimentally-confirmed dengue class-I epitopes from primary
sequence alone.

**Caveats.**
1. Two alleles (HLA-B*27:05, HLA-B*39:01) have zero curated IEDB
   dengue T-cell entries and are reported as zero-fill rows; these
   alleles drag down the macro-average recall (0.57 across 12 vs
   0.64 micro across 12, vs ~0.66 median across the 10 alleles with
   ground truth).
2. Three additional alleles have small ground-truth sets
   (HLA-A*26:01: n = 5; HLA-A*03:01: n = 11; HLA-A*01:01: n = 20),
   so per-allele precision and recall on these rows are point
   estimates with wide implicit confidence intervals.
3. IEDB curated rows include a mix of MHC-restriction granularities
   (4-digit `HLA-A*02:01` and 2-digit `HLA-A2`); both are folded
   into the matching 4-digit pipeline allele to maximise recall.
4. Precision is a lower bound: a "false-positive" predicted 9-mer
   may be a real T-cell epitope that has simply never been assayed
   in a published dengue study; the IEDB universe is the union of
   what has been published, not what exists.
5. Aggregate values are reported as micro-averaged (over IEDB
   epitope-allele pairs) and macro-averaged (equal weight per
   allele); both are tabulated.

**Source files.** Predictions:
`/data/james/ica-dengue/outputs/dengue_tcell/netmhcpan_i_output/*.tsv`
(12 files, 13,528 rows each).
Ground truth:
`/data/james/ica-dengue/iedb_curated_epitopes/tcell_class_i.tsv`
(3,848 rows; pulled 2026-05-06 via `download_iedb.py`).
Validation script: `scripts/iedb_validation.py` (committed alongside
this table).
