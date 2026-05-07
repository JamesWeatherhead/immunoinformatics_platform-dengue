#!/usr/bin/env python3
"""
S5 IEDB validation: precision/recall of NetMHCpan-I strong-binder predictions
(rank <= 2.0) against IEDB-curated experimental class-I dengue epitopes.

For each of the 12 HLA-I alleles in the pipeline output:
  - PREDICTED set: 9-mer peptides with rank <= 2.0
  - GROUND TRUTH set: IEDB epitopes restricted to that allele (4-digit exact
    OR 2-digit prefix match, e.g., HLA-A2 matches HLA-A*02:xx) with
    qualitative_measure starting with "Positive"
  - An IEDB epitope is "recovered" if any of its 9-mer windows is in PREDICTED
  - A predicted 9-mer is a "true positive" if it appears as a 9-mer window
    of any positive IEDB epitope for that allele
  - precision = TP_predictions / |PREDICTED|
  - recall    = recovered_IEDB / |GROUND_TRUTH|
  - F1        = harmonic mean

Outputs: TSV at /tmp/S5_iedb_validation.tsv with per-allele rows + micro/macro.
"""

from __future__ import annotations
import csv
import os
import re
from collections import defaultdict
from statistics import median

PRED_DIR = "/data/james/ica-dengue/outputs/dengue_tcell/netmhcpan_i_output"
IEDB_TSV = "/data/james/ica-dengue/iedb_curated_epitopes/tcell_class_i.tsv"
OUT_TSV  = "/tmp/S5_iedb_validation.tsv"
RANK_CUT = 2.0

# Map from prediction filename allele to canonical "HLA-A*01:01"
PRED_FILES = sorted(os.listdir(PRED_DIR))
PRED_FILES = [f for f in PRED_FILES if f.endswith(".tsv")]

def file_to_allele(fname: str) -> str:
    # netmhcpan_i_HLA_A_01_01_out.tsv -> HLA-A*01:01
    m = re.match(r"netmhcpan_i_HLA_([ABCabc])_(\d{2})_(\d{2})_out\.tsv$", fname)
    assert m, f"Unrecognized prediction filename: {fname}"
    return f"HLA-{m.group(1).upper()}*{m.group(2)}:{m.group(3)}"

def allele_two_digit(canonical: str) -> str:
    # HLA-A*01:01 -> HLA-A1   (drop leading zero per IEDB convention)
    m = re.match(r"HLA-([ABC])\*(\d{2}):\d{2}$", canonical)
    assert m
    locus, two = m.group(1), m.group(2)
    return f"HLA-{locus}{int(two)}"

def kmers(seq: str, k: int = 9) -> set[str]:
    seq = (seq or "").strip().upper()
    if len(seq) < k:
        return set()
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}

# 1. Load IEDB ground truth, indexed by canonical-or-2digit allele key
iedb_by_allele: dict[str, list[str]] = defaultdict(list)
with open(IEDB_TSV, newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        epitope = (row["epitope_sequence"] or "").strip().upper()
        host_hla = (row["host_HLA"] or "").strip()
        # Schema is shifted: the user-described "qualitative_measure" lives
        # in the n_subjects column, e.g., "Positive; ex_vivo".
        qm = (row["n_subjects"] or "").strip()
        if not epitope or not host_hla:
            continue
        if not qm.lower().startswith("positive"):
            continue
        # Only keep peptide-only sequences (drop entries with non-AA chars)
        if not re.fullmatch(r"[A-Z]+", epitope):
            continue
        iedb_by_allele[host_hla].append(epitope)

# Deduplicate
for k in list(iedb_by_allele):
    iedb_by_allele[k] = sorted(set(iedb_by_allele[k]))

# 2. For each pipeline allele, compute precision/recall
rows = []
all_tp = 0
all_pred = 0
all_recovered = 0
all_truth = 0
prec_list, rec_list, f1_list = [], [], []

for fname in PRED_FILES:
    canon = file_to_allele(fname)            # HLA-A*01:01
    two = allele_two_digit(canon)            # HLA-A1
    # Aggregate IEDB ground truth from BOTH the 4-digit and 2-digit keys
    gt_seqs: set[str] = set()
    for key in iedb_by_allele:
        if key == canon or key == two:
            gt_seqs.update(iedb_by_allele[key])
    n_truth = len(gt_seqs)

    # Build the union of all 9-mer windows across positive IEDB epitopes
    truth_kmers: set[str] = set()
    for s in gt_seqs:
        truth_kmers |= kmers(s, 9)

    # Load predictions, filter rank <= 2.0, take peptide column
    pred_peps: set[str] = set()
    with open(os.path.join(PRED_DIR, fname), newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for prow in reader:
            try:
                rank = float(prow["rank"])
            except (TypeError, ValueError):
                continue
            if rank <= RANK_CUT:
                pep = (prow["peptide"] or "").strip().upper()
                if pep:
                    pred_peps.add(pep)
    n_pred = len(pred_peps)

    # True positives = predicted 9-mers that appear as a 9-mer window of any
    # positive IEDB epitope for this allele
    tp_set = pred_peps & truth_kmers
    tp = len(tp_set)

    # Recovered IEDB epitopes = those for which ANY 9-mer window is in pred
    recovered = sum(1 for s in gt_seqs if kmers(s, 9) & pred_peps)

    precision = (tp / n_pred) if n_pred else 0.0
    recall = (recovered / n_truth) if n_truth else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0

    rows.append({
        "allele": canon,
        "n_predicted": n_pred,
        "n_iedb_known": n_truth,
        # n_overlap is reported at epitope granularity (matches recall
        # numerator): the number of IEDB-positive epitopes for which any
        # 9-mer window was predicted as a strong binder. The 9-mer-side
        # true-positive count (recall denominator-side) is tracked
        # internally for micro precision.
        "n_overlap": recovered,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
    })
    all_tp += tp
    all_pred += n_pred
    all_recovered += recovered
    all_truth += n_truth
    prec_list.append(precision)
    rec_list.append(recall)
    f1_list.append(f1)

# 3. Micro and macro aggregates
micro_p = (all_tp / all_pred) if all_pred else 0.0
micro_r = (all_recovered / all_truth) if all_truth else 0.0
micro_f1 = (2 * micro_p * micro_r / (micro_p + micro_r)) if (micro_p + micro_r) else 0.0
macro_p = sum(prec_list) / len(prec_list) if prec_list else 0.0
macro_r = sum(rec_list) / len(rec_list) if rec_list else 0.0
macro_f1 = sum(f1_list) / len(f1_list) if f1_list else 0.0

rows.append({
    "allele": "MICRO_AVG",
    "n_predicted": all_pred,
    "n_iedb_known": all_truth,
    "n_overlap": all_recovered,
    "precision": round(micro_p, 4),
    "recall": round(micro_r, 4),
    "f1": round(micro_f1, 4),
})
rows.append({
    "allele": "MACRO_AVG",
    "n_predicted": "",
    "n_iedb_known": "",
    "n_overlap": "",
    "precision": round(macro_p, 4),
    "recall": round(macro_r, 4),
    "f1": round(macro_f1, 4),
})

# 4. Write TSV
fields = ["allele", "n_predicted", "n_iedb_known", "n_overlap",
          "precision", "recall", "f1"]
with open(OUT_TSV, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                       lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow(r)

# 5. Print summary stats for downstream scripting
median_p = median(prec_list) if prec_list else 0.0
median_r = median(rec_list) if rec_list else 0.0
print(f"N_ALLELES\t{len(prec_list)}")
print(f"MEDIAN_PRECISION\t{median_p:.4f}")
print(f"MEDIAN_RECALL\t{median_r:.4f}")
print(f"MICRO_PRECISION\t{micro_p:.4f}")
print(f"MICRO_RECALL\t{micro_r:.4f}")
print(f"MACRO_PRECISION\t{macro_p:.4f}")
print(f"MACRO_RECALL\t{macro_r:.4f}")
# Sort by F1 to identify top performers
sorted_rows = [r for r in rows
               if r["allele"] not in ("MICRO_AVG", "MACRO_AVG")]
sorted_rows.sort(key=lambda r: r["f1"], reverse=True)
print("TOP_BY_F1")
for r in sorted_rows:
    print(f"  {r['allele']}\tP={r['precision']:.3f}\tR={r['recall']:.3f}\t"
          f"F1={r['f1']:.3f}\tn_iedb={r['n_iedb_known']}\tn_pred={r['n_predicted']}\t"
          f"n_overlap={r['n_overlap']}")
print("LOW_GROUND_TRUTH")
for r in sorted_rows:
    if r["n_iedb_known"] != "" and r["n_iedb_known"] < 10:
        print(f"  {r['allele']}\tn_iedb_known={r['n_iedb_known']}")
print(f"OUT\t{OUT_TSV}")
