#!/usr/bin/env python3
"""
Recompute IEDB Population Coverage (Bui 2006, PMID:16545123) from AFND
allele frequencies and netMHCpan strong-binder predictions, per region.

Algorithm:
  For region R and protein p:
    P_covered(R, p) = 1 - prod_{a in A_p} (1 - f_R(a))^2
  where A_p = alleles for which protein p has at least 1 strong binder
  (rank <= 2.0). Final per-region coverage = mean over the 4 DENV
  polyproteins (seq_num 1..4).

Resolution policy (IEDB-style "low resolution" mode):
  Predictions are 2-field (A*01:01). AFND records are mixed: 1-field
  (A*01), 2-field (A*01:01), 3-field, 4-field. We collapse everything
  to 1-field by summing 2/3/4-field daughter frequencies up into the
  parent and using max(reported_one_field, sum_of_daughters) per
  population to avoid double-counting where a population reported
  at multiple resolutions. We then compute coverage at 1-field, which
  is the pragmatic IEDB approach when prediction sets and frequency
  tables disagree on resolution.

Inputs:
  AFND TSVs:   /data/james/ica-dengue/afnd_tables/*.tsv
  netMHCpan:   /data/james/ica-dengue/outputs/dengue_tcell/netmhcpan_i_output/*.tsv

Outputs:
  TSV:    /data/james/ica-dengue/outputs/table4_dengue/hla_coverage_by_region_v2.tsv
  JSON:   /data/james/ica-dengue/outputs/table4_dengue/hla_coverage_by_region_v2.json
  LOG:    /data/james/ica-dengue/outputs/table4_dengue/hla_coverage_by_region_v2.log
"""
from __future__ import annotations
import csv
import glob
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

AFND_DIR = "/data/james/ica-dengue/afnd_tables"
NETMHC_DIR = "/data/james/ica-dengue/outputs/dengue_tcell/netmhcpan_i_output"
OUT_DIR = "/data/james/ica-dengue/outputs/table4_dengue"
RANK_THRESHOLD = 2.0

REGIONS_REPORT = [
    "mexico",
    "brazil_all",
    "brazil_north",
    "brazil_ne",
    "brazil_central",
    "brazil_se",
    "brazil_south",
    "brazil_other",
    "india",
    "colombia",
    "thailand",
    "vietnam",
    "peru",
    "tanzania",
    "senegal",
    "nigeria",
]

os.makedirs(OUT_DIR, exist_ok=True)
log_lines: list[str] = []
def log(msg: str) -> None:
    print(msg, flush=True)
    log_lines.append(msg)


# 1. Map netMHCpan filename allele to AFND-style allele key (e.g. "A*01:01").
ALLELE_FROM_FILENAME = re.compile(
    r"netmhcpan_i_HLA_([A-Z])_(\d{2})_(\d{2})_out\.tsv"
)


def filename_to_allele_key(name: str) -> str | None:
    m = ALLELE_FROM_FILENAME.match(name)
    if not m:
        return None
    locus, field1, _field2 = m.groups()
    # Collapse to 1-field for AFND interoperability
    return f"{locus}*{field1}"


def to_one_field(allele: str) -> str | None:
    """Collapse any resolution allele to 1-field (e.g., A*01:01:02 -> A*01).
    Strip suffixes like N, L, Q, S. Returns None if unparseable."""
    if not allele or "*" not in allele:
        return None
    locus, fields = allele.split("*", 1)
    first = fields.split(":", 1)[0]
    # Strip any non-digit suffix (e.g., 02:15N at 2-field, but here just 1-field)
    m = re.match(r"^(\d+)", first)
    if not m:
        return None
    return f"{locus}*{m.group(1).zfill(2)}"


# 2. Determine which (allele, polyprotein) pairs have strong binders.
# binders[allele][seq_num] = True iff at least one row with rank <= 2.0 exists.
binders: dict[str, set[int]] = defaultdict(set)
total_strong = 0
total_rows = 0
allele_files = sorted(glob.glob(os.path.join(NETMHC_DIR, "netmhcpan_i_HLA_*_out.tsv")))
for fp in allele_files:
    base = os.path.basename(fp)
    allele_key = filename_to_allele_key(base)
    if allele_key is None:
        log(f"WARN: cannot parse allele from filename {base}")
        continue
    with open(fp, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            total_rows += 1
            try:
                rank = float(row["rank"])
                seq_num = int(row["seq_num"])
            except (KeyError, ValueError):
                continue
            if rank <= RANK_THRESHOLD:
                binders[allele_key].add(seq_num)
                total_strong += 1
log(f"netMHCpan: parsed {len(allele_files)} alleles, {total_rows} rows, "
    f"{total_strong} strong binders (rank<={RANK_THRESHOLD})")
log("Per-allele protein coverage (seq_nums with at least 1 strong binder):")
for a in sorted(binders):
    log(f"  {a}: seq_nums={sorted(binders[a])}")

PROTEINS = sorted({sn for sns in binders.values() for sn in sns})
log(f"DENV polyproteins observed: {PROTEINS}")
if not PROTEINS:
    log("FATAL: no proteins / strong binders parsed")
    sys.exit(1)


# 3. Build per-region allele frequencies from AFND TSVs.
# AFND row provides: locus, region, country, allele (e.g. "A*01"), population,
#   pop_id, pct_individuals, allele_frequency, sample_size, source.
# Aggregate across populations within a region: weighted mean of
#   allele_frequency by sample_size, restricted to rows with usable freq+N.
# Then group A*01:01 vs A*01: AFND frequencies are typically reported at
# 2-field resolution (A*01:01) or 1-field (A*01). The 12 netMHCpan alleles
# are 2-field. For each predicted allele a (e.g. A*01:01), we look up the
# matching 2-field AFND record. If only 1-field rows exist for that locus
# for the region, we fall back to summing the 2-field daughter frequencies
# (here we conservatively skip; coverage is then 0 contribution from that
# allele -- this is the correct behaviour because we only have predictions
# for that 2-field allele).

def aggregate_region(region: str) -> dict:
    """Build per-population, 1-field allele frequency tables; then
    aggregate to a regional sample-size-weighted mean.

    Per-population resolution policy:
      For each (population, 1-field allele): collect all reported rows
      for that 1-field bucket. If the population reported at 2/3/4-field,
      sum daughter frequencies. If it ALSO reported at 1-field, take the
      max (so we never double-count or lose information). Then aggregate
      across populations weighted by sample_size.
    """
    files = sorted(glob.glob(os.path.join(AFND_DIR, f"{region}_*.tsv")))
    if not files:
        return {"allele_freqs": {}, "populations": 0, "n_alleles": 0,
                "total_sample_size": 0, "loci_present": []}
    # raw[pop_id][1field allele] = {"one_field_freq": float|None,
    #                                "daughter_sum": float,
    #                                "N": int}
    raw: dict[str, dict[str, dict]] = defaultdict(lambda: defaultdict(
        lambda: {"one_field_freq": None, "daughter_sum": 0.0, "N": 0}
    ))
    loci_present: set[str] = set()
    n_total_records = 0
    for fp in files:
        with open(fp, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for row in rd:
                n_total_records += 1
                locus = (row.get("locus") or "").strip()
                allele = (row.get("allele") or "").strip()
                pop_id = (row.get("pop_id") or "").strip()
                if not locus or not allele or not pop_id:
                    continue
                loci_present.add(locus)
                try:
                    f = float(row.get("allele_frequency") or "")
                except ValueError:
                    continue
                try:
                    n = int(float(row.get("sample_size") or 0))
                except ValueError:
                    n = 0
                if not (0.0 <= f <= 1.0) or n <= 0:
                    continue
                one_field = to_one_field(allele)
                if one_field is None:
                    continue
                slot = raw[pop_id][one_field]
                # update sample size to the max seen (same population may
                # have varying N in different rows; sample_size is per
                # locus-typing study, so should be consistent within
                # locus-population).
                if n > slot["N"]:
                    slot["N"] = n
                # determine if this row is 1-field or higher
                is_one_field = ":" not in allele
                if is_one_field:
                    cur = slot["one_field_freq"]
                    slot["one_field_freq"] = f if cur is None else max(cur, f)
                else:
                    slot["daughter_sum"] += f

    # Resolve per population: per-1-field allele freq = max(one_field_freq,
    # daughter_sum). Then aggregate across populations.
    bucket: dict[str, list[tuple[float, int]]] = defaultdict(list)
    pop_ids = set(raw.keys())
    for pop_id, alleles in raw.items():
        for allele_1f, slot in alleles.items():
            f1 = slot["one_field_freq"]
            ds = slot["daughter_sum"]
            f = max(f1 if f1 is not None else 0.0, ds)
            n = slot["N"]
            if f <= 0 or n <= 0:
                continue
            # cap at 1.0 (additive sums can rarely exceed due to typing
            # ambiguity / rounding)
            f = min(f, 1.0)
            bucket[allele_1f].append((f, n))

    allele_freqs: dict[str, float] = {}
    sample_sizes_by_pop: dict[str, int] = {}
    for pop_id, alleles in raw.items():
        # use max N across loci within population as an approximation
        if alleles:
            sample_sizes_by_pop[pop_id] = max(
                slot["N"] for slot in alleles.values()
            )
    total_N = sum(sample_sizes_by_pop.values())
    for allele, items in bucket.items():
        wsum = sum(f * n for f, n in items)
        nsum = sum(n for _, n in items)
        if nsum > 0:
            allele_freqs[allele] = wsum / nsum

    return {
        "allele_freqs": allele_freqs,
        "populations": len(pop_ids),
        "n_alleles": len(allele_freqs),
        "total_sample_size": int(total_N),
        "loci_present": sorted(loci_present),
        "n_records": n_total_records,
    }


# 4. Compute coverage.
def coverage_for_region(region_summary: dict) -> dict:
    af = region_summary["allele_freqs"]
    per_protein: dict[int, float] = {}
    per_protein_alleles: dict[int, list[str]] = {}
    for sn in PROTEINS:
        # alleles which (a) have predictions with strong binders for protein sn
        # AND (b) are present in this region's AFND frequencies.
        alleles = []
        prod = 1.0
        for a, sns in binders.items():
            if sn not in sns:
                continue
            f = af.get(a, 0.0)
            if f <= 0:
                continue
            alleles.append((a, f))
            prod *= (1.0 - f) ** 2
        per_protein[sn] = 1.0 - prod
        per_protein_alleles[sn] = sorted(a for a, _ in alleles)
    mean_cov = sum(per_protein.values()) / len(per_protein)
    return {
        "per_protein_coverage": per_protein,
        "per_protein_alleles_used": per_protein_alleles,
        "mean_coverage_dengue_polyproteins": mean_cov,
    }


# 5. Run for all reported regions plus also any region present on disk.
all_present = sorted({
    os.path.basename(fp).rsplit("_", 1)[0]
    for fp in glob.glob(os.path.join(AFND_DIR, "*.tsv"))
})
log(f"AFND regions on disk: {all_present}")

results = []
for region in REGIONS_REPORT:
    if region not in all_present:
        log(f"WARN: region {region} not found in AFND tables; skipping")
        continue
    summ = aggregate_region(region)
    cov = coverage_for_region(summ)
    results.append({
        "region": region,
        "n_populations_aggregated": summ["populations"],
        "n_alleles": summ["n_alleles"],
        "total_sample_size": summ["total_sample_size"],
        "loci_present": summ["loci_present"],
        "n_records": summ.get("n_records", 0),
        "mean_coverage_dengue_polyproteins": cov["mean_coverage_dengue_polyproteins"],
        "per_protein_coverage": cov["per_protein_coverage"],
        "predicted_alleles_with_freq_in_region": sorted(
            a for a in binders if a in summ["allele_freqs"]
        ),
    })
    log(
        f"{region}: pops={summ['populations']} alleles={summ['n_alleles']} "
        f"loci={summ['loci_present']} N={summ['total_sample_size']} "
        f"coverage={cov['mean_coverage_dengue_polyproteins']:.4f}"
    )


# 6. Write outputs.
out_tsv = os.path.join(OUT_DIR, "hla_coverage_by_region_v2.tsv")
with open(out_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow([
        "region",
        "n_populations_aggregated",
        "n_alleles",
        "total_sample_size",
        "mean_coverage_dengue_polyproteins",
    ])
    for r in results:
        w.writerow([
            r["region"],
            r["n_populations_aggregated"],
            r["n_alleles"],
            r["total_sample_size"],
            f"{r['mean_coverage_dengue_polyproteins']:.4f}",
        ])
log(f"WROTE {out_tsv}")

out_json = os.path.join(OUT_DIR, "hla_coverage_by_region_v2.json")
with open(out_json, "w") as fh:
    json.dump({
        "rank_threshold": RANK_THRESHOLD,
        "predicted_alleles": sorted(binders),
        "proteins_seq_nums": PROTEINS,
        "results": results,
    }, fh, indent=2, default=str)
log(f"WROTE {out_json}")

out_log = os.path.join(OUT_DIR, "hla_coverage_by_region_v2.log")
with open(out_log, "w") as fh:
    fh.write("\n".join(log_lines) + "\n")
print(f"WROTE {out_log}")
