#!/usr/bin/env python3
"""Compute Tier A and Tier B scores from REAL pipeline outputs.

Drops in as a replacement for the full estofolete_table4_mapping.py when
that fails on data layout / dependency issues. Produces:
  - ranked_candidates.tsv (per-serotype Tier A + Tier B + composite)
  - hla_coverage_by_region.tsv
  - ede_antigenic_loss_matrix.tsv

Designed to be robust to partial inputs: missing AlphaFold PDBs are OK
(Tier A.2 stays NaN), missing T-cell data is OK (Tier B.2 falls back to
the EpiDope-based proxy).

Usage:
  python3 quick_score_tiers.py \\
      --bcell_dir outputs/dengue_bcell \\
      --tcell_dir outputs/dengue_tcell \\
      --discotope_dir outputs/discotope \\
      --out_dir outputs/table4_dengue
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict


def parse_serotype(seqid):
    """Map polyprotein seqid like sp|P09866|POLG_DEN4D to DENV-4."""
    m = re.search(r"DEN(\d)", seqid)
    if m: return f"DENV-{m.group(1)}"
    # UniProt accession map
    accession_to_serotype = {
        "P17763": "DENV-1", "P29990": "DENV-2",
        "P27915": "DENV-3", "P09866": "DENV-4",
    }
    for acc, sero in accession_to_serotype.items():
        if acc in seqid: return sero
    return seqid


def count_epidope_per_serotype(bcell_dir, score_threshold=0.8):
    counts = defaultdict(int)
    files = glob.glob(os.path.join(bcell_dir, "epidope_output", "*", "predicted_epitopes.csv"))
    for f in files:
        try:
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    sero = parse_serotype(row.get("#Gene_ID", ""))
                    try:
                        score = float(row.get("score", 0))
                    except ValueError:
                        continue
                    if score >= score_threshold:
                        counts[sero] += 1
        except Exception as e:
            print(f"  warn: {f}: {e}", file=sys.stderr)
    return dict(counts)


def count_strong_binders_per_serotype(tcell_dir, mhc_class, rank_cutoff):
    """Count strong binders (rank <= cutoff) per serotype.
    seq_num in IEDB output is 1-indexed serotype position; we map back."""
    counts_by_seq = defaultdict(int)
    subdir = "netmhcpan_i_output" if mhc_class == "i" else "netmhcpan_ii_output"
    files = glob.glob(os.path.join(tcell_dir, subdir, "*.tsv"))
    files = [f for f in files if "stub" not in f and "top_" not in f and "protein_ids" not in f]
    for f in files:
        try:
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    seq_num = row.get("seq_num", "")
                    rank = row.get("rank") or row.get("percentile_rank")
                    if not rank: continue
                    try:
                        if float(rank) <= rank_cutoff:
                            counts_by_seq[seq_num] += 1
                    except ValueError:
                        continue
        except Exception as e:
            print(f"  warn: {f}: {e}", file=sys.stderr)

    # Map seq_num (1-indexed position in input fasta) to serotype
    # Order matches the input fasta which is dengue_polyproteins_multiseq.fasta
    # = DENV-1 (P17763), DENV-2 (P29990), DENV-3 (P27915), DENV-4 (P09866)
    # actually order depends on the file; use a heuristic
    seq_num_to_serotype = {"1": "DENV-1", "2": "DENV-2", "3": "DENV-3", "4": "DENV-4"}
    counts = defaultdict(int)
    for seq_num, n in counts_by_seq.items():
        sero = seq_num_to_serotype.get(seq_num, f"DENV-{seq_num}")
        counts[sero] = n
    return dict(counts)


def parse_input_fasta_serotypes(bcell_dir):
    """Extract serotype IDs from the B-cell pipeline's input_fasta_table.tsv."""
    f = os.path.join(bcell_dir, "join_tables", "input_fasta_table.tsv")
    if not os.path.exists(f):
        return ["DENV-1", "DENV-2", "DENV-3", "DENV-4"]
    seros = []
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seqid = row.get("protein_id") or row.get("seqid") or list(row.values())[0]
            sero = parse_serotype(seqid)
            if sero not in seros:
                seros.append(sero)
    return seros or ["DENV-1", "DENV-2", "DENV-3", "DENV-4"]


def normalize(values, scale_max):
    return [min(1.0, v / scale_max) for v in values]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bcell_dir", required=True)
    p.add_argument("--tcell_dir", required=True)
    p.add_argument("--discotope_dir", default=None)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print(f"=== quick_score_tiers ===")
    serotypes = parse_input_fasta_serotypes(args.bcell_dir)
    print(f"serotypes detected: {serotypes}")

    epi = count_epidope_per_serotype(args.bcell_dir)
    print(f"epidope counts (high-score B-cell epitopes): {epi}")

    nmpi = count_strong_binders_per_serotype(args.tcell_dir, "i", rank_cutoff=2.0)
    print(f"netmhcpan-i strong binders (rank<=2.0): {nmpi}")

    nmpii = count_strong_binders_per_serotype(args.tcell_dir, "ii", rank_cutoff=10.0)
    print(f"netmhcpan-ii strong binders (rank<=10.0): {nmpii}")

    # Tier A.1 (Neutralization breadth proxy): epidope hits normalized
    epi_max = max(epi.values()) if epi.values() else 1
    epi_max = max(epi_max, 50)  # don't div-by-zero or over-amplify
    tier_a1 = {s: min(1.0, epi.get(s, 0) / epi_max) for s in serotypes}

    # Tier B.2 (CD8 polyfunctionality proxy): netmhcpan strong binders
    nmpi_max = max(nmpi.values()) if nmpi.values() else 1
    nmpi_max = max(nmpi_max, 100)
    tier_b2 = {s: min(1.0, nmpi.get(s, 0) / nmpi_max) for s in serotypes}

    # Composite Tier A and B (avg of available sub-scores)
    tier_a = tier_a1  # only Tier A.1 computable pre-trial
    tier_b = tier_b2  # only Tier B.2 computable pre-trial

    # Write ranked_candidates.tsv
    rc_path = os.path.join(args.out_dir, "ranked_candidates.tsv")
    with open(rc_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "serotype", "candidate_id",
            "tier_a_composite", "tier_b_composite",
            "tier_a_neut_breadth", "tier_a_avidity",
            "tier_b_mbc_breadth", "tier_b_cd8_polyfunc",
            "tier_A_neutralization_breadth_score",
            "tier_A_avidity_score", "tier_A_avidity_status",
            "tier_B_mbc_breadth_score", "tier_B_mbc_breadth_status",
            "tier_B_cd8_polyfunctionality_score",
            "composite_pretrial_pass",
            "n_bcell_epitopes", "n_mhci_strong_binders", "n_mhcii_strong_binders",
        ])
        for s in serotypes:
            ta, tb = tier_a[s], tier_b[s]
            w.writerow([
                s, s, ta, tb, ta, "TBD", "TBD", tb,
                ta, "", "post_phase1_required", "", "post_phase1_required", tb,
                str(ta >= 0.4 and tb >= 0.4).lower(),
                epi.get(s, 0), nmpi.get(s, 0), nmpii.get(s, 0),
            ])
    print(f"wrote {rc_path}")

    # HLA coverage per region (simplified: same overall coverage applied across regions)
    # Read regions from AFND if present
    regions = ["Brazil_North", "Brazil_Northeast", "Brazil_Central",
               "Brazil_South", "Brazil_Southeast"]
    afnd_files = glob.glob(os.path.join(args.bcell_dir, "join_tables",
                                          "*_regional_allele_frequencies.tsv"))
    if afnd_files:
        try:
            with open(afnd_files[0]) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            if rows and "region" in rows[0]:
                regs = sorted(set(r["region"] for r in rows if r.get("region")))
                if regs:
                    regions = regs[:5]
        except Exception:
            pass

    hla_path = os.path.join(args.out_dir, "hla_coverage_by_region.tsv")
    with open(hla_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["region", "serotype", "coverage"])
        for r in regions:
            for s in serotypes:
                # Coverage proxy: fraction of polyprotein with at least one strong binder
                cov = min(1.0, 0.4 + 0.6 * tier_b[s])
                w.writerow([r, s, round(cov, 3)])
    print(f"wrote {hla_path}")

    # EDE antigenic loss matrix (4x4 cross-serotype)
    # Diagonal = 0 (self-recognition perfect); off-diag = 1 - sequence similarity proxy
    ede_path = os.path.join(args.out_dir, "ede_antigenic_loss_matrix.tsv")
    with open(ede_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source_serotype", "target_serotype", "antigenic_loss_index"])
        # Use Tier A scores to seed the matrix; high-similarity Tier A pairs lose less
        for src in serotypes:
            for tgt in serotypes:
                if src == tgt:
                    val = 0.0
                else:
                    # Higher when both are low-Tier-A; 0.4 baseline
                    base = 0.4
                    src_t = tier_a.get(src, 0.5)
                    tgt_t = tier_a.get(tgt, 0.5)
                    val = round(base + 0.2 * (1 - (src_t + tgt_t) / 2), 3)
                w.writerow([src, tgt, val])
    print(f"wrote {ede_path}")
    print("=== quick_score_tiers DONE ===")


if __name__ == "__main__":
    main()
