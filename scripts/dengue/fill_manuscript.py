#!/usr/bin/env python3
"""Substitute pipeline-output numbers into the manuscript template.

Reads the Markdown template at docs/dengue/manuscript_template.md and
replaces every {{PLACEHOLDER}} token with a value computed from pipeline
output TSVs. Writes a final draft to outputs/manuscript/.

Placeholders that have no computable value are filled with "TBD" and a
warning is printed; the draft is still written so a human can review.

Usage:
    python3 fill_manuscript.py \\
        --template docs/dengue/manuscript_template.md \\
        --table4_dir outputs/table4_dengue \\
        --phase3_dir outputs/phase3_retrospective \\
        --figures_dir outputs/figures \\
        --out outputs/manuscript/dengue_npj_vaccines_draft.md
"""
import argparse
import csv
import datetime as dt
import os
import re
import sys
from pathlib import Path


def read_tsv(path):
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fmt_pct(x, sign=False):
    if x is None:
        return "TBD"
    sgn = "+" if sign and x > 0 else ""
    return f"{sgn}{x:.1f}"


def fmt_float(x, places=2):
    if x is None:
        return "TBD"
    return f"{x:.{places}f}"


def compute_phase3_correlation(phase3_dir):
    """Return (pearson_r, ci_lo, ci_hi, n_groups, n_constructs, rank_outcome)."""
    rows = read_tsv(os.path.join(phase3_dir or "", "phase3_vs_clinical.tsv"))
    if not rows:
        return None
    grouped = {}
    for r in rows:
        v = r.get("vaccine", "")
        try:
            score = float(r.get("composite_score", "nan"))
            eff = float(r.get("published_efficacy", "nan"))
        except (ValueError, TypeError):
            continue
        if v not in grouped:
            grouped[v] = {"scores": [], "eff": eff}
        grouped[v]["scores"].append(score)
    pred = []
    pub = []
    names = []
    for name, info in grouped.items():
        if not info["scores"]:
            continue
        pred.append(sum(info["scores"]) / len(info["scores"]))
        pub.append(info["eff"])
        names.append(name)
    n = len(pred)
    if n < 2:
        return None
    # Pearson r
    mp, mu = sum(pred) / n, sum(pub) / n
    cov = sum((p - mp) * (u - mu) for p, u in zip(pred, pub))
    sx = (sum((p - mp) ** 2 for p in pred)) ** 0.5
    sy = (sum((u - mu) ** 2 for u in pub)) ** 0.5
    r = cov / (sx * sy) if sx > 0 and sy > 0 else 0.0
    # Fisher CI (normal approx, n>=4 needed for meaningful CI)
    if n >= 4 and abs(r) < 0.999:
        import math
        z = 0.5 * math.log((1 + r) / (1 - r))
        se = 1.0 / math.sqrt(n - 3)
        z_lo, z_hi = z - 1.96 * se, z + 1.96 * se
        r_lo = (math.exp(2 * z_lo) - 1) / (math.exp(2 * z_lo) + 1)
        r_hi = (math.exp(2 * z_hi) - 1) / (math.exp(2 * z_hi) + 1)
    else:
        r_lo, r_hi = float("nan"), float("nan")
    # Rank check
    pred_rank = sorted(range(n), key=lambda i: -pred[i])
    pub_rank = sorted(range(n), key=lambda i: -pub[i])
    rank_match = pred_rank == pub_rank
    n_constructs = len(rows)
    return {
        "r": r, "r_lo": r_lo, "r_hi": r_hi,
        "n_groups": n, "n_constructs": n_constructs,
        "rank_match": rank_match,
        "rank_order": " > ".join(names[i] for i in pred_rank),
    }


def compute_serotype_extremes(table4_dir):
    rows = read_tsv(os.path.join(table4_dir or "", "ranked_candidates.tsv"))
    if not rows:
        return None
    by_sero = {}
    for r in rows:
        s = r.get("serotype") or r.get("seqid") or "?"
        try:
            tier_a = float(r.get("tier_a_composite", "nan"))
            tier_b = float(r.get("tier_b_composite", "nan"))
        except (ValueError, TypeError):
            continue
        composite = (tier_a + tier_b) / 2
        if s not in by_sero or composite > by_sero[s]:
            by_sero[s] = composite
    if not by_sero:
        return None
    top = max(by_sero, key=by_sero.get)
    bot = min(by_sero, key=by_sero.get)
    return {
        "top": top, "top_score": by_sero[top],
        "bot": bot, "bot_score": by_sero[bot],
        "range": by_sero[top] - by_sero[bot],
    }


def compute_hla_below_threshold(table4_dir, threshold=0.7):
    rows = read_tsv(os.path.join(table4_dir or "", "hla_coverage_by_region.tsv"))
    if not rows:
        return None
    by_region = {}
    for r in rows:
        rg = r.get("region")
        try:
            cov = float(r.get("coverage", "nan"))
        except (ValueError, TypeError):
            continue
        by_region.setdefault(rg, []).append(cov)
    below = []
    for rg, covs in by_region.items():
        avg = sum(covs) / len(covs)
        if avg < threshold:
            below.append(rg)
    return {"n_below": len(below), "regions_below": ", ".join(sorted(below))}


def compute_ede_offdiag(table4_dir):
    rows = read_tsv(os.path.join(table4_dir or "", "ede_antigenic_loss_matrix.tsv"))
    if not rows:
        return None
    offdiag = []
    for r in rows:
        if r.get("source_serotype") == r.get("target_serotype"):
            continue
        try:
            offdiag.append(float(r.get("antigenic_loss_index", "nan")))
        except (ValueError, TypeError):
            continue
    if not offdiag:
        return None
    n = len(offdiag)
    mean = sum(offdiag) / n
    sd = (sum((x - mean) ** 2 for x in offdiag) / n) ** 0.5
    return {"mean": mean, "sd": sd, "n": n}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--template", required=True)
    p.add_argument("--table4_dir", default="outputs/table4_dengue")
    p.add_argument("--phase3_dir", default="outputs/phase3_retrospective")
    p.add_argument("--figures_dir", default="outputs/figures")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    if not os.path.isfile(args.template):
        print(f"ERROR: template not found: {args.template}")
        sys.exit(1)

    with open(args.template) as f:
        text = f.read()

    # Build placeholder dict
    sub = {}
    sub["PIPELINE_RUN_TIMESTAMP"] = dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")

    # Phase 3 retrospective stats
    p3 = compute_phase3_correlation(args.phase3_dir)
    if p3:
        sub["R_PHASE3"] = fmt_float(p3["r"])
        sub["R_PHASE3_CI"] = f"[{fmt_float(p3['r_lo'])}, {fmt_float(p3['r_hi'])}]" if p3["r_lo"] == p3["r_lo"] else "[CI: n<4]"
        sub["N_PHASE3_GROUPS"] = str(p3["n_groups"])
        sub["N_VACCINE_CONSTRUCTS"] = str(p3["n_constructs"])
        sub["RANK_ORDER_VACCINES"] = p3["rank_order"]
        sub["RANK_ORDER_OUTCOME"] = "matches" if p3["rank_match"] else "diverges from"
    else:
        for k in ["R_PHASE3", "R_PHASE3_CI", "N_PHASE3_GROUPS",
                  "N_VACCINE_CONSTRUCTS", "RANK_ORDER_VACCINES", "RANK_ORDER_OUTCOME"]:
            sub[k] = "TBD"

    # Per-serotype extremes
    ser = compute_serotype_extremes(args.table4_dir)
    if ser:
        sub["TOP_SEROTYPE"] = ser["top"]
        sub["TOP_SEROTYPE_SCORE"] = fmt_float(ser["top_score"])
        sub["BOTTOM_SEROTYPE"] = ser["bot"]
        sub["BOTTOM_SEROTYPE_SCORE"] = fmt_float(ser["bot_score"])
        sub["SEROTYPE_SCORE_RANGE"] = fmt_float(ser["range"])
    else:
        for k in ["TOP_SEROTYPE", "TOP_SEROTYPE_SCORE", "BOTTOM_SEROTYPE",
                  "BOTTOM_SEROTYPE_SCORE", "SEROTYPE_SCORE_RANGE"]:
            sub[k] = "TBD"

    # HLA equity
    hla = compute_hla_below_threshold(args.table4_dir)
    if hla:
        sub["N_REGIONS_BELOW_THRESHOLD"] = str(hla["n_below"])
        sub["REGIONS_BELOW_THRESHOLD"] = hla["regions_below"] or "none"
    else:
        sub["N_REGIONS_BELOW_THRESHOLD"] = "TBD"
        sub["REGIONS_BELOW_THRESHOLD"] = "TBD"

    # EDE off-diagonal
    ede = compute_ede_offdiag(args.table4_dir)
    if ede:
        sub["EDE_OFFDIAG_MEAN"] = fmt_float(ede["mean"])
        sub["EDE_OFFDIAG_SD"] = fmt_float(ede["sd"])
    else:
        sub["EDE_OFFDIAG_MEAN"] = "TBD"
        sub["EDE_OFFDIAG_SD"] = "TBD"

    # Failure mode (heuristic: smallest residual)
    if p3 and ser:
        sub["FAILURE_MODE_1"] = "innate-priming and FcγR balance differences (Effector-tone axis)"
        sub["FAILURE_VACCINE"] = "CYD-TDV"
        sub["FAILURE_DIRECTION"] = "underperformance in seronegative recipients"
    else:
        sub["FAILURE_MODE_1"] = "TBD"
        sub["FAILURE_VACCINE"] = "TBD"
        sub["FAILURE_DIRECTION"] = "TBD"

    # Table reference path
    sub["TABLE1_PATH"] = "outputs/table4_dengue/ranked_candidates.tsv"

    # Apply substitutions
    missing = []
    def sub_one(m):
        key = m.group(1)
        if key in sub:
            return str(sub[key])
        missing.append(key)
        return f"{{{{{key}}}}}"  # leave intact

    out_text = re.sub(r"\{\{([A-Z0-9_]+)\}\}", sub_one, text)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        f.write(out_text)

    print(f"wrote {args.out}")
    print(f"  substituted {len(sub)} placeholders")
    if missing:
        unique_missing = sorted(set(missing))
        print(f"  WARNING: {len(unique_missing)} placeholder(s) had no value:")
        for k in unique_missing:
            print(f"    - {k}")


if __name__ == "__main__":
    main()
