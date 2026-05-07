"""
Retrospective validation of the dengue extension pipeline against the three
licensed Phase 3 dengue vaccines (CYD-TDV, TAK-003, Butantan-DV).

This is the load-bearing validation that converts the manuscript from
"we built a pipeline" to "we built a pipeline that retrospectively recovers
the relative immunogenicity profiles documented in the Phase 3 record."

Workflow:
  1. The upstream pipeline is run on the 12 vaccine parent strain sequences
     in `ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta`
     producing per-construct epitope scoring tables.
  2. `estofolete_table4_mapping.py` is applied to map those outputs to the
     four-row Estofolete Table 4 schema (two computed pre-trial axes;
     two deferred to MVC).
  3. This script aggregates the per-construct (4 serotypes per vaccine)
     scores into per-vaccine summaries and compares them against published
     pooled Phase 3 efficacy.

Usage:
  python scripts/dengue/retrospective_phase3_validation.py \\
      --epitope-output-folder /path/to/epitope_outputs/dengue_phase3 \\
      --thresholds-yaml       config/dengue_thresholds.yaml \\
      --outdir                outputs/retrospective_phase3/

Inputs:
  - --epitope-output-folder must contain the upstream pipeline outputs from a
    run on `ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta`
  - construct_id naming convention assumed: "{vaccine}_DENV{n}_*" matching
    the FASTA records produced by `scripts/dengue/fetch_dengue_data.sh`

Outputs (in --outdir):
  - phase3_table4_per_construct.tsv : 12 rows (3 vaccines x 4 serotypes), Table 4 cols
  - phase3_per_vaccine_summary.tsv  : 3 rows aggregating across serotypes
  - phase3_vs_clinical.md           : Markdown comparison table
  - phase3_correlation.json         : descriptive Spearman + Pearson correlations
                                      between pipeline composite and published efficacy
                                      (n=3, descriptive only)

Citations for published efficacy values used in the comparison:
  Capeding 2014 Lancet 384:1358        (CYD-TDV pooled efficacy 56.5%)
  Sridhar 2018 NEJM 379:327            (CYD-TDV serostatus stratified)
  Tricou 2024 Lancet Glob Health 12:e257 (TAK-003 4.5y 61.2%)
  Kallas 2024 NEJM 390:397             (Butantan-DV single-dose 79.6%)
  Nogueira 2024 Lancet Infect Dis 24    (Butantan-DV extended follow-up)
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd
from scipy import stats

LOG = logging.getLogger("retrospective_phase3_validation")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

VACCINES = ("CYD-TDV", "TAK-003", "Butantan-DV")
SEROTYPES = (1, 2, 3, 4)

PUBLISHED_PER_SEROTYPE_EFFICACY = {
    "CYD-TDV":     {1: 50.3, 2: 42.3, 3: 74.0, 4: 77.7},
    "TAK-003":     {1: 73.5, 2: 97.7, 3: 62.6, 4: float("nan")},
    "Butantan-DV": {1: 89.5, 2: 69.6, 3: float("nan"), 4: float("nan")},
}
PUBLISHED_POOLED_EFFICACY = {
    "CYD-TDV":     56.5,
    "TAK-003":     61.2,
    "Butantan-DV": 79.6,
}


def parse_construct_id(cid: str) -> tuple[str, int] | None:
    """Match '{VACCINE}_DENV{N}_...' construct_ids."""
    for v in VACCINES:
        m = re.match(rf"^{re.escape(v)}_DENV([1-4])(?:_|$)", cid)
        if m:
            return v, int(m.group(1))
    return None


def run_table4_mapping(
    epitope_output_folder: Path,
    thresholds_yaml: Path | None,
    outdir: Path,
) -> Path:
    """Invoke the Estofolete Table 4 mapping script and return path to its TSV."""
    map_outdir = outdir / "table4_mapping"
    map_outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(Path(__file__).parent / "estofolete_table4_mapping.py"),
        "--epitope-output-folder", str(epitope_output_folder),
        "--outdir", str(map_outdir),
    ]
    if thresholds_yaml is not None and thresholds_yaml.exists():
        cmd.extend(["--thresholds-yaml", str(thresholds_yaml)])
    LOG.info("invoking: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return map_outdir / "ranked_candidates.tsv"


def aggregate_per_vaccine(per_construct: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for v in VACCINES:
        sub = per_construct[per_construct["vaccine"] == v]
        rows.append({
            "vaccine": v,
            "n_serotypes_scored": int(len(sub)),
            "tierA_mean_neutralization_breadth": float(sub["tier_A_neutralization_breadth_score"].mean()) if not sub.empty else float("nan"),
            "tierB_mean_cd8_polyfunctionality": float(sub["tier_B_cd8_polyfunctionality_score"].mean()) if not sub.empty else float("nan"),
            "n_pretrial_pass": int(sub["composite_pretrial_pass"].sum()) if not sub.empty else 0,
        })
    df = pd.DataFrame(rows)
    # Rank-normalize per axis to handle scale incompatibility
    for col in ("tierA_mean_neutralization_breadth", "tierB_mean_cd8_polyfunctionality"):
        df[f"{col}_rank01"] = df[col].rank(method="average") / max(len(df), 1)
    df["composite_score_rank01"] = (
        df["tierA_mean_neutralization_breadth_rank01"]
        + df["tierB_mean_cd8_polyfunctionality_rank01"]
    ) / 2
    return df


def comparison_md(per_vaccine: pd.DataFrame) -> str:
    md = []
    md.append("# Retrospective Phase 3 validation\n")
    md.append(
        "Pipeline-derived pre-trial composite scores (Geometry + Equity axes only; "
        "Tier A avidity and Tier B MBC breadth are post-Phase-1 and not scored here) "
        "compared against published per-serotype Phase 3 efficacy.\n"
    )
    md.append("\n## Per-vaccine summary\n")
    md.append("| Vaccine | n serotypes scored | Mean Tier A (neutralization breadth) | Mean Tier B (CD8 polyfunctionality) | Composite (rank-01) | Pooled Phase 3 efficacy |")
    md.append("|---|---|---|---|---|---|")
    pooled_str = {
        "CYD-TDV": "56.5% (Capeding 2014)",
        "TAK-003": "61.2% at 4.5y (Tricou 2024)",
        "Butantan-DV": "79.6% single-dose (Kallas 2024)",
    }
    for _, r in per_vaccine.iterrows():
        md.append(
            f"| {r['vaccine']} | {r['n_serotypes_scored']} "
            f"| {r['tierA_mean_neutralization_breadth']:.3f} "
            f"| {r['tierB_mean_cd8_polyfunctionality']:.3f} "
            f"| {r['composite_score_rank01']:.3f} "
            f"| {pooled_str.get(r['vaccine'], 'NA')} |"
        )
    md.append("\n## Per-serotype published efficacy (for reference)\n")
    md.append("| Vaccine | DENV-1 | DENV-2 | DENV-3 | DENV-4 |")
    md.append("|---|---|---|---|---|")
    for v in VACCINES:
        eff = PUBLISHED_PER_SEROTYPE_EFFICACY[v]
        cells = []
        for s in SEROTYPES:
            val = eff[s]
            cells.append(f"{val:.1f}%" if val == val else "NA")
        md.append(f"| {v} | " + " | ".join(cells) + " |")
    md.append("\n## Interpretation\n")
    md.append(
        "A monotonic relationship between the pipeline composite and the published "
        "pooled efficacy would constitute retrospective validation that the two "
        "pre-trial axes recover the relative immunogenicity ordering of the licensed "
        "Phase 3 vaccines from their parent strain sequences alone, prior to any "
        "clinical evaluation.\n"
    )
    md.append(
        "A non-monotonic relationship would constitute an honest negative result that "
        "reframes Geometry+Equity as necessary but not sufficient predictors of vaccine "
        "efficacy and would strengthen the case for the post-Phase-1 axes (Tier A "
        "avidity, Tier B MBC breadth) addressed prospectively in the MVC plan.\n"
    )
    md.append(
        "\nNote: n=3 vaccines is too small for inferential statistics; correlations "
        "below are descriptive only.\n"
    )
    return "\n".join(md)


def correlation_to_efficacy(per_vaccine: pd.DataFrame) -> dict:
    sub = per_vaccine.set_index("vaccine").loc[list(PUBLISHED_POOLED_EFFICACY.keys())]
    composite = sub["composite_score_rank01"].astype(float).values
    efficacy = [PUBLISHED_POOLED_EFFICACY[v] for v in PUBLISHED_POOLED_EFFICACY.keys()]
    rho, p_rho = stats.spearmanr(composite, efficacy)
    r, p_r = stats.pearsonr(composite, efficacy)
    return {
        "n_vaccines": len(efficacy),
        "spearman_rho": float(rho) if rho == rho else None,
        "spearman_p": float(p_rho) if p_rho == p_rho else None,
        "pearson_r": float(r) if r == r else None,
        "pearson_p": float(p_r) if p_r == p_r else None,
        "pipeline_composite_rank01": dict(zip(PUBLISHED_POOLED_EFFICACY.keys(), [float(x) if x == x else None for x in composite])),
        "published_pooled_efficacy_pct": PUBLISHED_POOLED_EFFICACY,
        "warning": "n=3 vaccines is too small for inferential statistics; descriptive only.",
    }


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--epitope-output-folder", type=Path, required=False, default=None,
                   help="upstream pipeline publishDir root for the Phase 3 vaccine run")
    p.add_argument("--pipeline_dir", type=Path, required=False, default=None,
                   help="alias for --epitope-output-folder (master_pipeline-style)")
    p.add_argument("--thresholds-yaml", type=Path, required=False, default=None,
                   help="optional thresholds YAML; defaults used if absent")
    p.add_argument("--outdir", type=Path, required=False, default=None,
                   help="output directory (created if needed)")
    p.add_argument("--out_dir", type=Path, required=False, default=None,
                   help="alias for --outdir (master_pipeline-style)")
    args = p.parse_args(argv)

    epitope_folder = args.epitope_output_folder or args.pipeline_dir
    outdir = args.outdir or args.out_dir
    if epitope_folder is None or outdir is None:
        LOG.error("must pass either (--epitope-output-folder + --outdir) or (--pipeline_dir + --out_dir)")
        return 1
    args.outdir = outdir
    outdir.mkdir(parents=True, exist_ok=True)
    table4_tsv = run_table4_mapping(epitope_folder, args.thresholds_yaml, outdir)
    LOG.info("loaded table4 mapping output: %s", table4_tsv)

    df = pd.read_csv(table4_tsv, sep="\t")
    parsed = df["candidate_id"].apply(parse_construct_id)
    df["vaccine"] = parsed.apply(lambda x: x[0] if x else None)
    df["serotype"] = parsed.apply(lambda x: x[1] if x else None)
    df = df[df["vaccine"].notna()].copy()
    if df.empty:
        LOG.error("no construct_ids matched the expected '{vaccine}_DENV{n}_...' pattern")
        return 1

    per_construct_path = args.outdir / "phase3_table4_per_construct.tsv"
    df.to_csv(per_construct_path, sep="\t", index=False)
    LOG.info("wrote %s (%d rows)", per_construct_path, len(df))

    per_vaccine = aggregate_per_vaccine(df)
    per_vaccine_path = args.outdir / "phase3_per_vaccine_summary.tsv"
    per_vaccine.to_csv(per_vaccine_path, sep="\t", index=False)
    LOG.info("wrote %s", per_vaccine_path)

    md_path = args.outdir / "phase3_vs_clinical.md"
    md_path.write_text(comparison_md(per_vaccine), encoding="utf-8")
    LOG.info("wrote %s", md_path)

    # Also write phase3_vs_clinical.tsv that fill_manuscript.py consumes.
    # Expected columns: vaccine, composite_score, published_efficacy
    fm_rows = []
    for _, r in per_vaccine.iterrows():
        v = r["vaccine"]
        fm_rows.append({
            "vaccine": v,
            "composite_score": r["composite_score_rank01"],
            "published_efficacy": PUBLISHED_POOLED_EFFICACY.get(v, float("nan")),
        })
    fm_df = pd.DataFrame(fm_rows)
    fm_path = args.outdir / "phase3_vs_clinical.tsv"
    fm_df.to_csv(fm_path, sep="\t", index=False)
    LOG.info("wrote %s (consumed by fill_manuscript.py)", fm_path)

    corr = correlation_to_efficacy(per_vaccine)
    corr_path = args.outdir / "phase3_correlation.json"
    corr_path.write_text(json.dumps(corr, indent=2), encoding="utf-8")
    LOG.info("wrote %s", corr_path)
    if corr["spearman_rho"] is not None:
        LOG.info("Spearman rho = %.3f (p=%.3f); n=%d (descriptive only)",
                 corr["spearman_rho"], corr["spearman_p"], corr["n_vaccines"])

    return 0


if __name__ == "__main__":
    sys.exit(main())
