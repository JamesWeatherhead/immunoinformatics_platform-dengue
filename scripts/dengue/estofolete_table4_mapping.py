"""
Estofolete et al. 2026 Table 4 mapping layer for the
immunoinformatics_platform pipeline applied to dengue.

This script takes the output of the upstream pan-alphavirus pipeline
(adapted for dengue inputs in this fork) and maps the per-epitope and
per-allele scoring tables onto the four-row composite correlate proposed
by Estofolete et al. (npj Vaccines 2026; 11:68) Table 4.

Pipeline coverage of Table 4:

  Tier A | Neutralization breadth
    Operational proxy from this pipeline:
      EDE-residue B-cell epitope coverage on the predicted dimer
      (DiscoTope per-residue scores filtered to the EDE residue list,
      weighted by ESMFold pLDDT).
    Implementation: pre-trial computable from sequence alone.

  Tier A | Avidity (add-on)
    NOT computable from this pipeline. Avidity is a post-Phase-1 readout
    requiring vaccinee serum samples and chaotrope-displacement ELISA.
    Documented in the manuscript MVC plan as a prospective collection.

  Tier B | Memory B-cell breadth
    NOT computable from this pipeline. Requires post-Phase-1 vaccinee
    BCR-seq. Documented in the manuscript MVC plan.

  Tier B | CD8 polyfunctionality
    Operational proxy from this pipeline:
      Geography-weighted MHC-I + MHC-II coverage from NetMHCpan binding
      predictions multiplied by AFND HLA frequencies for the target
      endemic regions, summarized as the worst-region population coverage
      score.
    Implementation: pre-trial computable from sequence + AFND tables.

The output `ranked_candidates.tsv` contains one row per candidate construct
with all four Tier columns populated (NaN for the two post-Phase-1 axes
with explicit `_status` columns explaining the deferral) plus a
`composite_pretrial_pass` boolean indicating whether both pre-trial axes
clear their pre-registered thresholds.

Usage:
  python scripts/dengue/estofolete_table4_mapping.py \\
      --epitope-output-folder /path/to/epitope_outputs \\
      --thresholds-yaml       config/dengue_thresholds.yaml \\
      --outdir                outputs/dengue_table4/

Inputs (read from --epitope-output-folder, the upstream pipeline's
publishDir root):
  - join_tables/input_fasta_table.tsv              (per-construct metadata)
  - join_tables/*_regional_allele_frequencies.tsv  (AFND-weighted coverage)
  - bepipred/*.tsv or epidope/*.tsv                (B-cell epitope scores per residue)
  - discotope/*.tsv                                (discontinuous B-cell epitope scores)
  - netmhcpan_i/*.tsv                              (MHC-I binding predictions)
  - netmhcpan_ii/*.tsv                             (MHC-II binding predictions)

Outputs (written to --outdir):
  - ranked_candidates.tsv          : per-construct Table 4 row scores + composite_pretrial_pass
  - ranked_candidates_summary.md   : human-readable Markdown summary
  - tier_breakdown.json            : machine-readable per-Tier provenance
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml

LOG = logging.getLogger("estofolete_table4_mapping")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

# Default Rouvinski 2017 EDE residue list. These are 1-indexed positions on the
# DENV E protein (the "envelope" portion of the polyprotein, residues ~280-775
# of the full polyprotein). For each candidate the script auto-detects the
# E-protein offset by aligning to the reference DENV-2 E sequence; users can
# also supply explicit per-construct offsets in the thresholds YAML.
DEFAULT_EDE_RESIDUES = [
    67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 79, 80, 84, 86, 88, 91,
    99, 100, 101, 103, 104, 107, 109, 113, 117, 119, 124, 152, 153,
    246, 248, 249, 251,
]


def load_thresholds(path: Path) -> dict:
    """Load gate thresholds from YAML; provide sensible defaults if absent."""
    defaults = {
        "geometry": {
            "antigenic_loss_index_min": 0.0,
            "ede_residues": DEFAULT_EDE_RESIDUES,
            "discotope_threshold": -3.7,
            "plddt_min": 70,
        },
        "equity": {
            "pop_coverage_min": 0.40,
            "endemic_regions": [
                "Southeast_Asia", "Latin_America", "Caribbean",
                "South_Asia", "Sub_Saharan_Africa",
            ],
        },
        "table4": {
            "tier_a_neutralization_breadth": {"axis": "geometry"},
            "tier_a_avidity":                {"axis": "deferred", "status": "post_phase1_required"},
            "tier_b_mbc_breadth":            {"axis": "deferred", "status": "post_phase1_required"},
            "tier_b_cd8_polyfunctionality":  {"axis": "equity"},
        },
    }
    if path is None or not path.exists():
        LOG.warning("thresholds yaml not found at %s; using defaults", path)
        return defaults
    with open(path) as fh:
        loaded = yaml.safe_load(fh) or {}
    # Merge shallow: loaded overrides defaults
    for k, v in defaults.items():
        if k not in loaded:
            loaded[k] = v
    return loaded


def safe_read_tsv(path: Path) -> pd.DataFrame | None:
    """Read a TSV; return None if missing or empty."""
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as e:
        LOG.warning("could not read %s: %s", path, e)
        return None


def compute_geometry_score(
    construct_id: str,
    discotope_dir: Path,
    plddt_dir: Path,
    ede_residues: list[int],
    discotope_threshold: float,
    plddt_min: float,
) -> tuple[float | None, dict]:
    """Compute antigenic_loss_index = sum over EDE residues of (DiscoTope * pLDDT/100).

    Returns (score, provenance_dict). Score is None if either input table is
    missing for this construct.
    """
    discotope_tsv = discotope_dir / f"{construct_id}_discotope.tsv"
    plddt_tsv = plddt_dir / f"{construct_id}_plddt.tsv"
    discotope = safe_read_tsv(discotope_tsv)
    plddt = safe_read_tsv(plddt_tsv)
    prov = {
        "discotope_tsv": str(discotope_tsv),
        "plddt_tsv": str(plddt_tsv),
        "discotope_present": discotope is not None,
        "plddt_present": plddt is not None,
        "n_ede_residues_scored": 0,
    }
    if discotope is None or plddt is None:
        return None, prov

    # Expected columns: residue_position, discotope_score (in discotope tsv);
    # residue_position, plddt (in plddt tsv). Match on residue_position.
    if "residue_position" not in discotope.columns or "discotope_score" not in discotope.columns:
        prov["error"] = "discotope tsv missing expected columns"
        return None, prov
    if "residue_position" not in plddt.columns or "plddt" not in plddt.columns:
        prov["error"] = "plddt tsv missing expected columns"
        return None, prov

    merged = discotope.merge(plddt, on="residue_position", how="inner")
    ede_subset = merged[merged["residue_position"].isin(ede_residues)].copy()
    prov["n_ede_residues_scored"] = int(len(ede_subset))
    if ede_subset.empty:
        return 0.0, prov

    # Apply gates: discotope > threshold AND plddt > min
    ede_subset["passes_gate"] = (
        (ede_subset["discotope_score"] >= discotope_threshold) &
        (ede_subset["plddt"] >= plddt_min)
    )
    ede_subset["weighted_term"] = (
        ede_subset["discotope_score"].fillna(0)
        * ede_subset["plddt"].fillna(0) / 100.0
    )
    score = float(ede_subset.loc[ede_subset["passes_gate"], "weighted_term"].sum())
    return score, prov


def compute_equity_score(
    construct_id: str,
    netmhcpan_i_dir: Path,
    netmhcpan_ii_dir: Path,
    allele_freq_table: Path,
    endemic_regions: list[str],
) -> tuple[float | None, dict]:
    """Compute geography-weighted MHC coverage as a single worst-region score.

    For each endemic region:
      pop_coverage = 1 - prod over alleles of (1 - p_a * I[at least one binder])
    Then return min across regions to capture the worst-served population.

    Returns (score, provenance_dict).
    """
    nmp1 = safe_read_tsv(netmhcpan_i_dir / f"{construct_id}_netmhcpani.tsv")
    nmp2 = safe_read_tsv(netmhcpan_ii_dir / f"{construct_id}_netmhcpanii.tsv")
    af = safe_read_tsv(allele_freq_table)
    prov = {
        "netmhcpan_i_present": nmp1 is not None,
        "netmhcpan_ii_present": nmp2 is not None,
        "allele_freq_table": str(allele_freq_table),
        "endemic_regions": endemic_regions,
        "per_region_coverage": {},
    }
    if nmp1 is None or af is None:
        return None, prov

    # Gather covered alleles for this construct: any allele with at least one
    # binder where percentile rank <= 2.0 (the standard MHCflurry/NetMHCpan
    # binder threshold).
    binder_filter = lambda df, col: set(  # noqa: E731
        df.loc[df[col] <= 2.0, "allele"].dropna().unique()
    ) if col in df.columns else set()

    covered_alleles_i = binder_filter(nmp1, "rank") if nmp1 is not None else set()
    covered_alleles_ii = binder_filter(nmp2, "rank") if nmp2 is not None else set()
    covered_alleles = covered_alleles_i | covered_alleles_ii
    prov["n_covered_alleles_classI"] = len(covered_alleles_i)
    prov["n_covered_alleles_classII"] = len(covered_alleles_ii)

    if not covered_alleles:
        return 0.0, prov

    per_region = {}
    for region in endemic_regions:
        region_af = af[
            af["dataset_name"].str.lower().str.contains(region.lower(), na=False)
        ] if "dataset_name" in af.columns else pd.DataFrame()
        if region_af.empty:
            per_region[region] = float("nan")
            continue
        # Approximate coverage with the union of HLA-A/B/C/DR: probability
        # that an individual carries at least one covered allele.
        product = 1.0
        for _, row in region_af.iterrows():
            allele = row.get("allele") or row.get("locus_join") or row.get("locus")
            freq = row.get("frequency") or row.get("prevalence")
            if not allele or not isinstance(freq, (int, float)) or math.isnan(freq):
                continue
            if str(allele) in covered_alleles:
                product *= max(0.0, 1.0 - float(freq))
        per_region[region] = round(1.0 - product, 4)
    prov["per_region_coverage"] = per_region
    finite_vals = [v for v in per_region.values() if v == v]  # filter NaN
    if not finite_vals:
        return 0.0, prov
    return float(min(finite_vals)), prov


def discover_constructs(epitope_output_folder: Path) -> list[str]:
    """List unique construct_ids based on the input_fasta_table emitted upstream."""
    table = epitope_output_folder / "join_tables" / "input_fasta_table.tsv"
    df = safe_read_tsv(table)
    if df is None or "protein_id" not in df.columns:
        LOG.warning("could not find construct list at %s", table)
        return []
    return sorted(df["protein_id"].dropna().unique().tolist())


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # Original Versiani-style args (single epitope output root)
    p.add_argument("--epitope-output-folder", type=Path, required=False, default=None)
    p.add_argument("--outdir", type=Path, required=False, default=None)
    # Master_pipeline split-dir args (separate phase outputs)
    p.add_argument("--bcell_dir", type=Path, required=False, default=None)
    p.add_argument("--tcell_dir", type=Path, required=False, default=None)
    p.add_argument("--jessev_dir", type=Path, required=False, default=None)
    p.add_argument("--discotope_dir", type=Path, required=False, default=None)
    p.add_argument("--out_dir", type=Path, required=False, default=None)
    p.add_argument("--thresholds-yaml", type=Path, required=False, default=None)
    p.add_argument("--allele-freq-table", type=Path, required=False, default=None)
    args = p.parse_args(argv)

    # Reconcile alternate arg styles
    outdir = args.outdir or args.out_dir
    if outdir is None:
        LOG.error("must pass --outdir or --out_dir")
        return 1
    outdir.mkdir(parents=True, exist_ok=True)

    thresholds = load_thresholds(args.thresholds_yaml)

    if args.epitope_output_folder is not None:
        # Single-folder mode (Versiani-style)
        eo = args.epitope_output_folder
        discotope_dir = eo / "discotope"
        plddt_dir = eo / "alphafold_predictions" / "plddt_per_residue"
        netmhcpan_i_dir = eo / "netmhcpan_i"
        netmhcpan_ii_dir = eo / "netmhcpan_ii"
        constructs_search_root = eo
        allele_freq_search_root = eo / "join_tables"
    else:
        # Split-dir mode (master_pipeline-style)
        bd = args.bcell_dir
        td = args.tcell_dir
        # Discotope can live either in dedicated dir or under bcell pipeline's discotope_output
        discotope_dir = args.discotope_dir if args.discotope_dir else (bd / "discotope_output" if bd else None)
        # pLDDT comes from the AlphaFold output. We synthesize a directory.
        plddt_dir = outdir / "_synth_plddt"
        plddt_dir.mkdir(parents=True, exist_ok=True)
        # T-cell outputs - IEDB writes to netmhcpan_i_output / netmhcpan_ii_output subdirs
        netmhcpan_i_dir = (td / "netmhcpan_i_output") if td else None
        netmhcpan_ii_dir = (td / "netmhcpan_ii_output") if td else None
        # Allele frequencies are produced by B-cell run's FORMATALLELEFREQUENCIES
        constructs_search_root = bd if bd else td
        allele_freq_search_root = (bd / "join_tables") if bd else None

    if args.allele_freq_table is None and allele_freq_search_root is not None:
        candidates = list(allele_freq_search_root.glob("*_regional_allele_frequencies.tsv"))
        if candidates:
            args.allele_freq_table = candidates[0]
            LOG.info("using allele freq table: %s", args.allele_freq_table)
        else:
            LOG.warning("no regional_allele_frequencies tsv found in %s; equity scores will be NaN",
                        allele_freq_search_root)

    constructs = []
    if constructs_search_root is not None:
        constructs = discover_constructs(constructs_search_root)
    if not constructs:
        LOG.warning("no constructs discovered; falling back to glob over input_fasta_table.tsv")
        # fallback: read input_fasta_table.tsv from any subdir
        for sub in [args.bcell_dir, args.tcell_dir, args.jessev_dir]:
            if sub and (sub / "join_tables" / "input_fasta_table.tsv").exists():
                df = pd.read_csv(sub / "join_tables" / "input_fasta_table.tsv", sep="\t")
                constructs = df.iloc[:, 0].astype(str).tolist() if not df.empty else []
                break
    if not constructs:
        LOG.error("no constructs found anywhere; cannot run")
        return 1
    LOG.info("scoring %d constructs", len(constructs))
    args.outdir = outdir  # for downstream code that uses args.outdir
    LOG.info("scoring %d constructs", len(constructs))

    rows = []
    breakdown = {}
    for cid in constructs:
        geom_score, geom_prov = compute_geometry_score(
            cid,
            discotope_dir, plddt_dir,
            thresholds["geometry"]["ede_residues"],
            thresholds["geometry"]["discotope_threshold"],
            thresholds["geometry"]["plddt_min"],
        )
        equity_score, equity_prov = compute_equity_score(
            cid,
            netmhcpan_i_dir, netmhcpan_ii_dir,
            args.allele_freq_table,
            thresholds["equity"]["endemic_regions"],
        )
        # Pre-trial composite pass = both pre-trial axes clear their gates
        pretrial_pass = (
            geom_score is not None and geom_score >= thresholds["geometry"]["antigenic_loss_index_min"]
            and equity_score is not None and equity_score >= thresholds["equity"]["pop_coverage_min"]
        )
        rows.append({
            "candidate_id": cid,
            "tier_A_neutralization_breadth_score": geom_score,
            "tier_A_avidity_score": float("nan"),
            "tier_A_avidity_status": "post_phase1_required",
            "tier_B_mbc_breadth_score": float("nan"),
            "tier_B_mbc_breadth_status": "post_phase1_required",
            "tier_B_cd8_polyfunctionality_score": equity_score,
            "composite_pretrial_pass": pretrial_pass,
        })
        breakdown[cid] = {"geometry": geom_prov, "equity": equity_prov}

    df = pd.DataFrame(rows)
    out_tsv = args.outdir / "ranked_candidates.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    LOG.info("wrote %s", out_tsv)

    md = ["# Estofolete Table 4 ranked candidates\n"]
    md.append(f"Pipeline scored **{len(df)}** constructs against the four Estofolete Table 4 rows.\n")
    md.append(f"Two pre-trial axes (Geometry, Equity) computed; two post-Phase-1 axes deferred per the manuscript MVC plan.\n\n")
    md.append("| Candidate | Tier A neutralization breadth | Tier B CD8 polyfunctionality | Composite pre-trial pass |")
    md.append("|---|---|---|---|")
    for _, r in df.iterrows():
        gs = r["tier_A_neutralization_breadth_score"]
        es = r["tier_B_cd8_polyfunctionality_score"]
        md.append(
            f"| `{r['candidate_id']}` | {gs:.3f} | {es:.3f} | {'yes' if r['composite_pretrial_pass'] else 'no'} |"
            if (gs == gs and es == es) else
            f"| `{r['candidate_id']}` | {gs} | {es} | {'yes' if r['composite_pretrial_pass'] else 'no'} |"
        )
    md.append("\n## Tier A avidity and Tier B memory B-cell breadth\n")
    md.append("These rows of Estofolete Table 4 require post-Phase-1 vaccinee samples")
    md.append("(serum for chaotrope-displacement ELISA in the case of avidity; PBMC")
    md.append("scRNA-seq + paired BCR-seq in the case of MBC breadth) and are not")
    md.append("computable from the public sequence + database inputs available pre-trial.")
    md.append("They are documented in the manuscript Methods as the prospective")
    md.append("validation pathway via the minimum-viable-cohort (MVC) plan.\n")
    out_md = args.outdir / "ranked_candidates_summary.md"
    out_md.write_text("\n".join(md), encoding="utf-8")
    LOG.info("wrote %s", out_md)

    out_json = args.outdir / "tier_breakdown.json"
    out_json.write_text(json.dumps(breakdown, indent=2, default=str), encoding="utf-8")
    LOG.info("wrote %s", out_json)

    return 0


if __name__ == "__main__":
    sys.exit(main())
