# Changelog

All notable changes to this fork are documented here. The upstream
[`pmccaffrey6/immunoinformatics_platform`](https://github.com/pmccaffrey6/immunoinformatics_platform)
maintains its own history.

## [v1.0-dengue-prep] - 2026-05-06

### Added
- DENV 1-4 reference polyproteins (UniProt P17763, P29990, P27915, P09866) in `ref_fastas/dengue/`
- 12 Phase 3 vaccine parent strain CDS sequences from GenBank in `ref_fastas/dengue_phase3_vaccines/` covering CYD-TDV / Dengvaxia, TAK-003 / Qdenga, and Butantan-DV / TV003
- `scripts/dengue/fetch_dengue_data.sh` to pull the above from UniProt + GenBank with full provenance documentation
- `scripts/dengue/estofolete_table4_mapping.py` (~330 LOC) mapping pipeline outputs to the four-row composite correlate of Estofolete et al. (npj Vaccines 2026 11:68) Table 4
- `scripts/dengue/retrospective_phase3_validation.py` (~200 LOC) orchestrating the retrospective scoring of the three Phase 3 vaccines and producing a comparison table against published clinical efficacy
- `docs/dengue/DENGUE_EXTENSION.md` documenting the extension methodology, attribution, and how to run

### Security
- Scrubbed Gurobi WLS credentials inherited from upstream `containers/jess_ev/src/gurobi.lic` and the `Dockerfile`'s `ENV` variables; replaced with placeholders and added runtime-injection documentation in `containers/jess_ev/README.md`. The upstream git history still contains the originals; those credentials should be considered compromised and rotated by their owner.

### Notes
- This fork preserves attribution to the upstream pipeline of Versiani et al. (Sci. Adv. 2026; 12:eaeb2066). All Nextflow workflows, container definitions, and reference data infrastructure inherit unchanged.
- Two of the four Estofolete Table 4 rows (Tier A avidity, Tier B MBC breadth) are explicitly NOT implemented in code because they require post-Phase-1 vaccinee samples that do not exist for never-administered candidates. They are documented in the manuscript Methods as the prospective MVC validation pathway.

## [2026-05-06 session-end notes]

### Verified on cluster
- Built `cd-hit-4-8-1.sif` and `cdhit-to-tsv.sif` from upstream Docker contexts
- CD-HIT process **actually ran** on the 4 dengue polyproteins via
  `-profile hyperplane_singularity` and produced the expected 4-cluster output
- Pete's licenses confirmed live on cluster:
  - Gurobi WLS credentials baked into `pmccaffrey6/jess_ev:latest` (env vars + /gurobi.lic)
  - NetMHCpan-4.1 + NetMHCIIpan-4.3 binaries committed in `containers/netmhcpan_*/src/`
- No new license registrations required for the dengue extension; the `register
  for licenses` advice in the earlier runbook is superseded.

### Known container compatibility issue (tractable)
- Pete's CDHITTOTSV container fails on `np.bool` (deprecated in numpy >=1.24,
  but pandas.util.testing still references it). Per-container fix: pin
  pandas to >=2.0 or pin numpy to <1.24 in the relevant Dockerfile, then
  rebuild that single SIF.
