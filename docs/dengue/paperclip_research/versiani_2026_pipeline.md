# Versiani, McCaffrey et al. 2026 — Immunoinformatics Platform: Pipeline Audit

## Bibliographic header

- Target citation supplied by user: Versiani AF, McCaffrey P et al. 2026. Sci. Adv. 12:eaeb2066. Title approx. "Immunoinformatics platform for epitope scoring."
- Repo-confirmed working title (GitHub release tag 1.0, published 2026-02-25): "Iterative Vaccine Design in an Era of Emerging Infectious Disease."
- Software: github.com/pmccaffrey6/immunoinformatics_platform (master branch; release 1.0; last push 2026-02-25). Nextflow DSL2; Python 96.8% / Nextflow 2.4%. Owner pmccaffrey6.
- Verification status: Paper NOT indexed in paperclip full-text or abstract corpora as of 2026-05-06 (searched: "Versiani McCaffrey immunoinformatics", "Versiani dengue", "aeb2066", "Vasilakis Versiani Nogueira", abstracts and full-text). DOI URL (sciadv.aeb2066) returned 403 to WebFetch (Cloudflare block, not 404). Repository release metadata is the strongest existence signal available; treat the Sci. Adv. citation as unverified pending direct journal-site retrieval.

## Verbatim list of pipeline processes (from `nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf`)

Container directory in `containers/` shows the canonical tool/version slugs the authors ship:

| Process (Nextflow) | Tool | Version (per container slug) |
|---|---|---|
| PROCESSINPUTFASTA | Biopython FASTA cleanup | n/a |
| FORMATALLELEFREQUENCIES | in-house allele-frequency table builder | host-data table |
| GETUNIPROTBYACCESSION | EBI Proteins API | live |
| BLASTFROMFILES | NCBI BLAST+ (`makeblastdb`, `blastp`) | n/a stated |
| FILTERPROTEINSBYBLAST | pident > 80 filter | n/a |
| CLUSTALOMEGAMSA | Clustal Omega | n/a stated |
| CONSERVEDSEQSFROMCLUSTAL / COMBINEDCLUSTALEPITOPEFASTAS | conservation parser | n/a |
| CDHIT / CDHITTOTSV | CD-HIT | container `cd-hit` (similarity 0.95 default) |
| EPIDOPE | EpiDope | container `epidope_03` (v0.3) |
| BEPIPRED / BEPIPREDTOTSV(2) | BepiPred (IEDB B-cell standalone) | container `iedb_b_cell_31` (v3.1) |
| MIGRATEAF2 | AlphaFold2 ranked_0 staging | container `alphafold2` |
| DISCOTOPE / DISCOTOPETOTSV | DiscoTope | container `discotope_11` (v1.1; invokes `/discotope-1.1/discotope`) |
| NETMHCPANI / NETMHCPANIIEDB | NetMHCpan | container `netmhcpan_i_4.1` (v4.1; `-l 9 -BA`) |
| NETMHCPANII / NETMHCPANIIIEDB | NetMHCIIpan | container `netmhcpan_ii_4.1` (v4.1; `-length 15 -BA -rankF 10`) |
| SPLITFASTAS | local FASTA splitter | n/a |
| CONSOLIDATEEPITOPES / GATHEREPITOPEFASTAS | pandas merge (epidope/bepipred/discotope/netmhcpan-I/-II) | n/a |
| TCRPMHC | TCRpMHCmodels (Galaxy PepDock template) | container `tcrpmhc-models`, `galaxy_pepdock` |
| PDBEPISA / PDBEPISATOTABLE | PISA (interface area / solvation energy) | host `run_mod` script |
| PREPAREDATAFORJESSEV / RUNJESSEV | JessEV ILP | container `jess_ev` (`pmccaffrey6/jess_ev:latest`, `design.py -e 3 -s 4`) |

Default switches (all "no" unless overridden): cdhit_input_proteins, bepipred, epidope, netmhcpani, netmhcpanii, dc_bcell, consolidate_epitopes, tcrpmhc, jessev, include_docking_in_immunogenicity. Default region: Brazil. Default MHC-I alleles: 12 HLA-A/B; MHC-II: 7 DRB1/3/4/5.

## Methods excerpt per process (1–2 sentences each, distilled from the workflow source)

- PROCESSINPUTFASTA: strips FASTA descriptors and emits a clean multi-FASTA plus a `protein_id`/`description` join table.
- BLAST + Clustal subworkflow: pulls user-supplied B-cell antigen template UniProt accessions (default `Q8QZ72,Q8JUX5`), BLASTs input proteomes (evalue 0.0005, retains pident>80), Clustal-Omega aligns hits, and exports conserved fragments ≥9 residues as a FASTA for downstream B-cell scoring.
- FORMATALLELEFREQUENCIES: filters a precomputed all-allele-frequencies TSV by region substrings, computes per-locus prevalence, and reformats keys to `HLA-A01:01` / `DRB1_0301` style.
- CDHIT: clusters input proteins or epitopes at user-set similarity (default 0.95) with `-l 5 -aS 0.9 -uS 0.1 -g 1`.
- EPIDOPE: runs EpiDope deep-learning B-cell predictor with 12 threads on either raw or BLAST-conserved FASTA.
- BEPIPRED: cleans non-canonical residues then runs IEDB `predict_antibody_epitope.py -m Bepipred`; downstream filter retains peptides 5–22 aa.
- MIGRATEAF2 + DISCOTOPE: pulls AlphaFold2 `ranked_0.pdb` predictions from `/tmp/alphafold/*` and runs DiscoTope-1.1 on chain A for conformational B-cell scoring.
- NETMHCPANI: runs netMHCpan-4.1 with `-BA -l 9` against 12 HLA-A/B alleles, exporting EL/BA scores and ranks.
- NETMHCPANII: runs netMHCIIpan-4.1 with `-BA -length 15 -rankF 10` against 7 DRB alleles.
- CONSOLIDATEEPITOPES: merges five tables (EpiDope, BepiPred, DiscoTope, NetMHCpan-I, NetMHCpan-II), filters BA rank ≤ user threshold (default 0.05), and writes per-tool feathers/TSVs/FASTAs.
- TCRPMHC + PDBEPISA: builds TCR-pMHC homology models (Galaxy PepDock / TCRpMHCmodels) for the top NetMHCpan-I peptides (capped at 100), runs PISA for interface SOLVENTAREA / SOLVATIONENERGY metrics, retains sequences <10 aa.
- PREPAREDATAFORJESSEV: weights NetMHCpan-I BA score by regional allele prevalence ("weighted_immunogen"), optionally multiplied by PISA docking output.
- RUNJESSEV: iterative ILP (`design.py -e 3 -s 4`) selects 3 epitopes with min spacer length 4 per round; previously selected epitopes are masked from the next round and outputs concatenated.

## Canonical citation form (recommended by authors)

GitHub release v1.0 names the publication as "Iterative Vaccine Design in an Era of Emerging Infectious Disease" (Sci. Adv., 2026). The repo carries no `CITATION.cff`, no `LICENSE`, and no DOI badge; the authors' preferred citation is therefore the journal article itself, with a software footnote of the form:

> McCaffrey P. et al. immunoinformatics_platform, v1.0. github.com/pmccaffrey6/immunoinformatics_platform (release 2026-02-25). Cited as the software companion to Versiani AF, McCaffrey P et al., Sci. Adv. 12:eaeb2066 (2026).

User-supplied DOI/eLocator (`sciadv.aeb2066`) could not be programmatically verified.

## Limitations the authors acknowledge

No `LIMITATIONS`, `CAVEATS`, `Discussion`, or methods-paper text is present in the public repo; limitations would live in the (unverified) Sci. Adv. manuscript. Implicit limitations visible from code:

- DiscoTope locked at v1.1 (DiscoTope-3 not used); BepiPred at IEDB v3.1.
- MHC panels hard-coded to 12 class-I + 7 class-II European-skewed alleles; allele-frequency weighting only as good as the upstream `all_allele_frequencies.tsv`.
- TCRpMHC step is capped at 100 peptides (`take(100)`) and PISA filtering retains sequences <10 aa, so longer or lower-ranked candidates are discarded.
- JessEV defaults to 3 epitopes per vaccine, spacer ≥4, MAX_NUM_PROTEIN_IDS=10 cap → vaccine breadth bounded by these constants.
- Hard-coded absolute paths under `/home/pathinformatics/...` and a Docker socket (`docker.from_env()`) inside RUNJESSEV → not portable to standard Nextflow executors without modification.

## Source files inspected

- `nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf` (1180 lines)
- `containers/` directory listing (17 tool subdirs)
- Repo metadata + release `1.0` (GitHub API)
- Paperclip full-text + abstract corpora (paper not indexed)
