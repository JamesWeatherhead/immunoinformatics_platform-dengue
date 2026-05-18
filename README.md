# immunoinformatics_platform-dengue

This is a fork of [`pmccaffrey6/immunoinformatics_platform`](https://github.com/pmccaffrey6/immunoinformatics_platform) (the published pan-alphavirus epitope-selection pipeline of Versiani et al. *Sci. Adv.* 2026; 12:eaeb2066) extended to dengue and adapted to compute pre-trial estimators of the four-row composite correlate of dengue vaccine efficacy proposed by Estofolete et al. *npj Vaccines* 2026; 11:68 Table 4.

<img width="2354" height="498" alt="CleanShot 2026-05-18 at 09 11 01@2x" src="https://github.com/user-attachments/assets/13fecbfc-075d-4440-93ba-1fb7abd26c2d" />

## What this fork adds

- **DENV 1-4 reference polyproteins** plus **12 Phase 3 vaccine parent strain sequences** (CYD-TDV, TAK-003, Butantan-DV) for retrospective validation.
- **Estofolete Table 4 mapping layer** that joins the upstream pipeline outputs (DiscoTope, ESMFold, NetMHCpan I + II, AFND-weighted coverage) into per-construct Tier A neutralization breadth and Tier B CD8 polyfunctionality scores.
- **Retrospective Phase 3 validation** that runs the pipeline on the three licensed dengue vaccines and compares pipeline-derived rankings to published clinical efficacy (Sridhar 2018 NEJM, Tricou 2024 Lancet Glob Health, Kallas 2024 NEJM).
- **Honest scope**: two of four Tier 4 rows are computed pre-trial (Geometry, Equity); two require post-Phase-1 vaccinee samples (Tier A avidity, Tier B MBC breadth) and are documented as the prospective minimum-viable-cohort (MVC) validation pathway in the companion manuscript Methods.
- **Security fix**: scrubbed Gurobi credentials that the upstream had committed to `containers/jess_ev/src/gurobi.lic` and the `Dockerfile`. See `containers/jess_ev/README.md` for runtime injection.

See [`docs/dengue/DENGUE_EXTENSION.md`](docs/dengue/DENGUE_EXTENSION.md) for full methodology and how-to.

## Quick start

```bash
# 1. Clone
git clone https://github.com/JamesWeatherhead/immunoinformatics_platform-dengue.git
cd immunoinformatics_platform-dengue

# 2. Pull dengue inputs (DENV references + Phase 3 vaccine parent strains)
bash scripts/dengue/fetch_dengue_data.sh

# 3. Register for NetMHCpan-4.1 + NetMHCIIpan-4.3 academic licenses
#    https://services.healthtech.dtu.dk/

# 4. Get a Gurobi academic license for the JessEV step
#    https://www.gurobi.com/academia/academic-program-and-licenses/

# 5. Run the upstream pipeline on dengue inputs (cluster recommended; ~24h)
nextflow run nextflow-workflows/full_epitope_scoring/full_epitope_scoring.nf \
    --protein_file ref_fastas/dengue/dengue_polyproteins_multiseq.fasta \
    --epitope_output_folder /path/to/epitope_outputs/dengue_main \
    --allele_target_region "Brazil,Thailand,Vietnam,Philippines"

# 6. Map outputs to Estofolete Table 4
python scripts/dengue/estofolete_table4_mapping.py \
    --epitope-output-folder /path/to/epitope_outputs/dengue_main \
    --outdir outputs/dengue_table4/
```

For the retrospective Phase 3 validation, see `scripts/dengue/retrospective_phase3_validation.py` and `docs/dengue/DENGUE_EXTENSION.md`.

## Citing

If you use this fork, please cite **both** the upstream pipeline and the dengue companion paper:

- Versiani AF, McCaffrey P, Ribeiro-Filho HV, Silva NIO, Lopes-de-Oliveira PS, Carrera J-P, Nogueira ML, Marques RE, Rossi SL, Vasilakis N. **Integrated reiterative pipeline for rapid epitope-based pan-alphavirus vaccines.** *Sci. Adv.* 2026; 12(11):eaeb2066. [doi:10.1126/sciadv.aeb2066](https://doi.org/10.1126/sciadv.aeb2066)
- Estofolete CF, Saivish MV, Nogueira ML, Vasilakis N. **From promise to pitfalls: immunological lessons from dengue vaccines and their implications.** *npj Vaccines* 2026; 11:68. [doi:10.1038/s41541-026-01400-4](https://doi.org/10.1038/s41541-026-01400-4)
- Saivish MV et al. (manuscript in revision; this fork is the implementation system the editor's note required)

## License

Inherits the upstream license. The dengue extension files added in this fork are released under the same terms.

---

The original upstream README is preserved as `README_upstream_pmccaffrey6.md`.
