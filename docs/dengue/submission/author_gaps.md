# Author and funding gaps checklist (npj Vaccines pre-submission)

This checklist enumerates remaining metadata required for the npj Vaccines submission portal. Cross-reference with `credit_author_contributions.md` and `competing_interests.md` (already finalised in this directory).

## 1. Author ORCIDs

ORCIDs were resolved against the ORCID public registry (pub.orcid.org v3.0). Disambiguation was performed against institutional affiliation, given/family-name match, and known publication history. Co-authors must each confirm their own ID at submission; the values below are the verified-most-likely matches and should be checked by each author before the corresponding-author cover sheet is finalised.

| # | Author | Affiliation | ORCID (verified candidate) |
|---|--------|-------------|----------------------------|
| 1 | Nikos Vasilakis (corresponding) | UTMB Galveston, Pathology | 0000-0002-0708-3289 |
| 2 | James Weatherhead (corresponding) | UTMB Galveston (MD/PhD/MPH candidate) | 0009-0003-5436-7641 |
| 3 | Peter McCaffrey | UTMB Galveston, Pathology Informatics | 0000-0002-4450-8143 |
| 4 | Daniela Weiskopf | La Jolla Institute for Immunology | 0000-0003-2968-7371 |
| 5 | Cassia F. Estofolete | FAMERP, Sao Jose do Rio Preto | 0000-0001-7324-1591 |
| 6 | Mauricio L. Nogueira | FAMERP, Sao Jose do Rio Preto | 0000-0003-1102-2419 |

Status: 6/6 ORCIDs identified. Each author to confirm at sign-off; if any candidate is incorrect, mark "to be supplied by author" on the portal.

## 2. Funding sources / NIH grant numbers

Active NIH grants on which Vasilakis (UTMB) is named PI/contact PI (queried via NIH RePORTER API, FY 2024-2025 active awards):

- **U01 AI151807** - Coordinating Research on Emerging Arboviral Threats Encompassing the Neotropics (CREATE-NEO), NIAID. UTMB.
- **D43 TW012246** - Training and Research on Arboviruses and Zoonoses In Nigeria and Sierra Leone (TRAIN), Fogarty International Center. UTMB.
- **R01 AI145918** - Trade-offs between Arbovirus Transmission and Clearance in Native and Novel Hosts, NIAID. (Subaward; New Mexico State University prime.)
- **U01 AI115577** - Mechanisms and Public Health Impact of Sylvatic Dengue Virus Emergence in Borneo (recently completed; reference for sylvatic-dengue methodology). UTMB.

McCaffrey (UTMB Pathology Informatics): grants supporting `immunoinformatics_platform` development - to be supplied by author.

Weiskopf (LJI): NIAID T-cell-epitope contracts (HHSN272201400045C, 75N93019C00065 historical; current funding) - to be supplied by author.

Estofolete and Nogueira (FAMERP): FAPESP and CNPq grant numbers underwriting the Estofolete et al. 2026 *npj Vaccines* composite-correlate work - to be supplied by authors.

## 3. AI-assisted-writing declaration (npj Vaccines policy, 2024-)

> The authors used Anthropic Claude Code (Opus 4.7) to assist with code development, citation verification via the paperclip CLI, and manuscript formatting. All scientific claims, data, and analyses were reviewed and verified by the authors. No generative-AI tool is listed as an author, and no AI-generated text was included without author review and revision.

## 4. Ethics statement

This work is entirely *in silico*. It uses publicly released amino-acid sequences (DENV reference polyproteins from GenBank/UniProt and the three Phase 3 parent strains from prior publications), publicly released HLA allele-frequency tables (AFND release 2024-12), and publicly released T-cell and B-cell epitope curations (IEDB). No human subjects were enrolled, no human biospecimens were analysed, no animal experiments were performed, and no individual-level patient data were used. UTMB IRB exemption is therefore not applicable; no ACUC protocol applies.

## 5. Suggested reviewers (3, no UTMB / LJI / FAMERP coauthorship)

1. **Aravinda M. de Silva, PhD.** Department of Microbiology and Immunology, University of North Carolina at Chapel Hill. Dengue B-cell immunology and EDE / quaternary epitope expert. No UTMB, LJI, or FAMERP coauthorship.
2. **Eva Harris, PhD.** Division of Infectious Diseases and Vaccinology, School of Public Health, University of California, Berkeley. Dengue clinical immunology and pediatric cohort lead (Nicaragua). No UTMB, LJI, or FAMERP coauthorship.
3. **Bjoern Peters, PhD.** La Jolla Institute for Immunology. IEDB lead and immunoinformatics. **Note: LJI affiliation overlap with Weiskopf - flag this conflict explicitly to handling editor; if the editor prefers no LJI overlap, alternates: Lauren Yauch (Stanford / IDEAYA, immunoinformatics) or Felix Rey (Institut Pasteur, EDE structural biology).**

## 6. Pre-print declaration

The manuscript has not been deposited as a preprint at the time of submission. Upon acceptance at npj Vaccines (or at the corresponding author's discretion at first revision), a version will be deposited at **bioRxiv** under DOI to be assigned, with the npj Vaccines submission ID cross-referenced. No earlier version of this manuscript has been posted to bioRxiv, medRxiv, arXiv, SSRN, or Research Square.
