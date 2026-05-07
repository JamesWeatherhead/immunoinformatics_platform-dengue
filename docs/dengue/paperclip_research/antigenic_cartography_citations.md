# Antigenic Cartography Reference Papers

Verified via paperclip CLI on 2026-05-06. Sources: PMC full-text and OpenAlex abstracts.

## 1. Smith et al. 2004 — Original antigenic cartography (influenza H3N2)

- **Citation:** Smith DJ, Lapedes AS, de Jong JC, Bestebroer TM, Rimmelzwaan GF, Osterhaus ADME, Fouchier RAM. Mapping the antigenic and genetic evolution of influenza virus. *Science* 305:371-376 (2004).
- **DOI:** 10.1126/science.1097211
- **PMID:** 15218094
- **Paperclip ID:** abstracts/1979158913
- **Methodology:** Modified multidimensional scaling (MDS) applied to a matrix of pairwise hemagglutination inhibition (HI) titers between H3N2 viruses (1968-2003) and post-infection ferret antisera, embedding strains in a low-dimensional Euclidean "antigenic map" where distance approximates serological dissimilarity.
- **Antigenic distance metric:** **Log2 fold-reduction in HI titer.** Specifically, distance d(v,s) between virus v and serum s = log2(max titer in row of s) - log2(measured titer for v against s); each unit on the map = one two-fold dilution of HI titer.

## 2. Katzelnick et al. 2015 — Dengue antigenic cartography (canonical)

- **Citation:** Katzelnick LC, Fonville JM, Gromowski GD, Bustos Arriaga J, Green A, James SL, Lau L, Montoya M, Wang C, VanBlargan LA, Russell CA, Thu HM, Pierson TC, Buchy P, Aaskov JG, Munoz-Jordan JL, Vasilakis N, Gibbons RV, Tesh RB, Osterhaus ADME, Fouchier RAM, Durbin A, Simmons CP, Holmes EC, Harris E, Whitehead SS, Smith DJ. Dengue viruses cluster antigenically but not as discrete serotypes. *Science* 349:1338-1343 (2015).
- **DOI:** 10.1126/science.aac5017
- **PMID:** 26383952
- **Paperclip ID:** abstracts/2205335187
- **Methodology:** Direct port of Smith 2004 MDS antigenic cartography to dengue, using a large matrix of plaque reduction neutralization test (PRNT) titers between 47 DENV1-4 strains and post-primary-infection nonhuman-primate (and human) sera; resulting map shows DENV types overlap rather than form four discrete serotype clusters.
- **Antigenic distance metric:** **Log2 fold-reduction in PRNT50 neutralization titer** (i.e., one map unit = one two-fold dilution drop in 50% plaque reduction neutralization titer relative to the homologous-serum maximum), analogous to Smith 2004 but on PRNT50 rather than HI.

## 3. Bedford et al. 2014 — Antigenic dynamics with WHO HI surveillance data

- **Citation:** Bedford T, Suchard MA, Lemey P, Dudas G, Gregory V, Hay AJ, McCauley JW, Russell CA, Smith DJ, Rambaut A. Integrating influenza antigenic dynamics with molecular evolution. *eLife* 3:e01914 (2014).
- **DOI:** 10.7554/eLife.01914
- **PMID:** 24497547
- **Paperclip IDs:** PMC3909918 (full text), abstracts/2158189516
- **Methodology:** Bayesian phylogenetic model that jointly infers genealogy and a continuous antigenic phenotype on the tree from HI-titer tables produced by WHO Collaborating Centres (Crick / Francis Crick Institute) for H3N2, H1N1, Vic and Yam, mapping antigenic effects to specific HA branches.
- **Antigenic distance metric:** **Log2 HI-titer drop**, decomposed as d_ij = sum of branch effects + virus-avidity term + serum-potency term; each unit corresponds to one two-fold HI dilution, matching the Smith 2004 convention.

## Cross-paper note on the metric

All three papers use the same underlying scale: **one antigenic unit = one two-fold dilution drop in serological titer** (HI for influenza, PRNT50 for dengue). This is why Smith 2004 is the methodological parent of both Katzelnick 2015 (dengue) and Bedford 2014 (Bayesian extension).
