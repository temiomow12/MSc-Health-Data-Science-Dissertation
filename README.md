# MSc-Health-Data-Science-Dissertation
Title: Tracing the Return of Psoriatic Lesional Skin Toward a Non‑Lesional State Through Gene Expression and Transcription-Factor Activity
This repository contains the R analysis code for the dissertation. It is a secondary analysis of the publicly available psoriasis RNA-seq dataset E-MTAB-14509 (Rider et al., 2026).

## Data availability
**Primary dataset (E-MTAB-14509)** — RNA-seq expression and sample metadata are available from ArrayExpress under accession
E-MTAB-14509:(https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14509).
**Inflammatome gene set** — the ranked inflammation signature is from Cort et al (2026)(https://github.com/isadpc/Inflammatome/tree/main/data/data/04_rank_agg_list.tsv); the file used here was stored locally and is not redistributed in this repository.

Analyses were run on the University of Birmingham BlueBEAR HPC service, in R 4.5.0. 
Key packages:
| biomaRt | 2.64.0 | Ensembl to HGNC identifier mapping |
| decoupleR | 2.14.0 (Bioconductor 3.21) | TF activity inference from CollecTRI regulons |
| glmnet | 4.1-9 | Penalised logistic regression (LASSO, elastic net) |
| pROC | 1.18.5 | ROC curves and AUC |
| limma | 3.61.4 | Expression data handling |
| ggplot2 | 3.5.2 | Visualisation |
| ggrepel | 0.9.8 | Non-overlapping plot labels |
| dplyr / tidyverse | 1.1.4 / 2.0.0 | Data manipulation |
| tibble | 3.2.1 | Data frames |
| readr | 2.1.5 | File reading |
| formatR | 1.14 | Code formatting |

## Scripts
Run in this order. Each script is self-contained (loads its own data); scripts 3 and 4 also re-derive the week-0/week-12 objects they need, so they do not depend on 1–2 having been run first in the same session.

- **`01_EDA.R`** — 	Exploratory data analysis, cohort summary, PASI trajectories, and unsupervised PCA (including PC1 loadings and inflammatome highlight)
  
- **`02_inflammation_score.R`** — Per-sample inflammatome score (mean log-normalised expression of the top-100 inflammatome genes) and its correlation with PASI

- **`03_LASSO_ElasticNet.R`** — 	LASSO + Elastic Net classifiers on gene expression; week-12 application; gene-signature to inflammatome overlap.

- **`04_CollecTRI_TF_activity.R`** — Transcription-factor activity inference (CollecTRI + decoupleR run_ulm); LASSO + EN on TF activity; week-12 application; TF ↔ gene-signature intersection

2. Figure/table to code map
Every numbered figure and table in the dissertation, and where it is produced.

Thesis item	Produced by:

Table 3.1 — cohort demographics	01_EDA.R (table() / summary() blocks)

Figure 3.1A — sample flow diagram	Drawn manually (not code-generated); counts come from 01_EDA.R

Figure 3.1B — PASI Week 0 vs 12	01_EDA.R (pasi_plot)

Figure 3.1C — inflammatome score vs PASI	02_inflammation_score.R

Figure 3.2 — Pearson r = 0.37 / Spearman 0.40	02_inflammatome_score.R (cor.test)

Figure 3.2A — scree plot	01_EDA.R

Figure 3.2B — PCA, all samples	01_EDA.R

Figure 3.2C/D — PCA, Week 0 / Week 12	01_EDA.R(pca_df0, pca_df12)

Figure 3.2E/F — PC1 loadings + inflammatome	01_EDA.R (pc1_genes, loadings_pca.tsv)

3.3.1 — PC1 ∩ inflammatome overlaps	01_EDA.R

Table 3.2 — EN gene signature	03_LASSO_ElasticNet.R (en_df)

Figure 3.3A/B — gene ROC (LASSO / EN)	03_LASSO_ElasticNet.R

Figure 3.3C–E — GJB2 / KYNU / KLK13	03_LASSO_ElasticNet.R

Figure 3.4.3 — signature ∩ inflammatome	03_LASSO_ElasticNet.R

Tables 3.3 / 3.4 — TF signatures (LASSO / EN)	04_CollecTRI_TF_activity.R (lasso_df, en_df)

Figure 3.4A/B — TF ROC (LASSO / EN)	04_CollecTRI_TF_activity.R

Figure 3.4C–F — STAT4 / BACH1 / CEBPE / MSX1	04_CollecTRI_TF_activity.R

Figure 3.5.2 — TF ∩ gene-signature (IL17A, PI3)	04_CollecTRI_TF_activity.R (links)

Figure 3.5A/B — gene ROC at Week 12	03_LASSO_ElasticNet.R

Figure 3.5C/D — TF ROC at Week 12	04_CollecTRI_TF_activity.R

Figure 3.5E — paired probabilities (TF)	04_CollecTRI_TF_activity.R

Figure 3.5F — paired probabilities (genes)	03_LASSO_ElasticNet.R

3.6.2 — Wilcoxon p (genes / TF)	03_LASSO_ElasticNet.R / 04_CollecTRI_TF_activity.R


**3. Input data**
Raw and patient-level data are not redistributed. The scripts expect the following files.

**File	Description	Source**
Skin_norm_counts_d.txt	Normalised, log-CPM expression matrix (genes × samples; ~15,735 genes)E-MTAB-14509, ArrayExpress
E-MTAB-14509.sdrf.txt	Sample and clinical metadata	E-MTAB-14509, ArrayExpress
Final_Annotation_List_BioMart.tsv	Ensembl ID - HGNC symbol mapping	Generated via biomaRt
04_rank_agg_list.tsv	Inflammatome rank-aggregated gene list	Cort et al. (2026) supplementary / repository
collectri_regulons.csv	CollecTRI signed regulon network	decoupleR / CollecTRI
The primary dataset is available from ArrayExpress under accession E-MTAB-14509.
I downloaded and added the E-MATB-14509 dataset and collecTRI to github.

**4. Reproducing the analysis**
Environment. Package versions are locked with renv. 
Analyses were run in R 4.5.0 on the University of Birmingham BlueBEAR HPC service.
Paths. File paths in the scripts are currently absolute to the BlueBEAR project directory (/rds/projects/b/bravol-inf-multi/Esther/...). To run elsewhere, set the base data directory once at the top of each script
