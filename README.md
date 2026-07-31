# MSc-Health-Data-Science-Dissertation
Title: Tracing the Return of Psoriatic Lesional Skin Toward a Non‑Lesional State Through Gene Expression and Transcription-Factor Activity
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

All scripts read input files from absolute paths (e.g. the BlueBEAR `/rds/projects/...` locations used here). **These paths must be edited at the top of each script to point to wherever the corresponding data sit on your own system** before the scripts will run. Run the scripts in numerical order.

## Scripts

Run in order; each builds on outputs from the previous.

- **`01_EDA.R`** — Loads and cleans the E-MTAB-14509 cohort. Removes whole-blood and week 1 samples, filters missing covariates, and characterises the dataset at patient and sample level (cohort tables, PASI distributions, PCA).

- **`02_inflammation_score.R`** — Computes a per-sample inflammation score from the inflammatome gene set and correlates it with clinical PASI scores (Pearson and Spearman).

- **`03_LASSO_ElasticNet.R`** — Trains penalised logistic-regression classifiers(LASSO and elastic net) to distinguish lesional from non-lesional skin at week 0, using gene expression as features. Includes cross-validated λ tuning, a 70/30 hold-out, ROC/AUC evaluation, and application of the frozen models to week 12 to assess transcriptomic convergence.

- **`04_CollecTRI_TF_activity.R`** — Infers transcription-factor activity from the expression data using CollecTRI regulons (via decoupleR), then repeats the classification and convergence analysis using TF activities as features. Relates the TF and gene signatures through shared regulator–target links.
