library(dplyr)
library(readr)
library(decoupleR)
library(tibble)
library(tidyr)
library(glmnet)
library(ggplot2)
library(pROC)
options(device = 'RStudioGD')
options(device = 'png')
install.packages("ragg")
options(device = ragg::agg_png)

#EDA
collectri<-read_csv('/rds/projects/b/bravol-inf-multi/Esther/Data/prior_knowledge/collectri_regulons.csv')
length(unique(collectri$source))
colnames(collectri) #1176
head(collectri)
length(collectri$target) #47643

# EDA: regulon size per TF 
collectri%>%count(source,name = 'regulon_size')%>%
  arrange(desc(regulon_size))%>%head(20)

#match gene names
expr<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/E-MTAB-14509/Skin_norm_counts_d.txt',header = TRUE,
                 check.names = FALSE)
expr<-expr[rowSums(is.na(expr))==0,]
expr$ensembl_id<-rownames(expr)

expr$ensembl_id <- sub(
  "\\..*",
  "",
  expr$ensembl_id
)
biomart<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/annotations/Final_Annotation_List_BioMart.tsv',header = TRUE,
                    check.names = FALSE)
expr<-merge(
  expr,
  biomart[,c('ENSG.ID','Gene.name')],
  by.x = 'ensembl_id',
  by.y = 'ENSG.ID',
  all.x = TRUE
)
expr <- expr %>%
  filter(!is.na(Gene.name), Gene.name != "")
expr <- expr %>%
  select(-ensembl_id) %>%
  group_by(Gene.name) %>%
  summarise(
    across(where(is.numeric), mean),
    .groups = "drop"
  )
expr<-as.data.frame(expr)
rownames(expr)<-expr$Gene.name
expr$Gene.name<-NULL

expr_genes<-rownames(expr)
cov<- collectri%>%
  mutate(in_data=target %in% expr_genes) %>%
  group_by(source)%>%
  summarise(n=n(), n_in=sum(in_data), frac=mean(in_data)) %>%
  arrange(frac)
cov %>%
  arrange(desc(frac)) %>%
  head(20)
#n = Total number of target genes for that source (e.g. transcription factor)
#n_in = Number of those target genes that are present in your expression dataset
#frac = Fraction (or proportion) of target genes present in your dataset, calculated as n_in / n

#filtering to usable network
ct_use<-collectri %>%
  filter(target %in% expr_genes) %>%  #keep measurable targets only
  add_count(source,name = 'regulon_size') %>%
  filter(regulon_size>=5) %>%
  select(source,target,mor)
#minsize >=5 is the standard decoupleR cutoff; TFs with fewer than five measured targets give unstable activity estimates and are excluded

#per-sample TF activity
acts<-run_ulm(mat=expr, net=ct_use,
              .source='source',.target = 'target',
              .mor = 'mor',minsize = 5)
acts$statistic<-NULL
#this reduces the feature space (the Univariate Linear model; for one sample and one TF, decoupleR fits a simple linear regression model)

gene_to_tf<-ct_use%>%
  group_by(target)%>%
  summarise(TFs=paste(sort(unique(source)),collapse = ';'),
            n_TF=n_distinct(source))
gene_to_tf%>%filter(target=='IL36G')

#EDA of new dataset
summary(acts$score)
length(unique(acts$source)) #694
unique(acts$source)

#prep table for modelling
acts<-acts%>%
  select(source,condition,score)%>%
  pivot_wider(
    names_from = source,
    values_from = score)
dataset <- read.delim("/rds/projects/b/bravol-inf-multi/Esther/Data/E-MTAB-14509/E-MTAB-14509.sdrf.txt",
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE,
                      check.names = FALSE)
psoriasis<-dataset[dataset$`Characteristics[organism part]`!='Whole Blood',]
psoriasis<-psoriasis[psoriasis$`Factor Value[time]`!='1',]
psoriasis0<- psoriasis[
  !(psoriasis$`Factor Value[time]` %in% c(1, 12)),]

expr_ensg<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/E-MTAB-14509/Skin_norm_counts_d.txt',header = TRUE,
                 check.names = FALSE)
expr_ensg<-expr_ensg[rowSums(is.na(expr_ensg))==0,]
expr_ensg$ensembl_id<-rownames(expr_ensg)

expr_ensg$ensembl_id <- sub(
  "\\..*",
  "",
  expr_ensg$ensembl_id
)

week0.samples<-psoriasis0$`Source Name`
week0.samples<-intersect(week0.samples,colnames(expr_ensg))
expr_week0<-expr_ensg[,week0.samples]
#transposed<-t(expr_week0)
acts_week0 <- acts %>%
  filter(condition %in% colnames(expr_week0))

acts_week0<-merge(
  acts_week0,
  psoriasis[,c("Source Name","Characteristics[organism part]")],
  by.x='condition',
  by.y='Source Name',
  all.x = TRUE)
acts$`Characteristics[organism part]`

is.na(acts_week0$`Characteristics[organism part]`)
acts_week0<-acts_week0[rowSums(is.na(acts_week0))==0,]
is.na(acts_week0$`Characteristics[organism part]`)
#formatting
y<-factor(acts_week0$`Characteristics[organism part]`) #converts labels (i.e.lesional/non-lesional) into factors
x <- acts_week0[, !(names(acts_week0) %in%
                     c("condition", 'Characteristics[organism part]'))]
x <- as.matrix(x)

#LASSO model
# formatting (unchanged up to here)
y <- factor(acts_week0$`Characteristics[organism part]`,
            levels = c("Lesional Skin", "Nonlesional Skin"))  # set reference explicitly
x <- acts_week0[, !(names(acts_week0) %in%
                      c("condition", "Characteristics[organism part]"))]
x <- as.matrix(x)

stopifnot(!anyNA(y))

# 70/30 random split
set.seed(42)
train_index_vec <- sample(seq_len(nrow(x)), size = floor(0.70 * nrow(x)))
train_index <- seq_len(nrow(x)) %in% train_index_vec
test_index  <- !train_index

x_train <- x[train_index, , drop = FALSE]
x_test  <- x[test_index,  , drop = FALSE]
y_train <- y[train_index]
y_test  <- y[test_index]

# --- LASSO on the 70% training set ---
set.seed(42)
cv_lasso <- cv.glmnet(x_train, y_train, family = "binomial",alpha = 1, type.measure = "auc",nfolds = 10, keep = TRUE)
max(cv_lasso$cvm)

# --- Apply to untouched 30% ---
test_prob_lasso <- as.numeric(predict(cv_lasso, newx = x_test,
                                      s = "lambda.1se", type = "response"))
test_pred_lasso <- factor(
  ifelse(test_prob_lasso >= 0.5, "Lesional Skin", "Nonlesional Skin"),
  levels = levels(y))

confusion_matrix_lasso <- table(Predicted = test_pred_lasso, Actual = y_test)
confusion_matrix_lasso
test_accuracy_lasso <- mean(test_pred_lasso == y_test)
test_accuracy_lasso

test_roc_lasso <- roc(y_test, test_prob_lasso,
                      levels = levels(y), direction = "<", quiet = TRUE)
auc(test_roc_lasso)
ci.auc(test_roc_lasso)
plot(test_roc_lasso,main=paste0('LASSO test-set ROC:AUC =',
                             round(as.numeric(auc(test_roc_lasso)),3)))

# --- Signature from the training model ---
lasso_coef  <- coef(cv_lasso, s = "lambda.1se")
lasso_TF <- setdiff(rownames(lasso_coef)[which(as.numeric(lasso_coef) != 0)],
                       "(Intercept)")
coef_values_lasso   <- as.numeric(lasso_coef)
feature_names_lasso <- rownames(lasso_coef)
keep_lasso <- coef_values_lasso != 0 & feature_names_lasso != "(Intercept)"
lasso_df <- data.frame(TF_name = feature_names_lasso[keep_lasso],
                       coefficient = coef_values_lasso[keep_lasso],
                       stringsAsFactors = FALSE)
lasso_df <- lasso_df[order(-abs(lasso_df$coefficient)), ]
write.table(lasso_df, '/rds/projects/b/bravol-inf-multi/Esther/Data/lasso_coef_TF',
            sep = '\t', row.names = TRUE, quote = FALSE)
write.table(lasso_df,'/rds/projects/b/bravol-inf-multi/Esther/Data/lasso_coef_TF',sep = '\t',row.names = T,quote = F)

#EN Model
set.seed(42)
train_index_vec_en <- sample(seq_len(nrow(x)), size = floor(0.70 * nrow(x)))
train_index_en <- seq_len(nrow(x)) %in% train_index_vec_en
test_index_en  <- !train_index_en

x_train <- x[train_index, , drop = FALSE]
x_test  <- x[test_index,  , drop = FALSE]
y_train <- y[train_index]
y_test  <- y[test_index]

#EN on the 80% training set
set.seed(42)
cv.en <- cv.glmnet(x_train, y_train, family = "binomial",alpha = 0.5, type.measure = "auc",nfolds = 10, keep = TRUE)
max(cv.en$cvm)

# --- Apply to untouched 20% ---
test_prob_en <- as.numeric(predict(cv.en, newx = x_test,
                                      s = "lambda.1se", type = "response"))
test_pred_en <- factor(
  ifelse(test_prob_en >= 0.5, "Lesional Skin", "Nonlesional Skin"),
  levels = levels(y))

confusion_matrix_en <- table(Predicted = test_pred_en, Actual = y_test)
confusion_matrix_en
test_accuracy_en <- mean(test_pred_en == y_test)
test_accuracy_en

test_roc_en <- roc(y_test, test_prob_en,
                      levels = levels(y), direction = "<", quiet = TRUE)
auc(test_roc_en)
ci.auc(test_roc_en)

plot(test_roc_en,main=paste0('Elastic Net test-set ROC:AUC =',
                          round(as.numeric(auc(test_roc_en)),3)))
# --- Signature from the training model ---
en_coef  <- coef(cv.en, s = "lambda.1se")
en_TF <- setdiff(rownames(en_coef)[which(as.numeric(en_coef) != 0)],
                       "(Intercept)")
coef_values_en   <- as.numeric(en_coef)
feature_names_en<- rownames(en_coef)
keep_en <- coef_values_en != 0 & feature_names_en != "(Intercept)"
en_df <- data.frame(TF_name = feature_names_en[keep_en],
                       coefficient = coef_values_en[keep_en],
                       stringsAsFactors = FALSE)
en_df <- en_df[order(-abs(en_df$coefficient)), ]
write.table(lasso_df, '/rds/projects/b/bravol-inf-multi/Esther/Data/lasso_coef_TF',sep = '\t', row.names = TRUE, quote = FALSE)
write.table(en_df,'/rds/projects/b/bravol-inf-multi/Esther/Data/en_coef_TF',sep = '\t',row.names = T,quote = F)

#how many TFs did each model keep?
length(lasso_TF);length(en_TF)
overlap<-intersect(en_TF,lasso_TF)
print(overlap)

##psoriasis dataset week 12
psoriasis12 <- psoriasis[
  !(psoriasis$`Factor Value[time]` %in% c(0,1)),
]
names(psoriasis12) <- make.unique(names(psoriasis12))
week12.samples<-psoriasis12$`Source Name`
week12.samples<-intersect(week12.samples,colnames(expr_ensg))
expr_week12<-expr_ensg[,week12.samples]
acts_week12 <- acts %>%
  filter(condition %in% colnames(expr_week12))

acts_week12<-merge(
  acts_week12,
  psoriasis[,c("Source Name","Characteristics[organism part]")],
  by.x='condition',
  by.y='Source Name',
  all.x = TRUE)
#acts$`Characteristics[organism part]`

is.na(acts_week12$`Characteristics[organism part]`)
acts_week12<-acts_week12[rowSums(is.na(acts_week12))==0,]
is.na(acts_week12$`Characteristics[organism part]`)
#formatting
y_12<-factor(acts_week12$`Characteristics[organism part]`) #converts labels (i.e.lesional/non-lesional) into factors
x_12 <- acts_week12[, !(names(acts_week12) %in% 
                          c("condition", "Characteristics[organism part]"))]

x_12 <- as.matrix(x_12)
x_12 <- x_12[, colnames(x)]          # same TFs, same order as the trained model
stopifnot(all(colnames(x_12) == colnames(x)))

#pred week 12 lasso
pred12_lasso<- as.numeric(predict(cv_lasso, newx = x_12,
                                   s = "lambda.1se", type = "response"))
#pred week 12 en
pred12_en<- as.numeric(predict(cv.en, newx = x_12,
                                   s = "lambda.1se", type = "response"))
#roc
roc12_lasso<- roc(y_12, pred12_lasso, quiet = TRUE)
roc12_en<- roc(y_12, pred12_en,    quiet = TRUE)

ci.auc(roc12_lasso)     # AUC + 95% CI, LASSO on week 12
ci.auc(roc12_en) #AUC +95%, EN week 12

pred0_lasso <- as.numeric(predict(cv_lasso, newx = x,
                                  s = "lambda.1se", type = "response"))
w0 <- data.frame(
  sample  = colnames(expr_week0),
  prob    = pred0_lasso,
  status  = y,
  week    = "Week 0"
)
w12 <- data.frame(
  sample  = colnames(expr_week12),
  prob    = pred12_lasso,
  status  = y_12,
  week    = "Week 12"
)
#aligns each sample to patient (for pairing)
w0$patient <- sub("-.*", "", w0$sample)   # keep everything before the first hyphen
w12$patient<-sub('-.*','',w12$sample)
combined<- rbind(w0,w12) #combines w0 and w12

#isolating lesional samples
les <- combined[combined$status == "Lesional Skin", ] 
paired <- merge(les[les$week == "Week 0",  c("patient","prob")], #matches each patient's week 0 lesional probability to their own week 12 probability, producing one row per patient with two columns : prob wk0, prob wk12
                les[les$week == "Week 12", c("patient","prob")],
                by = "patient", suffixes = c("_wk0","_wk12"))
paired$prob_wk0  <- 1 - paired$prob_wk0
paired$prob_wk12 <- 1 - paired$prob_wk12
#Wilcoxon test
convergence<-wilcox.test(paired$prob_wk0, paired$prob_wk12, paired=TRUE)
convergence


#ROC Figure discrimination @week 12 LASSO
roc_df <- data.frame(fpr = 1 - roc12_lasso$specificities,
                     tpr = roc12_lasso$sensitivities)
ggplot(roc_df, aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = "dashed") +
  geom_path(colour = "#534AB7", linewidth = 0.9) +
  coord_equal() +
  labs(title = sprintf("Week 12 ROC (frozen model): AUC %.3f",
                       as.numeric(auc(roc12_lasso))),
       x = "False positive rate", y = "True positive rate") +
  theme_minimal()

#ROC Figure discrimination @week 12 EN
roc_df_en <- data.frame(fpr_en = 1 - roc12_en$specificities,
                     tpr_en = roc12_en$sensitivities)
ggplot(roc_df_en, aes(fpr_en, tpr_en)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = "dashed") +
  geom_path(colour = "#534AB7", linewidth = 0.9) +
  coord_equal() +
  labs(title = sprintf("Week 12 ROC (frozen model) Elastic Net: AUC %.3f",
                       as.numeric(auc(roc12_en))),
       x = "False positive rate", y = "True positive rate") +
  theme_minimal()

# (b) the convergence figure: paired probabilities, wk0 vs wk12
ggplot(paired) +
  geom_segment(aes(x = "Week 0", xend = "Week 12",
                   y = prob_wk0, yend = prob_wk12, group = patient),
               colour = "grey80") +
  geom_point(aes(x = "Week 0",  y = prob_wk0),  colour = "#993C1D") +
  geom_point(aes(x = "Week 12", y = prob_wk12), colour = "#185FA5") +
  labs(title = "Predicted lesional probability, lesional samples",
       subtitle = "Frozen week-0 model applied to week 12; paired by patient",
       x = NULL, y = "Predicted probability of lesional") +
  theme_minimal()

#TF coefficients plotted 
# MSX1
msx1_df <- data.frame(expression = x[, "MSX1"], status = y)
ggplot(msx1_df, aes(x = status, y = expression)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  ggtitle("MSX1") +
  labs(x = "Lesion status", y = "Expression (log-normalised)")

#CEBPE
cebpe_df <- data.frame(expression = x[, "CEBPE"], status = y)
ggplot(cebpe_df, aes(x = status, y = expression)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  ggtitle("CEBPE") +
  labs(x = "Lesion status", y = "Expression (log-normalised)")
#BACH1
bach1_df <- data.frame(expression = x[, "BACH1"], status = y)
ggplot(bach1_df, aes(x = status, y = expression)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  ggtitle("BACH1") +
  labs(x = "Lesion status", y = "Expression (log-normalised)")

#STAT4
stat4_df <- data.frame(expression = x[, "STAT4"], status = y)
ggplot(stat4_df, aes(x = status, y = expression)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  ggtitle("STAT4") +
  labs(x = "Lesion status", y = "Expression (log-normalised)")


#SEPARATE
lasso_genes<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/Dataframes created/lasso_coef',header = TRUE,
                        check.names = FALSE)
en_genes<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/Dataframes created/en_coef',header = TRUE,
                     check.names = FALSE)
en_genes<-en_genes$HGNC
lasso_genes<-lasso_genes$HGNC
tf_set_lasso <- as.character(lasso_TF)# <- swap for en_tf or union(lasso_tf, en_tf)
tf_set_lasso <- unique(tf_set_lasso[!is.na(tf_set_lasso)])

tf_targets_lasso <- collectri %>%
  filter(source %in% tf_set_lasso) %>%
  distinct(source, target)

target_pool_lasso <- as.character(unique(tf_targets_lasso$target))
cat("TFs matched:", length(intersect(tf_set_lasso, collectri$source)),
    "| target genes:", length(target_pool_lasso), "\n")

# ---- 2. Overlap ----
lasso_overlap <- base::intersect(as.character(en_genes), as.character(target_pool_lasso))
cat("\nLASSO genes that are targets of selected TFs:\n"); print(lasso_overlap)

tf_set_en <- as.character(en_TF)# <- swap for en_tf or union(lasso_tf, en_tf)
tf_set_en <- unique(tf_set_en[!is.na(tf_set_en)])

tf_targets_en <- collectri %>%
  filter(source %in% tf_set_en) %>%
  distinct(source, target)

target_pool_en <- as.character(unique(tf_targets_en$target))
cat("TFs matched:", length(intersect(tf_set_en, collectri$source)),
    "| target genes:", length(target_pool_en), "\n")

# ---- 2. Overlap ----
en_overlap   <- base::intersect(as.character(en_genes),    as.character(target_pool_en))
cat("\nEN genes that are targets of selected TFs:\n");    print(en_overlap)

# ---- 3. Which TF regulates which overlapping gene ----
links <- tf_targets_en %>%
  filter(target %in% union(lasso_overlap, en_overlap)) %>%
  arrange(target, source)
print(links)

