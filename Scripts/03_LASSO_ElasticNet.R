library(formatR)
library(glmnet)
library(pROC)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
options(device = 'RStudioGD')
options(device = 'png')
install.packages("ragg")
options(device = ragg::agg_png)

#load and clean datasets
expr<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/E-MTAB-14509/Skin_norm_counts_d.txt',header = TRUE,
                 check.names = FALSE)
expr$ensembl_id<-rownames(expr)
expr$ensembl_id <- sub(
  "\\..*",
  "",
  expr$ensembl_id
)

biomart<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/annotations/Final_Annotation_List_BioMart.tsv',header = TRUE,
                    check.names = FALSE)
expr_mapped<-merge(
  expr,
  biomart[,c('ENSG.ID','Gene.name')],
  by.x = 'ensembl_id',
  by.y = 'ENSG.ID',
  all.x = TRUE
)
rownames(expr_mapped)<-expr_mapped$ensembl_id
expr_mapped$ensembl_id <- NULL
dataset<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/E-MTAB-14509/E-MTAB-14509.sdrf.txt',header = TRUE,
                    check.names = FALSE,)
psoriasis<-dataset[dataset$`Characteristics[organism part]`!='Whole Blood',]
psoriasis0 <- psoriasis[
  !(psoriasis$`Factor Value[time]` %in% c(1, 12)),
]

#week 0 samples only
week0.samples<-psoriasis0$`Source Name`
week0.samples<-intersect(week0.samples,colnames(expr_mapped))
expr_week0<-expr_mapped[,week0.samples]

#formatting the dataset for modelling
X<-t(as.matrix(expr_week0))
X<-X[,apply(X,2,var,na.rm=TRUE)>0] #keeps columns that variance is above 0
psoriasis0_index<-match(rownames(X),psoriasis0$'Source Name')
y<- psoriasis0$`Characteristics[organism part]`[match(rownames(X),psoriasis0$`Source Name`)] #for each sample in the rows of X, find the position of that same name in the metadata's Source name column- aligns by key
y<-factor(y) #converts labels (i.e.lesional/non-lesional) into factors

source_name<-psoriasis0$'Source Name'[psoriasis0_index]
#70/30 random split
set.seed(42)  #makes it reproducible for cross validation
train_index_vec<-sample(seq_len(nrow(X)),size = floor(0.70*nrow(X)))
train_index<-seq_len(nrow(X))%in%train_index_vec
test_index<-!train_index

X_train<-X[train_index,,drop=FALSE]
X_test<-X[test_index,,drop=FALSE]
y_train<-y[train_index]
y_test<-y[test_index]

#remove zero variance genes
keep_genes<-apply(X_train,2,var,na.rm=TRUE)>0
X_train<-X_train[,keep_genes,drop=FALSE]
X_test<-X_test[,keep_genes,drop=FALSE]

#LASSO model
cv.lasso<-cv.glmnet(x=X_train,y=y_train,family='binomial',alpha=1,type.measure='auc',nfolds=10,keep=TRUE) #'family=binomial' is for logistic regression, 'alpha=1' for LASSO specifically (L1 penalty), 'nfolds=10' 10 folds for cross validation
max(cv.lasso$cvm)

#Apply to untouched 30%
test_probability<-as.numeric(predict(cv.lasso,newx=X_test,s='lambda.1se',type='response'))
test_prediction<-factor(ifelse(test_probability>=0.5,'Lesional Skin','Nonlesional Skin'),levels = levels(y))

#confusion matrix
confusion_matrix<-table(Predicted=test_prediction, Actual=y_test)
confusion_matrix
test_accuracy<-mean(test_prediction==y_test)
test_accuracy

#test-set ROC and AUC
test_roc<-roc(response=y_test,predictor=test_probability, levels=levels(y),direction='<')
auc(test_roc)
ci.auc(test_roc)

plot(test_roc,main=paste0('LASSO test-set ROC:AUC =',
                          round(as.numeric(auc(test_roc)),3)))

#signature from training model
lasso.coef<-coef(cv.lasso,s='lambda.1se')
lasso.genes <- setdiff(rownames(lasso.coef)[which(as.numeric(lasso.coef) != 0)],"(Intercept)")
lasso.genes
coef_values_lasso <- as.numeric(lasso.coef)
feature_names_lasso <- rownames(lasso.coef)

keep_lasso <- coef_values_lasso != 0 & feature_names_lasso != "(Intercept)"
lasso_df <- data.frame(
  ensembl = feature_names_lasso[keep_lasso],
  coefficient = coef_values_lasso[keep_lasso],
  stringsAsFactors = FALSE)

lasso_df$HGNC <- expr_mapped$Gene.name[
  match(lasso_df$ensembl, rownames(expr_mapped))]
lasso_df <- lasso_df[
  order(-abs(lasso_df$coefficient)),
]
write.table(lasso_df,'/rds/projects/b/bravol-inf-multi/Esther/Data/lasso_coef',sep = '\t',row.names = T,quote = F)
# GJB2
gjb2_df <- data.frame(expression = X[, "ENSG00000165474"], status = y)
ggplot(gjb2_df, aes(x = status, y = expression)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  ggtitle("GJB2") +
  labs(x = "Lesion status", y = "Expression (log-normalised)")

# KYNU
kynu_df <- data.frame(expression = X[, "ENSG00000115919"], status = y)
ggplot(kynu_df, aes(x = status, y = expression)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  ggtitle("KYNU") +
  labs(x = "Lesion status", y = "Expression (log-normalised)")

# KLK13
klk13_df <- data.frame(expression = X[, "ENSG00000167759"], status = y)
ggplot(klk13_df, aes(x = status, y = expression)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  ggtitle("KLK13") +
  labs(x = "Lesion status", y = "Expression (log-normalised)")

#fit the model: EN

cv.en <- cv.glmnet(x = X_train, y = y_train, family = "binomial",
                   alpha = 0.5, type.measure = "auc", nfolds = 10, keep = TRUE)
max(cv.en$cvm)

# Apply to untouched 30%
test_probability_en <- as.numeric(
  predict(cv.en, newx = X_test, s = "lambda.1se", type = "response"))

test_prediction_en <- factor(
  ifelse(test_probability_en >= 0.5, "Lesional Skin", "Nonlesional Skin"),
  levels = levels(y))

# Confusion matrix + accuracy (EN)
confusion_matrix_en <- table(Predicted = test_prediction_en, Actual = y_test)
confusion_matrix_en
test_accuracy_en <- mean(test_prediction_en == y_test)
test_accuracy_en

# Test-set ROC + AUC (EN)
test_roc_en <- roc(response = y_test, predictor = test_probability_en,
                   levels = levels(y), direction = "<")
auc(test_roc_en)
ci.auc(test_roc_en)
plot(test_roc_en, main = paste0("Elastic Net test-set ROC: AUC = ",
                                round(as.numeric(auc(test_roc_en)), 3)))

# Signature from the EN training model
en.coef <- coef(cv.en, s = "lambda.1se")
en.genes <- setdiff(rownames(en.coef)[which(as.numeric(en.coef) != 0)],
                    "(Intercept)")

coef_values_en  <- as.numeric(en.coef)
feature_names_en <- rownames(en.coef)
keep_en <- coef_values_en != 0 & feature_names_en != "(Intercept)"

en_df <- data.frame(
  ensembl = feature_names_en[keep_en],
  coefficient = coef_values_en[keep_en],
  stringsAsFactors = FALSE)

en_df$HGNC <- expr_mapped$Gene.name[match(en_df$ensembl, rownames(expr_mapped))]
en_df <- en_df[order(-abs(en_df$coefficient)), ]
en_df
write.table(en_df,'/rds/projects/b/bravol-inf-multi/Esther/Data/en_coef',sep = '\t',row.names = T,quote = F)

#how many genes did each model keep?
length(lasso.genes);length(en.genes)

#linking to inflammatome
rankedlist<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Inflammatome/Inflammatome-main/data/04_rank_agg_list.tsv',header=TRUE,check.names= FALSE)
inflammatome100<-rankedlist$Gene.name[1:100]
inflammatome2000<-rankedlist$Gene.name[1:2000]

lasso.symbols<- expr_mapped$Gene.name[match(lasso.genes, rownames(expr_mapped))]
en.symbols<-expr_mapped$Gene.name[match(en.genes,rownames(expr_mapped))]
sum(lasso.symbols %in% inflammatome100)
sum(lasso.symbols %in% inflammatome2000)
sum(en.symbols %in% inflammatome100)
sum(en.symbols %in% inflammatome2000)

overlap_lasso<-intersect(lasso.symbols, inflammatome2000)
overlap_en<-intersect(en.symbols, inflammatome2000)

write.table(overlap_en,'/rds/projects/b/bravol-inf-multi/Esther/Data/overlap_en_genes',sep = '\t',row.names = T,quote = F)
write.table(overlap_lasso,'/rds/projects/b/bravol-inf-multi/Esther/Data/overlap_lasso_genes',sep = '\t',row.names = T,quote = F)


#psoriasis dataset week 12
psoriasis12 <- psoriasis[
  !(psoriasis$`Factor Value[time]` %in% c(0,1)),
]
names(psoriasis12) <- make.unique(names(psoriasis12))
week12_samples<-psoriasis12$`Source Name`
week12_samples<-intersect(week12_samples,colnames(expr_mapped))
expr_week12<-expr_mapped[,week12_samples]
x_12 <- t(as.matrix(expr_week12))
x_12<-x_12[,colnames(X_train)]
y_12<- psoriasis12$`Characteristics[organism part]`[match(rownames(x_12),psoriasis12$`Source Name`)]
y_12<-factor(y_12) #converts labels (i.e.lesional/non-lesional) into factors

stopifnot(all(colnames(x_12) == colnames(X)))

#pred week 12 lasso
pred12_lasso<- as.numeric(predict(cv.lasso, newx = x_12,
                                  s = "lambda.1se", type = "response"))
#pred week 12 en
pred12_en<- as.numeric(predict(cv.en, newx = x_12,
                               s = "lambda.1se", type = "response"))
#roc
roc12_lasso<- roc(y_12, pred12_lasso,levels=levels(y),direction='<',quiet = TRUE)
roc12_en<- roc(y_12, pred12_en, levels=levels(y),direction='<',quiet = TRUE)

ci.auc(roc12_lasso)     # AUC + 95% CI, LASSO on week 12
ci.auc(roc12_en) #AUC +95%, EN week 12

#use full week-0 X so every patient as a week0 score
pred0_lasso <- as.numeric(predict(cv.lasso, newx = X[,colnames(X_train)],
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
  labs(title = sprintf("Week 12 ROC (frozen model)Genes: AUC %.3f",
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
  labs(title = sprintf("Week 12 ROC (frozen model) Elastic Net Genes: AUC %.3f",
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
       subtitle = "Frozen week-0 model applied to week 12; paired by patient; Genes",
       x = NULL, y = "Predicted probability of lesional") +
  theme_minimal()

