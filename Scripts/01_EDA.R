library(dplyr)
library(tidyverse)
dataset <- read.delim("/rds/projects/b/bravol-inf-multi/Esther/Data/E-MTAB-14509/E-MTAB-14509.sdrf.txt",
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
#EDA
head(dataset)
dim(dataset)
colnames(dataset)
table(dataset$'Characteristics[organism part]')
table(dataset$'Characteristics[drug]')
summary(dataset$`Characteristics[psoriasis area severity index (pasi)]`)
summary(dataset$`Characteristics[disease onset type]`)
unique(dataset$`Characteristics[developmental stage]`)
length(unique(dataset$`Characteristics[developmental stage]`))
length(unique(dataset$`Characteristics[disease onset type]`))
unique(dataset$`Characteristics[disease onset type]`)
table(dataset$`Characteristics[hla_c0602 carrier]`)
table(dataset$`Factor Value[time]`)
psoriasis<-dataset[dataset$`Characteristics[organism part]`!='Whole Blood',]
psoriasis<-psoriasis[psoriasis$`Factor Value[time]`!='1',]
table(psoriasis$`Characteristics[organism part]`)
table(psoriasis$`Characteristics[drug]`)
table(psoriasis$`Characteristics[hla_c0602 carrier]`)
table(psoriasis$'Factor Value[time]')
library(ggplot2)
options(device = 'RStudioGD')
options(device = 'png')
install.packages("ragg")
options(device = ragg::agg_png)
names(psoriasis)[duplicated(names(psoriasis))]
names(psoriasis) <- make.unique(names(psoriasis))

psoriasis<-na.omit(psoriasis)

#psoriasis dataset week 0
psoriasis0 <- psoriasis[
  !(psoriasis$`Factor Value[time]` %in% c(1, 12)),
]
names(psoriasis0) <- make.unique(names(psoriasis0))
#psoriasis dataset week 12
psoriasis12 <- psoriasis[
  !(psoriasis$`Factor Value[time]` %in% c(0,1)),
]
names(psoriasis12) <- make.unique(names(psoriasis12))


#paired boxplot
pasi_plot <- psoriasis %>%
  mutate(
    patient_id = sub("-.*", "", `Source Name`),
    PASI = as.numeric(
      `Characteristics[psoriasis area severity index (pasi)]`
    ),
    timepoint = factor(`Factor Value[time]`)
)    %>%
  filter(timepoint %in% c("0", "12")) %>%
  distinct(patient_id, timepoint, .keep_all = TRUE) %>%
  filter(!is.na(PASI))
 
ggplot(pasi_plot,aes(x = timepoint,y = PASI,group = patient_id)) +
  geom_boxplot(aes(group = timepoint),fill = "skyblue",alpha = 0.5,outlier.shape = NA) +
  geom_line(alpha =0.15,colour = "grey40") +
  geom_point(aes(colour = timepoint),size = 2,alpha =0.1 ) +
  labs(title = "Change in PASI score between Week 0 and Week 12",x = "Timepoint",y = "PASI score")+
  theme_minimal()+ theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),axis.line = element_line(),           # keep axis lines
    axis.ticks = element_line()           # keep axis ticks
  )
#boxplot early vs late onset
ggplot(psoriasis, aes(
  x = `Characteristics[disease onset type]`,
  y = `Characteristics[psoriasis area severity index (pasi)]`
)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +ggtitle("PASI Score according to disease onset")+
  labs(
    x = "Psoriasis Onset Type (Early vs Late)",
    y = "PASI Score"
  )

table(psoriasis$`Derived Array Data File`)
expr1<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/E-MTAB-14509/Skin_norm_counts_d.txt',
                  header = TRUE,
                  check.names = FALSE)

expr1$ensembl_id<-rownames(expr1)
biomart<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/annotations/Final_Annotation_List_BioMart.tsv',
                    check.names = FALSE)

expr1$ensembl_id <- rownames(expr1)

expr1$ensembl_id_clean <- sub(
  "\\..*",
  "",
  expr1$ensembl_id)


expr1_map <- merge(
  expr1,
  biomart[, c("ENSG.ID", "Gene.name")],   
  by.x = "ensembl_id_clean",
  by.y = "ENSG.ID",
  all.x = TRUE
)
head(expr1_map)
colnames(expr1_map)
expr1_map <- expr1_map[, !colnames(expr1_map) %in% c("ensembl_id_clean")]

#PCA
expr1_pca<-expr1_map[,!(colnames(expr1_map) %in% c('Gene.name','ensembl_id'))]
expr1_pca<-as.data.frame(lapply(expr1_pca,as.numeric))
rownames(expr1_pca)<-expr1_map$ensembl_id
expr1_t<-t(expr1_pca)
pca<-prcomp(expr1_t,scale. = TRUE)

pca_df<-as.data.frame(pca$x)
pca_df$sample<-rownames(pca_df)
pca_df$sample_id<-rownames(pca_df)
pca_df$sample_id<-sub('^X','',pca_df$sample_id)
small<-psoriasis[,c('Source Name', 'Characteristics[organism part]')]
colnames(small) <- c("sample_id", "Group")
pca_df$sample_id <- as.character(pca_df$sample_id)
small$sample_id <- as.character(small$sample_id)
small$sample_id <- gsub("-", ".", small$sample_id)
pca_df<-merge(pca_df,small,by='sample_id', all.x = TRUE)
pca_df<-merge(pca_df,small,by='sample_id', all.x = TRUE)
pca_df[is.na(pca_df$Group.y),
       c("sample_id")]
pca_df$sample_id[is.na(pca_df$Group.y)] %in%
  small$sample_id

pca_df <- subset(pca_df,
                      !is.na(Group.y))

#plot PCA
ggplot(pca_df,aes(PC1,PC2, colour = Group.y)) + 
  geom_point(size=3) +
  theme_classic()

ggplot(pca_df,aes(PC2,PC3, colour=Group.y)) + 
  geom_point(size=3) +
  theme_classic()

ggplot(pca_df,aes(PC3,PC4, colour=Group.y)) + 
  geom_point(size=3) +
  theme_classic()

ggplot(pca_df,aes(PC4,PC5, colour=Group.y)) + 
  geom_point(size=3) +
  theme_classic()

#scree plot
variance <- pca$sdev^2 / sum(pca$sdev^2)

variancedf <- data.frame(PC = factor(1:10),variance = variance[1:10]*100)

ggplot(variancedf, aes(x = PC, y = variance)) +
  geom_line(group = 1) +
  geom_point(size = 4) +
  xlab("Principal Component") +
  ylab("Variance (%)") +
  ggtitle("Scree Plot") + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),axis.line = element_line(),           # keep axis lines
    axis.ticks = element_line()           # keep axis ticks
+ ylim(0, 100))

 #week 0
small0 <- psoriasis0[, c("Source Name",
                         "Characteristics[organism part]")]

colnames(small0) <- c("sample_id","Group")

small0$sample_id <- gsub("-", ".", small0$sample_id)

pca_df0 <- merge(
  pca_df,
  small0,
  by = "sample_id",
  all.x = TRUE
)

pca_df0 <- subset(
  pca_df0,
  !is.na(Group)
)
ggplot(
  pca_df0,
  aes(PC1, PC2, colour = Group)
) +
  geom_point(size = 3) +
  theme_classic()

#week 12
small12 <- psoriasis12[, c("Source Name",
                           "Characteristics[organism part]")]

colnames(small12) <- c("sample_id","Group")

small12$sample_id <- gsub("-", ".", small12$sample_id)
pca_df12 <- merge(
  pca_df,
  small12,
  by = "sample_id",
  all.x = TRUE
)

pca_df12 <- subset(
  pca_df12,
  !is.na(Group)
)
ggplot(
  pca_df12,
  aes(PC1, PC2, colour = Group)
) +
  geom_point(size = 3) +
  theme_classic()

#loadings PC1 general
loadings<-pca$rotation
head(loadings)
pc1.df<-data.frame(gene=rownames(loadings),loading=loadings[,1])
pc1.df<-pc1.df[order(abs(pc1.df$loading),decreasing=TRUE),]
pc1_genes <- merge(
  pc1.df,
  biomart[, c("ENSG.ID", "Gene.name")],   
  by.x = "gene",
  by.y = "ENSG.ID",
  all.x = TRUE
)
pc1_genes <- pc1_genes[order(abs(pc1_genes$loading), decreasing = TRUE), ]

rownames(pc1_genes) <- pc1_genes$gene
pc1_genes<-pc1_genes[,!colnames(pc1_genes)%in%c('gene')]
write.table(pc1_genes, "/rds/projects/b/bravol-inf-multi/Esther/Data/loadings_pca.tsv", sep = "\t", row.names = T, quote=F)

#pca variance
pca$sdev
var<-(pca$sdev^2)/sum(pca$sdev^2)*100
var[1:10]

ggplot(pca_df, aes(PC1, PC2, colour = Group.y)) +
  geom_point(size = 3) +
  labs(
    x = paste0("PC1 (", round(var[1], 1), "%)"),
    y = paste0("PC2 (", round(var[2], 1), "%)")
  ) +
  theme_classic()
rankedlist<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Inflammatome/Inflammatome-main/data/04_rank_agg_list.tsv',header = TRUE,
                       check.names = FALSE)
                      
inflammatome<-rankedlist$ENSG.ID[1:100]
loadings$gene<-rownames(loadings)
loadings$inflammatome<-ifelse(loadings$gene %in% inflammatome,'Inflammatome','Other')

#why is there a diff number of samples in psoriasis than expr?
setdiff(psoriasis$`Source Name`, colnames(expr1))
