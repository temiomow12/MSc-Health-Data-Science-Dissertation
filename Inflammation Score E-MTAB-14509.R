# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load packages-----------------------------------------------------------------
list.of.packages <- c("ggplot2","dplyr","tidyr","limma","biomaRt", "ggrepel")
library(biomaRt)
lapply(list.of.packages, library, character.only=TRUE)
library(ggrepel)
options(device = 'RStudioGD')
options(device = 'png')
install.packages("ragg")
options(device = ragg::agg_png)

# Read the ranked list obtained by aggregation  ----------------------
ranked.list <- read.csv("/rds/projects/b/bravol-inf-multi/Esther/Inflammatome/Inflammatome-main/data/04_rank_agg_list.tsv",sep="\t",header=TRUE)

# Preprocessing ----------------------------------------------------------------
dataset <- read.delim("/rds/projects/b/bravol-inf-multi/Esther/Data/E-MTAB-14509/E-MTAB-14509.sdrf.txt",
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE,
                      check.names = FALSE)
psoriasis<-dataset[dataset$`Characteristics[organism part]`!='Whole Blood',]
biomart<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/annotations/Final_Annotation_List_BioMart.tsv',
                    check.names = FALSE)
expr1<-read.delim('/rds/projects/b/bravol-inf-multi/Esther/Data/E-MTAB-14509/Skin_norm_counts_d.txt',
                  header = TRUE,
                  check.names = FALSE)
expr1$ensembl_id <- rownames(expr1)
expr1_map <- merge(
  expr1,
  biomart[, c("ENSG.ID", "Gene.name")],   
  by.x = "ensembl_id",
  by.y = "ENSG.ID",
  all.x = TRUE
)

# Read initial data
df.prot = expr1_map

# Extract sample info based on column names
coldata <- psoriasis[, c("Source Name","Characteristics[organism part]")]

colnames(coldata) <- c("sample","condition")

coldata$sample <- gsub("\\.", "-", coldata$sample)

coldata <- coldata[
  match(colnames(df.prot), coldata$sample),
]
coldata<- coldata[-1, ]   
coldata <- na.omit(coldata)  
df.prot<-na.omit(df.prot)
rownames(df.prot) <- df.prot$ensembl_id
df.prot<-df.prot[,!colnames(df.prot)%in%c('ensembl_id','Gene.name')]
match(colnames(df.prot),coldata$sample)
samplecolumns<-colnames(df.prot)[
  !(colnames(df.prot)%in% c('ensembl_id','Gene.name'))
]
match(samplecolumns,coldata$sample)

# Reorder such as to fit the metadata
expr_match <- df.prot[, !(colnames(df.prot) %in% c("ensembl_id", "Gene.name"))]


# added by me: defining expr.psoriasis.markers.100
ranked.list$ENSG.ID[1:100]
expr.psoriasis.100 <- expr_match[ranked.list$ENSG.ID[1:100],] 
expr.psoriasis.100.mean <- apply(expr.psoriasis.100, 2, mean)
expr.psoriasis.100.m = expr.psoriasis.100.mean; type="mean"
coldata <- psoriasis[, c("Source Name","Characteristics[organism part]","Characteristics[psoriasis area severity index (pasi)]")]
colnames(coldata) <- c(
  "sample",
  "condition",
  "PASI"
)

coldata <- coldata[
 match(samplecolumns, coldata$sample),]

score <- as.numeric(coldata$PASI)
# build the score frame FIRST
data <- data.frame( sample = names(expr.psoriasis.100.mean), inflammation_score = as.numeric(expr.psoriasis.100.mean) )

coldata <- coldata[ match(samplecolumns, coldata$sample),] 

# merge score and PASI on the shared sample key 
merged <- merge( data, coldata, by.x = "sample", by.y = "sample" )

pearson  <- cor.test(merged$inflammation_score,
                     merged$`PASI`,
                     method = "pearson")
spearman <- cor.test(merged$inflammation_score,
                     merged$`PASI`,
                     method = "spearman")

pearson; spearman

write.table(data, "/rds/projects/b/bravol-inf-multi/Esther/Data/inflammation_score_psoriasis.tsv", sep = "\t", row.names = T, quote=F)


ggplot(data, aes(y = expr.psoriasis.100.m, x = score)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", se = FALSE) + 
  labs(y = "Inflammation signature-based score", x = "PASI Score", title = "") +
  theme_classic() +
  # Add correlation coefficient as text
  annotate("text", y = max(expr.psoriasis.100.mean) + 0.3 , x = max(score), 
           label = paste("Pearson r =", round(pearson$estimate, 2)), 
           hjust = 1.5, vjust = 3, color = "red", size = 3.5)

