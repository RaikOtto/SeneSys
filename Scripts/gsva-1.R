# Dieses Skript dient zur Erprobung des GSVA packages.
# Functional Enrichment und 
# Molecular signature identification

setwd("~/GSVA")
library(GSVA)

# read in expression data as matrix object
gct_file <- as.matrix(read.table("GSE98588.gct", header = T, 
                                 skip = 2, row.names = 1))
gct_file <- gct_file[,-1]

# read in gene sets as list object
library(limma)
library(qusage)
gmt_file <- read.gmt("SeneSys_gene_sets.gmt")

# read in phenotype file (suvarness level)
cls_file <- read.table("SUVARness.cls", header = F, skip = 2)
cls_file <- unname(unlist(c(cls_file[1,])))

# read in phenotype file (dlbcl cluster)
cls_file_dlbcl <- read.table("clusters.cls", header = F, skip = 2)
cls_file_dlbcl <- unname(unlist(c(cls_file_dlbcl[1,])))

####### functional enrichment #######
# Because we want to use limma on the resulting GSVA enrichment scores, 
# we leave deliberately unchanged the default argument mx.diff=TRUE 
# to obtain approximately normally distributed ES
fe_es <- gsva(gct_file, gmt_file, min.sz=10, max.sz=500, verbose=TRUE)
# We test whether there is a difference between the GSVA enrichment scores 
# from each pair of phenotypes using a simple linear model and moderated 
# t-statistics computed by the limma package using an 
# empirical Bayes shrinkage method
adjPvalueCutoff <- 0.001

# Test SUVARness Level
design <- model.matrix(~ cls_file)
colnames(design) <- c("ALL", "HighvsLow")
fit <- lmFit(fe_es, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="HighvsLow", number=Inf)
DEgeneSets <- topTable(fit, coef="HighvsLow", number=Inf, 
                       p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)
#        ALL HighvsLow
#Down     8        25
#NotSig  16         6
#Up      13         6

# visualize in heatmap 
# library(ggplot2)
# library(tidyr)
# library(reshape2)
# 
# melted_es <- melt(fe_es)
# colnames(melted_es) <- c("pathway", "sample", "ES")
# 
# ggplot(data = melted_es, aes(x=sample, y=pathway, fill=ES)) + 
#   geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_y_discrete(label = function(x) stringr::str_trunc(x, 12)) 

# Annotation data frame
library(ComplexHeatmap)
fe_es_short <- fe_es
rownames(fe_es_short) <- sapply(rownames(fe_es_short), 
                                function(x) stringr::str_trunc(x, 12))

# Combine the heatmap and the annotation
ht = Heatmap(fe_es_short, name = "GSVA ES", column_split = cls_file, 
             row_names_gp = gpar(fontsize = 7), 
             column_names_gp = gpar(fontsize = 4))

draw(ht, column_title = "Samples", column_title_side = "bottom",
     row_title = "Pathway")



# test DLBCL Cluster
design2 <- model.matrix(~ factor(cls_file_dlbcl))
colnames(design2) <- c("ALL", "1vsRest", "2vsRest", "3vsRest", "4vsRest", "5vsRest")
fit2 <- lmFit(fe_es, design2)
fit2 <- eBayes(fit2)
allGeneSets2 <- topTable(fit2,  number=Inf)
DEgeneSets2 <- topTable(fit2, number=Inf, p.value=adjPvalueCutoff, adjust="BH")
res2 <- decideTests(fit2, p.value=adjPvalueCutoff)
summary(res2)
#        ALL 1vsRest 2vsRest 3vsRest 4vsRest 5vsRest
#Down     1       0       4       7       0       0
#NotSig  36      37      33      30      37      37
#Up       0       0       0       0       0       0

ht2 = Heatmap(fe_es_short, name = "GSVA ES", column_split = cls_file_dlbcl, 
             row_names_gp = gpar(fontsize = 7), 
             column_names_gp = gpar(fontsize = 4))

draw(ht2, column_title = "Samples", column_title_side = "bottom",
     row_title = "Pathway")

# Which pathways/gene sets show differential enrichment between phenotypes?
res_df <- as.data.frame(res)
down_res <- rownames(res_df[which(res_df[2] == -1),])
up_res <- rownames(res_df[which(res_df[2] == 1),])

#res_df2 <- as.data.frame(res2)
#down_res2 <- rownames(res_df2[which(res_df2[2] == -1),])
#up_res2 <- rownames(res_df2[which(res_df2[2] == 1),])

####### molecular signature identification #######
# employing GSVA to transform the gene expression measurements into
# enrichment scores for the gene sets, without taking the 
# sample subtype grouping into account
msi_es <- gsva(gct_file, gmt_file, min.sz=10, max.sz=500, verbose=TRUE, mx.diff = F)

# visualize to show the GSVA enrichment scores obtained for 
# the up-regulated gene sets across the samples of the subtypes.

msi_es_short <- msi_es
rownames(msi_es_short) <- sapply(rownames(msi_es_short), 
                                function(x) stringr::str_trunc(x, 12))

ht3 = Heatmap(msi_es_short, name = "GSVA ES", column_split = cls_file, 
              row_names_gp = gpar(fontsize = 7), 
              column_names_gp = gpar(fontsize = 4))

draw(ht3, column_title = "Samples", column_title_side = "bottom",
     row_title = "Pathway")

ht4 = Heatmap(msi_es_short, name = "GSVA ES", column_split = cls_file_dlbcl, 
              row_names_gp = gpar(fontsize = 7), 
              column_names_gp = gpar(fontsize = 4))

draw(ht4, column_title = "Samples", column_title_side = "bottom",
     row_title = "Pathway")

### export fe_es und msi_es ###

write.table(fe_es, file = "GSVA_ES_normal.gct", quote = F, sep = "\t", 
            row.names = T, col.names = T)

write.table(msi_es, file = "GSVA_ES_bimodal.gct", quote = F, sep = "\t", 
            row.names = T, col.names = T)

