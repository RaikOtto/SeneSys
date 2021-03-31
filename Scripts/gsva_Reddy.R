# Dieses Skript dient zur sampleweisen Exploration von SUVARness im Reddy Dataset
# nutze dafür die gct Tabelle von Reddy (Salmon_DESEQ4.gct)
# und die bereits existierende gmt Datei mit den gene sets (SeneSys_gene_sets.gmt)

setwd("/home/fattohim/GSVA")
library(GSVA)

# read in expression data as matrix object
gct_file <- as.matrix(read.table("/home/fattohim/Reddy_data/Salmon_DESEQ4.gct", header = T,  
                                 row.names = 1))
gct_file <- gct_file[,-1]

# read in gene sets as list object
library(limma)
library(qusage)
gmt_file <- read.gmt("SeneSys_gene_sets.gmt")

# gsva with normally distributed enrichment scores
reddy_fe_es <- gsva(gct_file, gmt_file, min.sz=10, max.sz=500, verbose=TRUE)

# gsva with bimodally distributed enrichment scores
reddy_msi_es <- gsva(gct_file, gmt_file, min.sz=10, max.sz=500, verbose=TRUE, mx.diff = F)


### export reddy_fe_es und reddy_msi_es ###

write.table(reddy_fe_es, file = "Reddy_GSVA_ES_normal.gct", quote = F, sep = "\t", 
            row.names = T, col.names = T)

write.table(reddy_msi_es, file = "Reddy_GSVA_ES_bimodal.gct", quote = F, sep = "\t", 
            row.names = T, col.names = T)
##########################GSE10846####################

install.packages("GSE10846")
data("GSE10846")

