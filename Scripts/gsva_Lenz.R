# Dieses Skript dient zur sampleweisen Exploration von SUVARness im Reddy Dataset
# nutze dafür die gct Tabelle von Reddy (Salmon_DESEQ4.gct)
# und die bereits existierende gmt Datei mit den gene sets (SeneSys_gene_sets.gmt)

setwd("/home/fattohim/GSVA")
library(GSVA)

# read in expression data as matrix object
expr_18 <- as.matrix(read.table("/home/fattohim/Lenz_data/GSE11318/GSE11318.tsv", header = T,  
                                 row.names = 1))
expr_46 <- as.matrix(read.table("/home/fattohim/Lenz_data/GSE10846/GSE10846.tsv", header = T,  
                                row.names = 1))

# turn ensemble ids into hgnc ids
ensembl_18 <- rownames(expr_18)
ensembl_46 <- rownames(expr_46)

library("biomaRt")
# prepare for getBM function 
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")

# turn ensembl ids into hgnc symbols
results18 <- getBM(attributes=c('hgnc_symbol', 'external_gene_name', 'ensembl_gene_id'), 
                 filters = 'ensembl_gene_id', 
                 values = ensembl_18, 
                 mart = mart)
# match the rownames of the expr matrix to the hgnc translation
match_idx <- match(rownames(expr_18), results18$ensembl_gene_id)
#na_idx <- which(is.na(match(rownames(expr_18), results18$ensembl_gene_id)))

# function that converts the match idx to the hgnc symbol or to NA 
# (which is later removed from matrix)
match_ens_hgnc <- function(match_vc){
  new_rownames <- c()
  for (i in 1:length(match_vc)) {
    if(!is.na(match_vc[i])){
      new_rownames[i] <- results18$external_gene_name[match_vc[i]]
    }
    else{
      new_rownames[i] <- NA
    }
  }
  return(new_rownames)
}

new_rownms <- match_ens_hgnc(match_idx)
rownames(expr_18) <- new_rownms
rownames(expr_46) <- new_rownms
expr_18 <- expr_18[-which(is.na(rownames(expr_18))),]
expr_46 <- expr_46[-which(is.na(rownames(expr_46))),]


# read in gene sets as list object
library(limma)
library(qusage)
gmt_file <- read.gmt("SeneSys_gene_sets.gmt")

# gsva with normally distributed enrichment scores
lenz_fe_es18 <- gsva(expr_18, gmt_file, min.sz=10, max.sz=500, verbose=TRUE)

# gsva with bimodally distributed enrichment scores
lenz_msi_18 <- gsva(expr_18, gmt_file, min.sz=10, max.sz=500, verbose=TRUE, mx.diff = F)


### export lenz_fe_es18 und lenz_msi_es18 ###

write.table(lenz_fe_es18, file = "Lenz18_GSVA_ES_normal.gct", quote = F, sep = "\t", 
            row.names = T, col.names = T)

write.table(lenz_msi_18, file = "Lenz18_GSVA_ES_bimodal.gct", quote = F, sep = "\t", 
            row.names = T, col.names = T)



# gsva with normally distributed enrichment scores
lenz_fe_es46 <- gsva(expr_46, gmt_file, min.sz=10, max.sz=500, verbose=TRUE)

# gsva with bimodally distributed enrichment scores
lenz_msi_46 <- gsva(expr_46, gmt_file, min.sz=10, max.sz=500, verbose=TRUE, mx.diff = F)


### export lenz_fe_es18 und lenz_msi_es18 ###

write.table(lenz_fe_es46, file = "Lenz46_GSVA_ES_normal.gct", quote = F, sep = "\t", 
            row.names = T, col.names = T)

write.table(lenz_msi_46, file = "Lenz46_GSVA_ES_bimodal.gct", quote = F, sep = "\t", 
            row.names = T, col.names = T)


