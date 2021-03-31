library(devtools)
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")

#meta_info = read.table("~/SeneSys///Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = read.table("~/Deko_Projekt//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$NEC_NET_PCA
res_scdc = as.data.frame(meta_info)

table(meta_info$ABC_GCB)

expr_raw = read.table("~/Deko_Projekt/Data/Human_differentiated_pancreatic_islet_cells_Bulk/GSE142720_rma_norm_log2_matrix.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = T)
#expr_raw = read.table("~/SeneSys/Data/GSE98588.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
col_names = colnames(expr_raw)
row_names = rownames(expr_raw)
expr_raw = as.data.frame(matrix(as.double(as.character(unlist(expr_raw))),ncol = length(col_names)))
colnames(expr_raw) = col_names
rownames(expr_raw) = row_names
expr_raw[1:5,1:5]

# prediction
#basis = readRDS("~/artdeco/inst/Models/bseqsc/PMBC.RDS")
basis = readRDS("~/artdeco/inst/Models/bseqsc/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.RDS")
basis = basis[[1]]
table(rownames(basis) %in% rownames(expr_raw))

basis_2 = rbind(basis,matrix(colnames(basis),nrow=1))
basis_2 = cbind(rownames(basis_2),basis_2)

expr_raw_2 = rbind(colnames(expr_raw),expr_raw)
expr_raw_2 = cbind(rownames(expr_raw_2),expr_raw_2)

results <- CIBERSORT(
    "~/Downloads/Basis.tsv",
    '~/Deko_Projekt/Data/Human_differentiated_pancreatic_islet_cells_Bulk/GSE142720_rma_norm_log2_matrix.HGNC.tsv',
    1000,
    TRUE,
    FALSE,
    'sig.score')

bseqsc_fit = bseqsc_proportions(
    expr_raw_2,
    reference = basis_2,
    verbose = FALSE,
    absolute = FALSE,
    log = FALSE,
    perm = 1000
)

results_tab = cbind ( t(bseqsc_fit$coefficients), bseqsc_fit$stats)
#write.table(results,"~/Deko_Projekt/Results/Bseq.7_Japan.tsv",sep ="\t", quote =F)
