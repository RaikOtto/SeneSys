library(devtools)
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")

meta_info = read.table("~/SeneSys///Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$NEC_NET_PCA
res_scdc = as.data.frame(meta_info)

table(meta_info$ABC_GCB)

expr_raw = read.table("~/SeneSys/Data/GSE98588.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
#expr_raw = read.table("~/SeneSys/Data/GSE98588.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = exp(expr_raw)
expr_raw[1:5,1:5]
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

# prediction
basis = readRDS("~/artdeco/inst/Models/bseqsc/PMBC.RDS")

bseqsc_fit = bseqsc::bseqsc_proportions(
    expr_raw,
    basis[[1]],
    verbose = FALSE,
    absolute = FALSE,
    log = FALSE,
    perm = 1000
)

results_tab = cbind ( t(bseqsc_fit$coefficients), bseqsc_fit$stats)
#write.table(results_tab,"~/SeneSys/Results/GSE98588.Cell_fraction_predictions.bseq-sc.tsv",sep ="\t", quote =F)
