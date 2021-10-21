library("qusage")
library("GSVA")
library("GSEABase")
library("GSVAdata")
library("Biobase")
library("stringr")
library("genefilter")
library("limma")
source("~/SeneSys/Scripts/Visualization_colors.R")

#meta_info = read.table("~/MAPTor_NET//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = meta_info[meta_info$Sample!="",]
meta_info = meta_info[!is.na(meta_info$Sample),]
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info = meta_info[meta_info$Study != "GSE11318",]
rownames(meta_info) = meta_info$Sample


#expr_raw = read.table("~/MAPTor_NET/BAMs/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t", as.is = T,header = T, row.names = 1, fill = T)
expr_raw = read.table("~/SeneSys/Data/Schmitz.HGNC.tsv",sep ="\t", as.is = T,header = T, row.names = 1, fill = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

expr_raw = expr_raw[ ,which(colnames(expr_raw) %in% meta_info$Sample ) ]
meta_data = meta_info[colnames(expr_raw),]

gmt_file = read.gmt("~/SeneSys/Misc/SeneSys_gene_sets.gmt")

row_var = apply(expr_raw, FUN = var, MARGIN = 1)
summary(row_var)
expr = expr_raw[row_var > median(row_var),]

fe_es = gsva(as.matrix(expr), gmt_file, min.sz=10, max.sz=500, verbose=TRUE)
#write.table(fe_es,"~/SeneSys/Results/Schmitz.GSVA.tsv",sep ="\t",quote = F,row.names = T)

#####

fe_es = read.table("~/MAPTor_NET/Results/GSVAR_senesys_S57_not_normalized_hgnc.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
colnames(fe_es) = str_replace(colnames(fe_es),pattern = "^X","")
fe_es[1:5,1:5]

vis_mat = fe_es
#vis_mat_var = names(sort(apply(vis_mat,MARGIN = 1, FUN =var),decreasing = T))[1:20]
vis_mat_var = order(vis_mat["SUVARNESS",])
#quantile(vis_mat_var,seq(0,1,by=.1))[10]
#vis_mat = vis_mat[vis_mat_var >= quantile(vis_mat_var,seq(0,1,by=.01))[100],]
vis_mat = vis_mat[,vis_mat_var]
#vis_mat = vis_mat[,order(vis_mat["E2F target genes",])]
#meta_data = meta_info[colnames(t),]
#vis_mat = vis_mat[,order(meta_data[,"B"])]
meta_data[colnames(vis_mat),"SUVARNESS"] = vis_mat["SUVARNESS",]
meta_data["ABC_GCB"] = str_replace_all(meta_data[,"ABC_GCB"],pattern = " ", "")

##

