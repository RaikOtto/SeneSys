library(qusage)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library("stringr")
library(genefilter)
library(limma)
source("~/SeneSys/Scripts/Visualization_colors.R")

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = meta_info[meta_info$Sample!="",]
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

expr_raw = read.table("~/SeneSys/Data/Reddy.tsv",sep ="\t", as.is = T,header = T, row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
meta_data = meta_info[colnames(expr_raw),]

gmt_file = read.gmt("~/SeneSys/Misc/SeneSys_gene_sets.gmt")

row_var = apply(expr_raw, FUN = var, MARGIN = 1)
summary(row_var)
expr = expr_raw[row_var > median(row_var),]

fe_es = gsva(as.matrix(expr), gmt_file, min.sz=10, max.sz=500, verbose=TRUE)
#write.table(fe_es,"~/SeneSys/Results/Reddy.gsva.tsv",sep ="\t",quote = F,row.names = F)

#####

#fe_es = read.table("~/SeneSys/Results/GSVA_ES_bimodal.csv",sep ="\t", header = T,as.is = T,row.names = 1)
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
#meta_data$Drup_response = meta_data$ABC_GCB

#pdf("~/Downloads/E2F_and_other_pathways_of_interest.pdf")
pheatmap::pheatmap(
    #pca,
    vis_mat,
    #annotation_col = meta_data[,c("Drup_response","Drup_response")],
    annotation_col = meta_data["ABC_GCB"],
    annotation_colors = aka3,
    show_rownames = T,
    show_colnames = F,
    treeheight_col = 0,
    legend = F,
    fontsize_col = 7,
    cluster_cols = T,
    cluster_rows = T,
    clustering_method = "average"
)
#dev.off()
rownames(vis_mat)
meta_data = meta_info[colnames(vis_mat),]

num_vec = as.double(vis_mat["SUVARNESS",] * -1)
vis_vec = num_vec
vis_vec[ num_vec> 0 ] = "high"
vis_vec[num_vec <= 0 ] = "low"

d = as.data.frame(cbind(vis_vec,meta_data[,c("Cluster")]))
#pheatmap::pheatmap(table(d))

aggregate(as.double(vis_mat["SUVARNESS",]), by = list(meta_data[,"Cluster"]), FUN = mean)

plot(agg_vec)
aggregate(num_vec, by = list(meta_data[,"ABC_GCB"]), FUN = mean)

meta_info[colnames(vis_mat),"SUVARNESS"] = vis_mat["SUVARNESS",]

#write.table(meta_info,"~/SeneSys/Misc/Meta_information.tsv",sep ="\t",quote = F,row.names = F)
