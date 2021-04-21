library("stringr")
library("ggplot2")
library("dplyr")
library("grid")

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

source("~/SeneSys/Scripts/Visualization_colors.R")
meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = meta_info %>% dplyr::filter(Study != "GSE11318")
meta_info = meta_info[meta_info$Sample!="",]
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.DeSeq2.HGNC.tsv"
#expr_raw = read.table("~/SeneSys/Results/Mouse_GSVA_ES_normal.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
#expr_raw = expr_raw[,!(colnames(expr_raw) %in% "GSM2601431")]
meta_data = meta_info[colnames(expr_raw),]
#meta_data$Cluster = as.factor(meta_data$Cluster)
expr_raw[1:5,1:5]

<<<<<<< HEAD
genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/Senescence geneset.tsv",sep ="\t", stringsAsFactors = F, header = F,fill = T,as.is = T)
#genes_of_interest_hgnc_t = read.table("~/SeneSys/Results/LM22_basis.tsv",sep ="\t", stringsAsFactors = F, header = T,fill = TRUE)
#genes_of_interest_hgnc_t = read.table("~/SeneSys/Results/Data_9461.Cell_fraction_predictions.bseq-sc.tsv",sep ="\t", stringsAsFactors = F, header = T,fill = TRUE)

genes_of_interest_hgnc_t[,1]

for ( i in 1:nrow(genes_of_interest_hgnc_t)){
    
    print(i)
    stem_path = "~/Downloads/Plots/"
    expr_raw = expr_raw[,colnames(expr_raw) != "GSM2601431"]
=======
genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.gmt",sep ="\t", stringsAsFactors = F, header = F)
#genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/Senescence geneset.tsv",sep ="\t", stringsAsFactors = F, header = T,fill = TRUE)

genes_of_interest_hgnc_t[,1]
#i = 5
for ( i in 1:nrow(genes_of_interest_hgnc_t)){

    print(i)
    stem_path = "~/Downloads/Plots/DeSEQ2/"
    #expr_raw = expr_raw[,colnames(expr_raw) != "GSM2601431"]
>>>>>>> b7e982528d8b656e7092e25f591c2a7eeb066cde
    meta_data = meta_info[colnames(expr_raw),]
    sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
    #sad_genes = rownames(genes_of_interest_hgnc_t)[genes_of_interest_hgnc_t[,3]>1000]
    
    #gene_matrix = expr_raw[,i]
    #percentiles = quantile(gene_matrix,probs = seq(0,1,by=.1))[10]
    #sad_genes  = rownames(expr_raw)[gene_matrix>percentiles]
    #sad_genes  = rownames(genes_of_interest_hgnc_t)[genes_of_interest_hgnc_t[,3]>1000]
    sad_genes   = sad_genes[ sad_genes != ""]
    #genes_of_interest_hgnc_t[i,1]
    #pathway_name = colnames(genes_of_interest_hgnc_t)[i]
    pathway_name = genes_of_interest_hgnc_t[i,1]
    pathway_name
    pca_name = paste0(c(stem_path,"/",pathway_name,".pca.pdf"),collapse ="")

    table(str_to_upper(sad_genes) %in% str_to_upper(rownames(expr_raw)))
    sad_genes[! (str_replace_all(sad_genes,pattern = "_","") %in% str_replace_all(rownames(expr_raw),pattern = "_",""))]
    
    expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
    variance_vec = apply(expr,FUN=var, MARGIN=1)
    expr = expr[variance_vec > 1,]
    
    correlation_matrix = cor(expr)
<<<<<<< HEAD
    
    try_return = try({pca = prcomp(t(correlation_matrix))})
    if(class(try_return) == "try-error")
      next
    #meta_data$RES_RP_NR = meta_data[,"ABC_GCB"]
    
    pdf(pca_name)
=======
    class_message =  try({pca = prcomp(t(correlation_matrix))})
    
    if (typeof(class_message)== "character")
      next()
    #meta_data$RES_RP_NR = meta_data[,"ABC_GCB"]
    
    pdf(pca_name)
    
    pheatmap::pheatmap(
      correlation_matrix,
      annotation_col = meta_data[c("Drug_Treatment")],
      annotation_colors = aka3,
      show_rownames = T,
      show_colnames = F,
      treeheight_col = 0,
      legend = F,
      cluster_rows = T,
      cluster_cols = T,
      fontsize_col = 7,
      clustering_method = "average"
    )
    
>>>>>>> b7e982528d8b656e7092e25f591c2a7eeb066cde
    p = ggbiplot::ggbiplot(
        pca,
        choices = c(1,2),
        obs.scale = 1,
        groups = meta_data[,"Drug_Treatment"],
        ellipse = TRUE,
        circle = TRUE,
        labels = rownames(meta_data),
        var.axes = F,
        labels.size = 1.5
    )
<<<<<<< HEAD
    
=======
    p
>>>>>>> b7e982528d8b656e7092e25f591c2a7eeb066cde
    print(p)
    dev.off()
    genes_of_interest_hgnc_t[i,1]
    ## Figure 1
}
    meta_data$Cluster_SUV = meta_data$Cluster
    #meta_data$Cluster_SUV[meta_data$Cluster_SUV %in% c(0,1,4)] = "SUV_high"
    #meta_data$Cluster_SUV[meta_data$Cluster_SUV %in% c(2,3,5)] = "SUV_low"
    order_vec = order(as.double(expr["SUVARNESS",]),decreasing = T)
    expr = expr[,order_vec]
    meta_data = meta_info[colnames(expr),]
    meta_data$Drup_response = meta_data$ABC_GCB
    meta_data$SUVARNESS = as.double(expr["SUVARNESS",])
    
    pheatmap::pheatmap(
        #pca,
        expr,
        annotation_col = meta_data[,c("SUVARNESS","Drup_response")],
        #annotation_col = meta_data[c("SUVARNESS")],
        annotation_colors = aka3,
        show_rownames = T,
        show_colnames = F,
        treeheight_col = 0,
        legend = F,
        cluster_rows = T,
        cluster_cols = F,
        fontsize_col = 7,
        clustering_method = "average"
    )
## Figure 1

table(meta_data$Cluster,meta_data$ABC_GCB)
pheatmap::pheatmap(table(meta_data$Cluster,meta_data$ABC_GCB))
    
