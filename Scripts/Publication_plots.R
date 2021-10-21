library("umap")
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

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = meta_info %>% dplyr::filter(Study != "GSE11318")
meta_info = meta_info[meta_info$Sample!="",]
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

#path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.DeSeq2.HGNC.tsv"
path_transcriptome_file = "~/SeneSys/Data/Schmitz.HGNC.DESeq2.tsv"
#path_transcriptome_file = "~/SeneSys/Data/GSE98588.DESeq2.tsv"

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

table(colnames(expr_raw) %in% meta_info$Sample)
meta_data = meta_info[colnames(expr_raw),]
#meta_data$Cluster = as.factor(meta_data$Cluster)

#fe_es = read.table("~/SeneSys/Results/Schmitz.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
fe_es = read.table("~/SeneSys/Results/Data_9461.DESeq2.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
#fe_es = read.table("~/SeneSys/Results/GSE98588.DESeq2.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)

selection = c("Left_Right","CARDness","SUVARNESS","TIS up")
for (selector in selection){
  
  vector_ori = as.double(fe_es[selector,colnames(expr)])
  
  quantiles = as.double(quantile(vector_ori, seq(0,1,0.01)))
  threshold = quantiles[51]
  #thresh_low = quantiles[34]
  #thresh_high = quantiles[67]
  vector = rep("high",length(vector_ori) )
  #vector[vector_ori <= thresh_high] = "medium"
  #vector[vector_ori <= thresh_low] = "low"
  vector[vector_ori <= threshold] = "low"
  meta_data[,selector] = vector
}

source("~/SeneSys/Scripts/Visualization_colors.R")
meta_data$Drup_response = meta_data$ABC_GCB

#

genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.gmt.tsv",sep ="\t", stringsAsFactors = F, header = F)

genes_of_interest_hgnc_t$V1

i = 59
genes_of_interest_hgnc_t$V1[i]
sad_genes = genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)] 
sad_genes = sad_genes[sad_genes != ""]

table(sad_genes %in% rownames(expr_raw))
expr = expr_raw[match(sad_genes,  rownames(expr_raw), nomatch = 0),]
#expr[1:5,1:5]
correlation_matrix = cor(expr);pcr = prcomp(t(correlation_matrix))

p = pheatmap::pheatmap(
  #expr,
  correlation_matrix,
  annotation_col = meta_data[,c("Drup_response","Predicted",selection)],
  #annotation_col = meta_data[c("ABC_GCB","Cluster",selection)],
  #annotation_col = meta_data[c("ABC_GCB","Predictions","Progression",selection)],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = F,
  treeheight_col = 0,
  treeheight_row = 0,
  legend = F,
  cluster_rows = T,
  cluster_cols = T,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)
p

index = p$tree_col$labels[ p$tree_col$order]
subtype_vector = rep("right", ncol(expr))
subtype_vector[1:(which(index == "AS_222823_LR_34436"))] = "left"
meta_info$Left_Right =rep("",nrow(meta_info))
matcher = match(index, meta_info$Sample)
meta_info[matcher,"Left_Right"] = subtype_vector

## Figure 1

table(meta_data$Left_Right,meta_data$Progression)
pheatmap::pheatmap(table(meta_data$Cluster,meta_data$ABC_GCB))

###

p = ggbiplot::ggbiplot(
  pcr,
  #labels = meta_data$Drup_response,
  groups = meta_data$Drup_response,
  var.axes = F,
  ellipse = TRUE,
  obs.scale = 1,
  var.scale = 1,
  labels.size = 2
)
p = p + scale_color_manual(values = c("darkgreen","orange","black"))
#p = p + xlim(-1.9,1.3) + ylim(-.85,1.2) # Fig 5
p

# umap

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
#custom.config$random_state = 995
custom.config$n_components=2

genes_of_interest_hgnc_t$V1
i = 59
sad_genes = genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)] 
sad_genes = sad_genes[sad_genes != ""]
expr = expr_raw[match(sad_genes,  rownames(expr_raw), nomatch = 0),]
#expr[1:5,1:5]
correlation_matrix = cor(expr);pcr = prcomp(t(correlation_matrix))

umap_result = umap::umap(
  correlation_matrix,
  colvec = meta_data$Drup_response,
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data$Predicted) ))
#umap_p = umap_p+geom_text(size= 2,aes(label=rownames(meta_data),color = as.character(meta_data$Progression)),hjust=0, vjust=0)
#umap_p = umap_p+geom_text(size= 2,aes(label=rownames(meta_data),color = as.character(meta_data$Predictions)),hjust=0, vjust=0)
#umap_p = umap_p + theme(legend.position = "none") + xlab("") + ylab("")
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$Drup_response), level=.5, type ="t", size=1.5)
#umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$Predictions), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("darkgreen","green","yellow","darkred","black")) ##33ACFF ##FF4C33
umap_p
genes_of_interest_hgnc_t$V1[i]
#custom.config$random_state


################






for ( i in 1:nrow(genes_of_interest_hgnc_t)){
  
  print(i)
  stem_path = "~/Downloads/Plots/"
  #expr_raw = expr_raw[,colnames(expr_raw) != "GSM2601431"]
  
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
  p
  
  print(p)
  dev.off()
  genes_of_interest_hgnc_t[i,1]
  ## Figure 1
}

#######################

#write.table(meta_info,"~/SeneSys/Misc/Meta_information.tsv",sep ="\t",quote =F,row.names = F)

cluster_left_right = as.data.frame(cbind(meta_data$Cluster,meta_data$Left_Right))
cluster_left_right$V2 = factor(cluster_left_right$V2, levels = c("high","medium","low"))
table(cluster_left_right)

meta_data$Cluster_2 = meta_data$Cluster
meta_data$Cluster_2[meta_data$Cluster_2 %in% c("3","4")] = "34"
meta_data$Cluster_2[meta_data$Cluster_2 %in% c("0","1","2","5")] = "0125"
