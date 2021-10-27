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
#path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.HGNC.tsv"
#path_transcriptome_file = "~/SeneSys/Data/Schmitz.HGNC.DESeq2.tsv"
#path_transcriptome_file = "~/SeneSys/Data/Schmitz.HGNC.tsv"
#path_transcriptome_file = "~/SeneSys/Data/GSE98588.DESeq2.tsv"
path_transcriptome_file = "~/SeneSys/Data/GSE98588_new.HGNC.tsv"


expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

table(colnames(expr_raw) %in% meta_info$Sample)
meta_data = meta_info[colnames(expr_raw),]

#expr_raw = expr_raw[,meta_data$Drug_Treatment %in% c("NR","RP")]
#meta_data = meta_info[colnames(expr_raw),]

#fe_es = read.table("~/SeneSys/Results/Schmitz.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
#fe_es = read.table("~/SeneSys/Results/Data_9461.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
fe_es = read.table("~/SeneSys/Results/GSE98588.new.HGNC.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
rownames(fe_es) = make.names(rownames(fe_es))

selection = c("ML_core","TIS.down","Left_Right_Hema_combined")
for (selector in selection){
  
  vector_ori = as.double(fe_es[selector,colnames(expr_raw)])
  
  quantiles = as.double(quantile(vector_ori, seq(0,1,0.01)))
  threshold = quantiles[67]
  #thresh_low = quantiles[34]
  #thresh_high = quantiles[67]
  vector = rep("high",length(vector_ori) )
  #vector[vector_ori <= thresh_high] = "medium"
  #vector[vector_ori <= thresh_low] = "low"
  vector[vector_ori <= threshold] = "low"
  meta_data[,selector] = vector
}

source("~/SeneSys/Scripts/Visualization_colors.R")
#

genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.gmt.tsv",sep ="\t", stringsAsFactors = F, header = F)

genes_of_interest_hgnc_t$V1

i = 6
genes_of_interest_hgnc_t$V1[i]
sad_genes = genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)] 
sad_genes = c(genes_of_interest_hgnc_t[6,3:ncol(genes_of_interest_hgnc_t)] ,genes_of_interest_hgnc_t[7,3:ncol(genes_of_interest_hgnc_t)] )
#sad_genes = sad_genes[1:50]
sad_genes = sad_genes[sad_genes != ""]

table(sad_genes %in% rownames(expr_raw))
expr = expr_raw[match(sad_genes,  rownames(expr_raw), nomatch = 0),]

colnames(expr) = meta_data$Name
rownames(meta_data) = meta_data$Name
expr[1:5,1:5]
correlation_matrix = cor(expr);pcr = prcomp(t(correlation_matrix))

p = pheatmap::pheatmap(
  #expr,
  correlation_matrix,
  #annotation_col = meta_data[,c("Drug_Treatment",selection)],
  annotation_col = meta_data[c("Left_Right_ML_core_50","ABC_GCB","Cluster",selection)],
  #annotation_col = meta_data[,c("ABC_GCB","Progression",selection)],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = FALSE,
  #treeheight_col = 0,
  treeheight_row = 0,
  legend = F,
  cluster_rows = T,
  cluster_cols = T,
  fontsize_col = 7,
  clustering_method = "ward.D"
)
p

index = p$tree_col$labels[ p$tree_col$order]
subtype_vector = rep("right", ncol(expr))
subtype_vector[1:(which(index == "GSM2601448"))] = "left" #GSM2601412
meta_info$Left_Right_ML_core_50 =rep("",nrow(meta_info))
matcher = match(index, meta_info$Sample)
meta_info[matcher,"Left_Right_ML_core_50"] = subtype_vector

## Figure 1

###

p = ggbiplot::ggbiplot(
  pcr,
  groups = meta_data$Drug_Treatment,
  var.axes = F,
  ellipse = TRUE,
  labels = meta_data$Name
)
p = p + geom_point(aes(colour=meta_data$Drug_Treatment), size = 4)
p = p + scale_color_manual(values = c("darkgreen","orange"))
p

# umap

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
#custom.config$random_state = 995
custom.config$n_components=2

genes_of_interest_hgnc_t$V1
i = 62
sad_genes = genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)] 
sad_genes = sad_genes[sad_genes != ""]
sad_genes = sad_genes[1:50]
expr = expr_raw[match(sad_genes,  rownames(expr_raw), nomatch = 0),]
expr = expr[,meta_data$ABC_GCB != "Unclassified"]
meta_data_plot = meta_info[colnames(expr),]
#expr[1:5,1:5]
correlation_matrix = cor(expr);pcr = prcomp(t(correlation_matrix))

umap_result = umap::umap(
  correlation_matrix,
  colvec = meta_data_plot$Drug_Treatment,
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point( aes( size = 4, color = as.character(meta_data_plot$Drug_Treatment) ))
#umap_p = umap_p+geom_text(size= 4,aes(label=rownames(meta_data),color = as.character(meta_data$Cluster)),hjust=0, vjust=0)
umap_p = umap_p +  xlab("") + ylab("")  
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data_plot$Drug_Treatment), level=.5, type ="t", size=1.5)
#umap_p = umap_p + scale_color_manual( values = c("darkgreen","orange"),name="Drug_treatment") ##33ACFF ##FF4C33
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

data_t = as.data.frame(cbind(meta_data$ABC_GCB,meta_data$TIS.down))
table(data_t)

data_t = as.data.frame(cbind(meta_data$Left_Right_ML_core_50,meta_data$ABC_GCB))
table(data_t)

meta_data$Cluster_2 = meta_data$Cluster
meta_data$Cluster_2[meta_data$Cluster_2 %in% c("3","4")] = "34"
meta_data$Cluster_2[meta_data$Cluster_2 %in% c("0","1","2","5")] = "0125"

####

meta_data = meta_info[colnames(fe_es),]
fe_es = fe_es[,colnames(expr_raw)]

NR_vec =c()
NP_vec = c()

for (i in 1:nrow(fe_es)){
  values = as.double(as.character(fe_es[i,]))
  cohort = list(meta_data$Drug_Treatment)
  result = aggregate(
    values,
    FUN = mean,
    by = cohort
  )
  NR_vec = c(NR_vec, result[1,2])
  NP_vec = c(NP_vec, result[2,2])
}

NR_vec
NP_vec

which.max(abs(NR_vec - NP_vec))

rownames(fe_es)[order(abs(NR_vec - NP_vec), decreasing = TRUE)[1:3]]

###
selection = c("TANG_SENESCENCE_TP53_TARGETS_DN","ML_core","TIS.down","Left_Right_Hema_combined")
for (selector in selection){
  
  vector_ori = as.double(fe_es[selector,colnames(expr_raw)])
  
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


data = as.data.frame(cbind(meta_data$Cluster,meta_data$ML_core))
table(data)
chisq.test(table(data))
hist(as.double(fe_es["ML_core",]))

###
data = as.data.frame(cbind(meta_data$Progression,meta_data$ML_core))
table(data)
