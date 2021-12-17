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
meta_info = meta_info[meta_info$Sample!="",]
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

#path_transcriptome_file = "~/Dropbox/testproject/SABGal_exp.csv"
path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.DeSeq2.HGNC.tsv"
#path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.HGNC.tsv"
#path_transcriptome_file = "~/SeneSys/Data/Schmitz.HGNC.DESeq2.tsv"
#path_transcriptome_file = "~/SeneSys/Data/Schmitz.HGNC.tsv"
path_transcriptome_file = "~/SeneSys/Data/GSE98588.DESeq2.tsv"
#path_transcriptome_file = "~/SeneSys/Data/GSE98588_new.HGNC.tsv"

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

table(colnames(expr_raw) %in% meta_info$Sample)
meta_data = meta_info[colnames(expr_raw),]

selection = c("ADR","ADROHT","PS")
for (selector in selection){
  
  vector_ori = as.double(meta_data[colnames(expr_raw),selector])
  
  quantiles = as.double(quantile(vector_ori, seq(0,1,0.01)))
  threshold = quantiles[67]
  vector = rep("high",length(vector_ori) )
  vector[vector_ori <= threshold] = "low"
  meta_data[,paste(selector,"binary",sep ="_")] = vector
}

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
  vector = rep("High",length(vector_ori) )
  vector[vector_ori <= threshold] = "How"
  meta_data[,selector] = vector
}
source("~/SeneSys/Scripts/Visualization_colors.R")
#

genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.gmt.tsv",sep ="\t", stringsAsFactors = F, header = F)

genes_of_interest_hgnc_t$V1

i = 60
genes_of_interest_hgnc_t$V1[i]
sad_genes = genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)] 
#sad_genes = c(genes_of_interest_hgnc_t[6,3:ncol(genes_of_interest_hgnc_t)] ,genes_of_interest_hgnc_t[7,3:ncol(genes_of_interest_hgnc_t)] )
#sad_genes = sad_genes[1:50]
sad_genes = sad_genes[sad_genes != ""]

table(sad_genes %in% rownames(expr_raw))
expr = expr_raw[match(sad_genes,  rownames(expr_raw), nomatch = 0),]
expr = props
meta_data = meta_info[colnames(expr),]
meta_data$bgal_binary[meta_data$bgal_binary == ""] = "Unknown"

colnames(expr) = meta_data$Sample
#rownames(meta_data) = meta_data$Name
expr[1:5,1:5]
correlation_matrix = cor(t(expr));
pcr = prcomp(t(correlation_matrix))

source("~/SeneSys/Scripts/Visualization_colors.R")
p = pheatmap::pheatmap(
  #expr,
  correlation_matrix,
  #fe_es,
  #annotation_col = meta_data[,c("Drug_Treatment","ADROHT_binary","ADR_binary","bgal_binary")],
  annotation_col = meta_data[c("ABC_GCB","Cluster","ADR_binary","PS_binary")],
  #annotation_col = meta_data[,c("ABC_GCB","Progression","ADR_binary","ADROHT_binary","PS_binary")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = FALSE,
  #treeheight_col = 0,
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
subtype_vector[1:(which(index == "GSM2601448"))] = "left" #GSM2601412
meta_info$Left_Right_ML_core_50 =rep("",nrow(meta_info))
matcher = match(index, meta_info$Sample)
meta_info[matcher,"Left_Right_ML_core_50"] = subtype_vector

## Figure 1

###

p = ggbiplot::ggbiplot(
  pcr,
  groups = meta_data$bgal_binary,
  var.axes = F,
  ellipse = TRUE,
  labels = meta_data$Name
)
p = p + geom_point(aes(colour=meta_data$bgal_binary), size = 4)
p = p + scale_color_manual(values = c("darkgreen","orange"))
p

# umap

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
#custom.config$random_state = 995
custom.config$n_components=2

genes_of_interest_hgnc_t$V1
i = 60
sad_genes = genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)] 
sad_genes = sad_genes[sad_genes != ""]
sad_genes = sad_genes[]
expr = expr_raw[match(sad_genes,  rownames(expr_raw), nomatch = 0),]
expr = expr[,meta_data$ABC_GCB != "Unclassified"]
meta_data_plot = meta_info[colnames(expr),]
#expr[1:5,1:5]
correlation_matrix = cor(expr);pcr = prcomp(t(correlation_matrix))

umap_result = umap::umap(
  correlation_matrix,
  colvec = meta_data$Drug_Treatment,
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
#umap_p = umap_p + geom_point( size =4, aes(color=meta_data$adr_binary) )#+geom_text(aes(label=meta_data$Sample),hjust=0, vjust=0)
#umap_p = umap_p + geom_point( size =4, aes(color=meta_data$bgal_binary) )#+geom_text(aes(label=meta_data$Sample),hjust=0, vjust=0)
umap_p = umap_p + geom_point( size =4, aes(color=meta_data$Drug_Treatment) )#+geom_text(aes(label=meta_data$Sample),hjust=0, vjust=0)
umap_p = umap_p +  xlab("") + ylab("") 
#umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$bgal_binary), level=.5, type ="t", size=1.5)
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$Drug_Treatment), level=.5, type ="t", size=1.5)
#umap_p = umap_p + scale_color_manual( values = c("darkgreen","orange"),name="Drug_treatment") ##33ACFF ##FF4C33
umap_p
#genes_of_interest_hgnc_t$V1[i]
#custom.config$random_state

meta_ic50 = meta_data[! is.na(meta_data$ADR_IC50),]
cor.test(meta_ic50$ADR_IC50,meta_ic50$ADR)
cor.test(meta_ic50$ADR_IC50,meta_ic50$ADROHT)
cor.test(meta_ic50$ADR_IC50,meta_ic50$ADR)
cor.test(meta_ic50$ADR_IC50,meta_ic50$ADROHT)

meta_Sa_Bgal = meta_data[! is.na(meta_data$Sa_Bgal),]
cor.test(meta_Sa_Bgal$Sa_Bgal,meta_Sa_Bgal$ADR)
cor.test(meta_Sa_Bgal$Sa_Bgal,meta_Sa_Bgal$ADROHT)

chisq.test(meta_Sa_Bgal$bgal_binary,meta_Sa_Bgal$adr_binary)
chisq.test(meta_Sa_Bgal$bgal_binary,meta_Sa_Bgal$ADROHT_Binary)

chisq.test(meta_Sa_Bgal$Drug_Treatment,meta_Sa_Bgal$adr_binary)
table(meta_Sa_Bgal$Drug_Treatment,meta_Sa_Bgal$adr_binary)

chisq.test(meta_Sa_Bgal$Drug_Treatment,meta_Sa_Bgal$ADROHT_Binary)

chisq.test(meta_data$Drug_Treatment,meta_data$ADROHT_Binary)
chisq.test(meta_data$Drug_Treatment,meta_data$ADR_binary)
table(meta_data$Drug_Treatment,meta_data$ADR_binary)
table(meta_data$Drug_Treatment,meta_data$ADROHT_Binary)

chisq.test(meta_data$Progression,meta_data$ADROHT_binary)

table(meta_data$Cluster,meta_data$PS_binary)
subvector = meta_data$Cluster
subvector[subvector %in% c(1,2,3)] = "123"
subvector[subvector %in% c(0,4,5)] = "045"
table(meta_data$PS_binary,meta_data$ABC_GCB)
