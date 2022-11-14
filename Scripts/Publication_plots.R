library("umap")
library("stringr")
library("ggplot2")
library("dplyr")
library("grid")

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
#meta_info = read.table("~/SeneSys/Misc/Mouse_cell_types/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = meta_info[meta_info$Sample!="",]
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

#path_transcriptome_file = "~/Dropbox/testproject/SABGal_exp.csv"
#path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.DeSeq2.HGNC.tsv"
#path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.HGNC.tsv"
#path_transcriptome_file = "~/SeneSys/Data/Schmitz.HGNC.DESeq2.tsv"
#path_transcriptome_file = "~/SeneSys/Data/Schmitz.clean.tsv"
#path_transcriptome_file = "~/SeneSys/Data/GSE98588.DESeq2.tsv"
#path_transcriptome_file = "~/SeneSys/Data/GSE98588_new.HGNC.tsv"
path_transcriptome_file = "~/SeneSys/Data/Reddy_hgnc_S255.tsv"

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors = F, header = T, as.is = TRUE,row.names = 1, dec= ".")
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

table(colnames(expr_raw) %in% meta_info$Sample)
meta_data = meta_info[colnames(expr_raw),]
rownames(meta_data) = meta_data$Sample

selection_vec = which(meta_data$Doro_annot != "")
expr_raw = expr_raw[,selection_vec]
dim(expr_raw)

#selection = c("ADR","ADROHT","PS")

#for (selector in selection){
  
#  vector_ori = as.double(meta_data[colnames(expr_raw),selector])
  
#  quantiles = as.double(quantile(vector_ori, seq(0,1,0.01)))
#  threshold = quantiles[51]
#  vector = rep("high",length(vector_ori) )
#  vector[vector_ori <= threshold] = "low"
#  meta_data[,paste(selector,"binary",sep ="_")] = vector
#}

#expr_raw = expr_raw[,meta_data$Drug_Treatment %in% c("NR","RP")]
#meta_data = meta_info[colnames(expr_raw),]

#fe_es = read.table("~/SeneSys/Results/Schmitz.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
#fe_es = read.table("~/SeneSys/Results/Data_9461.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
#fe_es = read.table("~/SeneSys/Results/GSE98588.new.HGNC.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
#fe_es = read.table("~/SeneSys/Results/Reddy.GSVA.tsv",sep ="\t", header = T,as.is = T,row.names = 1)
#rownames(fe_es) = make.names(rownames(fe_es))
#colnames(fe_es ) = str_replace_all(colnames(fe_es ), "X","")

#selection = c("ML_core","TIS.up","TIS.down")
#for (selector in selection){
  
#  vector_ori = as.double(fe_es[selector,colnames(expr_raw)])
  
#  quantiles = as.double(quantile(vector_ori, seq(0,1,0.01)))
#  threshold = quantiles[67]
#  vector = rep("High",length(vector_ori) )
#  vector[vector_ori <= threshold] = "Medium_Low"
#  meta_data[,selector] = vector
#}
source("~/SeneSys/Scripts/Visualization_colors.R")
#
#selection = c("S01","S02","S03","S04","S05")
#indeces = as.integer(apply(meta_data[,selection], FUN = function(vec){return(which.max(vec))}, MARGIN = 1))
#meta_data$Cell_state = selection[indeces]
#selection = c("L1","L2","L3","L4","L5","L6","L7","L8","L9")
#indeces = as.integer(apply(meta_data[,selection], FUN = function(vec){return(which.max(vec))}, MARGIN = 1))
#meta_data$Ecotype = selection[indeces]

#genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.gmt.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SAS_TF.gmt.tsv",sep ="\t", stringsAsFactors = F, header = F)

genes_of_interest_hgnc_t$V1

i = 38
genes_of_interest_hgnc_t$V1[i]
sad_genes = genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)] 
#sad_genes = c(genes_of_interest_hgnc_t[6,3:ncol(genes_of_interest_hgnc_t)] ,genes_of_interest_hgnc_t[7,3:ncol(genes_of_interest_hgnc_t)] )
#sad_genes = sad_genes[1:50]
sad_genes = sad_genes[sad_genes != ""]

table(sad_genes %in% rownames(expr_raw))
expr = expr_raw[match(sad_genes,  rownames(expr_raw), nomatch = 0),]

rownames(meta_data) = meta_data$Sample
expr[1:5,1:5]
expr = expr
#expr = expr[,meta_data[colnames(expr),"ABC_GCB"] != "Unclassified"]
source("~/SeneSys/Scripts/Visualization_colors.R")

p = pheatmap::pheatmap(
  t(expr),
  #correlation_matrix,
  #fe_es,
  #annotation_col = meta_data[,c("Doro_annot","Response_to.initial.therapy","CNS_Relapse","CNS_Involvement")],
  #annotation_col = meta_data[,c("ABC_GCB","TIS.up")],
  #annotation_col = meta_data[,c("ABC_GCB","Left_Right","TIS.down")],
  annotation_colors = aka3,
  show_rownames = TRUE,
  show_colnames = FALSE,
  #treeheight_col = 0,
  treeheight_row = 0,
  legend = F,
  cluster_rows = T,
  cluster_cols = T,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)

correlation_matrix = cor(t(expr))
correlation_matrix = expr
pcr = prcomp(t(correlation_matrix))

source("~/SeneSys/Scripts/Visualization_colors.R")
p = pheatmap::pheatmap(
  #scale(expr),
  cor(t(expr)),
  #fe_es,
  #annotation_col = meta_data["Drug_Treatment"],
  #annotation_col = meta_data[,c("ABC_GCB","Cell_state","TIS.down","TIS.up")],
  annotation_col = meta_data["Doro_annot"],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = FALSE,
  #treeheight_col = 0,
  treeheight_row = 0,
  legend = F,
  cluster_rows = T,
  cluster_cols = T,
  fontsize_col = 3,
  clustering_method = "ward.D2"
)

#svg(filename = "~/Downloads/Chapuy_3.svg", width = 10, height = 10)
p
#dev.off()

index = p$tree_col$labels[ p$tree_col$order]
subtype_vector = rep("TIS down", ncol(expr))
subtype_vector[1:(which(index == "GSM2601460"))] = "TIS mixed" #GSM2601460 #s236  #AS_22222824_LR_34436
subtype_vector[1:(which(index == "GSM2601468"))] = "TIS up" #GSM2601468 #s529
matcher = match(index, meta_data$Sample)
meta_data[matcher,"TIS_cluster"] = subtype_vector

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
meta_data_plot = meta_data[colnames(expr),]

expr = meta_data[,c("L1","L2","L3","L4","L5","L6","L7","L8","L9")]
expr[1:5,1:5]
expr = expr[,colnames(expr) != "GSM2601431"]
expr = expr[rownames(expr) != "GSM2601431",]
expr = expr[meta_data[rownames(expr),"ABC_GCB"] != "Unclassified",]
meta_data_plot = meta_data[rownames(expr),]

correlation_matrix = cor(t(expr))

umap_result = umap::umap(
  correlation_matrix,
  colvec = meta_data$Doro_annot,
  preserve.seed = TRUE#,
  #config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y)
)
#umap_p = umap_p + geom_point( size =4, aes(color=meta_data$adr_binary) )#+geom_text(aes(label=meta_data$Sample),hjust=0, vjust=0)
#umap_p = umap_p + geom_point( size =4, aes(color=meta_data$bgal_binary) )#+geom_text(aes(label=meta_data$Sample),hjust=0, vjust=0)
umap_p = umap_p + geom_point( size =3, aes(color = meta_data$Doro_annot) )#+geom_text(aes(label=meta_data$Name[meta_data$Drug_Treatment == "RES"]))
umap_p = umap_p +  xlab("") + ylab("")
#umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$bgal_binary), level=.5, type ="t", size=1.5)
#umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data_plot$Drug_Treatment), level=.5, type ="t", size=1.5)
#umap_p = umap_p + scale_color_manual( values = c("darkgreen","brown","black","gray")) ##33ACFF ##FF4C33
umap_p = umap_p +  theme(legend.position="top")
#umap_p = umap_p + scale_colour_manual(values = c("darkgreen","brown"))

#svg(filename = "~/Downloads/Chapuy_3b.svg", width = 10, height = 10)
umap_p  
dev.off()

#custom.config$random_state

meta_ic50 = meta_data[! is.na(meta_data$ADR_IC50),]
cor.test(meta_ic50$ADR_IC50,meta_ic50$ADR)
cor.test(meta_ic50$ADR_IC50,meta_ic50$ADROHT)
cor.test(meta_ic50$ADR_IC50,meta_ic50$ADR)
cor.test(meta_ic50$ADR_IC50,meta_ic50$ADROHT)

meta_Sa_Bgal = meta_data[! is.na(meta_data$Sa_Bgal),]
meta_Sa_Bgal = meta_Sa_Bgal[!is.na(meta_Sa_Bgal$TIS_cluster),]
aggregate(meta_Sa_Bgal$Sa_Bgal, by = list(meta_Sa_Bgal$TIS_cluster), FUN = mean)
down = meta_Sa_Bgal[meta_Sa_Bgal$TIS_cluster  == "TIS down","Sa_Bgal"]
up = meta_Sa_Bgal[meta_Sa_Bgal$TIS_cluster  == "TIS up","Sa_Bgal"]
t.test(down,up)

table(meta_data$Cluster,meta_data$PS_binary)
subvector = meta_data$Cluster
subvector[subvector %in% c(1,2,3)] = "123"
subvector[subvector %in% c(0,4,5)] = "045"
table(meta_data$PS_binary,meta_data$ABC_GCB)

props = props[meta_data$Sample,]
meta_data$p_value= p_value = props$P_value
meta_data$p_value[p_value < 0.05] = "Significant"
meta_data$p_value[meta_data$p_value != "Significant"] = "Not significant"
table(meta_data$PS == 0,meta_data$Progression)
vec = meta_data$Progression
vec[meta_data$Progression == "Progression"] = 1
vec[vec != 1] = 0

library(ROCR)
df <- data.frame(cbind(props$P_value,vec))
pred <- prediction(as.double(df$V1), as.integer(df$vec))
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)

####

expr_vis = expr
expr_vis$Drug_Treatment = meta_data$Drug_Treatment
vis_mat = reshape2::melt(expr_vis)
colnames(vis_mat) = c("Drug_Treatment","Characteristic","Value")
vis_mat_export = vis_mat %>% group_by(Drug_Treatment,Characteristic ) %>% dplyr::summarise("Mean_value" = mean(Value))

ggplot(vis_mat_export, aes(x = Characteristic, y = Mean_value, fill = Drug_Treatment)) + geom_bar(stat="identity", position=position_dodge())

###

#selection = c("S01","S02","S03","S04","S05")
#indeces = as.integer(apply(meta_data[,selection], FUN = function(vec){return(which.max(vec))}, MARGIN = 1))
#meta_data$Cell_state = selection[indeces]
#meta_info[meta_data$Sample,"Cell_state"]

table(meta_data$Cell_state,meta_data$TIS_cluster)
table(meta_data$Cell_state,meta_data$Drug_Treatment)
table(meta_data$Cell_state,meta_data$TIS.down)
table(meta_data$Cell_state,meta_data$TIS.up)
table(meta_data$Progression,meta_data$Cell_state)

dd = meta_data[meta_data$TIS_cluster %in% c("TIS up", "TIS down"),]
dd = meta_data[(meta_data$TIS.up %in% "High") | (meta_data$TIS.down %in% "High"),]
table(dd$Ecotype,dd$TIS_cluster) 


dd = meta_data
dd[dd$Cell_state %in% c("S01","S04","S05"),"Cell_state"] = "1_4_5"
dd[dd$Cell_state %in% c("S02","S03"),"Cell_state"] = "2_3"
dd$Cell_state = factor(dd$Cell_state, levels = c("1_4_5","2_3"))

dd = dd[!is.na(dd$TIS_cluster),]
vis_data = reshape2::melt(table(dd[,c("TIS_cluster","Cell_state")]))
colnames(vis_data)= c("TIS_cluster","Cell_state","Count")
vis_data[vis_data$TIS_cluster == "TIS up","Count"] = vis_data[vis_data$TIS_cluster == "TIS up","Count"] * -1
vis_data$Cell_state = factor(vis_data$Cell_state, levels = rev(unique(vis_data$Cell_state)))

p = ggplot( data = vis_data,aes(x = TIS_cluster, fill = Cell_state ,  y = Count))
p = p + geom_bar(stat="identity", position=position_dodge())

p = p + scale_fill_manual(values = c("red","blue"))
p = p + scale_fill_manual(values = (c("#9D33FF","#ee9f53","#b01511","blue","red")))
p = p + coord_flip()

#svg(filename = "~/Downloads/Chapuy_raw_cell_state_tis_cluster.svg", width = 10, height = 10)
p
dev.off()

ggbiplot::ggbiplot(
    pcr,
    groups = as.character(meta_data$Drug_Treatment[meta_data$Drug_Treatment == "RES"]),
    ellipse = TRUE,
    circle = TRUE,
    var.axes = F,
    labels = meta_data$Name[meta_data$Drug_Treatment == "RES"]
)  + scale_color_manual(name="Drug response", values=c("black"))
#dev.off()

###

table(meta_data[,c("Doro_annot","CNS_Involvement")])
table(meta_data[,c("Doro_annot","CNS_Relapse")])
table(meta_data[,c("Doro_annot","CNS_Relapse","CNS_Involvement")])
table(meta_data[,c("Doro_annot","Response_to.initial.therapy")])


meta_data[,c("Doro_annot","Response_to.initial.therapy","CNS_Relapse","CNS_Involvement")]