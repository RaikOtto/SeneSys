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
meta_info = meta_info[meta_info$Sample!="",]
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.DeSeq2.HGNC.tsv"
#path_transcriptome_file = "~/SeneSys/Data/GSE98588.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
#expr_raw = expr_raw[,!(colnames(expr_raw) %in% "GSM2601431")]
meta_data = meta_info[colnames(expr_raw),]
#meta_data$Cluster = as.factor(meta_data$Cluster)
expr_raw[1:5,1:5]

#genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.gmt",sep ="\t", stringsAsFactors = F, header = F)
#genes_of_interest_hgnc_t = read.table("~/SeneSys/Results/LM22_basis.tsv",sep ="\t", stringsAsFactors = F, header = T,fill = TRUE)
genes_of_interest_hgnc_t = read.table("~/SeneSys/Results/Data_9461.Cell_fraction_predictions.bseq-sc.tsv",sep ="\t", stringsAsFactors = F, header = T,fill = TRUE)

#genes_of_interest_hgnc_t[,1]

#for ( i in 1:nrow(genes_of_interest_hgnc_t)){
    
    stem_path = "~/Downloads/Plots/"
    expr_raw = expr_raw[,colnames(expr_raw) != "GSM2601431"]
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
    pca_name = paste0(c(stem_path,"/",pathway_name,".pca.pdf"),collapse ="")
    
    table(str_to_upper(sad_genes) %in% str_to_upper(rownames(expr_raw)))
    sad_genes[! (str_replace_all(sad_genes,pattern = "_","") %in% str_replace_all(rownames(expr_raw),pattern = "_",""))]
    sad_genes
    
    expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
    expr = genes_of_interest_hgnc_t[,1:6]
    correlation_matrix = cor(expr)
    pca = prcomp(t(correlation_matrix))
    meta_data$RES_RP_NR = meta_data[,"ABC_GCB"]
    
    #pdf(pca_name)
    p = ggbiplot::ggbiplot(
        expr,
        choices = c(1,2),
        obs.scale = 1,
        groups = meta_data[,"RES_RP_NR"],
        ellipse = TRUE,
        circle = TRUE,
        #labels = rownames(meta_data),
        var.axes = F#,
    )
    print(p)
    dev.off()
    genes_of_interest_hgnc_t[i,1]
    ## Figure 1
#}    
    pheatmap::pheatmap(
        t(expr),
        #correlation_matrix,
        annotation_col = meta_data[c("ABC_GCB")],
        #annotation_colors = aka3,
        show_rownames = F,
        show_colnames = T,
        treeheight_col = 0,
        legend = F,
        fontsize_col = 7,
        clustering_method = "average"
    )
## Figure 1

library("umap")
umap_plot = umap::umap((correlation_matrix))
vis_data = as.data.frame(umap_plot$layout)
colnames(vis_data) = c("x","y")
dist_mat = dist((vis_data))
p = ggplot2::qplot(
    x = vis_data$x,
    y = vis_data$y,
    #color = meta_data$Complex,
    color = meta_data$ABC_GCB,
    geom=c("point"),
    xlab = "Umap dim 1",
    ylab = "Umap dim 2"#,
    #shape = meta_data[colnames(expr),"IFN_I"]
)
p

aggregate( t(expr), by = list(meta_data$ABC_GCB), FUN = mean)

###

data_t = read.table("~/Deco/Results/Cell_fraction_predictions/Baron_Bseqsc_All_Datasets.tsv",sep="\t", stringsAsFactors =  F, header = T, as.is = F)

#vis_mat = reshape2::melt(meta_data[,c("ABC_GCB","NK","CD8","Cluster","Complex")])
vis_mat = reshape2::melt(meta_data[,c("Cluster","NK")])
colnames(vis_mat) = c("Subtype","Cell_type","Proportion")
vis_mat_agg = aggregate(vis_mat[,"Proportion"],by = list(vis_mat[,"Subtype"]),FUN = mean)
colnames(vis_mat_agg) = c("Subtype","Proportion")
#svg(filename = "~/Deco/Results/Images/Figure_3_Proportion_MKI67_versus_Grading.svg", width = 10, height = 10)

###

p = ggplot(
    data = vis_mat_agg,
    aes(
        x = Subtype,
        y = Proportion,
        fill = Subtype
    )
)
p = p + geom_bar(stat="identity", position=position_dodge(), color = "black")
p = p + scale_fill_manual(values = c("darkgreen", "black"))
p = p + ylab(label = "P-value nu-SVR regression models") + theme(legend.position="top") + xlab(label = "Grading")
p = p + geom_errorbar(aes(),  position = "dodge")
p = p + guides(fill=guide_legend(title="Deconvolution model")) 
p = p + scale_fill_manual(labels = c("endocrine only", "endocrine & exocrine"), values = c("darkgreen", "black"))
p

svg("~/Deco/Results/Images/Figure_2a.svg")
p
dev.off()
