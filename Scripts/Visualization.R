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
rownames(meta_info) = meta_info$Sample_ID
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

#source("~/Deko/Scripts/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t[,1]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[38,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]

path_transcriptome_file = "~/SeneSys/Data/GSE98588_series_matrix.prepared.tsv"
#visualization_data_path = str_replace(path_transcriptome_file,pattern  ="\\.tsv",".vis.tsv")

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

meta_data = meta_info[colnames(expr_raw),]
meta_data = meta_data[ meta_data$Include,]

expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
expr = expr[,colnames(expr) %in% rownames(meta_data)]
expr[1:5,1:5]
dim(expr)

correlation_matrix = cor(expr)
pcr = prcomp(t(correlation_matrix))

p = ggbiplot::ggbiplot(
    pcr,
    obs.scale = 1,
    #    groups = meta_data$Grading,
    ellipse = TRUE,
    circle = TRUE,
    #labels = rownames(meta_data)
    var.axes = F#,
)
p
## Figure 1

#meta_data$Location[!str_detect(meta_data$Location,pattern = "Primary")] = "Metastasis"
#meta_data$Grading[meta_data$Grading == ""] = "G0"
pheatmap::pheatmap(
    correlation_matrix,
    #annotation_col = meta_data[c("NEC_NET","Grading")],
    #annotation_colors = aka3,
    show_rownames = F,
    show_colnames = F,
    treeheight_col = 0,
    legend = F,
    fontsize_col = 7,
    clustering_method = "average"
)

