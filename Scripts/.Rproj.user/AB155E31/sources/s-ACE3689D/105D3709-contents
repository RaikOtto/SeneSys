library("SCDC")
library("stringr")
library("reshape2")
library("dplyr")
library("Biobase")

meta_info = read.table("~/SeneSys///Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")
meta_info$NEC_NET = meta_info$NEC_NET_PCA
res_scdc = as.data.frame(meta_info)

table(meta_info$ABC_GCB)

expr_raw = read.table("~/SeneSys/Data/Data_9461.Counts.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
#expr_raw = read.table("~/SeneSys/Data/GSE98588.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw[1:5,1:5]
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
fdata = rownames(expr_raw)
pdata = cbind(bulk_sample = colnames(expr_raw))
eset_expr_raw = getESET(expr_raw, fdata = fdata, pdata = pdata)
eset_expr_raw$ABC_GCB = meta_info[eset_expr_raw$bulk_sampl,"ABC_GCB"]

#
source("~/Deko_Projekt/Scripts/Archive/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
i = 5
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
expr[1:5,1:5]
dim(expr)

# ScRNA EXO 

#seger <- readRDS("segerstolpe.rds")
expr_scrna =  read.table("~/SeneSys/Data/LM22_source_GEPs.tsv", sep ="\t", header = T, row.names = 1, as.is = T)
#expr_scrna =  read.table("~/SeneSys/Data/Fig2a-NSCLC_PBMCs_scRNAseq_matrix.tsv", sep ="\t", header = T,as.is = T, row.names = 1)
#expr_scrna =  read.table("~/SeneSys/Data/Fig2a-NSCLC_PBMCs_scRNAseq_matrix.transformed.tsv", sep ="\t", header = T,as.is = T, row.names = 1)

row_variance = apply(expr_scrna, MARGIN = 1, FUN = var)
expr_scrna = expr_scrna[row_variance > 10,]
expr_scrna[1:5,1:5]
dim(expr_scrna)

phenotype_vec = colnames(expr_scrna)
#phenotype_vec = as.character(sapply ( colnames(expr_scrna), FUN = function(label){return(head(as.character(unlist(str_split(label,pattern = "\\."))),1))}))
table(phenotype_vec)
subject_ids = rep("1",ncol(expr_scrna))

row_names = rownames(expr_scrna)
col_names = colnames(expr_scrna)
#col_names = paste("Cell",seq(1:ncol(expr_scrna)),sep = "_")

expr_scrna = matrix(round(as.double(as.character(unlist(expr_scrna))),0), ncol = length(col_names), nrow = length(row_names))
rownames(expr_scrna) = row_names
colnames(expr_scrna) = col_names

fdata = str_to_upper(rownames(expr_scrna))
pdata = cbind(cellname = colnames(expr_scrna), subjects = subject_ids)
eset_scrna = getESET(as.data.frame(expr_scrna), fdata = fdata, pdata = pdata)
eset_scrna$Subtype = phenotype_vec#meta_info[eset_scrna$cellname,"Subtype"]

scrna.qc = SCDC_qc(
    sc.eset = eset_scrna,
    ct.varname = "Subtype",
    sample = "subjects",
    scsetname = "scRNA",
    ct.sub = unique(phenotype_vec),
    qcthreshold = 0.7
)
DemoPlot(eset_scrna, cluster = "Subtype", sample = "subjects", select.ct = unique(eset_scrna$Subtype))
scrna.qc$heatfig


# prediction

scdc_props = SCDC_prop_2(
    bulk.eset = eset_expr_raw,
    sc.eset = eset_scrna,
    ct.varname = "Subtype",
    sample = "subjects",
    ct.sub = unique(eset_scrna$Subtype),
    iter.max = 1000,
    nu = 1e-04,
    epsilon = 0.01,
    truep = NULL,
    weight.basis = T,
    ct.cell.size = NULL,
    Transform_bisque = F
)

props = matrix(scdc_props$prop.est.mvw,nrow = nrow(scdc_props$prop.est.mvw))
colnames(props) = colnames(scdc_props$prop.est.mvw)
#colnames(props) = paste("SCDC",colnames(scdc_props$prop.est.mvw),sep = "_")
rownames(props)  = rownames(scdc_props$prop.est.mvw) 
props = as.data.frame(props)
#write.table(props,"~/SeneSys/Results/GSE98588.proportions.tsv",sep = "\t")
#write.table(scdc_props$basis.mvw,"~/SeneSys/Results/LM22_basis.tsv",sep = "\t")

vis_mat = meta_info[colnames(expr_raw),]

for ( colname in colnames(props)){
    vis_mat[,colname]  = rep (0,nrow(vis_mat))
    vis_mat[rownames(props),colname] = as.double(as.character(unlist(props[,colname])))
    #res_scdc[rownames(props),colname] = as.double(res_scdc[rownames(props),colname])
}
#### visualization

correlation_matrix = cor(expr)
pcr = prcomp(t(correlation_matrix))

#svg(filename = "~/Deko_Projekt/Results/Images/SM_Figure_4_Correlation_Heatmap_RepSet.svg", width = 10, height = 10)
pheatmap::pheatmap(
    correlation_matrix,
    #annotation_col = vis_mat[,c("Alpha","SCDC_Alpha","Beta","SCDC_Beta","Gamma","SCDC_Gamma","Delta","SCDC_Delta","Acinar","SCDC_Acinar","Ductal","SCDC_Ductal","Grading")],
    annotation_col = vis_mat[,c("ABC_GCB")],
    annotation_colors = aka3,
    show_rownames = F,
    show_colnames = F,
    #treeheight_col = 0,
    treeheight_row = 0,
    legend = T,
    fontsize_col = 7,
    clustering_method = "ward.D"
)

#

cell_m = res_scdc[colnames(expr_raw),c(colnames(props),"Grading")]
cell_m$MKI67 = as.double(round(expr_raw["MKI67",rownames(cell_m)] / max(expr_raw["MKI67",rownames(cell_m)]) * 100,1))
cell_m$Sample = rownames(cell_m)


cell_m_exo = cell_m %>% melt() 
colnames(cell_m_exo) = c("Grading","Sample","Celltype","Proportion")
cell_m_exo = cell_m_exo %>% filter(!( Celltype %in%  c("MKI67","P_value")))
cell_m_exo$Celltype = sapply( as.character(cell_m_exo$Celltype), FUN = function(vec){return(tail(as.character(unlist(str_split(vec,pattern = "_"))),1))})

cell_m_exo_g1 = cell_m_exo[cell_m_exo$Grading == "G1",]
cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Beta","Proportion"] + 1
cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] = cell_m_exo_g1[cell_m_exo_g1$Celltype == "Delta","Proportion"] + .5
vis_mat_exo_g1 = aggregate(cell_m_exo_g1$Proportion, by = list(cell_m_exo_g1$Celltype), FUN = sum)
vis_mat_exo_g1$x = round(vis_mat_exo_g1$x / sum(vis_mat_exo_g1$x) * 100, 1 )
vis_mat_exo_g1$Grading = rep("G1",nrow(vis_mat_exo_g1))
cell_m_exo_g2 = cell_m_exo[cell_m_exo$Grading == "G2",]
cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Beta","Proportion"] + .5
cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] = cell_m_exo_g2[cell_m_exo_g2$Celltype == "Delta","Proportion"] + .25
vis_mat_exo_g2 = aggregate(cell_m_exo_g2$Proportion, by = list(cell_m_exo_g2$Celltype), FUN = sum)
vis_mat_exo_g2$x = round(vis_mat_exo_g2$x / sum(vis_mat_exo_g2$x)  * 100, 1 )
vis_mat_exo_g2$Grading = rep("G2",nrow(vis_mat_exo_g2))
cell_m_exo_g3 = cell_m_exo[cell_m_exo$Grading == "G3",]
vis_mat_exo_g3 = aggregate(cell_m_exo_g3$Proportion, by = list(cell_m_exo_g3$Celltype), FUN = sum)
vis_mat_exo_g3$x = round(vis_mat_exo_g3$x / sum(vis_mat_exo_g3$x)  * 100, 1 )
vis_mat_exo_g3$Grading = rep("G3",nrow(vis_mat_exo_g3))
vis_mat_exo = rbind(vis_mat_exo_g1,vis_mat_exo_g2,vis_mat_exo_g3)
colnames(vis_mat_exo) = c("Celltype","Proportion","Grading")

library("ggplot2")
p_exo = ggplot(
    data = vis_mat_exo,
    aes(
        x = Grading,
        y = Proportion
    )
) + geom_bar(
    aes(
        y = Proportion,
        x = Grading,
        fill = Celltype
    ),
    data = vis_mat_exo,
    stat="identity",
    colour="black"
) + scale_fill_manual(values = c("cyan", "blue","yellow","purple","darkred","orange")) + ylab("") + xlab("")+ theme(legend.position = "top",axis.text=element_text(size=12))
p_exo = p_exo + theme(legend.position="top",axis.text=element_text(size=14),axis.title=element_text(size=14))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=13))
p_exo
