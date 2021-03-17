library("devtools")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("SCDC")

meta_info = read.table("~/SeneSys///Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

table(meta_info$ABC_GCB)

#expr_raw = read.table("~/SeneSys/Data/Data_9461.Counts.HGNC.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw = read.table("~/SeneSys/Data/GSE98588.tsv",sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
expr_raw[1:5,1:5]
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
fdata = rownames(expr_raw)
pdata = cbind(bulk_sample = colnames(expr_raw))
eset_expr_raw = getESET(expr_raw, fdata = fdata, pdata = pdata)
eset_expr_raw$ABC_GCB = meta_info[eset_expr_raw$bulk_sampl,"ABC_GCB"]

# ScRNA EXO 

#seger <- readRDS("segerstolpe.rds")
#expr_scrna =  read.table("~/SeneSys/Data/LM22_source_GEPs.tsv", sep ="\t", header = T, row.names = 1, as.is = T)
expr_scrna =  read.table("~/SeneSys/Data/Fig2a-NSCLC_PBMCs_scRNAseq_matrix.tsv", sep ="\t", header = T,as.is = T, row.names = 1)
#expr_scrna =  read.table("~/SeneSys/Data/Fig2a-NSCLC_PBMCs_scRNAseq_matrix.transformed.tsv", sep ="\t", header = T,as.is = T, row.names = 1)

row_variance = apply(expr_scrna, MARGIN = 1, FUN = var)
expr_scrna = expr_scrna[row_variance > 10,]
expr_scrna[1:5,1:5]
dim(expr_scrna)

phenotype_vec = colnames(expr_scrna)
phenotype_vec = as.character(sapply ( colnames(expr_scrna), FUN = function(label){return(head(as.character(unlist(str_split(label,pattern = "\\."))),1))}))
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

scrna.qc = SCDC_qc_ONE(
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

scdc_props = SCDC_prop_ONE(
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
#write.table(props,"~/SeneSys/Results/GSE98588.proportions.SCDC.tsv",sep = "\t")
