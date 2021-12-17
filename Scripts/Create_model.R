library("devtools")
library("NMF")
load_all("~/artdeco/")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")

expr_raw = read.table("~/Downloads/senescence_mouse.tsv",sep ="\t", header = T, row.names = 1)

dd = colnames(expr_raw)
colnames(expr_raw) = str_replace(colnames(expr_raw),"^X","")
expr_raw[1:5,1:5]
dim(expr_raw)

meta_info = read.table("~/SeneSys//Misc/Meta_information.tsv", sep ="\t", header = T)
rownames(meta_info) = meta_info$Sample

adr_oht = which(str_count(colnames(expr_raw), "_") == 2)
subtype_vector = as.character(sapply(colnames(expr_raw), FUN = function(vec){return(head(as.character(unlist(str_split(vec,"_"))),1))}))
subtype_vector[adr_oht] = "ADR_OHT"
table(subtype_vector)

#candidates = which(subtype_vector %in% c("alpha","beta","gamma","delta","acinar-s","acinar-reg+","acinar-i","ductal","muc5b+ ductal"))
#candidates = which(subtype_vector %in% c("Alpha","Beta","Gamma","Delta"))
#candidates = which(subtype_vector %in% c("Alpha","Beta","Gamma","Delta","Acinar","Ductal"))
#expr_raw = expr_raw[,candidates]
#meta_data = meta_info[colnames(expr_raw),]
#subtype_vector = meta_data$Cluster
table(subtype_vector)

#subtype_vector = as.character(str_remove_all(colnames(expr_raw), pattern = "_[0-9]*"))
#table(subtype_vector_reduced)

amount_genes = 400
amount_samples = 300
model_name = "Senescence_400_genes_300_cells_3_classes"

selected_samples = c()

for ( cell_type in unique(subtype_vector)){
    coords = which(subtype_vector == cell_type )
    
    if (length(coords) >= amount_samples)
        coords = sample(coords, size = amount_samples)
    
    selected_samples = c(selected_samples, coords)
}
length(selected_samples)

expr = expr_raw[,selected_samples]
dim(expr)
#meta_data_reduced = meta_info[colnames(expr),]
subtype_vector_reduced = subtype_vector[selected_samples]
table(subtype_vector_reduced)

rownames(expr) = str_to_upper(rownames(expr))

add_deconvolution_training_model_bseqsc(
    transcriptome_data = expr, 
    model_name = model_name,
    subtype_vector =  subtype_vector_reduced,
    training_p_value_threshold = 0.05,
    training_nr_permutations = 1000,
    training_nr_marker_genes = amount_genes
)
