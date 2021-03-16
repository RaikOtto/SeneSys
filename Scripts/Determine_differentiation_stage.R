library("devtools")
load_all("~/artdeco/")

meta_info = read.table("~/SeneSys//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

### add models

scRNA_file_path = "~/SeneSys/Data/Fig2a-NSCLC_PBMCs_scRNAseq_matrix.tsv"
model_name = str_replace_all(scRNA_file_path,pattern = "\\.tsv","")
model_name = "PMBC"#tail(str_split(model_name,pattern = "/")[[1]],1)

t = read.table(scRNA_file_path, sep ="\t", header = T, row.names = 1)
colnames(t) = str_replace_all(colnames(t),pattern ="\\.","_")
colnames(t) = str_replace_all(colnames(t),pattern ="^X","")
rownames(t) = str_replace_all(rownames(t),pattern ="\\.","_")
rownames(t) = str_replace_all(rownames(t),pattern ="-","")
rownames(t) = str_replace_all(rownames(t),pattern ="_","")
dim(t)
t[1:5,1:5]
#meta_data = meta_info[colnames(t),]
#subtype_vector = str_to_lower(meta_data$Subtype)
subtype_vector = str_replace_all(colnames(t),pattern = "_([A-Za-z0-9])*$","")
subtype_vector[subtype_vector == "B"] = "B_cells"
subtype_vector[subtype_vector == "NK"] = "NK_cells"
subtype_vector[subtype_vector == "NKT"] = "NKT_cells"
subtype_vector[subtype_vector == "CD14_"] = "CD14_cells"
subtype_vector[subtype_vector == "CD4_T"] = "CD4_T_cells"
subtype_vector[subtype_vector == "CD8_T"] = "CD8_T_cells"
subtype_vector[subtype_vector == "CD14_cells"] = "CD14__Monocytes"
table(subtype_vector)

#transcriptome_data = read.table(scRNA_file_path,sep="\t",header  = T)

add_deconvolution_training_model_bseqsc(
    transcriptome_data = t,
    model_name = model_name,
    subtype_vector =  str_to_lower(subtype_vector),
    training_p_value_threshold = 0.05,
    training_nr_permutations = 0,
    training_nr_marker_genes = 400
)

add_deconvolution_training_model_music(
    transcriptome_data = transcriptome_data,
    model_name = model_name,
    subtype_vector
)

add_deconvolution_training_model_NMF(
    transcriptome_data = transcriptome_data,
    model_name = model_name,
    subtype_vector = subtype_vector,
    parallel_processes = 8
)
