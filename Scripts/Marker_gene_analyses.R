library("stringr")

schlesinger_t = read.table("~/Deko_Projekt/Misc//Schlesinger.tsv", header = T, row.names = 1, sep ="\t")
schlesinger_genes = rownames(schlesinger_t)

baron = readRDS("~/artdeco/inst/Models/bseqsc//Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.RDS")
marker_genes_baron = baron[[3]]
sad_genes = marker_genes_baron_ductal = marker_genes_baron$ductal

round((sum(marker_genes_baron_ductal %in% schlesinger_genes == TRUE) / length(marker_genes_baron_ductal)) * 100,0)
sad_genes = meta_plastic_genes = marker_genes_baron_ductal[marker_genes_baron_ductal %in% schlesinger_genes == TRUE]
non_meta_plastic_genes = marker_genes_baron_ductal[marker_genes_baron_ductal %in% schlesinger_genes == FALSE]

### manipulation

baron = readRDS("~/artdeco/inst/Models/bseqsc//Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.RDS")
marker_genes_baron = baron[[3]]

marker_genes_baron$ductal = meta_plastic_genes
baron[[3]] = marker_genes_baron
signature_matrix = baron[[1]]
dim(signature_matrix)
signature_matrix = signature_matrix[rownames(signature_matrix) %in% meta_plastic_genes,]
dim(signature_matrix)
baron[[1]] = signature_matrix
saveRDS( baron, "~/artdeco/inst/Models/bseqsc//Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron_metaplastic.RDS")


baron = readRDS("~/artdeco/inst/Models/bseqsc//Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.RDS")
marker_genes_baron = baron[[3]]
marker_genes_baron$ductal = non_meta_plastic_genes
baron[[3]] = marker_genes_baron
signature_matrix = baron[[1]]
dim(signature_matrix)
signature_matrix = signature_matrix[!(rownames(signature_matrix) %in% meta_plastic_genes),]
dim(signature_matrix)
baron[[1]] = signature_matrix
saveRDS( baron, "~/artdeco/inst/Models/bseqsc//Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron_non_metaplatic.RDS")

### Segerstolpe
    
segerstolpe = readRDS("~/Downloads/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe.RDS")
marker_genes_segerstolpe = segerstolpe[[3]]
marker_genes_segerstolpe_ductal = marker_genes_segerstolpe$ductal

(sum(marker_genes_segerstolpe_ductal %in% schlesinger_genes == TRUE) / length(marker_genes_segerstolpe_ductal)) * 100

lawlor = readRDS("~/Downloads/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor.RDS")
marker_genes_lawlor = lawlor[[3]]
marker_genes_lawlor_ductal = marker_genes_lawlor$ductal

(sum(marker_genes_lawlor_ductal %in% schlesinger_genes == TRUE) / length(marker_genes_lawlor_ductal)) * 100

tosti = readRDS("~/artdeco/inst/Models/bseqsc/Baron_Tosti_800_genes_100_samples_all_endocrine_all_exocrine_no_metaplastic.RDS")
marker_genes_tosti = tosti[[3]]
names(marker_genes_tosti)
marker_genes_tosti_ductal = marker_genes_tosti$ductal

(sum(marker_genes_tosti_ductal %in% schlesinger_genes == TRUE) / length(marker_genes_tosti_ductal)) * 100
length(marker_genes_tosti_ductal)

marker_genes_tosti_new = marker_genes_tosti$alpha
(sum(marker_genes_tosti_new %in% schlesinger_genes == TRUE) / length(marker_genes_tosti_new)) * 100

####

baron = readRDS("~/artdeco/inst/Models/bseqsc/Senescence_400_genes_300_cells_3_classes.RDS")
marker_genes_baron = baron[[3]]

sad_genes = marker_genes_baron$ADR_OHT[1:100]
sad_genes = marker_genes_baron$PS[1:50]


baron[[3]] = marker_genes_baron
signature_matrix = baron[[1]]
dim(signature_matrix)
signature_matrix = signature_matrix[rownames(signature_matrix) %in% meta_plastic_genes,]
dim(signature_matrix)
baron[[1]] = signature_matrix
saveRDS( baron, "~/artdeco/inst/Models/bseqsc//Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron_metaplastic.RDS")
