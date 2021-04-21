library(stringr)
library("survival")

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = meta_info[meta_info$Study != "GSE11318",]
rownames(meta_info) = meta_info$Sample

#meta_info_11318 = meta_info[ meta_info$Study == "GSE11318",]
#rownames(meta_info_11318) = meta_info_11318$Sample
#gsva_11318 = as.data.frame(read.table("~/SeneSys/Results/GSE11318.GSVA.tsv",sep = "\t",header = T,stringsAsFactors = F))
#gsva_11318[1:5,1:5]
#meta_data = meta_info_11318[colnames(gsva_11318),]

#meta_info = meta_info[(meta_info$Sample!="") & (!(is.na(meta_info$Sample))),]
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

expr_raw = read.table("~/SeneSys/Data/GSE136971.HGNC.tsv",sep ="\t", as.is = T,header = T, row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
meta_data = meta_info[colnames(expr_raw),]
meta_data$OS = as.double(str_replace_all(meta_data$OS_Tissue, pattern = ",", "."))
meta_data = meta_data %>% dplyr::filter( ! is.na(meta_data$OS) )
meta_data = meta_data %>% dplyr::filter( meta_data$OS != "na" )
#meta_data = meta_data %>% dplyr::filter( meta_data$ABC_GCB %in% c("ABC","GCB") )

#agg = aggregate(meta_data$SUVARNESS, FUN = mean, by = list(meta_data$Cluster))
#agg = aggregate(meta_data$SUVARNESS, FUN = mean, by = list(meta_data$ABC_GCB))
values = meta_data[,"SUVARNESS"]
#values = meta_data[,"FRIDMAN_SENESCENCE_UP"]
hist(values)
#values = as.double(gsva_11318["SUVARNESS",])

thresh_mid = median(values)
thresh_mid = mean(values)

classification_vec = values
classification_vec[  values <= thresh_mid] = "low"
classification_vec[  values > thresh_mid] = "high"

fit = survival::survfit(
    survival::Surv(
        as.double(meta_data$Study), as.integer(meta_data$Censor_OS)
    ) ~ classification_vec
)
surv_p_value = survminer::surv_pvalue(fit, data = meta_data)$pval
surv_p_value

#pdf(graphics_path_survival_hisc,onefile = FALSE)#,width="1024px",height="768px")
#svg(filename = "~/Deko_Projekt/Results/Images/Figure_5_survival_grading_three.svg", width = 10, height = 10)
#svg(filename = "~/Deko_Projekt/Results/Images/Figure_5_survival_ductal_three.svg", width = 10, height = 10)
    print(survminer::ggsurvplot(fit, data = meta_data, risk.table = F, pval = T, censor.size = 10))
#dev.off()
