library("stringr")

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
sort(table(meta_info$Sample),decreasing = TRUE)

#path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.HGNC.tsv"
#path_transcriptome_file = "~/SeneSys/Data/Schmitz.HGNC.tsv"
path_transcriptome_file = "~/SeneSys/Data/GSE98588.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

meta_data = meta_info[colnames(expr_raw),]
#expr_raw = expr_raw[,meta_data$Drug_Treatment %in% c("NR","RP")]
#meta_data = meta_info[colnames(expr_raw),]
### cls file

cohort_vec = meta_data$ABC_GCB
cohort_vec = str_replace(cohort_vec," ","_")
expr = expr_raw[,cohort_vec %in% c("ABC", "GCB" )]
meta_data = meta_info[colnames(expr),]
cohort_vec = meta_data$ABC_GCB
table(cohort_vec)

first_entity = unique(cohort_vec)[1]
second_entity = unique(cohort_vec)[2]
third_entity = unique(cohort_vec)[3]

first_line = as.character(c(ncol(expr),length(unique(cohort_vec)),1))
first_line[4:(length(cohort_vec))] = ""
second_line = as.character(c("#",first_entity,second_entity,third_entity))
second_line[5:(length(cohort_vec))] = ""

cls_file = rbind(first_line,second_line,cohort_vec)

write.table(cls_file, "~/SeneSys/GSEA/Schmitz.ABC_GCB_new_new.cls",row.names = F, quote = F, sep ="\t",col.names = FALSE)
write.table(cls_file, "~/SeneSys/GSEA/Schmitz.ABC_GCB_new_new.cls.tsv",row.names = F, quote = F, sep ="\t",col.names = FALSE)

### gct file

first_line = "#1.0"
first_line[2:(ncol(expr)+2)] = ""
second_line = as.character(c(nrow(expr),ncol(expr)))
second_line[3:(ncol(expr)+2)] = ""
third_line = c("Name","Description",colnames(expr))
gct_matrix = cbind(rownames(expr),rownames(expr),expr)

output_gct = rbind(first_line,second_line, third_line, gct_matrix)
output_gct[1:5,1:5]

write.table(output_gct, "~/SeneSys/GSEA/GSE98588_new_new.gct",row.names = F, quote = F, sep ="\t",col.names = FALSE)
write.table(output_gct, "~/SeneSys/GSEA/GSE98588_new_new.gct.tsv",row.names = F, quote = F, sep ="\t",col.names = FALSE)
