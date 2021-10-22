library("stringr")

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample

path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.HGNC.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

meta_data = meta_info[colnames(expr_raw),]
#expr_raw = expr_raw[,meta_data$Drug_Treatment %in% c("NR","RP")]
#meta_data = meta_info[colnames(expr_raw),]
### cls file

cohort_vec = meta_data$Drug_Treatment
cohort_vec = str_replace(cohort_vec," ","_")
#expr = expr_raw[,cohort_vec %in% c("ABC", "GCB" )]
#meta_data = meta_info[colnames(expr),]
#cohort_vec = meta_data$ABC_GCB
table(cohort_vec)

first_entity = unique(cohort_vec)[1]
second_entity = unique(cohort_vec)[2]

first_line = as.character(c(ncol(expr_raw),length(unique(cohort_vec)),1))
first_line[4:(length(cohort_vec))] = ""
second_line = as.character(c("#",first_entity,second_entity))
second_line[4:(length(cohort_vec))] = ""

cls_file = rbind(first_line,second_line,cohort_vec)

write.table(cls_file, "~/SeneSys/GSEA/Data_9461_NR_RP_only.cls",row.names = F, quote = F, sep ="\t",col.names = FALSE)
write.table(cls_file, "~/SeneSys/GSEA/Data_9461_NR_RP_only.cls.tsv",row.names = F, quote = F, sep ="\t",col.names = FALSE)

### gct file

first_line = "#1.0"
first_line[2:(ncol(expr_raw)+2)] = ""
second_line = as.character(c(nrow(expr_raw),ncol(expr_raw)))
second_line[3:(ncol(expr_raw)+2)] = ""
third_line = c("Name","Description",colnames(expr_raw))
gct_matrix = cbind(rownames(expr_raw),rownames(expr_raw),expr_raw)

output_gct = rbind(first_line,second_line, third_line, gct_matrix)
output_gct[1:5,1:5]

write.table(output_gct, "~/SeneSys/GSEA/_NR_RP_only.gct",row.names = F, quote = F, sep ="\t",col.names = FALSE)
write.table(output_gct, "~/SeneSys/GSEA/_NR_RP_only.gct.tsv",row.names = F, quote = F, sep ="\t",col.names = FALSE)
