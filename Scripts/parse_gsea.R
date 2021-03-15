library("stringr")

expr_raw = read.table("~/SeneSys/Data/GSE98588_series_matrix.prepared.tsv",sep ="\t",header = T)
expr_raw[1:5,1:5]
dim(expr_raw)

expr_raw = cbind(rownames(expr_raw),expr_raw)
write.table(expr_raw, "~/SeneSys/Data/GSE98588_series_matrix.prepared.tsv",sep ="\t", row.names = F,quote = F)