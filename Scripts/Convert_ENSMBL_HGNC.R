library("stringr")


hgnc_list_uni = as.character( unique(hgnc_list) )

length(hgnc_list)
length(rownames(expr_raw))

### translate gpl80

mapping_t = read.table("~/SeneSys/Data/Mappings/GPL80-30376.tsv",sep ="\t", header = T,stringsAsFactors = T)
rownames(mapping_t) = mapping_t$ID
mapping_t[1:5,]

### parsing GSE68895

expr_raw = read.table("~/SeneSys/Data/GSE68895.tsv",sep ="\t", header = T,comment.char = '!',row.names = 1,stringsAsFactors = F)
dim(expr_raw)
expr_raw[1:5,1:5]

library("biomaRt")
affyids = rownames(expr_raw)

datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

getBM(
    attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol'),
    filters = 'affy_hg_u133_plus_2', 
    values = affyids, 
    mart = ensembl
)

gene_set = read.table("~/SeneSys/Data/SeneSys_gene_sets.gmt",header = F, stringsAsFactors = F, sep ="\t")

hgnc_list = as.character(unlist(mapping_t[ as.character(rownames(expr_raw)) ,"HGNC"]))
hgnc_list = as.character(unlist(sapply(hgnc_list, FUN= function(vec){return(str_split(vec,pattern = "/")[1])})))
hgnc_list = str_trim(hgnc_list)

write.table(expr_raw,"~/SeneSys/Data/GSE68895_series_matrix.gct", sep ="\t", quote = F,row.names = T)

#### GSE98588

library(org.Hs.eg.db)

expr_raw = read.table("~/SeneSys/Data/GSE98588_series_matrix.tsv",sep ="\t", header = T,comment.char = '!',row.names = 1,stringsAsFactors = F)
expr_raw[1:5,1:5]
dim(expr_raw)

proper_names= str_replace_all(rownames(expr_raw), pattern = "_at", "")

library("biomaRt")
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#listMarts(ensembl)
#datasets <- listDatasets(ensembl)
mart <- useEnsembl(
    biomart = "ensembl", 
    dataset = "hsapiens_gene_ensembl"#, 
    #mirror = "useast"
)

mapping_t = getBM(
    attributes = c('hgnc_symbol','ensembl_gene_id'),
    filters = 'ensembl_gene_id', 
    values = sad_genes, 
    mart = mart
)

# hgnc_id
mapping_t = getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    filters = 'ensembl_gene_id', 
    values = proper_names, 
    mart = mart
)
mapping_t = as.data.frame(as.matrix(mapping_t))
indices = match(proper_names , mapping_t$ensembl_gene_id, nomatch = 0)

mapping_vec = rep("",length(proper_names))
mapping_vec[indices != 0] = as.character(unlist(mapping_t[indices,2]))
map_t = cbind(proper_names, mapping_vec)
#write.table(map_t, "~/SeneSys/Misc/GSE98588_mapping_table_ENSEMBL_HGNC.tsv",sep="\t",row.names = F, quote = F)

expr_raw = expr_raw[!is.na(indices),]
proper_names = proper_names[!is.na(indices)]

mapping_t = mapping_t[ mapping_t$ensembl_gene_id %in% proper_names,]
#which(mapping_t$ensembl_gene_id == "ENSG00000230417")
#mapping_t = mapping_t[-11059,]
rownames(mapping_t) = mapping_t$ensembl_gene_id
hgnc_list = as.character(mapping_t[ proper_names, 2])

dim(expr_raw)
length(hgnc_list)

# var selection

expr_raw[1:5,1:5]
dim(expr_raw)
write.table(expr_raw, "~/SeneSys/Data/GSE98588.tsv",sep="\t",row.names = T, quote = F)

# .cls creation

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv", sep = "\t", header = T, stringsAsFactors = F)
meta_info = meta_info[meta_info$Sample != "",]
rownames(meta_info) = meta_info$Sample
meta_data = meta_info[colnames(expr_raw),]
dim(meta_data)

export_data = meta_data[,c("Sample","Complex")]
dim(export_data)

#write.table(export_data, "~/SeneSys/Data/GSE98588.cluster.cls",sep="\t",row.names = F, quote = F)

#### GSE34171

library("hgu95av2.db")
mapping_t = as.data.frame(hgu95av2SYMBOL)
rownames(mapping_t) = mapping_t$probe_id
mapping_t[which(mapping_t$symbol == "DDR1"), ]

expr_raw = read.table("~/SeneSys/Data/GSE34171.cls",sep ="\t", header = T,row.names = 1,stringsAsFactors = F)
expr_raw[1:5,1:5]
hgnc_list = mapping_t[rownames(expr_raw),2]
sample_names = colnames(expr_raw)

library("biomaRt")
affyids = rownames(expr_raw)

datasets <- listDatasets(ensembl)

ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

getBM(
    attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol'),
    filters = 'affy_hg_u133_plus_2', 
    values = affyids, 
    mart = ensembl
)

#### GSE10846

library("hgu133plus2.db")
mapping_t = as.data.frame(hgu133plus2SYMBOL)
rownames(mapping_t) = mapping_t$probe_id

expr_raw = read.table("~/SeneSys/Data/GSE10846.tsv",sep ="\t", header = T,row.names = 1,stringsAsFactors = F)
expr_raw[1:5,1:5]
hgnc_list = mapping_t[rownames(expr_raw),2]
table(!is.na(hgnc_list))

#var selection

write.table(expr_raw, "~/SeneSys/Data/GSE10846.tsv",sep="\t",row.names = T, quote = F)

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv", sep = "\t", header = T, stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
meta_data = meta_info[colnames(expr_raw),]
export_data = meta_data[,c("Sample","ABC_GCB")]
dim(export_data)

write.table(export_data, "~/SeneSys/Data/GSE10846.cls",sep="\t",row.names = F, quote = F)

#### GSE34171

library("hgu133plus2.db")
mapping_t = as.data.frame(hgu133plus2SYMBOL)
rownames(mapping_t) = mapping_t$probe_id

expr_raw = read.table("~/SeneSys/Data/GSE34171.gct",sep ="\t", header = T,row.names = 1,stringsAsFactors = F)
expr_raw[1:5,1:5]
hgnc_list = mapping_t[rownames(expr_raw),2]
table(!is.na(hgnc_list))

write.table(expr_raw, "~/SeneSys/Data/GSE34171.gct",sep="\t",row.names = T, quote = F)
dim(expr_raw)

#variance selection

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv", sep = "\t", header = T, stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
meta_data = meta_info[colnames(expr_raw),]
export_data = meta_data[,c("Sample","ABC_GCB")]
dim(export_data)

write.table(export_data, "~/SeneSys/Data/GSE10846.cls",sep="\t",row.names = F, quote = F)


### 45630

library("hgu133plus2.db")
mapping_t = as.data.frame(hgu133plus2SYMBOL)
rownames(mapping_t) = mapping_t$probe_id

expr_raw = read.table("~/SeneSys/Data/GSE34171.gct",sep ="\t", header = T,row.names = 1,stringsAsFactors = F)
expr_raw[1:5,1:5]
hgnc_list = mapping_t[rownames(expr_raw),2]
table(!is.na(hgnc_list))

write.table(expr_raw, "~/SeneSys/Data/GSE34171.gct",sep="\t",row.names = T, quote = F)
dim(expr_raw)

#variance selection

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv", sep = "\t", header = T, stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
meta_data = meta_info[colnames(expr_raw),]
export_data = meta_data[,c("Sample","ABC_GCB")]
dim(export_data)

write.table(export_data, "~/SeneSys/Data/GSE10846.cls",sep="\t",row.names = F, quote = F)

### GES11318 

library("hgu133plus2.db")
mapping_t = as.data.frame(hgu133plus2SYMBOL)
rownames(mapping_t) = mapping_t$probe_id

expr_raw = read.table("~/SeneSys/Data/GSE11318.gct",sep ="\t", header = T,row.names = 1,stringsAsFactors = F)
expr_raw[1:5,1:5]
hgnc_list = mapping_t[rownames(expr_raw),2]
table(!is.na(hgnc_list))

write.table(expr_raw, "~/SeneSys/Data/GSE11318.gct",sep="\t",row.names = T, quote = F)
dim(expr_raw)

#### GSE45630

library("biomaRt")
expr_raw = read.table("~/SeneSys/Data/GSE45630.gct",sep ="\t", header = T,row.names = 1,stringsAsFactors = F)
rownames(expr_raw) = str_replace_all(rownames(expr_raw), pattern = "_at", "")
expr_raw[1:5,1:5]

affyids = rownames(expr_raw)
mart <- useEnsembl(
    biomart = "ensembl", 
    dataset = "hsapiens_gene_ensembl", 
    mirror = "useast"
)
datasets <- listDatasets(mart)

mapping_t = getBM(
    attributes = c('entrezgene_id', 'hgnc_symbol'),
    filters = 'entrezgene_id', 
    values = str_replace_all(rownames(expr_raw), pattern = "_at", ""),
    mart = mart
)


hgnc_list = mapping_t[rownames(expr_raw),2]
table(!is.na(hgnc_list))

write.table(expr_raw, "~/SeneSys/Data/GSE45630.gct",sep="\t",row.names = T, quote = F)
dim(expr_raw)

#### GSE31312

library("biomaRt")
expr_raw = read.table("~/SeneSys/Data/GSE31312.gct",sep ="\t", header = T,row.names = 1,stringsAsFactors = F)
rownames(expr_raw) = str_replace_all(rownames(expr_raw), pattern = "_at", "")
expr_raw[1:5,1:5]


library("hgu133plus2.db")
mapping_t = as.data.frame(hgu133plus2SYMBOL)
rownames(mapping_t) = mapping_t$probe_id

hgnc_list = mapping_t[rownames(expr_raw),2]
table(!is.na(hgnc_list))

write.table(expr_raw, "~/SeneSys/Data/GSE31312.gct",sep="\t",row.names = T, quote = F)
dim(expr_raw)

#### GSE136971

library("biomaRt")
expr_raw = read.table("~/SeneSys/Data/GSE136971.gct",sep ="\t", header = T,row.names = 1,stringsAsFactors = F)
rownames(expr_raw) = str_replace_all(rownames(expr_raw), pattern = "_at", "")
expr_raw[1:5,1:5]


library("hgu133plus2.db")
mapping_t = as.data.frame(hgu133plus2SYMBOL)
rownames(mapping_t) = mapping_t$probe_id

hgnc_list = mapping_t[rownames(expr_raw),2]
table(!is.na(hgnc_list))

write.table(expr_raw, "~/SeneSys/Data/GSE136971.gct",sep="\t",row.names = T, quote = F)
dim(expr_raw)

