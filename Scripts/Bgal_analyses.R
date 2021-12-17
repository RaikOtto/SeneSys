library("DESeq2")
library("limma")
library("edgeR")
library("umap")
library("stringr")
library("ggplot2")
library("dplyr")
library("grid")

draw_colnames_45 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
    return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

#path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.DeSeq2.HGNC.tsv"
path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.HGNC.tsv"

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, as.is = F,row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

table(colnames(expr_raw) %in% meta_info$Sample)
meta_data = meta_info[colnames(expr_raw),]

### Bgal

bgal_mat = expr_raw[,!is.na(meta_data$Sa_Bgal)]
meta_data_bgal = meta_info[colnames(bgal_mat),]

bgal = meta_data_bgal$Sa_Bgal
hist(bgal)

bgal_binary = rep("high",length(bgal))
bgal_binary[bgal < mean(bgal) ] = "low"
meta_data_bgal$bgal_binary = bgal_binary

meta_info$bgal_binary = rep("",nrow(meta_info))
meta_info[colnames(bgal_mat),"bgal_binary"] = meta_data_bgal$bgal_binary

#write.table(meta_info,"~/SeneSys/Misc/Meta_information.tsv",sep ="\t", row.names = FALSE, quote =FALSE)

## ADR IC50

adr_mat = expr_raw[,!is.na(meta_data$ADR_IC50)]
meta_data_adr = meta_info[colnames(adr_mat),]

adr = log(meta_data_adr$ADR_IC50)
hist(adr)

adr_binary = rep("high",length(adr))
adr_binary[adr < mean(adr) ] = "low"
meta_data_adr$adr_binary = adr_binary

meta_info$adr_binary = rep("",nrow(meta_info))
meta_info[colnames(adr_mat),"adr_binary"] = meta_data_adr$adr_binary

#write.table(meta_info,"~/SeneSys/Misc/Meta_information.tsv",sep ="\t", row.names = FALSE, quote =FALSE)


###

design <- model.matrix(~0 + as.factor(meta_data_bgal$bgal_binary))
design <- model.matrix(~0 + as.factor(meta_data_adr$adr_binary))

dds_raw = DESeq2::DESeqDataSetFromMatrix(
    countData = round(adr_mat,0),
    colData = meta_data_adr,
    design = design,
    tidy = F
)

dds_raw = estimateSizeFactors(dds_raw)
#expr_raw = assay( DESeq2::varianceStabilizingTransformation(dds_raw) )

dds <- DESeq(dds_raw)
res = results(dds,)
res = res[!is.na(res$pvalue),]
summary(res$log2FoldChange)
res = res[order(res$pvalue,decreasing = F),]

#write.table(res, "~/SeneSys/Results/adr_ic50_dif_exp.tsv",sep ="\t", quote = F)
#write.table(res, "~/SeneSys/Results/bgal_dif_exp.tsv",sep ="\t", quote = F)
