library("devtools")
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")
library("stringr")
library("reshape2")
library("bseqsc")
library("dplyr")

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = meta_info[meta_info$Sample!="",]
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

#

#i_filename = "~/SeneSys/Data/Data_9461.Counts.HGNC.tsv"
#i_filename = "~/SeneSys/Data/GSE98588_new.HGNC.tsv"
i_filename = "~/SeneSys/Data/Schmitz.HGNC.tsv"
expr_raw = read.table(i_filename,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

meta_data = meta_info[colnames(expr_raw),]
table(meta_data$ABC_GCB)

#show_models_bseqsc()
model_name = "Senescence_400_genes_300_cells_3_classes"

props = Deconvolve_transcriptome(
    transcriptome_data = expr_raw,
    deconvolution_algorithm = "bseqsc",
    models = model_name,
    Cibersort_absolute_mode = TRUE,
    nr_permutations = 1000,
    output_file = ""
)

matcher_ori = match(rownames(props),meta_info$Sample, nomatch = 0)
meta_info[matcher_ori[matcher_ori != 0],c("ADR","ADROHT","PS")] = props[matcher_ori != 0,c("ADR","ADR_OHT","PS")]

meta_data = meta_info[colnames(expr_raw)[matcher_ori != 0],]

adr = meta_data$ADR
adr[adr < mean(meta_data$ADR) ] = "Low"
#adr[adr < 1000 ] = "Low"
adr[adr != "Low"] = "High"
meta_data$ADR_binary = adr

ADROHT = meta_data$ADROHT
ADROHT[ADROHT < mean(meta_data$ADROHT) ] = "Low"
#ADROHT[ADROHT < 2700 ] = "Low"
ADROHT[ADROHT != "Low"] = "High"
meta_data$ADROHT_binary = ADROHT

PS = meta_data$PS
PS[PS < mean(meta_data$PS) ] = "Low"
PS[PS != "Low"] = "High"
meta_data$PS_binary = PS

meta_info[meta_data$Sample,c("ADR_binary","ADROHT_binary","PS_binary")] = meta_data[,c("ADR_binary","ADROHT_binary","PS_binary")]

#write.table(meta_info,"~/SeneSys/Misc/Meta_information.tsv",sep = "\t",row.names = FALSE,quote =FALSE)

table(meta_data$ADROHT_binary, meta_data$Progression)
table(meta_data$ADR_binary, meta_data$Progression)
table(meta_data$PS_binary, meta_data$Progression)

#write.table(props,"~/Downloads/Data_9461_deco.tsv",sep ="\t", row.names = TRUE, quote =F )

table(meta_data$Progression,meta_data$ADROHT_binary)

vec = meta_data$Progression
vec[meta_data$Progression == "No Progression"] = "0"
vec[vec != "0"] = "1"
plot(meta_data$PS, as.integer(vec))

cor_mat = cor(t(props[,]))
heatmap(cor_mat)

matcher = match(colnames(expr_raw), meta_info$Sample, nomatch = 0)
vis_mat = meta_info[matcher,]
aggregate(vis_mat$ADR, by = list(vis_mat$Progression), FUN = mean)

prog = vis_mat$PS[vis_mat$Progression == "Progression"] 
no_prog = vis_mat$PS[vis_mat$Progression == "No Progression"] 
t.test(prog,no_prog)

table(meta_data$Progression,meta_data$PS_binary)
