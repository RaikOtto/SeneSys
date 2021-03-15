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

meta_info = read.table("~/Deco//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

#source("~/Deko/Scripts/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/Deco/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t[,1]
genes_of_interest_hgnc_t$V1
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]

path_transcriptome_file = "~/Deco/Data/Bench_data/MAPTor_NET.S57.tsv"
visualization_data_path = str_replace(path_transcriptome_file,pattern  ="\\.tsv",".vis.tsv")

expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

meta_data = meta_info[colnames(expr_raw),]
#meta_data_2 = rbind(meta_data_2,meta_data)
#table(meta_data_2$Study)
#table(meta_data_2$Location)

table(meta_data$Grading)
#meta_data = meta_data[which(meta_data[,"Histology"] == "Pancreatic_NEN"),]
expr_raw = expr_raw[,rownames(meta_data)]
#meta_data[colnames(expr_raw),"mki_67"] = log(as.double(expr_raw["MKI67",])+1)
#meta_data$mki_67[ which(meta_data$mki_67 > mean(meta_data$mki_67))] = "high"
#meta_data$mki_67[meta_data$mki_67 != "high"] = "low"
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw));colnames(expr) = colnames(expr_raw);rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
expr[1:5,1:5]
dim(expr)

correlation_matrix = cor(expr)
pcr = prcomp(t(correlation_matrix))

p = ggbiplot::ggbiplot(
    pcr,
    obs.scale = 1,
#    groups = meta_data$Grading,
    ellipse = TRUE,
    circle = TRUE,
    #labels = rownames(meta_data)
    var.axes = F#,
)
p
## Figure 1

#meta_data$Location[!str_detect(meta_data$Location,pattern = "Primary")] = "Metastasis"
#meta_data$Grading[meta_data$Grading == ""] = "G0"
pheatmap::pheatmap(
    correlation_matrix,
    annotation_col = meta_data[c("NEC_NET","Grading")],
    annotation_colors = aka3,
    show_rownames = F,
    show_colnames = T,
    treeheight_col = 0,
    legend = F,
    fontsize_col = 7,
    clustering_method = "average"
)

## Figure 1

# Plot 2

## Figure 2

# Plot 1

data_t = read.table("~/Deco/Results/ROC_curves.tsv",sep ="\t", header = T)

# hisc

data_t_hisc = subset(data_t, Predictor == "HISC" )
vis_mat_hisc = aggregate(
    as.integer(subset(data_t_hisc$ROC, data_t_hisc$Type == "In-silico")),
    by = list(as.character(subset(data_t_hisc$Dataset, data_t_hisc$Type == "In-silico"))),
    FUN = mean
)
vis_mat_hisc$SD = aggregate(
    as.integer(subset(data_t_hisc$ROC, data_t_hisc$Type == "In-silico")),
    by = list(as.character(subset(data_t_hisc$Dataset, data_t_hisc$Type == "In-silico"))),
    FUN = sd
)
vis_mat_hisc$Type = rep("HISC", nrow(vis_mat_hisc))
vis_mat_hisc = as.data.frame(as.matrix(vis_mat_hisc))
vis_mat_hisc = vis_mat_hisc[,-which(colnames(vis_mat_hisc)  == "SD.Group.1")]
colnames(vis_mat_hisc) = c("Dataset","ROC","SD","Type")

# ductal

data_t_ductal = subset(data_t, Predictor == "Ductal" )
vis_mat_ductal = aggregate(
    as.integer(subset(data_t_ductal$ROC, data_t_ductal$Type == "In-silico")),
    by = list(as.character(subset(data_t_ductal$Dataset, data_t_ductal$Type == "In-silico"))),
    FUN = mean
)
vis_mat_ductal$SD = aggregate(
    as.integer(subset(data_t_ductal$ROC, data_t_ductal$Type == "In-silico")),
    by = list(as.character(subset(data_t_ductal$Dataset, data_t_ductal$Type == "In-silico"))),
    FUN = sd
)
vis_mat_ductal$Type = rep("Ductal", nrow(vis_mat_ductal))
vis_mat_ductal = as.data.frame(as.matrix(vis_mat_ductal))
vis_mat_ductal = vis_mat_ductal[,-which(colnames(vis_mat_ductal)  == "SD.Group.1")]
colnames(vis_mat_ductal) = c("Dataset","ROC","SD","Type")

# MKI-67

data_t_mki67 = subset(data_t, Predictor == "MKI-67" )
vis_mat_mki67 = aggregate(
    subset(data_t_mki67$ROC, data_t_mki67$Type == "MKI-67"),
    by = list(subset(data_t_mki67$Dataset, data_t_mki67$Type == "MKI-67")),
    FUN = mean
)
vis_mat_mki67$SD = rep(0, nrow(vis_mat_mki67))
vis_mat_mki67$Type = rep("MKI_67", nrow(vis_mat_mki67))
colnames(vis_mat_mki67) = c("Dataset","ROC","SD","Type")

# merge

vis_mat = rbind(as.matrix(vis_mat_ductal),as.matrix(vis_mat_hisc),as.matrix(vis_mat_mki67))
vis_mat = as.data.frame(vis_mat)
vis_mat$ROC = as.double(as.character(vis_mat$ROC))
vis_mat$SD = as.double(as.character(vis_mat$SD))

### plots

p = ggplot( data = vis_mat,aes( x = Dataset, y = as.integer(ROC), min = ROC-SD, max = ROC+SD, fill =Type) )
p = p + geom_bar(stat="identity", position=position_dodge())
#p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#p = p + scale_fill_gradientn(colours = c("white","yellow","red"), breaks = c(0.0,.5,1.0))
#p = p + annotate("text", x=1:57,y = 5.5,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
#p = p + xlab("") + ylab("MEN1 expression in log TPM and MEN1 mutation allele frequency") + theme(legend.position = "top")
p + geom_errorbar(aes(),  position = "dodge")

# Figure 4 <- here be changes for supervised versus unsupervised

data_t = read.table("~/Deco/Results/SM_Table_2_Supervised_vs_Unsupervised_Grading_prediction.tsv",sep="\t",header = T)
data_t_sd = data_t[ !( data_t$X %in% c("Sensitivity","Specificity","Accuracy","PPV","Kappa")),]
data_t_sd = as.double(as.character(unlist((data_t_sd[,c(2,3)]))))
data_t = data_t[data_t$X %in% c("Sensitivity","Specificity","Accuracy","PPV","Kappa"),]
sd_min = as.double(as.character(unlist((data_t[,c(2,3)])))) - data_t_sd
sd_max = as.double(as.character(unlist((data_t[,c(2,3)])))) + data_t_sd
sd_max[sd_max > 100] = 100
data_t = reshape2::melt(data_t)  %>% dplyr::rename( 'Parameter' = 'X' ) %>% dplyr::rename( 'Type' = 'variable' )
data_t$Parameter = factor(data_t$Parameter, levels = c("Sensitivity","Specificity","Accuracy","PPV","Kappa"))
#data_t$SD = data_t_sd

p = ggplot( 
    data = data_t,
    aes( 
        x = Parameter,
        y = as.double(value),
        fill = Type,
        min = sd_min,
        max = sd_max
    )
)
p = p + geom_bar(stat="identity", position=position_dodge())
p = p + theme(axis.text.x = element_text(angle = 45, vjust = .5))
p = p + xlab("Performance paramter") + ylab("Performance in percent") + theme(legend.position = "top")
p = p + scale_fill_manual(values = c("red","blue"))
p = p+ geom_errorbar(aes(size = .5),  position = "dodge",size = 1)

#svg(filename = "~/Deco/Results/Images/Figure_3_classification_efficiency.svg", width = 10, height = 10)
p
#dev.off()
### survival plots ###

surv_cor = apply(expr_raw, MARGIN = 1, FUN = function(vec){return(cor(vec, as.double(meta_data$OS_Tissue)))})

subtype = as.double(expr_raw[ which( rownames(expr_raw) == "POU5F1"),df_map])
cut_off = quantile(sort(subtype), probs = seq(0,1,.1) )[9]
subtype[subtype > as.double(cut_off) ] = "Above"
subtype[subtype != "Above" ] = "Below"

#cell_type = meta_data$
cell_type[cell_type != "Not_sig"] = "Sig"
#col_vec = 
fit = survival::survfit( survival::Surv( meta_data$OS_Tissue ) ~ cell_type)

survminer::ggsurvplot(fit, data = meta_data, risk.table = T, pval = T)
# Visualize with survminer

### prediction NEC NET

mki_67 = deconvolution_results[rownames(meta_data),"MKI67"]
ductal = deconvolution_results[rownames(meta_data),"ductal"]
hisc = deconvolution_results[rownames(meta_data),"hisc"]
nec_net = meta_data$NEC_NET
target_vector = nec_net
target_vector[target_vector=="NEC"] = 0
target_vector[target_vector != 0] = 1

t_data = data.frame(
    "mki_67" = mki67,
    "nec_net" = nec_net,
    "ductal" = ductal,
    "hisc" = hisc
)

rf_fit <- glm(
    nec_net ~ hisc, data = t_data, family=binomial(link="logit")
)
predicted <- plogis(predict(rf_fit, t_data))  # predicted scores
optCutOff <- optimalCutoff(target_vector, predicted)[1] 
sensitivity = round(InformationValue::sensitivity(actuals = as.double(target_vector),predicted, threshold = optCutOff),2)
specificity = round(InformationValue::specificity(actuals = as.double(target_vector),predicted, threshold = optCutOff),2)

## Figure 1 # ESTIMATE

estimate_t = read.table("~/Deko/Results/MAPTor-NET.estimate.tsv",skip = 2,header = T)
colnames(estimate_t) = str_replace_all(colnames(estimate_t),pattern = "^X","")
rownames(estimate_t) = estimate_t$NAME
estimate_t = estimate_t[,c(-1,-2)]
estimate_t = estimate_t[-4,]

pheatmap::pheatmap(
    estimate_t,
    annotation_col = meta_data[c("Study")],
    annotation_colors = aka3,
    show_rownames = T,
    show_colnames = T,
    #treeheight_col = 0,
    #legend = F,
    #fontsize_col = 7,
    clustering_method = "average"
)

## Figure 3 Segerstolpe Heatmap

#deconvolution_results = readRDS("~/Deco/Results//TMP.RDS")

d_6 = read.table("~/Deco/Results/Segerstolpe_RepSet_6.tsv",sep =",", header = T, stringsAsFactors = F)
rownames(d_6) = d_6[,1]
d_6 = d_6[,-1]
meta_data = meta_info[rownames(d_6),]
grading = meta_data$Grading
d_6$Grading = grading

selector = c("Grading","ductal","acinar","delta","gamma","beta","alpha")
vis_mat_6 = d_6[,selector]
#vis_mat$MKI67 = round(vis_mat$MKI67 / max(vis_mat$MIK67) * 100,1)* 3
vis_mat_6$Grading = str_replace_all(vis_mat_6$Grading,pattern = "^G","")
vis_mat_6 = matrix(as.double(unlist(vis_mat_6)), ncol = length(selector))
colnames(vis_mat_6) = c("Grading","Ductal","Acinar","Delta","Gamma","Beta","Alpha")
vis_mat_6 = as.data.frame(vis_mat_6)
vis_mat_6$Grading = as.factor(vis_mat_6$Grading)

vis_mat_6 = reshape2::melt(vis_mat_6)
vis_mat_6$variable = factor(vis_mat_6$variable,levels = c("Grading","Ductal","Acinar","Delta","Gamma","Beta","Alpha"))

melt_mat_6 = vis_mat_6 %>% dplyr::group_by(Grading,variable)
melt_mat_6 = melt_mat_6 %>% dplyr::summarise( mean(value) ) %>% dplyr::rename( 'Mean_Proportion' = 'mean(value)' )
melt_mat_6 = melt_mat_6 %>% dplyr::group_by(variable) %>% dplyr::rename( 'Cell_Type_Proportion' = 'variable' )

col_vec = as.character(melt_mat_6$Cell_Type_Proportion)
col_vec[ col_vec =="Alpha"] = rep("blue", sum(col_vec =="Alpha"))
col_vec[ col_vec =="Beta"] = rep("darkgreen", sum(col_vec =="Beta"))
col_vec[ col_vec =="Gamma"] = rep("organge", sum(col_vec =="Gamma"))
col_vec[ col_vec =="Delta"] = rep("purple", sum(col_vec =="Delta"))
col_vec[ col_vec =="Ductal"] = rep("brown", sum(col_vec =="Ductal"))
col_vec[ col_vec =="Acinar"] = rep("lightblue", sum(col_vec =="Acinar"))

melt_mat_6$Cell_Type_Proportion = factor(melt_mat_6$Cell_Type_Proportion, levels = c("Acinar","Delta","Gamma","Beta","Alpha","Ductal"))

###

d_4 = read.table("~/Deco/Results/Segerstolpe_RepSet_4.tsv",sep ="\t", header = T, stringsAsFactors = F)
rownames(d_4) = d_4[,1]
d_4 = d_4[,-1]
#d_4 = read.table("~/Deco/Results/Baron_RepSet_4.tsv",sep ="\t", header = T, stringsAsFactors = F)
meta_data = meta_info[rownames(d_4),]
grading = meta_data$Grading
d_4$Grading = grading

selector = c("Grading","delta","gamma","beta","alpha")
vis_mat_4 = d_4[,selector]
vis_mat_4$Grading = str_replace_all(vis_mat_4$Grading,pattern = "^G","")
vis_mat_4 = matrix(as.double(unlist(vis_mat_4)), ncol = length(selector))
colnames(vis_mat_4) = c("Grading","Delta","Gamma","Beta","Alpha")
vis_mat_4 = as.data.frame(vis_mat_4)
vis_mat_4$Grading = as.factor(vis_mat_4$Grading)

vis_mat_4 = reshape2::melt(vis_mat_4)

melt_mat_4 = vis_mat_4 %>% dplyr::group_by(Grading,variable)
melt_mat_4 = melt_mat_4 %>% dplyr::summarise( mean(value) ) %>% dplyr::rename( 'Mean_Proportion' = 'mean(value)' )
melt_mat_4 = melt_mat_4 %>% dplyr::group_by(variable) %>% dplyr::rename( 'Cell_Type_Proportion' = 'variable' )
#mean_mat_4 = melt_mat_4 %>% dplyr::summarise( sd(value) ) %>% dplyr::rename( 'SD' = 'sd(value)' ) %>% dplyr::right_join(mean_mat) %>% dplyr::group_by(variable) %>% dplyr::rename( 'Cell_Type_Proportion' = 'variable' )

### merge g_4 and g_6

melt_mat_4$Type = rep("Endocrine",dim(melt_mat_4)[1])
melt_mat_6$Type = rep("Endocrine/Exocrine",dim(melt_mat_6)[1])
melt_mat = rbind(melt_mat_4,melt_mat_6)
melt_mat$Cell_Type_Proportion = factor(met_mat$Cell_Type_Proportion,levels = c("Alpha","Acinar","Beta","Ductal","Gamma","Delta"))

col_vec = as.character(melt_mat$Cell_Type_Proportion)
col_vec[ col_vec =="Alpha"] = rep("blue", sum(col_vec =="Alpha"))
col_vec[ col_vec =="Beta"] = rep("darkgreen", sum(col_vec =="Beta"))
col_vec[ col_vec =="Gamma"] = rep("organge", sum(col_vec =="Gamma"))
col_vec[ col_vec =="Delta"] = rep("purple", sum(col_vec =="Delta"))
col_vec[ col_vec =="Ductal"] = rep("brown", sum(col_vec =="Ductal"))
col_vec[ col_vec =="Acinar"] = rep("lightblue", sum(col_vec =="Acinar"))

#svg(filename = "~/Deco/Results/Images/Figure_3_Proportion_MKI67_versus_Grading.svg", width = 10, height = 10)

p = ggplot(
    data = melt_mat,
    aes(
        x = Grading,
        y = as.double(Mean_Proportion)
    )
)
p = p + geom_bar(
    aes(
        y = Mean_Proportion,
        x = Grading,
        fill = Cell_Type_Proportion
    ),
    data = melt_mat,
    stat="identity",
    colour="black"
)
p = p + ylab(label = "Averaged Cell-type fraction prediction per grading") + theme(legend.position="top") 
p = p + theme(axis.line = element_line(size=1, colour = "black"),
              panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank())  +
    theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          text=element_text(family="Tahoma"),
          axis.text.x=element_text(colour="black", size = 10),
          axis.text.y=element_text(colour="black", size = 10)
    )
p = p + labs(fill = "Cell-type fraction") + scale_fill_discrete(name = "Dose", labels = c("Ductal", "Acinar", "Beta","Gamma","Alpha","delta"))
p = p + scale_fill_manual(values = c("lightblue", "blue","darkgreen","purple","red","orange"))
p = p + facet_wrap( ~Type) #+ scale_fill_manual(values=c("orange", "black", "darkgreen"))
p = p + scale_x_discrete(labels=c("1" = "G1", "2" = "G2","3" = "G3"))
p

#svg(filename = "~/Deco/Results/Images/Figure_2_Cell_Type_fractions.svg", width = 10, height = 10)
p
#dev.off()

meta_data =meta_info[rownames(fractions),]
fractions = fractions[meta_data$Grading %in% c("G1","G2","G3"),]
meta_data =meta_info[rownames(fractions),]
fractions$Grading = meta_data$Grading
#write.table(fractions,"~/Deco/Results/fractions.tsv",row.names = T, quote =F , sep ="\t")

library(ggplot2)
ggplot(
    data = fractions,
    aes(
        x = Grading,
        y = ductal
    )
) + geom_point()
# Change the point size, and shape
ggplot(mtcars, aes(x=wt, y=mpg)) +
    geom_point(size=2, shape=23)

mat = fractions[,c("alpha","beta","gamma","delta","acinar","ductal")]
mat = t(mat)
correlation_matrix = cor(mat)
pcr = prcomp(t(correlation_matrix))

p = ggbiplot::ggbiplot(
    pcr,
    obs.scale = 1,
    groups = meta_data$Grading,
    ellipse = TRUE,
    circle = TRUE,
    #labels = rownames(meta_data)
    var.axes = F#,
)
p


library("umap")
umap_plot = umap::umap((correlation_matrix))
vis_data = as.data.frame(umap_plot$layout)
colnames(vis_data) = c("x","y")
dist_mat = dist((vis_data))
p = ggplot2::qplot(
    x = vis_data$x,
    y = vis_data$y,
    color = meta_data$Grading,
    geom=c("point"),
    xlab = "Umap dim 1",
    ylab = "Umap dim 2"#,
    #shape = meta_data[colnames(expr),"IFN_I"]
)
p
