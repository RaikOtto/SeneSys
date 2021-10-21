library("stringr")
library("caret")
library("skimr")
library("randomForest")

opt.cut = function( perf, pred ){
    cut.ind = mapply( FUN = function( x, y, p ){
        d = (x - 0)^2 + ( y - 1 )^2
        ind = which(d == min(d))
        c(
            sensitivity = y[[ ind ]],
            specificity = 1 - x[[ ind ]], 
            cutoff = p[[ ind ]])
    }, perf@x.values, perf@y.values, pred@cutoffs)
}

# split into training and testing
set.seed(23489)

draw_colnames_45 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
    return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

meta_info = read.table("~/SeneSys//Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

###

model = readRDS("~/Downloads/rpar.RDS")
model$finalModel$xNames = str_remove_all(model$finalModel$xNames, pattern = "`")
#rownames(model$finalModel$beta) = str_remove_all(rownames(model$finalModel$beta), pattern = "`")
#gene_names = rownames(model$finalModel$beta)
#colnames(model$finalModel$W) = str_remove_all(colnames(model$finalModel$W), pattern = "`")
#gene_names = model$finalModel$xNames
#gene_names = str_remove_all(gene_names, pattern = "`")

path_transcriptome_file = "~/SeneSys/Data/Schmitz.HGNC.DESeq2.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

meta_data = meta_info[colnames(expr_raw),]

x = t(expr_raw)

training_data = matrix(as.double(x), ncol = ncol(x), nrow=nrow(x))
rownames(training_data) = rownames(x)
colnames(training_data) = colnames(x)

#training_data = training_data[,colnames(training_data) %in% gene_names]
training_data[1:5,1:5]
#length(gene_names)
dim(training_data)

#model$finalModel$W = model$finalModel$W[,gene_names %in% colnames(training_data)]
#model$trainingData = model$trainingData[ ,colnames(model$trainingData) %in%colnames(training_data) ]

model$finalModel$xNames = model$finalModel$xNames[model$finalModel$xNames %in% colnames(training_data)]
#model$finalModel$beta = model$finalModel$beta[rownames(model$finalModel$beta) %in% colnames(training_data),]
model$finalModel$dim[1] = length(model$finalModel$xNames)
training_data = training_data[,model$finalModel$xNames]

fitted <- predict( model, data.frame(training_data))
fitted <- predict( model$finalModel, newx = training_data, 'prob')
meta_data$Predictions = fitted$predictions

table(meta_data$Predictions, meta_data$Progression)
