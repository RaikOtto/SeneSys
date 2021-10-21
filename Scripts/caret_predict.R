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

path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.DeSeq2.HGNC.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]

meta_data = meta_info[colnames(expr_raw),]

row_var = apply(expr_raw, FUN = var, MARGIN = 1)
quantiles = as.double(quantile(row_var,seq(0,1,0.01)))
x = expr_raw[row_var > quantiles[96],]
dim(x)

y = as.factor(meta_data$ABC_GCB)
x = t(x)

preProcess_missingdata_model <- preProcess(x, method='knnImpute')
preProcess_missingdata_model

trainData <- predict(preProcess_missingdata_model, newdata = x)

preProcess_range_model <- preProcess(trainData, method='range')
trainData <- predict(preProcess_range_model, newdata = trainData)

# Append the Y variable

apply(trainData[, 1:10], 2, FUN=function(x){c('min'=min(x), 'max'=max(x))})
#trainData$Purchase = y

featurePlot(
    x = trainData[, 1:9], 
    y = as.factor(y), 
    plot = "box",
    strip=strip.custom(par.strip.text=list(cex=.7)),
    scales = list(
        x = list(relation="free"), 
        y = list(relation="free")
    )
)

set.seed(100)
options(warn=-1)

subsets <- c(1:5, 10, 15, 18)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

#lmProfile <- rfe(
#    x=trainData[, 1:18], y=as.factor(y),
#    sizes = subsets,
#    rfeControl = ctrl)

#lmProfile
#regLogistic
training_data = matrix(as.double(x), ncol = ncol(x), nrow=nrow(x))
rownames(training_data) = rownames(x)
colnames(training_data) = colnames(x)

training_data = as.data.frame(training_data[,which(!(is.na(training_data[1,])))])

training_data[1:5,1:5]
dim(training_data)
#training_data$y = factor(y)

model_mars = train(factor(y) ~ ., data=training_data, method='regLogistic')
gene_names = model_mars$finalModel$xNames
gene_names = str_remove_all(gene_names, pattern = "`")
fitted <- predict( model_mars$finalModel, newx = training_data)

varimp_mars <- varImp(model_mars)
plot(varimp_mars, main="Variable Importance with MARS")

confusionMatrix(reference = as.factor(y), data = fitted, mode='everything', positive='MM')

# Define the training control
# multiClassSummary
fitControl <- trainControl(
    method = 'cv',                   # k-fold cross validation
    number = 5,                      # number of folds
    savePredictions = 'final',       # saves predictions for optimal tuning parameter
    classProbs = T,                  # should class probabilities be returned
    summaryFunction=multiClassSummary  # results summary function
) 

#saveRDS(model_mars,"~/Downloads/mars_large.RDS")

