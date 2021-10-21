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

path_transcriptome_file = "~/SeneSys/Data/Data_9461.Counts.HGNC.tsv"
expr_raw = read.table(path_transcriptome_file,sep="\t", stringsAsFactors =  F, header = T, row.names = 1,as.is = F)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
allowed_genes = unique(read.table("~/SeneSys/Data/Allowed_genes.tsv",sep = "\t", header = F)[,1])
expr_raw = expr_raw[rownames(expr_raw) %in% allowed_genes,]
expr_raw[1:5,1:5]
dim(expr_raw)

meta_data = meta_info[colnames(expr_raw),]
expr_raw = expr_raw[,meta_data$Drug_Treatment %in% c("NR","RP")]
meta_data = meta_info[colnames(expr_raw),]

row_var = apply(expr_raw, FUN = var, MARGIN = 1)
quantiles = as.double(quantile(row_var,seq(0,1,0.01)))
x = expr_raw[row_var > quantiles[90],]
dim(x)
x = t(x)

#lmProfile
#regLogistic
training_data = matrix(as.double(x), ncol = ncol(x), nrow=nrow(x))
rownames(training_data) = rownames(x)
colnames(training_data) = colnames(x)

training_data = as.data.frame(training_data[,which(!(is.na(training_data[1,])))])
#genes = c("RPL13","LGALS9","CALM1","CCT3","TUBB5","CD79A","RPS11","RPL8","PSAP","HDGF","HNRNPA2B1","EIF4G2","RPS9","CRIP1","CDK4","RPLP1","UBC","RPS18","RPL10","RPS15A","RPS25","RPL15","RPS5","HSPA8","ANP32E","SLC25A5","RPL19","SERINC3","PFN1","LSP1","UBB","MKNK2","CYFIP2","RACK1","RPS27A","DDX5","ATP2A3","PRPF8","TOP2A","HSP90AA1","SYK","HNRNPK","TKT","LCP1","CCT5","PABPC1","MYC","MYH9")
#training_data = training_data[,colnames(training_data) %in% genes]

y = as.character(meta_data$Drug_Treatment)
y = make.names(y)
training_data$outcome = as.factor(y)
colnames(training_data) = make.names(str_remove_all(colnames(training_data),"`"))

training_data[1:5,1:5]
dim(training_data)

#model_mars = train(y ~ ., data=training_data, method='regLogistic')
#gene_names = model_mars$finalModel$xNames
#gene_names = str_remove_all(gene_names, pattern = "`")
#fitted <- predict( model_mars$finalModel, newx = training_data)

#multiClassSummary
trControl <- trainControl(
    method = 'repeatedcv',
    number = 10,
    repeats =  3,
    search = 'grid',
    savePredictions = 'final',       
    classProbs = T,                  
    summaryFunction="twoClassSummary",
    allowParallel = TRUE)

fit.LR = caret::train(
    outcome~.,
    data = training_data,
    method="rf",
    metric= "Accuracy",
    preProc=c("center", "scale"),
    trControl=trControl,
    tuneLength = 3
)
print(fit.LR)
fit.LR$finalModel

varimp_mars <- varImp(fit.LR)
ordered_imp = data.frame(varimp_mars$importance)
imp_genes = str_replace_all(rownames(ordered_imp)[order(ordered_imp$Overall,decreasing = T)],"`","")
#write.table(imp_genes[1:50],"~/Downloads/imp_genes.tsv",sep ="\t", row.names = F, quote =F)
plot(fit.LR, main="Variable Importance with MARS",top = 20)

fitted = as.character(predict(fit.LR, as.matrix(training_data[,colnames(training_data) != "outcome"])))
confusionMatrix(reference = as.factor(y), data = as.factor(fitted), positive='NR')

meta_data$Predicted = fitted

#saveRDS(fit.LR,"~/Downloads/RandomForest.RDS")


#write.table(expr_raw,"~/SeneSys/Data/Data_9461.Counts.DeSeq2.HGNC.tsv",sep ="\t", row.names = T, quote =F)
