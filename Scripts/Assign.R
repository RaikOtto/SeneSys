library("ASSIGN")

# data

expr_raw = read.table("~/SeneSys/Data/GSE98588.tsv",sep ="\t", header = T,comment.char = '!',row.names = 1,stringsAsFactors = F)
dim(expr_raw)
expr_raw[1:5,1:5]

# meta data

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv", sep = "\t", header = T, stringsAsFactors = F)
meta_info = meta_info[meta_info$Sample != "",]
rownames(meta_info) = meta_info$Sample
meta_data = meta_info[colnames(expr_raw),]
expr_raw = expr_raw[,which(meta_data$Complex != "na")]
meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)
dim(meta_data)

cohort = meta_data$Complex

trainingLabel = list(
    control = c(),
    Complex = which(cohort != "clean")
)
trainingLabel$control$Complex = which(cohort == "clean")

# gene sets

genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t[,1]

geneList = list()
for (i in 1:nrow(genes_of_interest_hgnc_t)){
    hgnc_genes = genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]
    hgnc_genes = hgnc_genes[hgnc_genes != ""]
    hgnc_genes = hgnc_genes[ hgnc_genes %in% rownames(expr_raw)]
    geneList [genes_of_interest_hgnc_t[i,1] ] = list(hgnc_genes)
}
####

###

assign.wrapper(
    trainingData = NULL,
    testData = testData1,
    trainingLabel = NULL,
    #testLabel = trainingLabel,
    geneList = geneList1,
    n_sigGene = NULL,
    adaptive_B = TRUE,
    adaptive_S = TRUE,
    mixture_beta = TRUE,
    outputDir = "~/SeneSys/Scripts/tempdir/GSE98588",
    iter = 2000,
    burn_in = 1000
)
dir.create("tempdir")
tempdir <- "tempdir"

data(trainingData1)
data(testData1)
data(geneList1)

trainingLabel1 <- list(control = list(
    bcat=1:10, e2f3=1:10,
    myc=1:10, ras=1:10, src=1:10),
    bcat = 11:19, e2f3 = 20:28, myc= 29:38,
    ras = 39:48, src = 49:55
)
testLabel1 <- rep(c("Adeno", "Squamous"), c(53,58))

dir.create(file.path(tempdir,"wrapper_example1"))
assign.wrapper(trainingData=trainingData1, testData=testData1,
               trainingLabel=trainingLabel1, testLabel=testLabel1,
               geneList=NULL, n_sigGene=rep(200,5), adaptive_B=TRUE,
               adaptive_S=FALSE, mixture_beta=TRUE,
               outputDir=file.path(tempdir,"wrapper_example1"),
               iter=2000, burn_in=1000)

dir.create(file.path(tempdir,"wrapper_example2"))
assign.wrapper(trainingData=trainingData1, testData=testData1,
               trainingLabel=trainingLabel1, testLabel=NULL,
               geneList=geneList1, n_sigGene=NULL, adaptive_B=TRUE,
               adaptive_S=FALSE, mixture_beta=TRUE,
               outputDir=file.path(tempdir,"wrapper_example2"),
               iter=2000, burn_in=1000)

dir.create(file.path(tempdir,"wrapper_example3"))
assign.wrapper(
    trainingData = NULL,
    testData = trainingData1,
    trainingLabel = NULL,
    testLabel = trainingLabel1,
    geneList = geneList1,
    n_sigGene = NULL,
    adaptive_B = FALSE,
    adaptive_S = FALSE,
    mixture_beta = FALSE,
    outputDir= "/home/ottoraik/Downloads",
    iter = 2000,
    burn_in=1000
)
