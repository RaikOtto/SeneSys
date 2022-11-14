library("matrixStats")
library("circlize")
library("ComplexHeatmap")
library("data.table")
library("umap")
library("stringr")
library("ggplot2")
library("dplyr")
library("grid")

meta_info = read.table("~/SeneSys/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_info = meta_info[meta_info$Sample!="",]
rownames(meta_info) = meta_info$Sample
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
    
    # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
    
    # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)
        
        # Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            # pos: match (within the gene set)
            # neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos
            
            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
            
            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
            
            step_cdf_diff = step_cdf_pos - step_cdf_neg
            
            # Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes
            
            # Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })
    
    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
    
    # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))
    
    # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}

#####

expr_raw = read.table("~/SeneSys/Data/Reddy_hgnc_S255.tsv",sep = "\t", row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")

table(colnames(expr_raw) %in% meta_info$Sample)
meta_data = meta_info[colnames(expr_raw),]
rownames(meta_data) = meta_data$Sample

data = as.matrix(expr_raw)
data[1:5,1:5]

gene_set = read.table("~/SeneSys/GSEA/SAS_TF.csv",sep ="\t",header = TRUE)
head(gene_set)

gene_sets = as.list(as.data.frame(gene_set))

system.time(assign('res', ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)))
res1 = t(res)
head(res1)

#zscore the ssgsea output for comparative analysis
mat = (res - rowMeans(res))/(rowSds(as.matrix(res)))[row(res)]
result_t = t(mat)

res_mat = matrix(as.character(),nrow = nrow(meta_info),ncol = ncol(result_t))
colnames(res_mat) = colnames(result_t)
res_mat = data.frame(res_mat)
rownames(res_mat) = NA

matcher = match(rownames(result_t),meta_info$Sample,nomatch = 0)
dim(res_mat)
dim(meta_info)
res_mat[matcher,] = result_t
meta_info = cbind(meta_info,res_mat)

#write.table(meta_info,"~/SeneSys/Misc/Meta_information.tsv", sep = "\t", quote =F,row.names = FALSE)
