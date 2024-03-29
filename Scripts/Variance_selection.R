## hgnc
#library("org.Mm.eg.db")
library("hgu133plus2.db")
x <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])
rownames(xx) = xx$probe_id
hgnc_list = as.character(xx[rownames(expr_raw),2])
### hgnc

#expr_raw = bam_counts
#gene_ids = as.character(row.names(expr_raw))

hgnc_list = str_to_upper(as.character(  mapIds(
    org.Mm.eg.db,
    keys = gene_ids,
    column="SYMBOL",
    keytype="ENSEMBL",
    multiVals="first"
) ))

#hgnc_list = as.character(  mapIds(
#    org.Hs.eg.db,
#    keys = gene_ids,
#    column="SYMBOL",
#    keytype="ENSEMBL",
#    multiVals="first"
#) )

expr_raw = expr_raw[ ! is.na(hgnc_list),]
expr_raw = expr_raw[ hgnc_list != "",]
hgnc_list = hgnc_list[!is.na(hgnc_list)]
hgnc_list = hgnc_list[hgnc_list != ""]
hgnc_list = str_to_upper(hgnc_list)
hgnc_list_uni = as.character( unique(hgnc_list) )

dim(expr_raw)
length(hgnc_list)
length(hgnc_list_uni)

max_list = as.integer( sapply( hgnc_list_uni, FUN = function(gene){
    
    var_match = which( hgnc_list == gene )

    if (length(var_match) > 1){
        row_var_max = which.max(
            apply(
                as.matrix(
                    expr_raw[
                        var_match,
                    ]
                ),
                FUN = var,
                MARGIN = 1
            )
        )
        #if( length(which.max( apply( as.matrix(expr_raw[var_match,]), FUN = var, MARGIN = 1) )) == 0){
        #    print(gene)
        #}
    } else{
        row_var_max = 1
    }
    
    if( length(var_match[row_var_max])==0){
        return(0)
    } else {
        return( var_match[row_var_max] )    
    }
}))

which(is.na(hgnc_list[max_list] ))

expr_raw = expr_raw[max_list,]
rownames(expr_raw) = hgnc_list[max_list]
dim(expr_raw)

colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "" )
rownames(expr_raw) = str_to_upper(rownames(expr_raw))
expr_raw[1:5,1:5]

which(sapply(expr_raw[,1], FUN = is.na))
expr_raw = expr_raw[!(sapply(expr_raw[,1], FUN = is.na)),]
dim(expr_raw)

#write.table(expr_raw,"~/MAPTor_NET/BAMs_new/BON-1.illumina.tsv",sep="\t",quote =F)

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#listEnsembl()
#datasets <- listDatasets(ensembl)
head(datasets)
searchDatasets(mart = ensembl, pattern = "hsapiens")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

ensembl <- useEnsembl(
    biomart = "ensembl", 
    dataset = "hsapiens_gene_ensembl", 
    mirror = "useast")

genes_of_interest_hgnc_t = read.table("~/SeneSys/Misc/SeneSys_gene_sets.tsv",sep ="\t", stringsAsFactors = F, header = F)
converted_t = genes_of_interest_hgnc_t

for (i  in 1:nrow(genes_of_interest_hgnc_t)){

    print(i)    
    query_input = genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]
    query_input = query_input[query_input != ""]
    query_input = unique(query_input)
    conversion_result = getBM(attributes = c(
        'hgnc_symbol',
        'entrezgene_id'),
        filters = 'hgnc_symbol', 
        values = query_input, 
        mart = ensembl)
    
    res_vec = conversion_result[,2]
    res_vec = res_vec[!is.na(res_vec)]
    length_res_vec = length(res_vec)
    
    converted_t[i,3:ncol(genes_of_interest_hgnc_t)] = rep("",ncol(genes_of_interest_hgnc_t)-2)
    converted_t[i,3:(length_res_vec+2)] = as.character(res_vec)
}

converted_t[1:5,1:5]
write.table(converted_t,"~/Downloads/Senescence_genesets_entrez_gene_id.tsv",sep ="\t",quote = F, row.names = F,col.names = F
        )
