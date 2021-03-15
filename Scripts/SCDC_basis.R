SCDC_basis_2 = function (x, ct.sub = NULL, ct.varname, sample, ct.cell.size = NULL) 
{
    if (is.null(ct.sub)) {
        ct.sub <- unique(x@phenoData@data[, ct.varname])
    }
    ct.sub <- ct.sub[!is.na(ct.sub)]
    x.sub <- x[, x@phenoData@data[, ct.varname] %in% ct.sub]
    x.sub <- x.sub[rowSums(exprs(x.sub)) > 0, ]
    countmat <- exprs(x.sub)
    ct.id <- droplevels(as.factor(x.sub@phenoData@data[, ct.varname]))
    sample.id <- as.character(x.sub@phenoData@data[, sample])
    ct_sample.id <- paste(ct.id, sample.id, sep = "%")
    mean.mat <- sapply(unique(ct_sample.id), function(id) {
        y = as.matrix(countmat[, ct_sample.id %in% id])
        apply(y, 1, sum, na.rm = TRUE)/sum(y)
    })
    mean.id <- do.call("rbind", strsplit(unique(ct_sample.id), 
                                         split = "%"))
    sigma <- sapply(unique(mean.id[, 1]), function(id) {
        y = mean.mat[, mean.id[, 1] %in% id]
        if (is.null(dim(y))) {
            res = rep(0, length(y))
            message("Warning: the cell type [", id, "] is only available in at most 1 subject!")
        }
        else {
            res = apply(y, 1, var, na.rm = TRUE)
        }
        return(res)
    })
    sum.mat2 <- sapply(unique(sample.id), function(sid) {
        sapply(unique(ct.id), function(id) {
            y = as.matrix(countmat[, ct.id %in% id & sample.id %in% 
                                       sid])
            sum(y)/ncol(y)
        })
    })
    rownames(sum.mat2) <- unique(ct.id)
    colnames(sum.mat2) <- unique(sample.id)
    if (is.null(ct.cell.size)) {
        sum.mat <- rowMeans(sum.mat2, na.rm = T)
    }
    else {
        if (is.null(names(ct.cell.size))) {
            message("Cell size factor vector requires cell type names...")
            break
        }
        else {
            sum.mat <- ct.cell.size
        }
    }
    basis <- sapply(unique(mean.id[, 1]), function(id) {
        z <- sum.mat[mean.id[, 1]]
        mean.mat.z <- t(t(mean.mat) * z)
        y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id])
        apply(y, 1, mean, na.rm = TRUE)
    })
    my.max <- function(x, ...) {
        y <- apply(x, 1, max, na.rm = TRUE)
        y/median(y, na.rm = T)
    }
    var.adj <- sapply(unique(sample.id), function(sid) {
        my.max(sapply(unique(ct.id), function(id) {
            y = countmat[, ct.id %in% id & sample.id %in% sid, 
                         drop = FALSE]
            apply(y, 1, var, na.rm = T)
        }), na.rm = T)
    })
    colnames(var.adj) <- unique(sample.id)
    q15 <- apply(var.adj, 2, function(zz) {
        z1 = min(zz[zz > 0])
        z2 = quantile(zz, 0.15, na.rm = T)
        return(max(z1, z2))
    })
    q85 <- apply(var.adj, 2, quantile, probs = 0.85, na.rm = T)
    var.adj.q <- t(apply(var.adj, 1, function(y) {
        y[y < q15] <- q15[y < q15]
        y[y > q85] <- q85[y > q85]
        return(y)
    }))
    message("Creating Basis Matrix adjusted for maximal variance weight")
    mean.mat.mvw <- sapply(unique(ct_sample.id), function(id) {
        sid = unlist(strsplit(id, "%"))[2]
        y = as.matrix(countmat[, ct_sample.id %in% id])
        browser()
        yy = sweep(y, 1, sqrt(var.adj.q[1,sid ]), "/")
        apply(yy, 1, sum, na.rm = TRUE)/sum(yy)
    })
    basis.mvw <- sapply(unique(mean.id[, 1]), function(id) {
        z <- sum.mat[mean.id[, 1]]
        mean.mat.z <- t(t(mean.mat.mvw) * z)
        y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id])
        apply(y, 1, mean, na.rm = TRUE)
    })
    basis.mvw <- basis.mvw[, ct.sub]
    sigma <- sigma[, ct.sub]
    basis <- basis[, ct.sub]
    sum.mat <- sum.mat[ct.sub]
    return(list(basis = basis, sum.mat = sum.mat, sigma = sigma, 
                basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}