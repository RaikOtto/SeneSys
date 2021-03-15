SCDC_qc  = function (sc.eset, ct.varname, sample, scsetname = "Single Cell", 
          ct.sub, iter.max = 1000, nu = 1e-04, epsilon = 0.01, arow = NULL, 
          qcthreshold = 0.7, generate.figure = T, ct.cell.size = NULL, 
          cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), ...) 
{
    sc.basis = SCDC_basis_2(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, 
                          sample = sample, ct.cell.size = ct.cell.size)
    M.S <- sc.basis$sum.mat[ct.sub]
    xsc <- getCPM0(exprs(sc.eset)[rownames(sc.basis$basis.mvw), 
    ])
    N.sc <- ncol(xsc)
    m.basis <- sc.basis$basis.mvw[, ct.sub]
    sigma <- sc.basis$sigma[, ct.sub]
    valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(m.basis)) == 
                                                    0) & (!is.na(M.S))
    if (sum(valid.ct) <= 1) {
        stop("Not enough valid cell type!")
    }
    message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    m.basis <- m.basis[, valid.ct]
    M.S <- M.S[valid.ct]
    sigma <- sigma[, valid.ct]
    prop.qc <- NULL
    for (i in 1:N.sc) {
        message("Begin iterative weighted estimation...")
        basis.temp <- m.basis
        xsc.temp <- xsc[, i]
        sigma.temp <- sigma
        lm.qc <- nnls::nnls(A = basis.temp, b = xsc.temp)
        delta <- lm.qc$residuals
        wt.gene <- 1/(nu + delta^2 + colSums((lm.qc$x)^2 * t(sigma.temp)))
        x.wt <- xsc.temp * sqrt(wt.gene)
        b.wt <- sweep(basis.temp, 1, sqrt(wt.gene), "*")
        lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
        prop.wt <- lm.wt$x/sum(lm.wt$x)
        delta <- lm.wt$residuals
        for (iter in 1:iter.max) {
            wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x)^2 * 
                                                     t(sigma.temp)))
            x.wt <- xsc.temp * sqrt(wt.gene)
            b.wt <- sweep(basis.temp, 1, sqrt(wt.gene), "*")
            lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
            delta.new <- lm.wt$residuals
            prop.wt.new <- lm.wt$x/sum(lm.wt$x)
            if (sum(abs(prop.wt - prop.wt.new) < epsilon)) {
                prop.wt <- prop.wt.new
                delta <- delta.new
                message("Converged at iteration ", iter)
                break
            }
            prop.wt <- prop.wt.new
            delta <- delta.new
        }
        prop.qc <- rbind(prop.qc, prop.wt)
    }
    colnames(prop.qc) <- colnames(m.basis)
    rownames(prop.qc) <- colnames(xsc)
    if (!is.null(arow)) {
        df.arow <- data.frame(sc.eset@phenoData@data[, arow])
        rownames(df.arow) <- rownames(sc.eset@phenoData@data)
        colnames(df.arow) <- arow
    }
    else {
        df.arow <- NULL
    }
    if (generate.figure) {
        heat.anno <- pheatmap::pheatmap(prop.qc, annotation_row = df.arow, 
                                        annotation_names_row = FALSE, show_rownames = F, 
                                        annotation_names_col = FALSE, cutree_rows = length(ct.sub), 
                                        color = cbPalette[2:4], cluster_rows = T, cluster_cols = F)
    }
    else {
        heat.anno <- NULL
    }
    prop.qc.keep <- rowSums(prop.qc > qcthreshold) == 1
    sc.eset.qc <- sc.eset[, prop.qc.keep]
    return(list(prop.qc = prop.qc, sc.eset.qc = sc.eset.qc, heatfig = heat.anno))
}