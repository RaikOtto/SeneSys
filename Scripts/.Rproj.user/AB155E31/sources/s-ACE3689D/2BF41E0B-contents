SCDC_prop_2 = function (
    bulk.eset,
    sc.eset,
    ct.varname,
    sample,
    ct.sub,
    iter.max = 1000, 
    nu = 1e-04,
    epsilon = 0.01,
    truep = NULL,
    weight.basis = T, 
    ct.cell.size = NULL,
    Transform_bisque = F, ...) 
{
    bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset)) > 0, , drop = FALSE]
    ct.sub <- intersect(ct.sub, unique(sc.eset@phenoData@data[, 
                                                              ct.varname]))
    sc.basis <- SCDC_basis_2(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, 
                           sample = sample, ct.cell.size = ct.cell.size)
    commongenes <- intersect(rownames(sc.basis$basis.mvw), rownames(bulk.eset))
    if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])) {
        stop("Too few common genes!")
    }
    message(paste("Used", length(commongenes), "common genes..."))
    if (weight.basis) {
        basis.mvw <- sc.basis$basis.mvw[commongenes, ct.sub]
    }
    else {
        basis.mvw <- sc.basis$basis[commongenes, ct.sub]
    }
    if (Transform_bisque) {
        GenerateSCReference <- function(sc.eset, ct.sub) {
            cell.labels <- base::factor(sc.eset[[ct.sub]])
            all.cell.types <- base::levels(cell.labels)
            aggr.fn <- function(ct.sub) {
                base::rowMeans(Biobase::exprs(sc.eset)[, cell.labels == 
                                                           ct.sub, drop = F])
            }
            template <- base::numeric(base::nrow(sc.eset))
            sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
            return(sc.ref)
        }
        sc.ref <- GenerateSCReference(sc.eset, ct.sub)[genes, 
                                                       , drop = F]
        ncount <- table(sc.eset@phenoData@data[, sample], sc.eset@phenoData@data[, 
                                                                                 ct.varname])
        true.prop <- ncount/rowSums(ncount, na.rm = T)
        sc.props <- round(true.prop[complete.cases(true.prop), 
        ], 2)
        Y.train <- sc.ref %*% t(sc.props[, colnames(sc.ref)])
        dim(Y.train)
        X.pred <- exprs(bulk.eset)[commongenes, ]
        sample.names <- base::colnames(Biobase::exprs(bulk.eset))
        template <- base::numeric(base::length(sample.names))
        base::names(template) <- sample.names
        SemisupervisedTransformBulk <- function(gene, Y.train, 
                                                X.pred) {
            Y.train.scaled <- base::scale(Y.train[gene, , drop = T])
            Y.center <- base::attr(Y.train.scaled, "scaled:center")
            Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
            n <- base::length(Y.train.scaled)
            shrink.scale <- base::sqrt(base::sum((Y.train[gene, 
                                                          , drop = T] - Y.center)^2)/n + 1)
            X.pred.scaled <- base::scale(X.pred[gene, , drop = T])
            Y.pred <- base::matrix((X.pred.scaled * shrink.scale) + 
                                       Y.center, dimnames = base::list(base::colnames(X.pred), 
                                                                       gene))
            return(Y.pred)
        }
        Y.pred <- base::matrix(base::vapply(X = commongenes, 
                                            FUN = SemisupervisedTransformBulk, FUN.VALUE = template, 
                                            Y.train, X.pred, USE.NAMES = TRUE), nrow = base::length(sample.names))
        indices <- base::apply(Y.pred, MARGIN = 2, FUN = function(column) {
            base::anyNA(column)
        })
        if (base::any(indices)) {
            if (sum(!indices) == 0) {
                base::stop("Zero genes left for decomposition.")
            }
            Y.pred <- Y.pred[, !indices, drop = F]
            sc.ref <- sc.ref[!indices, , drop = F]
        }
        results <- base::as.matrix(base::apply(Y.pred, 1, function(b) {
            sol <- lsei::pnnls(sc.ref, b, sum = 1)
            return(sol$x)
        }))
        prop.est.mvw <- t(results)
        colnames(prop.est.mvw) <- colnames(sc.ref)
        rownames(prop.est.mvw) <- colnames(bulk.eset)
        yhat <- sc.ref %*% results
        colnames(yhat) <- colnames(bulk.eset)
        yobs <- exprs(bulk.eset)
        yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
        peval <- NULL
        if (!is.null(truep)) {
            peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw, 
                                pest.names = c("SCDC"), select.ct = ct.sub)
        }
    }
    else {
        #browser()
        xbulk <- getCPM0(exprs(bulk.eset)[commongenes, ])
        rownames(sc.basis$sigma) = str_to_upper(rownames(sc.basis$basis))
        sigma <- sc.basis$sigma[str_to_upper(commongenes), ct.sub]
        ALS.S <- sc.basis$sum.mat[ct.sub]
        N.bulk <- ncol(bulk.eset)
        valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(basis.mvw)) == 
                                                        0) & (!is.na(ALS.S))
        if (sum(valid.ct) <= 1) {
            stop("Not enough valid cell type!")
        }
        message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
        basis.mvw <- basis.mvw[, valid.ct]
        ALS.S <- ALS.S[valid.ct]
        sigma <- sigma[, valid.ct]
        prop.est.mvw <- NULL
        yhat <- NULL
        yhatgene.temp <- rownames(basis.mvw)
        for (i in 1:N.bulk) {
            basis.mvw.temp <- basis.mvw
            xbulk.temp <- xbulk[, i] * 100
            sigma.temp <- sigma
            message(paste(colnames(xbulk)[i], "has common genes", 
                          sum(xbulk[, i] != 0), "..."))
            lm <- nnls::nnls(A = basis.mvw.temp, b = xbulk.temp)
            delta <- lm$residuals
            wt.gene <- 1/(nu + delta^2 + colSums((lm$x * ALS.S)^2 * 
                                                     t(sigma.temp)))
            x.wt <- xbulk.temp * sqrt(wt.gene)
            b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
            lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
            prop.wt <- lm.wt$x/sum(lm.wt$x)
            delta <- lm.wt$residuals
            for (iter in 1:iter.max) {
                wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x * 
                                                          ALS.S)^2 * t(sigma.temp)))
                x.wt <- xbulk.temp * sqrt(wt.gene)
                b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), 
                              "*")
                lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
                delta.new <- lm.wt$residuals
                prop.wt.new <- lm.wt$x/sum(lm.wt$x)
                if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
                    prop.wt <- prop.wt.new
                    delta <- delta.new
                    R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% 
                                      as.matrix(lm.wt$x))/var(xbulk.temp)
                    message("WNNLS Converged at iteration ", iter)
                    break
                }
                prop.wt <- prop.wt.new
                delta <- delta.new
            }
            R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% as.matrix(lm.wt$x))/var(xbulk.temp)
            prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
            yhat.temp <- basis.mvw.temp %*% as.matrix(lm.wt$x)
            yhatgene.temp <- intersect(rownames(yhat.temp), yhatgene.temp)
            yhat <- cbind(yhat[yhatgene.temp, ], yhat.temp[yhatgene.temp, 
            ])
        }
        colnames(prop.est.mvw) <- colnames(basis.mvw)
        rownames(prop.est.mvw) <- colnames(xbulk)
        colnames(yhat) <- colnames(xbulk)
        yobs <- exprs(bulk.eset)
        yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
        peval <- NULL
        if (!is.null(truep)) {
            peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw, 
                                pest.names = c("SCDC"), select.ct = ct.sub)
        }
    }
    return(list(prop.est.mvw = prop.est.mvw, basis.mvw = basis.mvw, 
                yhat = yhat, yeval = yeval, peval = peval))
}