#' Batch correction with NPmatch
#'
#' @param X Normalized data matrix: features/genes in rows; samples in columns
#' @param y Factor vector indicating batch for each sample
#' @param dist.method Distance metric to use for matching ('cor' or 'euclidean')
#' @param center.x Logical for whether to center gene expression by row means
#' @param center.m Logical for whether to center expression by batch means
#' @param knn k nearest neighbor to be determined for each sample in each group
#' @param sdtop Number of top variable genes to use
#' @param return.B Logical for whether returing the full set of pairs along with the corrected data matrix
#' @param use.design Logical for whether using Limma's design matrix
#' @param use.cov Logical for whether using a scaled model matrix of the full pairs
#' 
#' @examples# Load the example dataset provided and run:
#' nX <- normalize.log2CPM(X)
#' nX <- limma::normalizeQuantiles(nX)
#' cX <- NPmatch(nX)
#'
#' @export
NPmatch <- function(X,
                    y,
                    dist.method = "cor",
                    center.x = TRUE,
                    center.m = TRUE,
                    knn = 1,
                    sdtop = 2000,
                    return.B = FALSE,
                    use.design = TRUE,
                    use.cov = FALSE) {

    ## Nearest-neighbour matching for batch correction.
    ## Creates a fully paired dataset with nearest
    ## matching neighbours when pairs are missing.

    ## Guard against unstable k-NN matching on very small sample sizes
    knn <- ifelse(ncol(X) <= 3, 1, knn)

    ## Compute distance matrix for NNM-pairing
    y1 <- paste0("y=", y)
    dX <- X

    ## Reduce for speed
    if(sdtop > nrow(dX)) sdtop <- nrow(dX)
    dX <- dX[order(-apply(dX, 1, sd, na.rm = TRUE)),][1:sdtop, ]

    if (center.x) {
        dX <- dX - rowMeans(dX, na.rm = TRUE)
    }

    if (center.m) {
        ## Center per condition group (takes out pheno differences)
        mX <- tapply(1:ncol(dX), y1, function(i) rowMeans(dX[, i, drop = FALSE], na.rm = TRUE))
        mX <- do.call(cbind, mX)
        dX <- dX - mX[, y1]
    }

    if (dist.method == "cor") {
        message("[NPmatch] computing correlation matrix D...")
        ## D <- 1 - crossprod(scale(dX)) / (nrow(dX) - 1) ## faster
        D <- 1 - cor(dX)
    } else {
        message("[NPmatch] computing distance matrix D...\n")
        D <- as.matrix(stats::dist(t(dX)))
    }

    D[is.na(D)] <- 0

    ## Find neighbours
    B <- matrix(0, 0, 0)
    if (knn > 1) {
        message(paste0("[NPmatch] finding ", knn, "-nearest neighbours..."))
        bb <- apply(D, 1, function(r) tapply(r, y1, function(s) head(names(sort(s)), knn)))
        B <- do.call(rbind, lapply(bb, function(x) unlist(x)))
        colnames(B) <- unlist(mapply(rep, names(bb[[1]]), sapply(bb[[1]], length)), use.names = FALSE)
    }
    if (knn == 1 || nrow(B) != ncol(X)) {
        message("[NPmatch] finding nearest neighbours...")
        B <- t(apply(D, 1, function(r) tapply(r, y1, function(s) names(which.min(s)))))
    }

    ## Sanity check. Bail out
    if (nrow(B) != ncol(X)) {
        message("[NPmatch] WARNING. FATAL ERROR. returning uncorrected X.")
        return(X)
    }

    ## Ensure sample is always present in own group
    rownames(B) <- colnames(X)
    idx <- cbind(1:nrow(B), match(y1, colnames(B)))
    B[idx] <- rownames(B)

    ## Imputing full paired data set
    kk <- match(as.vector(B), rownames(B))
    full.y <- y1[kk]
    full.pairs <- rep(rownames(B), ncol(B))
    full.X <- X[, kk]
    dim(full.X)

    ## Remove pairing effect
    message("[NPmatch] correcting for pairing effects...")
    design <- stats::model.matrix(~full.y)
    if (use.cov == FALSE) {
        if (!use.design)
            design <- matrix(1, ncol(full.X), 1)
        full.X <- limma::removeBatchEffect(full.X, batch = full.pairs, design = design)
    } else {
        V <- model.matrix(~ 0 + full.pairs)
        if (!use.design)
            design <- matrix(1, ncol(full.X), 1)
        full.X <- limma::removeBatchEffect(full.X, covariates = scale(V), design = design)
    }

    ## Contract to original samples
    message("[NPmatch] matching result...")
    full.idx <- rownames(B)[kk]
    cX <- do.call(cbind, tapply(1:ncol(full.X), full.idx,
                                function(i) rowMeans(full.X[, i, drop = FALSE])))
    cX <- cX[, colnames(X)]

    ## Retain original row means
    cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)
    res <- cX
    if (return.B)
        res <- list(X = cX, pairings = B)

    return(res)
}
