##--------------------------------------------##
##---SOURCE CODE FOR ANALYSES IN MANUSCRIPT---##
##--------------------------------------------##
## These functions were used in the NPM manuscript.
## All these functions are part of our playbase R package.
## To be able to run with all their functionalities and
## avoid unexpected issues, please install playbase and playdata.
## https://github.com/bigomics/playbase/tree/main

## https://github.com/bigomics/playbase/blob/main/R/pgx-stats.R
stats.numsig <- function(
                         X,
                         y,
                         lfc = 1,
                         q = 0.05,
                         set.na = NULL,
                         trend = TRUE,
                         verbose = TRUE,
                         gs.method = "rankcor"
                         ) {

  if (!is.null(set.na)) {
    y[y == set.na] <- NA
  }

  ## select non-missing
  sel <- !is.na(y)
  y <- y[sel]
  X <- X[, sel, drop = FALSE]

  ## Genes
  res <- gx.limmaF(X, y, fdr = 1, lfc = 0, trend = trend, verbose = 0)
  res <- res[order(res$P.Value), ]
  avx <- res[, grep("AveExpr", colnames(res))]
  ldiff <- apply(avx, 1, function(x) diff(range(x)))
  sig <- (abs(ldiff) > lfc & res$adj.P.Val < q)
  sig.genes <- rownames(res)[which(sig)]
  fc0 <- array(ldiff, dimnames = list(rownames(res)))
  names(fc0) <- toupper(names(fc0))

  ## Gene sets
  is.symbol <- (mean(names(fc0) %in% colnames(playdata::GSETxGENE)) > 0.2)
  sig.gsets <- c()
  gsa <- list(rho = NULL, p.value = NULL, q.value = NULL)
  if (is.symbol) {
    sel <- grep("^go|pathway|geo_human", rownames(playdata::GSETxGENE),
      value = TRUE, ignore.case = TRUE
    )
    gmt <- playdata::GSETxGENE[sel, ]

    if (gs.method == "rankcor") {
      gsa <- gset.rankcor(cbind(fc0), Matrix::t(gmt), compute.p = TRUE)
      gsa <- data.frame(rho = gsa$rho[, 1], p.value = gsa$p.value[, 1], q.value = gsa$q.value[, 1])
    } else if (gs.method == "fgsea") {
      gsets <- mat2gmt(t(gmt))
      gsa <- fgsea::fgsea(gsets, fc0)
      gsa <- data.frame(rho = gsa$NES, p.value = gsa$pval, q.value = gsa$padj)
    } else {
      stop("unknown geneset method: ", gs.method)
    }
    gsa <- gsa[order(gsa$p.value), ]
    sig.gs <- (gsa$q.value < q)
    sig.gsets <- rownames(gsa)[which(sig.gs)]
  }

  if (verbose) {
    cat("nsig.genes = ", length(sig.genes), "\n")
    cat("nsig.gsets = ", length(sig.gsets), "\n")
  }

  list(
    genes = sig.genes,
    gsets = sig.gsets,
    fc = fc0,
    p.genes = res$P.Value,
    p.gsets = gsa$p.value
  )
}

## https://github.com/bigomics/playbase/blob/main/R/pgx-correct.R
## From playbase::runBatchCorrectionMethods
runBatchCorrectionMethods <- function (X,
                                       batch,
                                       y,
                                       controls = NULL,
                                       ntop = 2000,
                                       sc = FALSE, 
                                       prefix = "",
                                       methods = NULL,
                                       remove.failed = TRUE
                                       ) {
  mod <- model.matrix(~factor(y))
  nlevel <- length(unique(y[!is.na(y)]))

  if (ntop < Inf) {
    X <- head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], ntop)
  }

  if (is.null(methods)) {
    methods <- c("uncorrected", "normalized_to_control", 
      "ComBat", "limma", "ComBat.no_mod", "limma.no_mod", 
      "superBC", "PCA", "RUV", "SVA", "NPM", "MNN", "Harmony")
  }

  xlist <- list()

  if ("uncorrected" %in% methods) {
    xlist[["uncorrected"]] <- X
  }

  ## --------------------------------------------------------------
  ## SUPERVISED METHODS (need batch parameter)
  ## --------------------------------------------------------------

  ## normalize to control
  if (!is.null(controls) && "normalized_to_control" %in% methods) {
    nX <- normalizeToControls(X, batch, y, controls)
    xlist[["normalized_to_control"]] <- nX
  }

  ## limma
  if ("limma" %in% methods && is.null(batch)) {
    xlist[["limma"]] <- X
  }
  if ("limma" %in% methods && !is.null(batch)) {
    cX <- try(limmaCorrect(X, batch, y = y))
    if ("try-error" %in% class(cX)) {
      ## if fails with model, try without
      cX <- try(limmaCorrect(X, batch, y = NULL))
    }
    xlist[["limma"]] <- cX
  }
  if ("limma.no_mod" %in% methods && is.null(batch)) {
    xlist[["limma.no_mod"]] <- X
  }
  if ("limma.no_mod" %in% methods && !is.null(batch)) {
    cX <- try(limmaCorrect(X, batch, y = NULL))
    xlist[["limma.no_mod"]] <- cX
  }

  ## ComBat
  if ("ComBat" %in% methods && is.null(batch)) {
    xlist[["ComBat"]] <- X
  }
  if ("ComBat" %in% methods && !is.null(batch)) {
    if (max(table(batch), na.rm = TRUE) > 1) {
      bX <- try(combatCorrect(X, batch, y = y))
      if ("try-error" %in% class(bX)) {
        ## if not successful, try without model
        bX <- try(combatCorrect(X, batch, y = NULL))
      }
      xlist[["ComBat"]] <- bX
    }
  }
  if ("ComBat.no_mod" %in% methods && !is.null(batch)) {
    bX <- try(combatCorrect(X, batch, y = NULL))
    xlist[["ComBat.no_mod"]] <- bX
  }
  if ("ComBat.no_mod" %in% methods && is.null(batch)) {
    xlist[["ComBat.no_mod"]] <- X
  }

  ## superbatchcorrect
  if ("superBC" %in% methods) {
    df <- data.frame(y = y)
    if (!is.null(batch)) 
      df <- cbind(df, batch = batch)
    xlist[["superBC"]] <- pgx.superBatchCorrect(X, df, model.par = "y", 
      batch.par = "*")$X
  }

  ## --------------------------------------------------------------
  ## UNSUPERVISED METHODS (need pheno vector)
  ## --------------------------------------------------------------

  ## PCA
  if ("PCA" %in% methods) {
    xlist[["PCA"]] <- try(pcaCorrect(X, y = y, p.notsig = 0.2))
  }

  ## RUV
  if ("RUV" %in% methods) {
    xlist[["RUV"]] <- try(ruvCorrect(X, y, k = NULL, type = "III"))
  }

  ## SVA
  if ("SVA" %in% methods) {
    xlist[["SVA"]] <- try(svaCorrect(X, y))
  }

  ## NPM
  if ("NPM" %in% methods) {
    xlist[["NPM"]] <- nnmCorrect(X, y, use.design = TRUE)
  }

  ## single-cell data specific
  if (sc) {
    if ("MNN" %in% methods) {
      xlist[["MNN"]] <- try(MNNcorrect(X, batch))
      if (!is.null(controls)) {
        kk <- (y %in% controls)
        xlist[["rMNN"]] <- try(MNNcorrect(X, batch, controls = kk))
      }
    }
    if ("Harmony" %in% methods) {
      res <- try(runHarmony(X, batch = batch))
      if (!"try-error" %in% class(res)) {
        xlist[["Harmony"]] <- as.matrix(res$corrected)
      }
    }
  }

  if (remove.failed) {
    is.error <- sapply(xlist, function(x) ("try-error" %in% class(x)))
    is.nullrow <- sapply(sapply(xlist, nrow), is.null)
    is.xnull <- sapply(xlist, is.null)
    xlist <- xlist[which(!is.xnull & !is.nullrow & !is.error)]
  }

  names(xlist) <- paste0(prefix, names(xlist))

  return(xlist)

}

## https://github.com/bigomics/playbase/blob/main/R/pgx-correct.R
## From playbase::bc.evaluateResults
bc.evaluateResults <- function(xlist,
                               pheno,
                               lfc = 0.2,
                               q = 0.2,
                               pos = NULL,
                               add.sil = TRUE,
                               plot = TRUE,
                               trend = TRUE,
                               ref = "uncorrected",
                               clust = "tsne"
                               ) {

  if (!ref %in% names(xlist)) {
    ref <- names(xlist)[1]
  }

  ## compute and make table
  message("[bc.evaluateResults] computing statistics...")
  numsig <- lapply(xlist, stats.numsig,
    y = pheno, lfc = lfc, q = q,
    trend = trend, verbose = FALSE
  )

  res <- t(sapply(numsig, function(r) {
    c(sapply(r[1:2], length), avg.fc = mean(abs(r[[3]]), na.rm = TRUE))
  }))

  sdx <- sapply(xlist, function(x) mean(matrixStats::rowSds(x, na.rm = TRUE)))
  snr <- res[, "avg.fc"] / sdx
  res <- cbind(res, avg.sd = sdx, SNR = snr)

  ## compute relative genes/geneset overlap
  message("[bc.evaluateResults] computing overlap...")
  g1 <- numsig[[ref]]$genes
  n1 <- sapply(numsig, function(s) length(intersect(s$genes, g1)))
  n2 <- sapply(numsig, function(s) length(g1))
  r.genes <- n1 / (1e-3 + n2)
  res <- cbind(res, r.genes)

  any.gsets <- any(sapply(numsig, function(s) length(s$gsets) > 0))
  if (any.gsets) {
    s1 <- numsig[[ref]]$gsets
    m1 <- sapply(numsig, function(s) length(intersect(s$gsets, s1)))
    m2 <- sapply(numsig, function(s) length(s1))
    r.gsets <- m1 / (1e-3 + m2)
    res <- cbind(res, r.gsets)
  }

  ## centered top
  xlist1 <- lapply(xlist, function(x) {
    x <- head(x[order(-matrixStats::rowSds(x, na.rm = TRUE)), ], 1000)
    x <- as.matrix(x)
    (x - rowMeans(x, na.rm = TRUE))
  })

  message("[bc.evaluateResults] computing silhouette scores...")
  silhouette <- rep(1, nrow(res))
  if (add.sil) {
    if (is.null(pos)) {
      if (clust == "tsne") {
        nb <- max(1, min(30, round(ncol(xlist[[1]]) / 5)))
        CLUSTFUN <- function(x) {Rtsne::Rtsne(scale(t(x)),
            check_duplicates = FALSE, perplexity = nb)$Y
        }
      } else {
        CLUSTFUN <- function(x) svd(scale(t(x), scale = FALSE))$u[, 1:2]
      }
      pos <- lapply(xlist1, function(x) CLUSTFUN(x))
    }
    pheno0 <- as.character(pheno)
    pheno0[is.na(pheno0)] <- "NA"
    silhouette <- sapply(pos, function(p) {
      score <- cluster::silhouette(as.integer(factor(pheno0)), stats::dist(p))
      mean(score[, "sil_width"], na.rm = TRUE)
    })
    silhouette <- pmax(silhouette, 1e-4)

    ## PCA score
    nu <- max(2, min(10, dim(xlist[[1]]) / 4))
    pca10 <- lapply(xlist1, function(x) {
      svd(t(x), nu = nu, nv = 0)$u
    })
    Y <- model.matrix(~pheno)[, -1]
    rho <- lapply(pca10, function(x) cor(x, Y))
    rho <- lapply(rho, function(x) rowMeans(abs(x), na.rm = TRUE))
    pc1.ratio <- sapply(rho, function(r) abs(r[1]) / sum(abs(r)))
    res <- cbind(res, silhouette, pc1.ratio)
  }

  sel <- c("genes", "SNR", "silhouette")
  sel <- intersect(sel, colnames(res))

  ##  score <- res.score * (silhouette / silhouette[1])**1
  overall.score <- t(t(1e-4 + res[, sel]) / (1e-4 + res[ref, sel]))
  overall.score[, "silhouette"] <- overall.score[, "silhouette"]**2 ## give more weight
  overall.score <- exp(rowMeans(log(overall.score), na.rm = TRUE)) ## geometric mean

  res1 <- cbind(score = overall.score, res)
  res1 <- res1[order(-res1[, "score"]), ]
  pos <- pos[rownames(res1)]

  if (plot) {
    nc <- ceiling(1.2 * sqrt(length(pos)))
    nr <- ceiling(length(pos) / nc)
    i <- 1
    xdim <- nrow(pos[[1]])
    cex1 <- cut(xdim, breaks = c(0, 20, 100, 400, 1000, 999999), c(1.8, 1.5, 1.2, 0.9, 0.6))
    cex1 <- as.numeric(as.character(cex1))

    par(mfrow = c(nr, nc))
    for (i in 1:length(pos)) {
      plot(pos[[i]][, 1:2],
        col = factor(pheno), pch = 20, cex = cex1,
        main = names(pos)[i], cex.main = 1.6
      )
      tt <- paste("score = ", round(res1[i, "score"], 3))
      legend("topright", legend = tt, cex = 1.1)
    }
  }

  p.genes <- lapply(numsig, function(s) s$p.genes)
  p.gsets <- lapply(numsig, function(s) s$p.gsets)

  return(list(scores = res1, pos = pos, p.genes = p.genes, p.gsets = p.gsets))

}

## Batch correction with NPM
## https://github.com/bigomics/NPM
NPM <- function(X,
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

  ## Compute distance matrix for NNM-pairing
  y1 <- paste0("y=", y)
  dX <- X

  ## Reduce for speed
  if(sdtop > nrow(dX)) sdtop <- nrow(dX)
  dX <- dX[order(-apply(dX, 1, sd)),][1:sdtop, ]

  if (center.x) {
    dX <- dX - rowMeans(dX, na.rm = TRUE)
  }
  
  if (center.m) {
    ## Center per condition group (takes out pheno differences)
    mX <- tapply(1:ncol(dX), y1, function(i) rowMeans(dX[, i, drop = FALSE]))
    mX <- do.call(cbind, mX)
    dX <- dX - mX[, y1]
  }

  if (dist.method == "cor") {
    message("[NPM] computing correlation matrix D...")
    ## D <- 1 - crossprod(scale(dX)) / (nrow(dX) - 1) ## faster
    D <- 1 - cor(dX)
  } else {
    message("[NPM] computing distance matrix D...\n")
    D <- as.matrix(stats::dist(t(dX)))
  }

  D[is.na(D)] <- 0

  ## Find neighbours
  if (knn > 1) {
    message(paste0("[NPM] finding ", knn, "-nearest neighbours..."))
    bb <- apply(D, 1, function(r) tapply(r, y1, function(s) head(names(sort(s)), knn)))
    B <- do.call(rbind, lapply(bb, function(x) unlist(x)))
    colnames(B) <- unlist(mapply(rep, names(bb[[1]]), sapply(bb[[1]], length)), use.names = FALSE)
  } else {
    message("[NPM] finding nearest neighbours...")
    B <- t(apply(D, 1, function(r) tapply(r, y1, function(s) names(which.min(s)))))
  }
  rownames(B) <- colnames(X)

  ## Ensure sample is always present in own group
  idx <- cbind(1:nrow(B), match(y1, colnames(B)))
  B[idx] <- rownames(B)

  ## Imputing full paired data set
  kk <- match(as.vector(B), rownames(B))
  full.y <- y1[kk]
  full.pairs <- rep(rownames(B), ncol(B))
  full.X <- X[, kk]
  dim(full.X)

  ## Remove pairing effect
  message("[NPM] correcting for pairing effects...")
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
  message("[NPM] matching result...")
  full.idx <- rownames(B)[kk]
  cX <- do.call(cbind, tapply(1:ncol(full.X), full.idx,
    function(i) rowMeans(full.X[, i, drop = FALSE])))
  cX <- cX[, colnames(X)]

  ## Retain original row means
  cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)
  res <- cX
  if (return.B) {
    res <- list(X = cX, pairings = B)
  }

  return(res)
}



##------------------------------##
##-----------RUN----------------##
##------------------------------##
DIR <- "/path/to/datasets/"

## Each dataset folder has:
## 1. counts.csv (gene expression matrix): features in rows, samples in columns
## 2. samples.csv (metadata matrix): samples in rows, metadata in columns.
## ps: samples must have at least 2 columns: pheno.var (phenotype) and batch.var (batch variable)
## ps: colnames(counts) must be same as rownames(sample).

IDS <- c(
  "GSE120099",
  "GSE82177",
  "GSE162760",
  "GSE171343", 
  "GSE153380",
  "GSE163214",
  "GSE182440",
  "GSE163857",
  "GSE117970_BC",
  "GSE173078",
  "GSE10846"
)

Methods <- c("uncorrected", "NPM", "SVA", "RUV", "PCA", "limma", "ComBat")

RES <- list() ## results all go here

i = 1;
for(i in 1:length(IDS)) {

  message("Processing ", IDS[i], "...")

  setwd(paste0(DIR,IDS[i]))

  files <- playbase::read_files('.')
  counts <- files$counts 
  samples <- files$samples

  ## Normalize data
  X <- playbase::logCPM(counts[rowSums(counts) != 0, ])
  X <- limma::normalizeQuantiles(X)

  ## Check counts and samples are aligned
  cc <- all.equal(colnames(X), rownames(samples))
  if(!cc) {
    stop("samples (columns) in X and metadata matrix do not match")
  }

  ## Take phenotype and batch vectors
  pheno.var <- "Status"
  batch.var <- "Batch"
  y <- samples[, pheno.var]
  Batch <- samples[, batch.var]

  ## Perform batch correction
  xlist <- runBatchCorrectionMethods(
    X = X,
    batch = Batch,
    y = y,
    methods = Methods,
    ntop = nrow(X)
  )

  ## tSNE
  ## nb <- max(1, min(30, round(ncol(X) / 5)))
  ## pos.tsne <- lapply(xlist, function(x) Rtsne::Rtsne(t(x), perplexity = nb)$Y)
  ## names(pos.tsne) <-  paste0(IDS[i], "--", names(pos.tsne))
  
  ## Evaluate batch correction; compute score
  res <- bc.evaluateResults(
    xlist = xlist,
    pheno = y,
    lfc = 0.5,
    q = 0.05,
    clust = "tsne",
    plot = FALSE
  )

  RES[[IDS[i]]] <- res$scores
  
}
##------------------------------##
##------------------------------##
##------------------------------##
