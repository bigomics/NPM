#' Normalize Data
#'
#' @param X raw data matrix (i.e., gene expression dataset): rows are features/genes; columns are samples.
#' @param total total count within each sample
#' @param prior offset integer added prior log transformation
#'
#' @return A matrix of log2-transformed counts per million (CPM) normalized expression data.
#'   The returned matrix has the same dimensions as the input matrix X, with genes as rows
#'   and samples as columns. Values are log2(CPM + prior) where CPM is calculated by
#'   normalizing each sample to the target total count (default 1e6).
#' @export
#' @import stats
#' @import utils
#'
#' @examples
#' # Load the example dataset provided and run:
#' data(GSE10846.Expression, package = "NPM")
#' nX <- normalize.log2CPM(GSE10846.Expression)
normalize.log2CPM <- function(
    X,
    total = 1e6,
    prior = 1
) {

    ## If total counts per sample is < 1e6 uses average count, else use 1e6.
    total0 <- mean(Matrix::colSums(X, na.rm = TRUE))
    total <- ifelse(total0 < 1e6, total0, 1e6)

    message("[log2CPM] setting column sums to = ", round(total, 2))
    totcounts <- Matrix::colSums(X, na.rm = TRUE)
    cpm <- sweep(X, 2, totcounts, FUN = "/") * total
    nX <- log2(prior + cpm)
    return(nX)

}
