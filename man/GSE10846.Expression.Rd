\name{GSE10846.Expression}
\docType{data}
\alias{GSE10846.Expression}
\title{Gene Expression Data from GSE10846 DLBCL Study}
\description{
    Gene expression profiling of DLBCL patient samples was performed to investigate, whether molecular gene expression signatures retain their prognostic significance in patients treated with chemotherapy plus Rituximab. The lymphnode, germinal center signature and a new angiogenesis signature were combined to a final multivariate model which defined quartile groups among Rituximab-CHOP-treated patients with distinct 3-year overall survival rates.
}
\usage{data(GSE10846.Expression)}
\format{
  A matrix with genes on the rows and samples on the columns.
  Row names are gene symbols and column names are sample identifiers (GSM IDs).
}
\source{
    Gene Expression Omnibus (GEO) accession GSE10846.
    \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846}
}
\references{
    Lenz G, Wright G, Dave SS, Xiao W et al. Stromal gene signatures in large-B-cell lymphomas. N Engl J Med 2008 Nov 27;359(22):2313-23
}
\seealso{
    \code{\link{GSE10846.Metadata}}
}
\examples{
    data(GSE10846.Expression)
    # Examine the structure of the expression data
    dim(GSE10846.Expression)
    head(GSE10846.Expression[, 1:5])
}
\keyword{datasets} 