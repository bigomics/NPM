\name{GSE10846.Metadata}
\docType{data}
\alias{GSE10846.Metadata}
\title{Clinical Metadata for GSE10846 DLBCL Study}
\description{
    Metadata: Gene expression profiling of DLBCL patient samples was performed to investigate, whether molecular gene expression signatures retain their prognostic significance in patients treated with chemotherapy plus Rituximab. The lymphnode, germinal center signature and a new angiogenesis signature were combined to a final multivariate model which defined quartile groups among Rituximab-CHOP-treated patients with distinct 3-year overall survival rates.
}
\usage{data(GSE10846.Metadata)}
\format{
  A data frame with samples on the rows and metadata on the columns.
  \describe{
    \item{dlbcl.type}{DLBCL subtype: "GCB" (germinal center B-cell-like) or "ABC" (activated B-cell-like)}
    \item{gender}{Patient gender: "male", "female", or NA}
    \item{age}{Patient age in years}
    \item{LDH.ratio}{Lactate dehydrogenase ratio}
    \item{Chemotherapy}{Treatment regimen: "CHOP" or "R-CHOP"}
    \item{OS.status}{Overall survival status: "ALIVE" or "DEAD"}
    \item{OS.years}{Overall survival time in years}
    \item{cluster}{Molecular cluster assignment: "C1", "C2", or "C3"}
    \item{group}{Combined grouping variable incorporating multiple clinical factors}
  }
  Row names are sample identifiers (GSM IDs) matching those in GSE10846.Expression.
}
\source{
    Gene Expression Omnibus (GEO) accession GSE10846.
    \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846}
}
\references{
    Lenz G, Wright G, Dave SS, Xiao W et al. Stromal gene signatures in large-B-cell lymphomas. N Engl J Med 2008 Nov 27;359(22):2313-23
}
\seealso{
    \code{\link{GSE10846.Expression}}
}
\examples{
    data(GSE10846.Metadata)
    # Examine the structure of the metadata
    str(GSE10846.Metadata)
    # Summary of DLBCL subtypes
    table(GSE10846.Metadata$dlbcl.type)
    # Summary of treatment regimens
    table(GSE10846.Metadata$Chemotherapy)
}
\keyword{datasets} 