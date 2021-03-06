% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{anscombe_inv}
\alias{anscombe_inv}
\title{inverse anscombe transform of counts}
\usage{
anscombe_inv(x)
}
\arguments{
\item{x}{vector of anscomb(counrs)}
}
\value{
return anscomb transform of x
}
\description{
inverse anscombe transform of counts
}
\examples{
\dontrun{
x = anscombe_inv(anscombe(scd$getgene("gene_a")))
}
}
