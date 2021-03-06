% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality_control.R
\name{QC_classification}
\alias{QC_classification}
\title{helper function for QC for scRNASeq data}
\usage{
QC_classification(scd, is_blank = scd$getfeature("cell_number") == 0,
  rt_result = F, quant = "")
}
\arguments{
\item{scd}{scdata object}

\item{is_blank}{boolean vector indicating lines where there is not cells}
}
\value{
return scdata object with a "QC_good" feature
}
\description{
helper function for QC for scRNASeq data
}
\examples{
\dontrun{
QC_classification(scd)
}
}
