% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_expression.R
\name{DEA_paraload_parameters}
\alias{DEA_paraload_parameters}
\title{helper function for DEA for scRNASeq data}
\usage{
DEA_paraload_parameters(paraload_file, scd, job_DEA_number = 5,
  formula_null, formula_full, b_cells, zi_threshold = 0.9, cpus,
  folder_name)
}
\arguments{
\item{paraload_file}{paraload table file}

\item{scd}{object to analyses}

\item{job_DEA_number}{number of DEA per paraload jobs}
}
\value{
write a paraload ready table with the parameters (as line) on which
to run every DEA_paraload() commands
}
\description{
helper function for DEA for scRNASeq data
}
\examples{
\dontrun{
DEA_paraload_parameters('results/cell_type/DEA/paraload_file.txt')
}
}
