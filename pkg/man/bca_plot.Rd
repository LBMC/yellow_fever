% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{bca_plot}
\alias{bca_plot}
\title{plot the 2 first axis of a BCA analysis}
\usage{
bca_plot(scd, by, color, ncomp = 2, top, n_groups, n_cells, norm_by,
  file, main = "", genes_list, xlimit, ylimit)
}
\arguments{
\item{scd}{is a scRNASeq data object}

\item{file}{output file where to save the graph (optional)}

\item{gene}{gene name to display}

\item{condition}{condition to compare (must be features of scd)}

\item{transform}{(bool) use anscomb transform or raw counts}
}
\value{
a list() with the number of time each cells is classied as 'good' (
non is_blank looking) or 'bad' (blank looking)
}
\description{
plot the 2 first axis of a BCA analysis
}
\examples{
\dontrun{
check_gene(scd, "genes_a", "sexe")
}
}
