#!/usr/bin/Rscript
library(scRNAtools)
devtools::load_all("pkg/", reset = T)
scRNAtools::QC_pbs(
  scd_file = "results/abundance.Rdata",
  QC_folder = "results/QC/QC_paraload/abundance/"
)
