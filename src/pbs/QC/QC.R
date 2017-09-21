#!/usr/bin/Rscript
library(scRNAtools)
QC_pbs(
  scd_file = "results/raw_counts.Rdata",
  QC_folder = "results/QC_paraload/"
)
