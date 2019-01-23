#!/usr/bin/Rscript
library(scRNAtools)
devtools::load_all("pkg/", reset = T)
scRNAtools::DEA_pbs()
