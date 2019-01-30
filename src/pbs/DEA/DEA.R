#!/usr/bin/Rscript
library(scRNAtools)
devtools::load_all("../scRNAtools/", reset = T)
scRNAtools::DEA_pbs()
