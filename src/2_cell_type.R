setwd("~/projects/yellow_fever/")
devtools::load_all("../scRNAtools/", reset = T)
load("results/QC/counts_QC.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T &
  scd$getfeature("clonality") %in% c(47, 152, 6, 162, 73, 177)
scd_test <- scdata$new(
  infos = scd$select(b_cells = b_cells)$getfeatures,
  counts = scd$select(b_cells = b_cells)$getcounts[ ,2100:2120],
  v = T
)



devtools::load_all("../scRNAtools/", reset = T)

system("mkdir -p results/test_DEA")

dea_test <- DEA(
  scd = scd_test,
  formula_null = "y ~ batch",
  formula_full = "y ~ batch + (1|clonality)",
  b_cells = scd_test$getfeature("QC_good") %in% T,
  cpus = 12,
  v = T,
  folder_name = "results/test_DEA"
)
dea_test
