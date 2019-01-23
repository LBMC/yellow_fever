rm(list=ls())
setwd("~/projects/yellow_fever/")
devtools::load_all("pkg/", reset = T)
load("results/QC/counts_QC.Rdata")
scd$setfeature(
  "experiment",
  gsub("(P\\d+)_\\d+", "\\1", scd$getfeature("id"), perl = T)
)
experiment <- scd$getfeature("experiment")
experiment[scd$getfeature("day") %in% "NONE"] <- "training"
scd$setfeature("experiment", experiment)

experiment <- "training"
b_cells <- scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1

QC_score <- ifelse(b_cells,
  0.5,
  scd$getfeature("QC_score")
)
to_QC <- ifelse(b_cells,
  TRUE,
  scd$getfeature("to_QC")
)
table(b_cells)
length(scd$getfeature("to_QC"))
length(QC_score)
scd$setfeature("QC_score",
  QC_score
)
scd$setfeature("to_QC",
  to_QC
)
scRNAtools::QC_classification(
  scd = scd,
  is_blank = scd$getfeature("QC_score") < 0.5
)
save(scd, file = paste0("results/QC/counts_QC_", experiment, ".Rdata"))
print(experiment)
print(summary(scd$select(b_cells = b_cells)$getfeature("QC_good")))
print(table(scd$getfeature('sex'), scd$getfeature("QC_good")))

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "QC_good", color_name = "antigen",
  tmp_file = paste0("results/tmp/pca_counts_", experiment, "_tmp.Rdata"),
  main = "training"
)

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells & scd$getfeature("QC_good") %in% T),
  color = "batch", color_name = "clonality",
  tmp_file = paste0("results/tmp/pca_counts_QC_", experiment, "_tmp.Rdata"),
  main = "training"
)

save(scd, file = "results/QC/counts_QC_training.Rdata")

# cells effect normalization
load("results/QC/counts_QC_training.Rdata")
devtools::load_all("pkg/", reset = T)

experiment <- "training"
b_cells <- scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1 &
    scd$getfeature("QC_good") %in% T
bad_cells <- scd$
  select(b_cells = b_cells, genes = ERCC(scd, minus = T))$getcounts
bad_cells <- rowSums(bad_cells > 5) < 1000
QC_good <- scd$getfeature("QC_good")
QC_good[bad_cells] <- FALSE
scd$setfeature("QC_good", QC_good)

system(paste0("rm results/tmp/normalization_", experiment, "_tmp.Rdata"))
b_cells <- scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1 &
    scd$getfeature("QC_good") %in% T
scd <- normalize(
  scd = scd,
  b_cells = b_cells,
  method = "SCnorm",
  cpus = 3,
  tmp_file = paste0("results/tmp/normalization_", experiment, "_tmp.Rdata")
)
save(scd, file = "results/QC/cells_counts_QC_training.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "QC_good", color_name = "antigen",
  tmp_file = paste0("results/tmp/pca_cells_counts_", experiment, "_tmp.Rdata"),
  main = "training"
)

