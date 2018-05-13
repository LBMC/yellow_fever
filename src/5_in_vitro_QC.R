setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("../scRNAtools/", reset = T)
load("results/QC/counts_QC.Rdata")
scd$setfeature(
  "experiment",
  gsub("(P\\d+)_\\d+", "\\1", scd$getfeature("id"), perl = T)
)

day <- "InVitro"
for (experiment in c("P1902", "P3128")) {
  b_cells <- scd$getfeature("day") %in% day &
    scd$getfeature("experiment") %in% experiment &
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
    main = "all day"
  )

  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("QC_good") %in% T),
    color = "batch", color_name = "clonality",
    tmp_file = paste0("results/tmp/pca_counts_QC_", experiment, "_tmp.Rdata"),
    main = "all day"
  )
}

save(scd, file = "results/QC/counts_QC_in_vitro_P1902_P3128.Rdata")

# cells effect normalization
load("results/QC/counts_QC_in_vitro_P1902_P3128.Rdata")
devtools::load_all("../scRNAtools/", reset = T)

day <- "InVitro"
for (experiment in c("P1902", "P3128")) {
  b_cells <- scd$getfeature("day") %in% day &
    scd$getfeature("experiment") %in% experiment &
    scd$getfeature("cell_number") %in% 1
    scd$getfeature("QC_good") %in% T
  bad_cells <- scd$
    select(b_cells = b_cells, genes = ERCC(scd, minus = T))$getcounts
  bad_cells <- rowSums(bad_cells > 5) < 1000
  QC_good <- scd$getfeature("QC_good")
  QC_good[bad_cells] <- FALSE
  scd$setfeature("QC_good", QC_good)
}

for (experiment in c("P1902", "P3128")) {
  system(paste0("rm results/tmp/normalization_", experiment, "_tmp.Rdata"))
  b_cells <- scd$getfeature("day") %in% day &
    scd$getfeature("experiment") %in% experiment &
    scd$getfeature("cell_number") %in% 1 &
    scd$getfeature("QC_good") %in% T
  scd <- normalize(
    scd = scd,
    b_cells = b_cells,
    method = "SCnorm",
    cpus = 3,
    tmp_file = paste0("results/tmp/normalization_", experiment, "_tmp.Rdata")
  )
}
save(scd, file = "results/QC/cells_counts_QC_in_vitro_P1902_P3128.Rdata")

load("results/QC/cells_counts_QC_in_vitro_P1902_P3128.Rdata")

for (experiment in c("P1902", "P3128")) {
  system(paste0(
      "rm results/tmp/normalization_cells_combat_,", experiment, "_tmp.Rdata"
    ))
  b_cells <- scd$getfeature("day") %in% day &
    scd$getfeature("experiment") %in% experiment &
    scd$getfeature("cell_number") %in% 1
    scd$getfeature("QC_good") %in% T
  scd <- normalize(
    scd = scd,
    b_cells = b_cells,
    method = "ComBat",
    cpus = 5,
    tmp_file = paste0(
      "results/tmp/normalization_cells_combat_,", experiment, "_tmp.Rdata"
    )
  )
}
save(scd, file = "results/QC/CB_counts_QC_in_vitro_P1902_P3128.Rdata")
load("results/QC/CB_counts_QC_in_vitro_P1902_P3128.Rdata")

summary(scd$select(b_cells = b_cells)$getgene("CCR7"))
