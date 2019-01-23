rm(list = ls())
setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("pkg/", reset = T)
load("results/QC/counts_QC.Rdata")
bad_F_cells <- paste0("P1292_", 1097:1192)
scd <- scd$select(b_cells = !( scd$getfeature("id") %in% bad_F_cells ))
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
rm(list = ls())
load("results/QC/counts_QC_in_vitro_P1902_P3128.Rdata")
devtools::load_all("pkg/", reset = T)

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
  QC_good[b_cells][bad_cells] <- FALSE
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
    tmp_file = paste0("results/tmp/normalization_",
                      experiment, "_tmp.Rdata"),
    FilterCellNum = 20
  )
}
save(scd, file = "results/QC/cells_counts_QC_in_vitro_P1902_P3128.Rdata")
system("~/scripts/sms.sh \"SCnorm done\"")

load("results/QC/cells_counts_QC_in_vitro_P1902_P3128.Rdata")
scd_invitro <- scd
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
bad_F_cells <- paste0("P1292_", 1097:1192)
scd <- scd$select(b_cells = !( scd$getfeature("id") %in% bad_F_cells ))
scd$setfeature(
  "experiment",
  gsub("(P\\d+)_\\d+", "\\1", scd$getfeature("id"), perl = T)
)
day <- "InVitro"
experiments <- c("P1902", "P3128")
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiments &
  scd$getfeature("cell_number") %in% 1
b_cells_invitro <- scd_invitro$getfeature("day") %in% day &
  scd_invitro$getfeature("experiment") %in% experiments &
  scd_invitro$getfeature("cell_number") %in% 1
infos <- scd$getfeatures
infos[b_cells, names(scd$getfeatures) %in% names(scd_invitro$getfeatures)] <- scd_invitro$select(b_cells = b_cells_invitro)$getfeatures
counts <- scd$getcounts
counts[b_cells, ] <- scd_invitro$select(b_cells = b_cells_invitro)$getcounts
scd <- scdata$new(
  infos = infos,
  counts = counts
)
save(scd, file = "results/QC/cells_counts_QC_in_vitro_P1902_P3128.Rdata")

genes_to_rm <- read.table("data/Genes_exclude.csv", h = T)
scd <- scd$select(genes = scd$getgenes[!scd$getgenes %in% genes_to_rm])
day <- "InVitro"

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
  write.csv(
    scd$select(b_cells = b_cells)$getcounts,
    file = paste0("results/cell_type/CB_counts_",
                  day ,"_", experiment, ".csv")
  )
  system("~/scripts/sms.sh \"Norm done\"")
}
save(scd, file = "results/QC/CB_counts_QC_in_vitro_P1902_P3128.Rdata")
system("src/dump_dropbox.sh")

load("results/QC/CB_counts_QC_in_vitro_P1902_P3128.Rdata")

summary(scd$select(b_cells = b_cells)$getgene("CCR7"))
