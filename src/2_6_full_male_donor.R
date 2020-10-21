rm(list=ls())
library(tidyverse)
library(readxl)
setwd("~/projects/mold/yellow_fever")
devtools::load_all("pkg/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")

training_cells <- read_xlsx("data/2019_08_27_Laurent_Cell_ID.xlsx", sheet=1) %>%
  pull("cell tag")

###############################################################################
# DEA PLS for the training data
b_cells_training <- scd$getfeature("id") %in% training_cells 
gene_expressed <- rowSums(scd$select(b_cells = b_cells_training)$getcounts > 0)
training_cells <- names(gene_expressed)[gene_expressed > 1000]

b_cells_training <- scd$getfeature("id") %in% training_cells
scd <- normalize(
  scd = scd,
  b_cells = b_cells_training,
  method = "SCnorm",
  cpus = 4,
  tmp_file = paste0("results/tmp/normalization_full_M_tmp.Rdata")
)
save(scd, file = "results/QC/cells_counts_QC_full_M.Rdata")

load(file = "results/QC/cells_counts_QC_full_M.Rdata")

b_cells_training <- scd$getfeature("id") %in% training_cells
table(b_cells_training)

b_cells <-  scd$getfeature("id") %in% training_cells

table(scd$getfeature("id") %in% training_cells)
table(training_cells %in% scd$getfeature("id"))

phenotype_surface_marker <- scd$getfeature("phenotype_surface_marker")
phenotype_surface_marker[phenotype_surface_marker == ""] <- NA
levels(phenotype_surface_marker) <- c("", "MEM", "MEM", "EFF", "EFF", "EFF",
  "MEM", NA, "EFF", "EFF", "EFF", "MEM", "MEM")
phenotype_surface_marker <- as.factor(as.vector(phenotype_surface_marker))
scd$setfeature("phenotype_surface_cell_type", phenotype_surface_marker)
scd$setfeature("surface_cell_type", scd$getfeature("phenotype_surface_cell_type"))
table(b_cells)
table(!is.na(scd$select(b_cells = b_cells)$getfeature("surface_cell_type")))
b_cells = b_cells & !(scd$getfeature("surface_cell_type") %in% "Naive")

dim(scd$select(b_cells = b_cells)$getcounts)
summary(as.factor(scd$select(b_cells = b_cells)$getfeature("surface_cell_type")))

system("cp results/cell_type/DEA_cell_types_force_training_lsplsstab.Rdata results/cell_type/DEA_cell_types_force_M_full_training_lsplsstab.Rdata")
system("cp results/cell_type/DEA_cell_types_force_classification_lplscv.Rdata results/cell_type/DEA_cell_types_force_M_full_classification_lplscv.Rdata")
system("rm results/cell_type/DEA_cell_types_force_M_full_weighting.Rdata")
system("rm results/cell_type/DEA_cell_types_force_M_full_classification_lpls.Rdata")

load("results/cell_type/DEA_cell_types_force_splsstab.Rdata")
PLS_genes <- DEA_cell_type_classification$classification$fit_spls$fit$selected

devtools::load_all("pkg/", reset = T)
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "surface_cell_type",
  features = c(),
  genes = PLS_genes,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_force_M_full",
  force = PLS_genes
)

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_force_splsstab_M_full.Rdata"
)

load("results/cell_type/DEA_cell_types_force_splsstab_M_full.Rdata")

load(file = "results/cell_type/cells_counts_QC_DEA_cell_type_M.Rdata")
b_cells_b <- scd$select(b_cells = b_cells)$getfeature("id") %in% training_cells
infos_M <- scd$getfeatures

DEA_cell_type_classification$classification$fit_spls$fit$selected
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells_training] <- DEA_cell_type_classification$groups[b_cells_b]
scd$setfeature("DEA_cell_type", cell_type_groups)
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells_training] <- DEA_cell_type_classification$pgroups[b_cells_b]
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
infos_M[b_cells, ] <- scd$select(b_cells = b_cells)$getfeatures

scd <- scdata$new(
  infos = infos_M,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type_training.Rdata")

training_cells_infos <- scd$select(b_cells = b_cells_training)$getfeatures
write.csv(training_cells_infos, file = "results/cell_type/training_DEA_cell_type.csv")
