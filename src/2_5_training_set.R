rm(list=ls())
setwd("~/projects/yellow_fever")
devtools::load_all("pkg/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
load("results/QC/counts_QC_M.Rdata")

b_cells <- scd$getfeature("day") %in% "NONE" & scd$getfeature("cell_number") == 1
b_bcells <- scd$getfeature("day") %in% "NONE" & scd$getfeature("cell_number") > 1

DEG_list <- c()
for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
  b_genes <- !is.na(mbatch_pDEA_cell_type_DEA$padj) &
    mbatch_pDEA_cell_type_DEA$padj < 0.05
  DEG <- mbatch_pDEA_cell_type_DEA$gene[b_genes]
  DEG_list <- c(DEG_list, DEG, mbatch_pDEA_cell_type_DEA$gene[b_genes])
}

phenotype_surface_marker <- scd$getfeature("phenotype_surface_marker")
levels(phenotype_surface_marker)
phenotype_surface_marker[phenotype_surface_marker == ""] <- NA
levels(phenotype_surface_marker) <- c("", "CM", "CM", "EFF", "EM", "EM", "MEM",
                                      "Naive", "Temra", "Temra", "Temra",
                                      "TSCM", "TSCM")
phenotype_surface_marker <- as.factor(as.vector(phenotype_surface_marker))
scd$setfeature("phenotype_surface_cell_type", phenotype_surface_marker)

cell_type_color <- c("#6BB6D6", "#2B336B", "#5BA87A", "#AC3828", "#A47395")
names(cell_type_color) <- c("CM", "EM", "Naive", "Temra", "TSCM")

DEG_list_a <- expressed(scd$select(b_cells = b_cells, genes = DEG_list), zi_threshold = 0.8)
DEG_list_b <- expressed(scd$select(b_cells = b_bcells, genes = DEG_list), zi_threshold = 0.8)
pca_plot_a_space(
  data_a = scale(ascb(scd$select(b_cells = b_cells, genes = DEG_list_a)$getcounts)),
  data_b = scale(ascb(scd$select(b_cells = b_bcells, genes = DEG_list_b)$getcounts)),
  color_a = scd$select(b_cells = b_cells)$getfeature("phenotype_surface_cell_type"),
  color_b = scd$select(b_cells = b_bcells)$getfeature("phenotype_surface_cell_type"),
  main = "training set",
  name_a = "1 cell",
  name_b = "50 cells",
  file = "results/cell_type/pca_training_dataset",
  axes = c(1, 2), size = 4,
  color_pal = cell_type_color
)
