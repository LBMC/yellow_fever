rm(list=ls())
setwd("~/projects/mold/yellow_fever")
devtools::load_all("pkg/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
load("results/QC/counts_QC_M.Rdata")

###############################################################################
# DEA PLS for the training data
training_cells <- paste0("P1299_", 1001:1088)
load(file = "results/cell_type/cells_counts_QC_DEA_cell_type_M.Rdata")
b_cells_training <- scd$getfeature("id") %in% training_cells
gene_expressed <- rowSums(scd$select(b_cells = b_cells_training)$getcounts > 0)
training_cells <- names(gene_expressed)[gene_expressed > 1000]

b_cells_training <- scd$getfeature("id") %in% training_cells
scd <- normalize(
  scd = scd,
  b_cells = b_cells_training,
  method = "SCnorm",
  cpus = 4,
  tmp_file = paste0("results/tmp/normalization_training_tmp.Rdata")
)
save(scd, file = "results/QC/cells_counts_QC_training.Rdata")

load("results/QC/cells_counts_QC_training.Rdata")
b_cells_training <- scd$getfeature("id") %in% training_cells
table(b_cells_training)

b_cells <- !is.na(scd$getfeature("phenotype_surface_cell_type")) &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("sex") %in% "M"
b_cells <- b_cells | scd$getfeature("id") %in% training_cells
table(b_cells)
table(!is.na(scd$select(b_cells = b_cells)$getfeature("surface_cell_type")))

dim(scd$select(b_cells = b_cells)$getcounts)
summary(as.factor(scd$select(b_cells = b_cells)$getfeature("surface_cell_type")))

system("cp results/cell_type/DEA_cell_types_force_training_lsplsstab.Rdata results/cell_type/DEA_cell_types_force_Mtraining_training_lsplsstab.Rdata")
system("cp results/cell_type/DEA_cell_types_force_classification_lplscv.Rdata results/cell_type/DEA_cell_types_force_Mtraining_classification_lplscv.Rdata")
system("rm results/cell_type/DEA_cell_types_force_Mtraining_weighting.Rdata")
system("rm results/cell_type/DEA_cell_types_force_Mtraining_classification_lpls.Rdata")

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
  output_file = "results/cell_type/DEA_cell_types_force_Mtraining",
  force = PLS_genes
)

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_force_splsstab_Mtraining.Rdata"
)

load("results/cell_type/DEA_cell_types_force_splsstab_Mtraining.Rdata")

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
table(training_cells_infos$DEA_cell_type, training_cells_infos$phenotype_surface_marker)

###############################################################################
load("results/cell_type/cells_counts_QC_DEA_cell_type_training.Rdata")

training_cells <- paste0("P1299_", 1001:1088)
b_cells_training <- scd$getfeature("id") %in% training_cells
gene_expressed <- rowSums(scd$select(b_cells = b_cells_training)$getcounts > 0)
gene_expressed[order(gene_expressed)]
training_cells <- names(gene_expressed)[gene_expressed > 3000]

b_cells <- scd$getfeature("day") %in% "NONE" & scd$getfeature("cell_number") == 1
b_cells <- scd$getfeature("id") %in% training_cells
b_bcells <- scd$getfeature("day") %in% "NONE" & scd$getfeature("cell_number") > 1

DEG_list <- c()
for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"), v=T)
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  DEG <- mbatch_DEA_cell_type_DEA$gene[b_genes]
  DEG_list <- c(DEG_list, DEG, mbatch_DEA_cell_type_DEA$gene[b_genes])
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

devtools::load_all("pkg/", reset = T)
loading <- pca_plot_a_space(
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
str(loading)
gene_loading <- loading$v
rownames(gene_loading) <- colnames(scd$select(b_cells = b_cells, genes = DEG_list_a)$getcounts)
write.csv(gene_loading, file="results/cell_type/pca_training_dataset_gene_loading.csv")
cell_loading <- loading$u
rownames(cell_loading) <- rownames(scd$select(b_cells = b_cells, genes = DEG_list_a)$getcounts)
write.csv(cell_loading, file="results/cell_type/pca_training_dataset_cell_loading.csv")

################################################################################
