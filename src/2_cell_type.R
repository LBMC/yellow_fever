setwd("~/projects/yellow_fever/")
devtools::load_all("../scRNAtools/", reset = T)
load("results/QC/CB_counts_QC.Rdata")
system("mkdir -p results/cell_type")

# csv table to select feature to classify the cell-type one
feature_to_select = c("tcr_found", "norm_fsc", "all_events_fsc_h_mean",
  "norm_ssc", "all_events_ssc_h_mean", "cd57", "fas", "ptprc_cd45ra", "cd4",
  "il7ra", "dextramer", "cd3e", "cd8a", "ccr7", "lineage_neg", "itga6.cd49f",
  "pdcd1", "cd27", "cfse", "quality", "phenotype_surface_marker")
to_select <- data.frame(features = c(feature_to_select, rep(NA, scd$getngenes - length(feature_to_select))),
  genes = scd$getgenes)
write.csv(to_select, file = "results/cell_type/feature_to_select.csv")

# load selection off genes and makers to classify on
genes_PLS <- read.csv("data/genes_PLS.csv")
surface_marker <- c()
genes_marker <- c()
for (marker_type in colnames(genes_PLS)) {
  for (marker in genes_PLS[[marker_type]]) {
    if (marker %in% scd$getgenes) {
      genes_marker <- c(genes_marker, marker)
    }
    if (marker %in% colnames(scd$getfeatures)) {
      surface_marker <- c(surface_marker, marker)
    }
  }
}

# build cell_type factor
phenotype_surface_marker <- scd$getfeature("phenotype_surface_marker")
phenotype_surface_marker[phenotype_surface_marker == ""] <- NA
levels(phenotype_surface_marker) <- c("", "MEM", "MEM", "EFF", "EFF", "EFF",
  "MEM", "Naive", "EFF", "EFF", "EFF", "MEM", "MEM")
phenotype_surface_marker <- as.factor(as.vector(phenotype_surface_marker))
scd$setfeature("surface_cell_type", phenotype_surface_marker)
b_cells <- scd$getfeature("QC_good") %in% T


################################################################################
# classification on surface_marker

surface_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "surface_cell_type",
  features = surface_marker,
  genes = genes_marker,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/surface_cell_types"
)
save(
  surface_cell_type_classification,
  file = "results/cell_type/surface_cell_types_all_smplscv.Rdata"
)
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- cell_type_classification$groups
scd$setfeature("surface_cell_type", cell_type_groups)

save(scd, file = "results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
scd_norm <- scd
load("results/QC/cells_counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/cells_counts_QC_surface_cell_type.Rdata")

load("results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
ggplot(data = data.frame(
  ccr7 = scd_norm$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd_norm$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd_norm$select(b_cells = b_cells)$getfeature("surface_cell_type")
), aes(x = ccr7, y = il7ra, color = cell_type)) +
  geom_point() +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_surface_cell_type.pdf")


system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells),
  color = "surface_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_counts_QC_surface_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "surface_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_counts_QC_surface_cell_type_", day, ".pdf"
  ))
}

system("mkdir -p results/cell_type/pcmf")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells),
  color = "surface_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/cell_type/pcmf/pcmf_counts_QC_surface_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "surface_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pcmf_counts_QC_surface_cell_type_", day, ".pdf"
  ))
}




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
