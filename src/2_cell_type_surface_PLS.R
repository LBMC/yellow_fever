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
genes_PLS <- read.csv("~/data/yellow_fever/2017_11_28_List_Laurent_Genes_PLS.csv")
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
scd$setfeature("phenotype_surface_cell_type", phenotype_surface_marker)
b_cells <- scd$getfeature("QC_good") %in% T


################################################################################
# classification on surface_marker
load("results/cell_type/CB_counts_QC_surface_cell_type_weighted.Rdata")
devtools::load_all("../scRNAtools/", reset = T)
b_cells <- scd$getfeature("QC_good") %in% T
data_gplot <- data.frame(
  ccr7 = scd$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd$select(b_cells = b_cells)$getfeature("phenotype_surface_cell_type")
)
cell_type_color <- scRNAtools::cell_type_palette(levels(data_gplot$cell_type))
data_gplot$cell_type <- as.vector(data_gplot$cell_type)
data_gplot$cell_type[is.na(data_gplot$cell_type)] <- "unknown"
data_gplot$cell_type <- as.factor(data_gplot$cell_type)
data_gplot$ccr7 <- as.numeric(as.vector(data_gplot$ccr7))
data_gplot$il7ra <- as.numeric(as.vector(data_gplot$il7ra))
summary(data_gplot)
cell_type_color <- c(cell_type_color, "gray")
names(cell_type_color) <- c(names(cell_type_color)[-3], "unknown")
ggplot() +
  scale_fill_manual( values = cell_type_color ) +
  scale_color_manual( values = cell_type_color ) +
  geom_point(data = data_gplot[data_gplot$cell_type == "unknown", ],
    aes(x = ccr7, y = il7ra, color = cell_type)) +
  geom_point(data = data_gplot[data_gplot$cell_type != "unknown", ],
    aes(x = ccr7, y = il7ra, color = cell_type)) +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_phenotype_surface_cell_type.pdf")
data_gplot <- data.frame(
  ccr7 = scd$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd$select(b_cells = b_cells)$getfeature("surface_cell_type"),
  pcell_type = as.numeric(as.vector(scd$select(b_cells = b_cells)$getfeature("psurface_cell_type")))
)
data_gplot$ccr7 <- as.numeric(as.vector(data_gplot$ccr7))
data_gplot$il7ra <- as.numeric(as.vector(data_gplot$il7ra))
ggplot(data = data_gplot, aes(x = ccr7, y = il7ra, color = cell_type)) +
  geom_point() +
  scale_fill_manual(
      values = scRNAtools::cell_type_palette(levels(data_gplot$cell_type))
  ) +
  scale_color_manual(
    values = scRNAtools::cell_type_palette(levels(data_gplot$cell_type))
  ) +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_surface_cell_type.pdf")

ggplot(data = data_gplot, aes(x = pcell_type)) +
  geom_histogram() +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_surface_cell_type_histogram.pdf")



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
