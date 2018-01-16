################################################################################
# classification on DEA genes for surface_cell_type

setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
load("results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata")

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
load(
  file = "results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata",
  v = T
)
b_genes <- !is.na(mbatch_day_surface_cell_type_DEA$padj) &
  mbatch_day_surface_cell_type_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_DEA$gene[b_genes]
b_cells <- scd$getfeature("QC_good") %in% T

b_cells <- scd$getfeature("QC_good") %in% T
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "phenotype_surface_cell_type",
  features = surface_marker,
  genes = c(genes_marker, DEA_genes),
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types",
  v = T
)
save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_all_splsstab.Rdata"
)

devtools::load_all("../scRNAtools/", reset = T)
b_cells <- scd$getfeature("QC_good") %in% T
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "phenotype_surface_cell_type",
  features = surface_marker,
  genes = c(genes_marker, DEA_genes),
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_force",
  force  = c(surface_marker, genes_marker),
  v = T
)
save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_force_splsstab.Rdata"
)

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type"))
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "surface_cell_type",
  features = surface_marker,
  genes = c(genes_marker, DEA_genes),
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_full"
)
save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_full_splsstab.Rdata"
)

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type"))
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "surface_cell_type",
  features = surface_marker,
  genes = c(genes_marker, DEA_genes),
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_full_force",
  force  = c(surface_marker, genes_marker),
  v = T
)
save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_full_force_splsstab.Rdata"
)

################################################################################
load("results/cell_type/DEA_cell_types_all_splsstab.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
str(DEA_cell_type_classification$classification$fit_spls$fit$selected)
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
scd$setfeature("DEA_cell_type", cell_type_groups)
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type"))

save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
scd_norm <- scd
load("results/QC/cells_counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")

table(scd$getfeature("surface_cell_type"))
table(scd$getfeature("DEA_cell_type"))
table(scd$getfeature("surface_cell_type"), scd$getfeature("DEA_cell_type"))

data_gplot <- data.frame(
  ccr7 = scd$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"),
  pcell_type = as.numeric(as.vector(scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")))
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
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type.pdf")

ggplot(data = data_gplot, aes(x = pcell_type)) +
  geom_histogram() +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type_histogram.pdf")


system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_counts_QC_DEA_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_counts_QC_DEA_cell_type_", day, ".pdf"
  ))
}

system("mkdir -p results/cell_type/pcmf")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type_", day, ".pdf"
  ))
}

system("mkdir -p results/cell_type/heatmap/")
DEA_cell_type_palette <- cell_type_palette
devtools::load_all("../scRNAtools/", reset = T)
hm <- heatmap_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("DEA_cell_type", "day", "pDEA_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")
    ))
  ),
  genes_order = scRNAtools::order_2_groups(
    scd = scd,
    b_cells = b_cells,
    genes = DEA_genes,
    by = scd$select(b_cells = b_cells)$getfeature("DEA_cell_type")
  ),
  title = "DE genes between DEA_cell_type",
  factor = c(T, T, F),
  file = "results/cell_type/heatmap/hm_CB_counts_QC_DEA_DEA_cell_type.pdf"
)
print(hm)
hm_corr <- heatmap_corr_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("DEA_cell_type", "day", "pDEA_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")
    ))
  ),
  title = "corr DE genes between DEA_cell_type",
  factor = c(T, T, F),
  file = "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_DEA_cell_type.pdf"
)
print(hm_corr)

for (day in c("D15", "D136", "D593")) {
  hm <- heatmap_genes(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    features = c("DEA_cell_type", "day", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("pDEA_cell_type")
      ))
    ),
    genes_order = order_2_groups(
      scd = scd,
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes,
      by = scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
        getfeature("DEA_cell_type")
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_CB_counts_QC_DEA_DEA_cell_type_",
      day, ".pdf"
    )
  )
  print(hm)
  hm_corr <- heatmap_corr_genes(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    features = c("DEA_cell_type", "day", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_DEA_cell_type_",
      day, ".pdf"
    )
  )
  print(hm_corr)
}

genes_list <- c("NELL2", "LTB", "GZMK", "IL7R", "GZMM", "CCL5", "NKG7", "CST7",
  "GZMA", "GZMH", "GZMB", "PRF1", "IFNG", "CCL3", "CCL4", "CCL4L1", "FASLG",
  "GNLY", "CCR7", "SELL", "ITGB1", "CXCR3", "CXCR4", "S1PR1", "S1PR5", "GPR56",
  "CXCR6", "CX3CR1", "BACH2", "TCF7", "LEF1", "TSC22D3", "JUNB", "ID2",
  "PRDM1", "ZNF683", "TBX21", "ZEB2", "HOPX", "CD69", "KLRB1", "KLRG1",
  "PDCD1", "CTLA4", "TIGIT", "TIMD4", "HAVCR2", "KLRD1", "CD2", "CD3E", "CD3D",
  "CD3G", "CD8A", "CD8B", "CD4")
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c("ccr7", "il7ra", "pDEA_cell_type"),
  order_by = "pDEA_cell_type",
  color_by = "DEA_cell_type",
  file = "results/cell_type/per_genes_barplot_CB_counts_QC_DEA_DEA_cell_type.pdf",
  main = "DEA DEA_cell_type"
)

################################################################################
load("results/cell_type/DEA_cell_types_force_splsstab.Rdata", v = T)
b_cells <- scd$getfeature("QC_good") %in% T
str(DEA_cell_type_classification$classification$fit_spls$fit$selected)
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
scd$setfeature("DEA_cell_type", cell_type_groups)
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
table(scd$getfeature("surface_cell_type"))
table(scd$getfeature("DEA_cell_type"))
table(scd$getfeature("surface_cell_type"), scd$getfeature("DEA_cell_type"))
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type"))

data_gplot <- data.frame(
  ccr7 = scd$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"),
  pcell_type = as.numeric(as.vector(scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")))
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
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type_force.pdf")

ggplot(data = data_gplot, aes(x = pcell_type)) +
  geom_histogram() +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type_histogram_force.pdf")

system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_counts_QC_DEA_cell_type_force.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_counts_QC_DEA_cell_type_", day, "_force.pdf"
  ))
}

system("mkdir -p results/cell_type/pcmf")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type_force.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type_", day, "_force.pdf"
  ))
}

system("mkdir -p results/cell_type/heatmap/")
DEA_cell_type_palette <- cell_type_palette
devtools::load_all("../scRNAtools/", reset = T)
hm <- heatmap_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("DEA_cell_type", "day", "pDEA_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")
    ))
  ),
  genes_order = scRNAtools::order_2_groups(
    scd = scd,
    b_cells = b_cells,
    genes = DEA_genes,
    by = scd$select(b_cells = b_cells)$getfeature("DEA_cell_type")
  ),
  title = "DE genes between DEA_cell_type",
  factor = c(T, T, F),
  file = "results/cell_type/heatmap/hm_CB_counts_QC_DEA_DEA_cell_type_force.pdf"
)
print(hm)
hm_corr <- heatmap_corr_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("DEA_cell_type", "day", "pDEA_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")
    ))
  ),
  title = "corr DE genes between DEA_cell_type",
  factor = c(T, T, F),
  file = "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_DEA_cell_type_force.pdf"
)
print(hm_corr)

for (day in c("D15", "D136", "D593")) {
  hm <- heatmap_genes(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    features = c("DEA_cell_type", "day", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("pDEA_cell_type")
      ))
    ),
    genes_order = order_2_groups(
      scd = scd,
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes,
      by = scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
        getfeature("DEA_cell_type")
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_CB_counts_QC_DEA_DEA_cell_type_",
      day, "_force.pdf"
    )
  )
  print(hm)
  hm_corr <- heatmap_corr_genes(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    features = c("DEA_cell_type", "day", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_DEA_cell_type_",
      day, "_force.pdf"
    )
  )
  print(hm_corr)
}

genes_list <- c("NELL2", "LTB", "GZMK", "IL7R", "GZMM", "CCL5", "NKG7", "CST7",
  "GZMA", "GZMH", "GZMB", "PRF1", "IFNG", "CCL3", "CCL4", "CCL4L1", "FASLG",
  "GNLY", "CCR7", "SELL", "ITGB1", "CXCR3", "CXCR4", "S1PR1", "S1PR5", "GPR56",
  "CXCR6", "CX3CR1", "BACH2", "TCF7", "LEF1", "TSC22D3", "JUNB", "ID2",
  "PRDM1", "ZNF683", "TBX21", "ZEB2", "HOPX", "CD69", "KLRB1", "KLRG1",
  "PDCD1", "CTLA4", "TIGIT", "TIMD4", "HAVCR2", "KLRD1", "CD2", "CD3E", "CD3D",
  "CD3G", "CD8A", "CD8B", "CD4")
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c("ccr7", "il7ra", "pDEA_cell_type"),
  order_by = "pDEA_cell_type",
  color_by = "DEA_cell_type",
  file = "results/cell_type/per_genes_barplot_CB_counts_QC_DEA_DEA_cell_type_force.pdf",
  main = "DEA_DEA_cell_type_force"
)

################################################################################
load("results/cell_type/DEA_cell_types_full_splsstab.Rdata", v = T)
b_cells <- scd$getfeature("QC_good") %in% T
str(DEA_cell_type_classification$classification$fit_spls$fit$selected)
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
scd$setfeature("DEA_cell_type", cell_type_groups)
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
table(scd$getfeature("surface_cell_type"))
table(scd$getfeature("DEA_cell_type"))
table(scd$getfeature("surface_cell_type"), scd$getfeature("DEA_cell_type"))
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type"))

data_gplot <- data.frame(
  ccr7 = scd$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"),
  pcell_type = as.numeric(as.vector(scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")))
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
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type_full.pdf")

ggplot(data = data_gplot, aes(x = pcell_type)) +
  geom_histogram() +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type_histogram_full.pdf")

system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_counts_QC_DEA_cell_type_full.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_counts_QC_DEA_cell_type_", day, "_full.pdf"
  ))
}

system("mkdir -p results/cell_type/pcmf")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type_full.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type_", day, "_full.pdf"
  ))
}

system("mkdir -p results/cell_type/heatmap/")
DEA_cell_type_palette <- cell_type_palette
devtools::load_all("../scRNAtools/", reset = T)
hm <- heatmap_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("DEA_cell_type", "day", "pDEA_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")
    ))
  ),
  genes_order = scRNAtools::order_2_groups(
    scd = scd,
    b_cells = b_cells,
    genes = DEA_genes,
    by = scd$select(b_cells = b_cells)$getfeature("DEA_cell_type")
  ),
  title = "DE genes between DEA_cell_type",
  factor = c(T, T, F),
  file = "results/cell_type/heatmap/hm_CB_counts_QC_DEA_DEA_cell_type_full.pdf"
)
print(hm)
hm_corr <- heatmap_corr_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("DEA_cell_type", "day", "pDEA_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")
    ))
  ),
  title = "corr DE genes between DEA_cell_type",
  factor = c(T, T, F),
  file = "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_DEA_cell_type_full.pdf"
)
print(hm_corr)

for (day in c("D15", "D136", "D593")) {
  hm <- heatmap_genes(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    features = c("DEA_cell_type", "day", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("pDEA_cell_type")
      ))
    ),
    genes_order = order_2_groups(
      scd = scd,
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes,
      by = scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
        getfeature("DEA_cell_type")
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_CB_counts_QC_DEA_DEA_cell_type_",
      day, "_full.pdf"
    )
  )
  print(hm)
  hm_corr <- heatmap_corr_genes(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    features = c("DEA_cell_type", "day", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_DEA_cell_type_",
      day, "_full.pdf"
    )
  )
  print(hm_corr)
}

genes_list <- c("NELL2", "LTB", "GZMK", "IL7R", "GZMM", "CCL5", "NKG7", "CST7",
  "GZMA", "GZMH", "GZMB", "PRF1", "IFNG", "CCL3", "CCL4", "CCL4L1", "FASLG",
  "GNLY", "CCR7", "SELL", "ITGB1", "CXCR3", "CXCR4", "S1PR1", "S1PR5", "GPR56",
  "CXCR6", "CX3CR1", "BACH2", "TCF7", "LEF1", "TSC22D3", "JUNB", "ID2",
  "PRDM1", "ZNF683", "TBX21", "ZEB2", "HOPX", "CD69", "KLRB1", "KLRG1",
  "PDCD1", "CTLA4", "TIGIT", "TIMD4", "HAVCR2", "KLRD1", "CD2", "CD3E", "CD3D",
  "CD3G", "CD8A", "CD8B", "CD4")
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c("ccr7", "il7ra", "pDEA_cell_type"),
  order_by = "pDEA_cell_type",
  color_by = "DEA_cell_type",
  file = "results/cell_type/per_genes_barplot_CB_counts_QC_DEA_DEA_cell_type_full.pdf",
  main = "DEA DEA_cell_type full"
)

################################################################################
load("results/cell_type/DEA_cell_types_full_force_splsstab.Rdata", v = T)
b_cells <- scd$getfeature("QC_good") %in% T
str(DEA_cell_type_classification$classification$fit_spls$fit$selected)
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
scd$setfeature("DEA_cell_type", cell_type_groups)
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
table(scd$getfeature("surface_cell_type"))
table(scd$getfeature("DEA_cell_type"))
table(scd$getfeature("surface_cell_type"), scd$getfeature("DEA_cell_type"))
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type"))

data_gplot <- data.frame(
  ccr7 = scd$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"),
  pcell_type = as.numeric(as.vector(scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")))
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
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type_full_force.pdf")

ggplot(data = data_gplot, aes(x = pcell_type)) +
  geom_histogram() +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type_histogram_full_force.pdf")

system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_counts_QC_DEA_cell_type_full_force.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_counts_QC_DEA_cell_type_", day, "_full_force.pdf"
  ))
}

system("mkdir -p results/cell_type/pcmf")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type_full_force.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type_", day, "_full_force.pdf"
  ))
}

system("mkdir -p results/cell_type/heatmap/")
DEA_cell_type_palette <- cell_type_palette
devtools::load_all("../scRNAtools/", reset = T)
hm <- heatmap_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("DEA_cell_type", "day", "pDEA_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")
    ))
  ),
  genes_order = scRNAtools::order_2_groups(
    scd = scd,
    b_cells = b_cells,
    genes = DEA_genes,
    by = scd$select(b_cells = b_cells)$getfeature("DEA_cell_type")
  ),
  title = "DE genes between DEA_cell_type",
  factor = c(T, T, F),
  file = "results/cell_type/heatmap/hm_CB_counts_QC_DEA_DEA_cell_type_full_force.pdf"
)
print(hm)
hm_corr <- heatmap_corr_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("DEA_cell_type", "day", "pDEA_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")
    ))
  ),
  title = "corr DE genes between DEA_cell_type",
  factor = c(T, T, F),
  file = "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_DEA_cell_type_full_force.pdf"
)
print(hm_corr)

for (day in c("D15", "D136", "D593")) {
  hm <- heatmap_genes(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    features = c("DEA_cell_type", "day", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("pDEA_cell_type")
      ))
    ),
    genes_order = order_2_groups(
      scd = scd,
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes,
      by = scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
        getfeature("DEA_cell_type")
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_CB_counts_QC_DEA_DEA_cell_type_",
      day, "_full_force.pdf"
    )
  )
  print(hm)
  hm_corr <- heatmap_corr_genes(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    features = c("DEA_cell_type", "day", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_DEA_cell_type_",
      day, "_full_force.pdf"
    )
  )
  print(hm_corr)
}

genes_list <- c("NELL2", "LTB", "GZMK", "IL7R", "GZMM", "CCL5", "NKG7", "CST7",
  "GZMA", "GZMH", "GZMB", "PRF1", "IFNG", "CCL3", "CCL4", "CCL4L1", "FASLG",
  "GNLY", "CCR7", "SELL", "ITGB1", "CXCR3", "CXCR4", "S1PR1", "S1PR5", "GPR56",
  "CXCR6", "CX3CR1", "BACH2", "TCF7", "LEF1", "TSC22D3", "JUNB", "ID2",
  "PRDM1", "ZNF683", "TBX21", "ZEB2", "HOPX", "CD69", "KLRB1", "KLRG1",
  "PDCD1", "CTLA4", "TIGIT", "TIMD4", "HAVCR2", "KLRD1", "CD2", "CD3E", "CD3D",
  "CD3G", "CD8A", "CD8B", "CD4")

genes_list <- c("GZMB", "CX3CR1", "CCL4", "GNLY", "GZMH", "KLRD1", "GZMG",
  "PRF1", "HOPX", "CCL5", "GZMK", "SELL", "IL7R", "LEF1", "TCF7", "LTB",
  "NELL2", "CCR7")
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c("ccr7", "pDEA_cell_type"),
  order_by = "pDEA_cell_type",
  color_by = "DEA_cell_type",
  file = "results/cell_type/per_genes_barplot_CB_counts_QC_DEA_DEA_cell_type_full_force.pdf",
  main = "DEA DEA_cell_type full force"
)

################################################################################
# classification on DEA genes for old_surface_cell_type

setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/CB_counts_QC_old_surface_cell_type.Rdata")
load("results/cell_type/mbatch_day_old_surface_cell_type_DEA.Rdata")
b_genes <- !is.na(mbatch_day_old_surface_cell_type_DEA$padj) &
  mbatch_day_old_surface_cell_type_DEA$padj < 0.05
marker_list = c("cd57", "fas", "ptprc_cd45ra", "cd4", "il7ra", "dextramer",
"cd3e", "cd8a", "ccr7", "itga6.cd49f", "pdcd1", "cd27")
genes_list <- c("CCR7", "IL7R", "GNLY", "GZMB", "GZMH", "GZMK", "GZMM", "LTB")
DEA_genes <- mbatch_day_old_surface_cell_type_DEA$gene[b_genes]
b_cells <- scd$getfeature("QC_good") %in% T

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("old_surface_cell_type"))
system("rm -R results/cell_type/DEA_old_cell_types_full_force")
system("rm -R results/cell_type/DEA_old_cell_types_full_force_classification_lpls*")
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "old_surface_cell_type",
  features = marker_list,
  genes = unique(c(marker_list, genes_list, DEA_genes)),
  ncores = 16,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_old_cell_types_full_force",
  force  = genes_list,
  v = T
)
save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_old_cell_types_full_force_splsstab.Rdata"
)
system("~/scripts/sms.sh \"PLS done\"")

load("results/cell_type/DEA_old_cell_types_full_force_splsstab.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
str(DEA_cell_type_classification$classification$fit_spls$fit$selected)
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
scd$setfeature("DEA_old_cell_type", cell_type_groups)
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
scd$setfeature("pDEA_old_cell_type", cell_type_pgroups)
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("old_surface_cell_type"))

genes_list <- c("NELL2", "LTB", "GZMK", "IL7R", "GZMM", "CCL5", "NKG7", "CST7",
  "GZMA", "GZMH", "GZMB", "PRF1", "IFNG", "CCL3", "CCL4", "CCL4L1", "FASLG",
  "GNLY", "CCR7", "SELL", "ITGB1", "CXCR3", "CXCR4", "S1PR1", "S1PR5", "GPR56",
  "CXCR6", "CX3CR1", "BACH2", "TCF7", "LEF1", "TSC22D3", "JUNB", "ID2",
  "PRDM1", "ZNF683", "TBX21", "ZEB2", "HOPX", "CD69", "KLRB1", "KLRG1",
  "PDCD1", "CTLA4", "TIGIT", "TIMD4", "HAVCR2", "KLRD1", "CD2", "CD3E", "CD3D",
  "CD3G", "CD8A", "CD8B", "CD4")

genes_list <- c("GZMB", "CX3CR1", "CCL4", "GNLY", "GZMH", "KLRD1", "GZMG",
  "PRF1", "HOPX", "CCL5", "GZMK", "SELL", "IL7R", "LEF1", "TCF7", "LTB",
  "NELL2", "CCR7")
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c("ccr7", "il7ra", "pDEA_old_cell_type"),
  order_by = "pDEA_old_cell_type",
  color_by = "DEA_old_cell_type",
  file = "results/cell_type/per_genes_barplot_CB_counts_QC_DEA_DEA_old_cell_type_full_force.pdf",
  main = "DEA DEA_old_cell_type full force"
)

################################################################################
# classification on DEA genes for raw_old_surface_cell_type

setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_raw_old_surface_cell_type.Rdata")
load("results/cell_type/mbatch_day_raw_old_surface_cell_type_DEA.Rdata")
b_genes <- !is.na(mbatch_day_raw_old_surface_cell_type_DEA$padj) &
  mbatch_day_raw_old_surface_cell_type_DEA$padj < 0.05
marker_list = c("cd57", "fas", "ptprc_cd45ra", "cd4", "il7ra", "dextramer",
"cd3e", "cd8a", "ccr7", "itga6.cd49f", "pdcd1", "cd27")
genes_list <- c("CCR7", "IL7R", "GNLY", "GZMB", "GZMH", "GZMK", "GZMM", "LTB")
DEA_genes <- mbatch_day_raw_old_surface_cell_type_DEA$gene[b_genes]
b_cells <- scd$getfeature("QC_good") %in% T

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("raw_old_surface_cell_type"))
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "raw_old_surface_cell_type",
  features = marker_list,
  genes = unique(c(marker_list, genes_list, DEA_genes)),
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_raw_old_cell_types_full_force",
  force  = genes_list,
  v = T
)
save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_raw_old_cell_types_full_force_splsstab.Rdata"
)
system("~/scripts/sms.sh \"PLS done\"")

load("results/cell_type/DEA_raw_old_cell_types_full_force_splsstab.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
str(DEA_cell_type_classification$classification$fit_spls$fit$selected)
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
scd$setfeature("DEA_raw_old_cell_type", cell_type_groups)
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
scd$setfeature("pDEA_raw_old_cell_type", cell_type_pgroups)
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("raw_old_surface_cell_type"))

genes_list <- c("NELL2", "LTB", "GZMK", "IL7R", "GZMM", "CCL5", "NKG7", "CST7",
  "GZMA", "GZMH", "GZMB", "PRF1", "IFNG", "CCL3", "CCL4", "CCL4L1", "FASLG",
  "GNLY", "CCR7", "SELL", "ITGB1", "CXCR3", "CXCR4", "S1PR1", "S1PR5", "GPR56",
  "CXCR6", "CX3CR1", "BACH2", "TCF7", "LEF1", "TSC22D3", "JUNB", "ID2",
  "PRDM1", "ZNF683", "TBX21", "ZEB2", "HOPX", "CD69", "KLRB1", "KLRG1",
  "PDCD1", "CTLA4", "TIGIT", "TIMD4", "HAVCR2", "KLRD1", "CD2", "CD3E", "CD3D",
  "CD3G", "CD8A", "CD8B", "CD4")

genes_list <- c("GZMB", "CX3CR1", "CCL4", "GNLY", "GZMH", "KLRD1", "GZMG",
  "PRF1", "HOPX", "CCL5", "GZMK", "SELL", "IL7R", "LEF1", "TCF7", "LTB",
  "NELL2", "CCR7")
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c("ccr7", "pDEA_raw_old_cell_type"),
  order_by = "pDEA_raw_old_cell_type",
  color_by = "DEA_raw_old_cell_type",
  file = "results/cell_type/per_genes_barplot_counts_QC_DEA_DEA_raw_sold_cell_type_full_force.pdf",
  main = "DEA DEA_raw_old_cell_type full force"
)

################################################################################

setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
scd_norm <- scd
load("results/QC/counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
load("results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata")

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
load(
  file = "results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata",
  v = T
)
b_genes <- !is.na(mbatch_day_surface_cell_type_DEA$padj) &
  mbatch_day_surface_cell_type_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_DEA$gene[b_genes]
b_cells <- scd$getfeature("QC_good") %in% T

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type"))
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "surface_cell_type",
  features = surface_marker,
  genes = c(genes_marker, DEA_genes),
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_full_force_raw",
  force  = c(surface_marker, genes_marker),
  v = T
)
save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_full_force_raw_splsstab.Rdata"
)
