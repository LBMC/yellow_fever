################################################################################
# classification on DEA genes for surface_cell_type

rm(list = ls())
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)

# with weights
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_surface_cell_type.Rdata")
load("results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata")
b_genes <- !is.na(mbatch_day_surface_cell_type_DEA$padj) &
  mbatch_day_surface_cell_type_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_DEA$gene[b_genes]
b_cells <- scd$getfeature("QC_good") %in% T & scd$getfeature("sex") %in% "M"
length(DEA_genes)

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

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type")) &
  scd$getfeature("sex") %in% "M"
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "phenotype_surface_cell_type",
  features = c(),
  genes = c(genes_marker, DEA_genes),
  ncores = 16,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_force",
  force = genes_marker,
)

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_force_splsstab.Rdata"
)
system("~/scripts/sms.sh \"PLS done\"")

load("results/cell_type/DEA_cell_types_force_splsstab.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type")) &
  scd$getfeature("sex") %in% "M"

factor_weight <- data.frame(
  factor = DEA_cell_type_classification$classification$fit_spls$fit$selected,
  weight = DEA_cell_type_classification$classification$model$X.weight)
factor_weight
write.csv(
  factor_weight,
  file = "results/cell_type/DEA_cell_types_force_factor_weight.csv"
)
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
cell_type_groups[cell_type_pgroups > 0.5] <- "MEM"
cell_type_groups[cell_type_pgroups < 0.5] <- "EFF"
scd$setfeature("DEA_cell_type", cell_type_groups)
scd$setfeature("pDEA_cell_type", cell_type_pgroups)

save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type_M.Rdata")
save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
load(file = "results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
scd_norm <- scd
load("results/QC/CB_counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type_M.Rdata")
save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")

load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")

################################################################################
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
PLS_type <- "DEA_cell_types_splsstab"
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("surface_cell_type"))
genes_list <- c("GZMB", "CX3CR1", "CCL4", "GNLY", "GZMH", "KLRD1", "GZMG",
  "PRF1", "HOPX", "CCL5", "GZMK", "SELL", "IL7R", "LEF1", "TCF7", "LTB",
  "NELL2", "CCR7")
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
cell_type_groups[cell_type_pgroups > 0.5] <- "MEM"
cell_type_groups[cell_type_pgroups < 0.5] <- "EFF"
scd$setfeature("DEA_cell_type", cell_type_groups)
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c("ccr7", "pDEA_cell_type"),
  order_by = "pDEA_cell_type",
  color_by = "DEA_cell_type",
  file = paste0(
    "results/cell_type/per_genes_barplot_cells_counts_QC_DEA_",
    PLS_type,
    ".pdf"),
  main = paste0("DEA DEA_cell_type ", PLS_type)
)

hm <- heatmap_genes(
  scd = scd$select(b_cells = b_cells, genes = genes_list),
  features = c("antigen", "day", "pDEA_cell_type"),
  cells_order = order(
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type")
    ))
  ),
  title = "PLS genes for DEA_cell_type",
  factor = c(T, T, F),
  file = paste0(
    "results/cell_type/per_genes_barplot_cells_counts_QC_DEA_",
    PLS_type,
    "_hm.pdf")
)
print(hm)

genes_list <- DEA_cell_type_classification$classification$fit_spls$fit$selected
features_list <- names(scd$getfeatures)[names(scd$getfeatures) %in% genes_list]
genes_list <- scd$getgenes[scd$getgenes %in% genes_list]
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c(features_list, "pDEA_cell_type"),
  order_by = "pDEA_cell_type",
  color_by = "DEA_cell_type",
  file = paste0(
    "results/cell_type/per_genes_barplot_cells_counts_QC_DEA_",
    PLS_type,
    "_selected.pdf"),
    main = paste0("DEA DEA_cell_type ", PLS_type, "_selected")
)

################################################################################
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
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type.pdf")

ggplot(data = data_gplot, aes(x = pcell_type)) +
  geom_histogram() +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type_histogram")

system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_cells_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_counts_QC_DEA_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_cells_counts_", day, "QC_good_tmp.Rdata"),
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
  tmp_file = "results/tmp/pCMF_cells_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_cells_counts_", day, "QC_good_tmp.Rdata"),
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

################################################################################
# DEA PLS for the F data
rm(list = ls())
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)

load("results/cell_type/cells_counts_QC_surface_cell_type.Rdata")
load("results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata")
b_genes <- !is.na(mbatch_day_surface_cell_type_DEA$padj) &
  mbatch_day_surface_cell_type_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_DEA$gene[b_genes]
b_cells <- scd$getfeature("QC_good") %in% T
length(DEA_genes)

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

load(file = "results/cell_type/cells_counts_QC_DEA_cell_type_M.Rdata")
scd_PLS <- scd
infos_M <- scd_PLS$getfeatures
load("results/QC/cells_counts_QC_F.Rdata")
scd$setfeature("phenotype_surface_cell_type", scd_PLS$getfeature("phenotype_surface_cell_type"))
scd$setfeature("surface_cell_type", scd_PLS$getfeature("surface_cell_type"))
scd$setfeature("psurface_cell_type", scd_PLS$getfeature("psurface_cell_type"))
scd$setfeature("DEA_cell_type", scd_PLS$getfeature("DEA_cell_type"))
scd$setfeature("pDEA_cell_type", scd_PLS$getfeature("pDEA_cell_type"))
b_cells_F <- scd$getfeature("sex") %in% "F"
infos_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getfeatures
scd <- scdata$new(
  infos = infos_M,
  counts = scd$getcounts
)

b_cells <- scd$getfeature("QC_good") %in% T
table(b_cells)
summary(as.factor(scd$getfeature("DEA_cell_type")))


# with only F data passing the QC

system("cp results/cell_type/DEA_cell_types_force_training_lsplsstab.Rdata results/cell_type/DEA_cell_types_force_MF_training_lsplsstab.Rdata")
system("cp results/cell_type/DEA_cell_types_force_classification_lplscv.Rdata results/cell_type/DEA_cell_types_force_MF_classification_lplscv.Rdata")
system("rm results/cell_type/DEA_cell_types_force_MF_classification_lpls.Rdata")

load("results/cell_type/DEA_cell_types_force_splsstab.Rdata")
PLS_genes <- DEA_cell_type_classification$classification$fit_spls$fit$selected
PLS_genes <- PLS_genes[!( PLS_genes %in% c("ZEB2") )]

devtools::load_all("../scRNAtools/", reset = T)
b_cells <- scd$getfeature("QC_good") %in% T
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "surface_cell_type",
  features = c(),
  genes = PLS_genes,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_force_MF",
  force = PLS_genes
)

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_force_splsstab_MF.Rdata"
)

load("results/cell_type/DEA_cell_types_force_splsstab_MF.Rdata")

load(file = "results/cell_type/cells_counts_QC_DEA_cell_type_M.Rdata")
b_cells_b <- scd$getfeature("QC_good") %in% T
b_cells_b <- (scd$getfeature("sex") %in% "F")[b_cells_b]
b_cells <- scd$getfeature("sex") %in% "F" & scd$getfeature("QC_good") %in% T
infos_M <- scd$getfeatures

DEA_cell_type_classification$classification$fit_spls$fit$selected
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups[b_cells_b]
scd$setfeature("DEA_cell_type", cell_type_groups)
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups[b_cells_b]
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
infos_M[b_cells, ] <- scd$select(b_cells = b_cells)$getfeatures

scd <- scdata$new(
  infos = infos_M,
  counts = scd$getcounts
)

b_cells <- scd$getfeature("QC_good") %in% T

save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type_F.Rdata")
save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
load(file = "results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
scd_norm <- scd
load("results/QC/CB_counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type_F.Rdata")
save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")

load(file = "results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
b_cells <- scd$getfeature("sex") %in% "F" & scd$getfeature("QC_good") %in% T
genes_list <- c("GZMB", "CX3CR1", "CCL4", "GNLY", "GZMH", "KLRD1", "GZMG",
  "PRF1", "HOPX", "CCL5", "GZMK", "SELL", "IL7R", "LEF1", "TCF7", "LTB",
  "NELL2", "CCR7")
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c("ccr7", "pDEA_cell_type"),
  order_by = "pDEA_cell_type",
  color_by = "DEA_cell_type",
  file = paste0(
    "results/cell_type/per_genes_barplot_cells_counts_QC_DEA_F.pdf"),
  main = paste0("DEA DEA_cell_type F")
)

genes_list <- DEA_cell_type_classification$classification$fit_spls$fit$selected
features_list <- names(scd$getfeatures)[names(scd$getfeatures) %in% genes_list]
genes_list <- scd$getgenes[scd$getgenes %in% genes_list]
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c(features_list, "pDEA_cell_type"),
  order_by = "pDEA_cell_type",
  color_by = "DEA_cell_type",
  file = paste0(
    "results/cell_type/per_genes_barplot_cells_counts_QC_DEA_F_selected.pdf"),
    main = paste0("DEA DEA_cell_type F_selected")
)
