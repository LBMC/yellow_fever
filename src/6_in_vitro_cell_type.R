################################################################################
# DEA PLS for the InVitro P1902 data
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)

load("results/cycling/CB_counts_QC_cycling_invitro_P1902_P3128.Rdata")
load("results/cell_type/mbatch_day_surface_cell_type_weighted_DEA.Rdata")
b_genes <- !is.na(mbatch_day_surface_cell_type_weighted_DEA$padj) &
  mbatch_day_surface_cell_type_weighted_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_weighted_DEA$gene[b_genes]
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

day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1

system("cp results/cell_type/DEA_cell_types_weighted_force_training_lsplsstab.Rdata results/cell_type/DEA_cell_types_weighted_force_invitro_P1902_training_lsplsstab.Rdata")
system("cp results/cell_type/DEA_cell_types_weighted_force_classification_lplscv.Rdata results/cell_type/DEA_cell_types_weighted_force_invitro_P1902_classification_lplscv.Rdata")
system("rm results/cell_type/DEA_cell_types_weighted_force_invitro_P1902_classification_lpls.Rdata")

load("results/cell_type/DEA_cell_types_weighted_force_splsstab.Rdata")

DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "founder_phenotype",
  features = c(),
  genes = DEA_cell_type_classification$classification$fit_spls$fit$selected,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_weighted_force_invitro_P1902",
  force = DEA_cell_type_classification$classification$fit_spls$fit$selected
)

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_weighted_force_splsstab_invitro_P1902.Rdata"
)

load("results/cell_type/DEA_cell_types_weighted_force_splsstab_invitro_P1902.Rdata")

founder_phenotype <- scd$getfeature("founder_phenotype")
founder_phenotype[b_cells & scd$getfeature("clonality") %in% "A7"] <- NA
scd$setfeature("founder_phenotype", founder_phenotype)

system("rm results/cell_type/DEA_cell_types_weighted_force_invitro_P1902_denovo*")
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "founder_phenotype",
  features = c(),
  genes = c(DEA_genes, genes_marker, "KLRD1", "SELL"),
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_weighted_force_invitro_P1902_denovo",
  force = unique(c(genes_marker, "KLRD1", "SELL"))
)

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_weighted_force_splsstab_invitro_P1902_denovo.Rdata"
)

load("results/cell_type/DEA_cell_types_weighted_force_splsstab_invitro_P1902_denovo.Rdata")

day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1

DEA_cell_type_classification$classification$fit_spls$fit$selected
cell_type_groups <- scd$getfeature("DEA_cell_type")
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
scd$setfeature("DEA_cell_type", cell_type_groups)
cell_type_pgroups <- scd$getfeature("pDEA_cell_type")
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P1902.Rdata")

x11()
load("results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
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
    "results/cell_type/per_genes_barplot_CB_counts_QC_DEA_invitro_P1902.pdf"),
  main = paste0("DEA DEA_cell_type InVitro P1902")
)

tmp_infos <- scd$select(b_cells = b_cells)$getfeatures
tmp_infos$clonality <- as.factor( as.vector(tmp_infos$clonality) )
tmp_infos$clonality <- factor(
  tmp_infos$clonality,
  levels = c("A7", "A8", "G6", "G8", "H9", "F3", "E4", "H2", "B4"))
tmp_infos$pDEA_clone_cell_type <- unlist(as.list(by(
  tmp_infos$pDEA_cell_type, tmp_infos$clonality, median
))[tmp_infos$clonality])
tmp_infos$DEA_clone_cell_type <- ifelse(
  tmp_infos$pDEA_clone_cell_type > 0.4,
  "MEM",
  "EFF"
)
g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log(cycling_score),
      color = DEA_clone_cell_type)) +
  geom_jitter() +
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$DEA_cell_type)))) +
  labs(x = "clones",
       y =  "cell-cycle score",
       color = "DEA_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P1902_cycling_vs_DEA_cell_type.pdf"
))

save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
load(file = "results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
scd_norm <- scd
load("results/QC/cells_counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")


day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1

# DEA analysis InVitro P1902

load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
system(
  paste0("mkdir -p results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA")
)
mbatch_DEA_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + DEA_cell_type",
  b_cells = b_cells,
  cpus = 16,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA")
)
save(
  mbatch_DEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
write.csv(
  mbatch_DEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA.csv")
)

load(
  paste0("results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA.Rdata")
)


for (alt in c("lesser", "greater")) {
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  DEA_genes <- mbatch_DEA_cell_type_DEA$gene
  DEA_genes_exp <- expressed(
    scd = scd$select(
      b_cells = b_cells,
      genes = DEA_genes
    ),
    zi_threshold = 0.50,
    alt = alt
  )
  scd_norm <- scd$select(
      b_cells = b_cells,
      genes = DEA_genes_exp
    )
  clonality <- scd_norm$getfeature("clonality")
  clonality <- factor(
    clonality,
    levels = c("A7", "A8", "G6", "G8", "H9", "F3", "E4", "H2", "B4")
  )
  scd_norm$setfeature("clone", clonality)

  clone_palette <- function(clonality, other_set = TRUE){
    clonality_color <- RColorBrewer::brewer.pal(
      10, "RdYlBu"
    )[-6]
    names(clonality_color) <- levels(clonality)
    return(clonality_color)
  }
  hm_title <- paste0("top ", length(DEA_genes_exp),
    "DE genes between cell-type
    ", day, "_", experiment, "(zi <= 0.5)")
  if (alt == "greater") {
    hm_title <- paste0("top ", length(DEA_genes_exp),
      "DE genes between cell-type
      ", day, "_", experiment, "(zi >= 0.5)")
  }
  gene_order <- order_by_factor(
    scd = scd_norm,
    score = scd_norm$getfeature("pDEA_cell_type"),
    tmp_file = paste0("results/cell_type/gene_cov_",
      day, "_", experiment, "_zi0.5_", alt, ".Rdata"),
    top = min(100, length(DEA_genes_exp))
  )
  gene_order <- order_by_factor(
    scd = scd_norm,
    score = scd_norm$getfeature("pDEA_cell_type"),
    tmp_file = paste0("results/cell_type/gene_cov_",
      day, "_", experiment, "_zi0.5_", alt, "_top100.Rdata"),
    top = min(100, length(DEA_genes_exp))
  )
  hm <- heatmap_genes(
    scd = scd_norm,
    features = c("clone", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    genes_order = gene_order,
    title = hm_title,
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_CB_counts_QC_DEA_cell_type_",
      day, "_", experiment, "_zi0.5_", alt, "_cov.pdf"
    ),
    gene_size = 3
  )
  print(hm)

  hm_corr <- heatmap_corr_genes(
    scd = scd_norm,
    features = c("clone", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    title = hm_title,
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, "_", experiment, "_zi0.5_", alt, "_cov_manhattan.pdf"
    ),
    dist_name = "manhattan",
    pca = F,
    pCMF = F,
    ncomp = 5,
    cpus = 10
  )
  print(hm_corr)
}

################################################################################
# DEA PLS for the InVitro P3128 data
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)

load("results/cycling/CB_counts_QC_cycling_invitro_P1902_P3128.Rdata")
load("results/cell_type/mbatch_day_surface_cell_type_weighted_DEA.Rdata")
b_genes <- !is.na(mbatch_day_surface_cell_type_weighted_DEA$padj) &
  mbatch_day_surface_cell_type_weighted_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_weighted_DEA$gene[b_genes]
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

day <- "InVitro"
experiment <- "P3128"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1

system("cp results/cell_type/DEA_cell_types_weighted_force_training_lsplsstab.Rdata results/cell_type/DEA_cell_types_weighted_force_invitro_P3128_training_lsplsstab.Rdata")
system("cp results/cell_type/DEA_cell_types_weighted_force_classification_lplscv.Rdata results/cell_type/DEA_cell_types_weighted_force_invitro_P3128_classification_lplscv.Rdata")
system("rm results/cell_type/DEA_cell_types_weighted_force_invitro_P3128_classification_lpls.Rdata")

load("results/cell_type/DEA_cell_types_weighted_force_splsstab.Rdata")

founder_phenotype <- scd$getfeature("founder_phenotype")
summary(as.factor(founder_phenotype[b_cells]))
founder_phenotype[
  b_cells & founder_phenotype %in% "UNK"] <- NA
scd$setfeature("founder_phenotype", founder_phenotype)

DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "founder_phenotype",
  features = c(),
  genes = DEA_cell_type_classification$classification$fit_spls$fit$selected,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_weighted_force_invitro_P3128",
  force = DEA_cell_type_classification$classification$fit_spls$fit$selected
)

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_weighted_force_splsstab_invitro_P3128.Rdata"
)

load("results/cell_type/DEA_cell_types_weighted_force_splsstab_invitro_P3128.Rdata")

day <- "InVitro"
experiment <- "P3128"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1

DEA_cell_type_classification$classification$fit_spls$fit$selected
cell_type_groups <- scd$getfeature("DEA_cell_type")
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
scd$setfeature("DEA_cell_type", cell_type_groups)
cell_type_pgroups <- scd$getfeature("pDEA_cell_type")
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P3128.Rdata")

load("results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P3128.Rdata")
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
    "results/cell_type/per_genes_barplot_CB_counts_QC_DEA_invitro_P3128.pdf"),
  main = paste0("DEA DEA_cell_type InVitro P3128")
)

tmp_infos <- scd$select(b_cells = b_cells)$getfeatures
tmp_infos$clonality <- as.factor( as.vector(tmp_infos$clonality) )
tmp_infos$pDEA_clone_cell_type <- unlist(as.list(by(
  tmp_infos$pDEA_cell_type, tmp_infos$clonality, median
))[tmp_infos$clonality])
tmp_infos$DEA_clone_cell_type <- ifelse(
  tmp_infos$pDEA_clone_cell_type > 0.4,
  "MEM",
  "EFF"
)
g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log(cycling_score),
      color = DEA_clone_cell_type)) +
  geom_jitter() +
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$DEA_cell_type)))) +
  labs(x = "clones",
       y =  "cell-cycle score",
       color = "DEA_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P3128_cycling_vs_DEA_cell_type.pdf"
))

save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P3128.Rdata")
load(file = "results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P3128.Rdata")
scd_norm <- scd
load("results/QC/cells_counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P3128.Rdata")


day <- "InVitro"
experiment <- "P3128"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1

# DEA analysis InVitro P3128

load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P3128.Rdata")
system(
  paste0("mkdir -p results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA")
)
mbatch_DEA_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + DEA_cell_type",
  b_cells = b_cells,
  cpus = 16,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA")
)
save(
  mbatch_DEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
write.csv(
  mbatch_DEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA.csv")
)


load(
  paste0("results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA.Rdata")
)

for (alt in c("lesser", "greater")) {
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  DEA_genes <- mbatch_DEA_cell_type_DEA$gene
  DEA_genes_exp <- expressed(
    scd = scd$select(
      b_cells = b_cells,
      genes = DEA_genes
    ),
    zi_threshold = 0.50,
    alt = alt
  )
  scd_norm <- scd$select(
      b_cells = b_cells,
      genes = DEA_genes_exp
    )
  clonality <- scd_norm$getfeature("clonality")
  clonality <- as.factor(as.vector(clonality))
  scd_norm$setfeature("clonality", clonality)

  clone_palette <- function(clonality, other_set = TRUE){
    clonality_color <- RColorBrewer::brewer.pal(
      length(levels(clonality)) + 1, "Set1"
    )[-(round(length(levels(clonality))/2) +1)]
    names(clonality_color) <- levels(clonality)
    return(clonality_color)
  }

  hm_title <- paste0("top ", length(DEA_genes_exp),
    "DE genes between cell-type
    ", day, "_", experiment, "(zi <= 0.5)")
  if (alt == "greater") {
    hm_title <- paste0("top ", length(DEA_genes_exp),
      "DE genes between cell-type
      ", day, "_", experiment, "(zi >= 0.5)")
  }
  gene_order <- order_by_factor(
    scd = scd_norm,
    score = scd_norm$getfeature("pDEA_cell_type"),
    tmp_file = paste0("results/cell_type/gene_cov_",
      day, "_", experiment, "_zi0.5_", alt, "_top100.Rdata"),
    top = min(100, length(DEA_genes_exp))
  )
  infos_norm <- scd_norm$getfeatures
  infos_norm$clonality <- as.factor(as.vector(scd_norm$getfeature("clonality")))
  scd_norm <- scdata$new(
    infos = infos_norm,
    counts = scd_norm$getcounts
  )
  table(scd_norm$getfeature("clonality"))

  hm <- heatmap_genes(
    scd = scd_norm,
    features = c("clonality", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    genes_order = gene_order,
    title = hm_title,
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_CB_counts_QC_DEA_cell_type_",
      day, "_", experiment, "_zi0.5_", alt, "_cov.pdf"
    )
  )
  print(hm)

  hm_corr <- heatmap_corr_genes(
    scd = scd_norm,
    features = c("clonality", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    title = hm_title,
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, "_", experiment, "_zi0.5_", alt, "_cov_manhattan.pdf"
    ),
    dist_name = "manhattan",
    pca = F,
    pCMF = F,
    ncomp = 5,
    cpus = 10
  )
  print(hm_corr)
}
