################################################################################
# DEA PLS for the InVitro P1902 data
Sys.setenv("DISPLAY"=":1")
rm(list=ls())
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)

load("results/cycling/cells_counts_QC_cycling_invitro_P1902_P3128.Rdata")
b_cells <- scd$getfeature("day") %in% "InVitro" &
  scd$getfeature("experiment") %in% c( "P1902", "P3128" ) &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1
infos_M <- scd$select(b_cells = b_cells)$getfeatures
counts_M <- scd$select(b_cells = b_cells)$getcounts
infos_M <- t(infos_M)
counts_M <- t(counts_M)
dim(infos_M)
dim(counts_M)
normalized_counts <- rbind(infos_M, counts_M)
dim(normalized_counts)
write.csv(
  normalized_counts,
  file = paste0("results/cell_type/cell_type_cells_counts_in_vitro.csv")
)

genes_to_rm <- read.table("data/Genes_exclude.csv", h = T)

rm(list = ls())
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata", v = T)
b_genes <- !is.na(mbatch_day_surface_cell_type_DEA$padj) &
  mbatch_day_surface_cell_type_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_DEA$gene[b_genes]
length(DEA_genes)


# DEA analysis InVitro P1902 at the founder level
load("results/cycling/cells_counts_QC_cycling_invitro_P1902_P3128.Rdata", v = T)
scd <- scd$select(genes = scd$getgenes[!scd$getgenes %in% genes_to_rm])
day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

system(
  paste0("rm -R results/cell_type/mbatch_", day, "_", experiment, "_founder_phenotype_DEA")
)
system(
  paste0("mkdir -p results/cell_type/mbatch_", day, "_", experiment, "_founder_phenotype_DEA")
)
mbatch_founder_phenotype_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + founder_phenotype",
  b_cells = b_cells,
  cpus = 8,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", day, "_", experiment, "_founder_phenotype_DEA")
)
save(
  mbatch_founder_phenotype_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_founder_phenotype_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_founder_phenotype_DEA$padj)))
print(table(mbatch_founder_phenotype_DEA$padj < 0.05))
write.csv(
  mbatch_founder_phenotype_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_founder_phenotype_DEA.csv")
)

scd_norm <- scd$select(
    b_cells = b_cells,
  )
clonality <- scd_norm$getfeature("clonality")
clonality <- factor(
  clonality,
  levels = c("A7", "A8", "G6", "G8", "H9", "F3", "E4", "H2", "B4")
)
scd_norm$setfeature("clone", clonality)

b_genes <- !is.na(mbatch_founder_phenotype_DEA$padj) &
  mbatch_founder_phenotype_DEA$padj < 0.05
DEA_genes <- mbatch_founder_phenotype_DEA$gene[b_genes]

cells_order <- order(scd_norm$getfeature("founder_phenotype"),
                    scd_norm$getfeature("clone"))
genes_order <- order_TMP(
  scd = scd_norm$select(genes = DEA_genes),
  by = scd_norm$getfeature("founder_phenotype"),
)
clone_palette <- function(clonality, other_set = TRUE){
  clonality_color <- RColorBrewer::brewer.pal(
    10, "RdYlBu"
  )[-6]
  names(clonality_color) <- levels(clonality)
  return(clonality_color)
}
founder_phenotype_palette <- cell_type_palette


devtools::load_all("../scRNAtools/", reset = T)
split_heatmap(
  scd = scd_norm,
  genes = DEA_genes,
  features = c("clone", "founder_phenotype"),
  title = "hm_CB_counts_QC_founder_phenotype_",
  factor = c(T, T),
  genes_order = genes_order,
  cells_order = cells_order,
  file = paste0(
    "results/cell_type/heatmap/hm_CB_counts_QC_", day, "_",
    experiment, "_founder_phenotype"
  )
)


###############################################################################

# 1st PLS based on founder_phenotype DEA and annotation
load("results/cycling/cells_counts_QC_cycling_invitro_P1902_P3128.Rdata", v = T)
day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1
founder_phenotype <- scd$getfeature("founder_phenotype")
founder_phenotype[b_cells & scd$getfeature("clonality") %in% "A7"] <- NA
scd$setfeature("founder_phenotype", founder_phenotype)
scd <- scd$select(genes = scd$getgenes[!scd$getgenes %in% genes_to_rm])

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

PLS_genes <- scd$getgenes[scd$getgenes %in%
                          unique(c(DEA_genes, genes_marker, "KLRD1", "SELL"))]
# filter out PLS genes with small counts
good_pls <- apply(scd$select(b_cells = b_cells,
                             genes = PLS_genes)$getcounts,
                  2, FUN = function(x){
                    if(max(x) == 0) {
                      return(FALSE)
                    }
                    return(median(x[x > 0]) > 10)
                  })
PLS_genes <- PLS_genes[good_pls]
table(good_pls)

summary(as.factor(scd$select(b_cells = b_cells)$getfeature("founder_phenotype")))

system("rm results/cell_type/founder_cell_types_invitro_P1902_denovo*")
founder_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "founder_phenotype",
  features = c(),
  genes = PLS_genes,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/founder_cell_types_invitro_P1902_denovo",
  force = unique(c(genes_marker, "KLRD1", "SELL"))[
            unique(c(genes_marker, "KLRD1", "SELL")) %in% PLS_genes
          ]
)

save(
  founder_cell_type_classification,
  file = "results/cell_type/founder_cell_types_splsstab_invitro_P1902_denovo.Rdata"
)
system("~/scripts/sms.sh \"PLS done\"")


load("results/cell_type/founder_cell_types_splsstab_invitro_P1902_denovo.Rdata")
day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

cell_type_pgroups <- scd$getfeature("pfounder_cell_type")
cell_type_pgroups[b_cells] <- founder_cell_type_classification$pgroups
scd$setfeature("pfounder_cell_type", cell_type_pgroups)
cell_type_groups <- scd$getfeature("founder_cell_type")
cell_type_groups[b_cells] <- ifelse(cell_type_pgroups[b_cells] > 0.5, "MEM", "EFF")
scd$setfeature("founder_cell_type", cell_type_groups)
save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")


load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
genes_list <- c("GZMB", "CX3CR1", "CCL4", "GNLY", "GZMH", "KLRD1", "GZMG",
  "PRF1", "HOPX", "CCL5", "GZMK", "SELL", "IL7R", "LEF1", "TCF7", "LTB",
  "NELL2", "CCR7")
per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = genes_list,
  features = c("ccr7", "pfounder_cell_type"),
  order_by = "pfounder_cell_type",
  color_by = "founder_cell_type",
  file = paste0(
    "results/cell_type/per_genes_barplot_CB_counts_QC_founder_invitro_P1902.pdf"),
  main = paste0("founder_cell_type InVitro P1902")
)

tmp_infos <- scd$select(b_cells = b_cells)$getfeatures
tmp_infos$clonality <- as.factor( as.vector(tmp_infos$clonality) )
tmp_infos$clonality <- factor(
  tmp_infos$clonality,
  levels = c("A7", "A8", "G6", "G8", "H9", "F3", "E4", "H2", "B4"))
tmp_infos$pfounder_clone_cell_type <- unlist(as.list(by(
  tmp_infos$pfounder_cell_type, tmp_infos$clonality, median
))[tmp_infos$clonality])
tmp_infos$founder_clone_cell_type <- ifelse(
  tmp_infos$pfounder_clone_cell_type > 0.5,
  "MEM",
  "EFF"
)
tmp_infos$KLRD1 <- scd$select(b_cells = b_cells)$getgene("KLRD1")
tmp_infos$SELL <- scd$select(b_cells = b_cells)$getgene("SELL")
tmp_infos$MKI67 <- scd$select(b_cells = b_cells)$getgene("MKI67")

table(tmp_infos$founder_clone_cell_type)

g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log(cycling_score),
      color = founder_clone_cell_type)) +
  geom_jitter() +
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$founder_cell_type)))) +
  labs(x = "clones",
       y =  "cell-cycle score",
       color = "founder_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P1902_cycling_vs_founder_cell_type.pdf"
))

g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log10(KLRD1 + 1),
      color = founder_clone_cell_type)) +
  geom_jitter() +
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$founder_cell_type)))) +
  labs(x = "clones",
       color = "foudner_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P1902_KLRD1_vs_founder_cell_type.pdf"
))

g <- ggplot(tmp_infos,
   aes(x = pfounder_cell_type,
      y = log10(KLRD1 + 1),
      color = founder_clone_cell_type)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$founder_cell_type)))) +
  labs(x = "clones",
       color = "foudner_cell_type",
       title = day)
print(g)

g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log10(SELL + 1),
      color = founder_clone_cell_type)) +
  geom_jitter() +
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$founder_cell_type)))) +
  labs(x = "clones",
       color = "founder_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P1902_SELL_vs_founder_cell_type.pdf"
))

g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log10(MKI67 + 1),
      color = founder_clone_cell_type)) +
  geom_jitter() +
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$founder_cell_type)))) +
  labs(x = "clones",
       color = "founder_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P1902_MKI67_vs_founder_cell_type.pdf"
))

# DEA analysis InVitro P1902
load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

system(
  paste0("rm -R results/cell_type/mbatch_", day, "_", experiment, "_pfounder_cell_type_DEA")
)
system(
  paste0("mkdir -p results/cell_type/mbatch_", day, "_", experiment, "_pfounder_cell_type_DEA")
)
mbatch_pfounder_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + pfounder_cell_type",
  b_cells = b_cells,
  continuous = "pfounder_cell_type",
  cpus = 4,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", day, "_", experiment, "_pfounder_cell_type_DEA")
)
save(
  mbatch_pfounder_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_pfounder_cell_type_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_pfounder_cell_type_DEA$padj)))
print(table(mbatch_pfounder_cell_type_DEA$padj < 0.05))
write.csv(
  mbatch_pfounder_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_pfounder_cell_type_DEA.csv")
)

# 2nd PLS based on DEA genes
day <- "InVitro"
experiment <- "P1903"
load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1903.Rdata")
load(paste1("results/cell_type/mbatch_", day, "_", experiment, "_pfounder_cell_type_DEA.Rdata"))
b_genes <- !is.na(mbatch_pfounder_cell_type_DEA$padj) &
  mbatch_pfounder_cell_type_DEA$padj < 1.05
DEA_genes <- mbatch_pfounder_cell_type_DEA$gene[b_genes]
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 2
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

system("rm results/cell_type/DEA_cell_types_invitro_P1902_denovo*")
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "founder_cell_type",
  features = c(),
  genes = c(genes_marker, DEA_genes),
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_invitro_P1902_denovo",
  force =  genes_marker
)
system("~/scripts/sms.sh \"PLS done\"")

summary( scd$select(b_cells = b_cells, genes = c(genes_marker, DEA_genes))$getcounts )

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_splsstab_invitro_P1902_denovo.Rdata"
)

load("results/cell_type/DEA_cell_types_splsstab_invitro_P1902_denovo.Rdata")

factor_weight <- data.frame(
  factor = DEA_cell_type_classification$classification$fit_spls$fit$selected,
  weight = DEA_cell_type_classification$classification$model$X.weight)
factor_weight
write.csv(
  factor_weight,
  file = "results/cell_type/DEA_cell_types_P1902_factor_weight.csv"
)

system("rm results/cell_type/DEA_cell_types_invitro_P1902_denovo_noforce*")
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "founder_cell_type",
  features = c(),
  genes = c(genes_marker, DEA_genes),
  ncores = 6,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_invitro_P1902_denovo_noforce",
  force = c()
)
system("~/scripts/sms.sh \"PLS done\"")

summary( scd$select(b_cells = b_cells, genes = c(genes_marker, DEA_genes))$getcounts )

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_splsstab_invitro_P1902_denovo_noforce.Rdata"
)

load("results/cell_type/DEA_cell_types_splsstab_invitro_P1902_denovo.Rdata")

factor_weight <- data.frame(
  factor = DEA_cell_type_classification$classification$fit_spls$fit$selected,
  weight = DEA_cell_type_classification$classification$model$X.weight)
factor_weight
write.csv(
  factor_weight,
  file = "results/cell_type/DEA_cell_types_P1902_noforce_factor_weight.csv"
)


day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

DEA_cell_type_classification$classification$fit_spls$fit$selected
cell_type_groups <- scd$getfeature("DEA_cell_type")
cell_type_groups[b_cells] <- DEA_cell_type_classification$groups
scd$setfeature("DEA_cell_type", cell_type_groups)
cell_type_pgroups <- scd$getfeature("pDEA_cell_type")
cell_type_pgroups[b_cells] <- DEA_cell_type_classification$pgroups
scd$setfeature("pDEA_cell_type", cell_type_pgroups)
save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")

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
  main = paste0("founder_cell_type InVitro P1902")
)

per_genes_barplot(
  scd = scd$select(b_cells = b_cells),
  genes = DEA_cell_type_classification$classification$fit_spls$fit$selected,
  features = c("ccr7", "pDEA_cell_type"),
  order_by = "pDEA_cell_type",
  color_by = "DEA_cell_type",
  file = paste0(
    "results/cell_type/per_genes_barplot_CB_counts_QC_DEA_invitro_P1902_selected.pdf"),
  main = paste0("founder_cell_type InVitro P1902")
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
  tmp_infos$pDEA_clone_cell_type > 0.5,
  "MEM",
  "EFF"
)
tmp_infos$KLRD1 <- scd$select(b_cells = b_cells)$getgene("KLRD1")
tmp_infos$SELL <- scd$select(b_cells = b_cells)$getgene("SELL")
tmp_infos$MKI67 <- scd$select(b_cells = b_cells)$getgene("MKI67")

g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log(cycling_score),
      color = DEA_clone_cell_type)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(aes(color = DEA_cell_type)) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$founder_cell_type)))) +
  labs(x = "clones",
       y =  "cell-cycle score",
       color = "founder_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P1902_cycling_vs_DEA_cell_type.pdf"
))

g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log10(SELL+1),
      color = DEA_clone_cell_type)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(aes(color = DEA_cell_type)) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$founder_cell_type)))) +
  labs(x = "clones",
       y =  "log10(SELL+1)",
       color = "founder_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P1902_SELL_vs_DEA_cell_type.pdf"
))

# 2nd DEA analysis InVitro P1902
load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

system(
  paste0("rm -R results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA")
)
system(
  paste0("mkdir -p results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA")
)
mbatch_pDEA_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + pDEA_cell_type",
  b_cells = b_cells,
  continuous = "pDEA_cell_type",
  cpus = 8,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA")
)
save(
  mbatch_pDEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_pDEA_cell_type_DEA$padj)))
print(table(mbatch_pDEA_cell_type_DEA$padj < 0.05))
write.csv(
  mbatch_pDEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA.csv")
)


system(
  paste0("rm -R results/cell_type/mbatch_", day, "_", experiment, "_noforce_pDEA_cell_type_DEA")
)

system(
  paste0("mkdir -p results/cell_type/mbatch_", day, "_", experiment, "_noforce_pDEA_cell_type_DEA")
)
mbatch_pDEA_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + pDEA_cell_type",
  b_cells = b_cells,
  continuous = "pDEA_cell_type",
  cpus = 10,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", day, "_", experiment, "_noforce_pDEA_cell_type_DEA")
)
save(
  mbatch_pDEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_noforce_pDEA_cell_type_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_pDEA_cell_type_DEA$padj)))
print(table(mbatch_pDEA_cell_type_DEA$padj < 0.05))
write.csv(
  mbatch_pDEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_noforce_pDEA_cell_type_DEA.csv")
)
############################

load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1
load(paste0("results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA.Rdata"))

for (alt in c("lesser", "greater")) {
  b_genes <- !is.na(mbatch_pDEA_cell_type_DEA$padj) &
    mbatch_pDEA_cell_type_DEA$padj < 0.05
  DEA_genes <- mbatch_pDEA_cell_type_DEA$gene[b_genes]
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
  system( paste0("rm results/cell_type/gene_cov_",
                 day, "_", experiment, "_zi0.5_", alt, "_top100*") )
  gene_order <- order_by_factor(
    scd = scd_norm,
    score = scd_norm$getfeature("pDEA_cell_type"),
    tmp_file = paste0("results/cell_type/gene_cov_",
      day, "_", experiment, "_zi0.5_", alt, "_top100"),
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
      "results/cell_type/heatmap/hm_CB_counts_QC_pDEA_cell_type_",
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
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_pDEA_cell_type_",
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
rm(list=ls())
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)

day <- "InVitro"
experiment <- "P1902"
load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
load(paste0("results/cell_type/mbatch_", day, "_", experiment, "_pfounder_cell_type_DEA.Rdata"))
b_genes <- !is.na(mbatch_pfounder_cell_type_DEA$padj) &
  mbatch_pfounder_cell_type_DEA$padj < 0.05
DEA_genes <- mbatch_pfounder_cell_type_DEA$gene[b_genes]
day <- "InVitro"
experiment <- "P3128"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1
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
PLS_genes <- scd$getgenes[scd$getgenes %in%
                          unique(c(DEA_genes, genes_marker, "KLRD1", "SELL"))]
good_pls <- apply(scd$select(b_cells = b_cells,
                             genes = PLS_genes)$getcounts,
                  2, FUN = function(x){
                    if(max(x) == 0) {
                      return(FALSE)
                    }
                    return(median(x[x > 0]) > 10)
                  })
PLS_genes <- PLS_genes[good_pls]

founder_phenotype <- scd$getfeature("founder_phenotype")
founder_phenotype[
  b_cells & founder_phenotype %in% "UNK"] <- NA
scd$setfeature("founder_phenotype", founder_phenotype)
summary(as.factor(founder_phenotype[b_cells]))

system("rm results/cell_type/DEA_cell_types_invitro_P3128_denovo*")
system("cp results/cell_type/DEA_cell_types_invitro_P1902_denovo_training_lsplsstab.Rdata results/cell_type/DEA_cell_types_invitro_P3128_denovo_training_lsplsstab.Rdata")
system("cp results/cell_type/DEA_cell_types_invitro_P1902_denovo_classification_lplscv.Rdata results/cell_type/DEA_cell_types_invitro_P3128_denovo_classification_lplscv.Rdata")
system("rm results/cell_type/DEA_cell_types_invitro_P3128_denovo_classification_lpls.Rdata")
load("results/cell_type/DEA_cell_types_splsstab_invitro_P1902_denovo.Rdata")

day <- "InVitro"
experiment <- c("P1902", "P3128")
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "founder_cell_type",
  features = c(),
  genes = c(genes_marker, DEA_genes),
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_invitro_P3128_denovo",
  force =  genes_marker
)
system("~/scripts/sms.sh \"PLS done\"")
traceback()


save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_splsstab_invitro_P3128_denovo.Rdata"
)

devtools::load_all("../scRNAtools/", reset = T)
b_cells <- scd$getfeature("QC_good") %in% T
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "founder_phenotype",
  features = c(),
  genes = DEA_cell_type_classification$classification$fit_spls$fit$selected,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_weighted_force_invitro_P3128_denovo",
  force = DEA_cell_type_classification$classification$fit_spls$fit$selected
)
system("~/scripts/sms.sh \"PLS done\"")

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_splsstab_invitro_P3128_denovo.Rdata"
)

load("results/cell_type/DEA_cell_types_splsstab_invitro_P3128_denovo.Rdata")

load("results/cell_type/DEA_cell_types_splsstab_invitro_P1902_denovo_noforce.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "founder_phenotype",
  features = c(),
  genes = DEA_cell_type_classification$classification$fit_spls$fit$selected,
  ncores = 6,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types_weighted_force_invitro_P3128_denovo",
  force = c()
)
system("~/scripts/sms.sh \"PLS done\"")

save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_splsstab_invitro_P3128_denovo_noforce.Rdata"
)

day <- "InVitro"
experiment <- c("P1902", "P3128")
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
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
day <- "InVitro"
experiment <- "P3128"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1
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
  tmp_infos$pDEA_clone_cell_type > 0.5,
  "MEM",
  "EFF"
)

tmp_infos$KLRD1 <- scd$select(b_cells = b_cells)$getgene("KLRD1")
tmp_infos$KLRD1_clone <- unlist(as.list(by(
  log10(tmp_infos$KLRD1+1), tmp_infos$clonality, median
))[tmp_infos$clonality])
tmp_infos$SELL <- scd$select(b_cells = b_cells)$getgene("SELL")
tmp_infos$SELL_clone <- unlist(as.list(by(
  log10(tmp_infos$SELL+1), tmp_infos$clonality, median
))[tmp_infos$clonality])
tmp_infos$MKI67 <- scd$select(b_cells = b_cells)$getgene("MKI67")

table(tmp_infos$DEA_clone_cell_type)

g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log(cycling_score),
      color = DEA_clone_cell_type)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(aes(color = DEA_cell_type)) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$DEA_cell_type)))) +
  labs(x = "clones",
       y =  "cell-cycle score",
       color = "founder_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P3128_cycling_vs_DEA_cell_type.pdf"
))

clone_name <- as.vector(tmp_infos$clonality)[order(tmp_infos$SELL_clone)]
tmp_infos$clonality <- factor(tmp_infos$clonality, levels = unique(clone_name))
g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log10(SELL+1),
      color = DEA_clone_cell_type)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(aes(color = DEA_cell_type)) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$DEA_cell_type)))) +
  labs(x = "clones",
       y =  "log10(SELL+1)",
       color = "founder_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P3128_SELL_vs_DEA_cell_type.pdf"
))

g <- ggplot(tmp_infos,
   aes(x = clonality,
      y = log10(KLRD1+1),
      color = DEA_clone_cell_type)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(aes(color = DEA_cell_type)) +
  theme_bw() +
  scale_color_manual(
    values = cell_type_palette(
      levels(factorize(tmp_infos$DEA_cell_type)))) +
  labs(x = "clones",
       y =  "log10(SELL+1)",
       color = "founder_cell_type",
       title = day)
print(g)
ggsave(file = paste0(
  "results/cycling/violing_invitro_P3128_KLRD1_vs_DEA_cell_type.pdf"
))

infos_M <- scd$getfeatures
write.csv(
  infos_M,
  file = paste0("results/cell_type/cell_type_infos_invivo_invitro.csv")
)


# DEA analysis InVitro P3128

day <- "InVitro"
experiment <- "P3128"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P3128.Rdata")
system(
  paste0("rm -R results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA")
)
system(
  paste0("mkdir -p results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA")
)
mbatch_pDEA_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + pDEA_cell_type",
  b_cells = b_cells,
  cpus = 10,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA")
)
save(
  mbatch_pDEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_pDEA_cell_type_DEA$padj)))
print(table(mbatch_pDEA_cell_type_DEA$padj < 0.05))
write.csv(
  mbatch_pDEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_pDEA_cell_type_DEA.csv")
)


load(
  paste0("results/cell_type/mbatch_", day, "_", experiment, "_DEA_cell_type_DEA.Rdata")
)

load("results/cycling/CB_counts_QC_cycling_invitro_P1902_P3128.Rdata")
day <- "InVitro"
experiment <- "P3128"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

for (alt in c("lesser", "greater")) {
  b_genes <- !is.na(mbatch_pDEA_cell_type_DEA$padj) &
    mbatch_pDEA_cell_type_DEA$padj < 0.05
  DEA_genes <- mbatch_pDEA_cell_type_DEA$gene
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
      day, "_", experiment, "_zi0.5_", alt),
    top = min(100, length(DEA_genes_exp))
  )

  devtools::load_all("../scRNAtools/", reset = T)
  system(paste0("rm results/cell_type/gene_cov_",
      day, "_", experiment, "_zi0.5_", alt, "_top100*"))
  gene_order <- order_by_factor(
    scd = scd_norm,
    score = scd_norm$getfeature("pDEA_cell_type"),
    tmp_file = paste0("results/cell_type/gene_cov_",
      day, "_", experiment, "_zi0.5_", alt, "_top100"),
    top = min(100, length(DEA_genes_exp))
  )
  traceback()

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
      "results/cell_type/heatmap/hm_CB_counts_QC_pDEA_cell_type_",
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
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_pDEA_cell_type_",
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

load("results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

scd_data <- scd$select(b_cells = b_cells)$getfeatures
data <- data.frame(
  ccr7 = as.numeric(as.vector(scd_data$ccr7)),
  il7ra = as.numeric(as.vector(scd_data$il7ra)),
  pMEM = scd_data$pDEA_cell_type,
  day = factor(scd_data$day, levels = c("D15", "D90")),
  classification = "PLS from DEA"
)

ggplot(data = data,
       aes(x = classification, y = ccr7, color = pMEM)) +
  scale_y_log10() +
  geom_jitter(height = 0) +
  facet_wrap(~day) +
  theme_bw()
ggsave("results/cell_type/P1902_in_vitro_ccr7_pMEM.pdf",
       width = 11,
       height = 10)

ggplot(data = data,
       aes(x = classification, y = il7ra, color = pMEM)) +
  scale_y_log10() +
  geom_jitter(height = 0) +
  facet_wrap(~day) +
  theme_bw()
ggsave("results/cell_type/P1902_in_vitro_il7ra_pMEM.pdf",
       width = 11,
       height = 10)

###############################################################################
############################ Clonality analysis################################

# P1902
load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
genes_to_rm <- read.table("data/Genes_exclude.csv", h = T)
scd <- scd$select(genes = scd$getgenes[!scd$getgenes %in% genes_to_rm])
day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1
table(scd$select(b_cells = b_cells)$getfeature("batch"))
table(as.factor(as.vector(scd$select(b_cells = b_cells)$getfeature("clonality"))))

system(
  paste0("rm -R results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA")
)
system(
  paste0("mkdir -p results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA")
)
mbatch_clonality_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ batch",
  formula_full = "y ~ batch + (1|clonality)",
  b_cells = b_cells,
  cpus = 8,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA")
)
save(
  mbatch_clonality_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_clonality_DEA$padj)))
print(table(mbatch_clonality_DEA$padj < 0.05))
write.csv(
  mbatch_clonality_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA.csv")
)

# P3128
load("results/cell_type/CB_counts_QC_DEA_cell_type_invitro_P3128.Rdata")
scd <- scd$select(genes = scd$getgenes[!scd$getgenes %in% genes_to_rm])
day <- "InVitro"
experiment <- "P3128"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1
table(scd$select(b_cells = b_cells)$getfeature("batch"))
table(as.factor(as.vector(scd$select(b_cells = b_cells)$getfeature("clonality"))))

system(
  paste0("rm -R results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA")
)
system(
  paste0("mkdir -p results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA")
)
mbatch_clonality_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + (1|clonality)",
  b_cells = b_cells,
  cpus = 8,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA")
)
save(
  mbatch_clonality_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_clonality_DEA$padj)))
print(table(mbatch_clonality_DEA$padj < 0.05))
write.csv(
  mbatch_clonality_DEA,
  file = paste0("results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA.csv")
)

##################################### Normalization ###########################

rm(list = ls())
library(tidyverse)
library(BiocParallel)
library(zinbwave)
devtools::load_all("../scRNAtools/", reset = T)

scd_sumarizedexp <- function(scd, b_cells, genes = scd$getgenes) {
  assay_data <- scd$select(b_cells = b_cells, genes = genes)$getcounts %>%
    as.matrix() %>%
    round() %>%
    t()
  col_data <- scd$select(b_cells = b_cells)$getfeatures %>%
    dplyr::mutate(batch = as.factor(as.vector(batch)),
          cycling = as.factor(as.vector(cycling)),
          DEA_cell_type = as.factor(as.vector(DEA_cell_type)),
          clonality = as.factor(as.vector(clonality))
    )
  data <- SummarizedExperiment(
    assays = list(counts = assay_data),
    colData = col_data
  )
  filter <- rowSums(assay(data)>5)>5
  data <- data[filter, ]
  return(data)
}

norm_sequence <- function(data,
                          norm_by = c("batch", "cycling_score", "pDEA_cell_type"),
                          K = 0) {
  norm_data <- list()
  X_factors <- "~ "
  for (factor in norm_by) {
    X_factors <- paste0(X_factors, factor)
    message(paste0("normalizing for ", X_factors))
    zinb_res <- zinbwave(data,
                        K = K,
                        X = X_factors,
                        normalizedValues=TRUE,
                        residuals = TRUE,
                        BPPARAM=MulticoreParam(12))
    norm_data[[factor]] <- assay(zinb_res, "normalizedValues")
    norm_data[[paste0(factor, "_model")]] <- zinb_res
    X_factors <- paste0(X_factors, " + ")
  }
  return(norm_data)
}

# p1902
load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
genes_to_rm <- read.table("data/Genes_exclude.csv", h = T)
scd <- scd$select(genes = scd$getgenes[!scd$getgenes %in% genes_to_rm])
day <- "InVitro"
experiment <- "P1902"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

data <- scd_sumarizedexp(scd, b_cells)
norm_data <- norm_sequence(data)
save(norm_data,
     file = "results/cell_type/counts_invitro_P1902_zinbwave_norm.Rdata"
)

# pca
data.frame(prcomp(t(norm_data[["batch_cycling_score_pDEA_cell_type"]]))$x[, 1:2],
           cell_type=colData(data)$DEA_cell_type,
           pcell_type=colData(data)$pDEA_cell_type,
           cycling=colData(data)$cycling,
           cycling_score=colData(data)$cycling_score,
           clonality=colData(data)$clonality) %>%
    ggplot(aes(PC1, PC2, colour=pcell_type, shape=cycling)) +
    geom_point() +
    theme_classic()

# p1902 DEG
load(
  paste0("results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA.Rdata")
)
b_genes <- !is.na(mbatch_clonality_DEA$padj) &
  mbatch_clonality_DEA$padj < 0.05
DEA_genes <- mbatch_clonality_DEA$gene

data <- scd_sumarizedexp(scd, b_cells, DEA_genes)
norm_data <- norm_sequence(data, c("batch", "cycling_score"))
save(norm_data,
     file = "results/cell_type/counts_invitro_P1902_zinbwave_norm_clonality_DEG.Rdata"
)
for (norm_type in c("batch", "cycling_score")) {
  message(norm_type)
  write.csv(
    norm_data[[norm_type]],
    file = paste0("results/cell_type/norm_counts_",
                  day ,"_", experiment, "_", norm_type, ".csv")
  )
}

# pca
x11()
data.frame(prcomp(t(norm_data[["batch"]]))$x[, 1:2],
           batch=colData(data)$batch,
           clonality=colData(data)$clonality,
           cycling=colData(data)$cycling,
           cycling_score=colData(data)$cycling_score,
           clonality=colData(data)$clonality) %>%
    ggplot(aes(PC1, PC2, colour=as.factor(batch), shape=as.factor(batch))) +
    geom_point() +
    theme_classic()

# P3128
load("results/cell_type/cells_counts_QC_DEA_cell_type_invitro_P1902.Rdata")
genes_to_rm <- read.table("data/Genes_exclude.csv", h = T)
scd <- scd$select(genes = scd$getgenes[!scd$getgenes %in% genes_to_rm])
day <- "InVitro"
experiment <- "P3128"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1

data <- scd_sumarizedexp(scd, b_cells, DEA_genes)
norm_data <- norm_sequence(data, c("batch", "cycling_score", "pDEA_cell_type"))
save(norm_data,
     file = "results/cell_type/counts_invitro_P3128_zinbwave_norm.Rdata"
)

# pca
data.frame(prcomp(t(norm_data[["batch_cycling_score_pDEA_cell_type"]]))$x[, 1:2],
           cell_type=colData(data)$DEA_cell_type,
           pcell_type=colData(data)$pDEA_cell_type,
           cycling=colData(data)$cycling,
           cycling_score=colData(data)$cycling_score,
           clonality=colData(data)$clonality) %>%
    ggplot(aes(PC1, PC2, colour=pcell_type, shape=cycling)) +
    geom_point() +
    theme_classic()

# p3128 DEG
load(
  paste0("results/cell_type/mbatch_", day, "_", experiment, "_clonality_DEA.Rdata")
)
b_genes <- !is.na(mbatch_clonality_DEA$padj) &
  mbatch_clonality_DEA$padj < 0.05
DEA_genes <- mbatch_clonality_DEA$gene

data <- scd_sumarizedexp(scd, b_cells, DEA_genes)
norm_data <- norm_sequence(data, c("batch", "cycling_score"))
save(norm_data,
     file = "results/cell_type/counts_invitro_P3128_zinbwave_norm_clonality_DEG.Rdata"
)
for (norm_type in c("batch", "cycling_score")) {
  message(norm_type)
  write.csv(
    norm_data[[norm_type]],
    file = paste0("results/cell_type/norm_counts_",
                  day ,"_", experiment, "_", norm_type, ".csv")
  )
}

scd$select(b_cells = b_cells)$getgene("RP9") %>%
  as.tibble() %>%
  dplyr::rename(RP9 = value) %>%
  dplyr::bind_cols( norm_data[["batch"]] %>%
                    t() %>%
                    as.tibble() %>%
                    dplyr::select(RP9) %>%
                    dplyr::rename(RP9_batch = RP9)
  ) %>%
  dplyr::bind_cols( norm_data[["cycling_score"]] %>%
                    t() %>%
                    as.tibble() %>%
                    dplyr::select(RP9) %>%
                    dplyr::rename(RP9_batch_cycling = RP9)
  ) %>%
  ggplot() +
  # geom_histogram(aes(x = RP9, fill = "red", alpha = 0.5), bins = 100) +
  geom_histogram(aes(x = RP9_batch, fill = "blue", alpha = 0.5), bins = 100) +
  geom_histogram(aes(x = RP9_batch_cycling, fill = "green", alpha = 0.5), bins = 100) +
  theme_classic()

scd$select(b_cells = b_cells)$getgene("RP9") %>%
  as.tibble() %>%
  dplyr::rename(RP9 = value) %>%
  dplyr::bind_cols( norm_data[["batch"]] %>%
                    t() %>%
                    as.tibble() %>%
                    dplyr::select(RP9) %>%
                    dplyr::rename(RP9_batch = RP9)
  ) %>%
  dplyr::bind_cols( norm_data[["cycling_score"]] %>%
                    t() %>%
                    as.tibble() %>%
                    dplyr::select(RP9) %>%
                    dplyr::rename(RP9_batch_cycling = RP9)
  ) %>%
  dplyr::mutate( RP9_batch = exp(RP9_batch + abs(min(RP9_batch))),
                 RP9_batch_cycling = exp(RP9_batch_cycling + abs(min(RP9_batch_cycling)))
  ) %>%
  ggplot() +
  geom_histogram(aes(x = RP9, fill = "red", alpha = 0.5), bins = 100) +
  geom_histogram(aes(x = RP9_batch, fill = "blue", alpha = 0.5), bins = 100) +
  geom_histogram(aes(x = RP9_batch_cycling, fill = "green", alpha = 0.5), bins = 100) +
  theme_classic()

