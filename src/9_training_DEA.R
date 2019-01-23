rm(list=ls())
setwd("~/projects/yellow_fever/")
devtools::load_all("pkg/", reset = T)
load("results/QC/cells_counts_QC_training.Rdata")

experiment <- "training"
b_cells <- scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1 &
  scd$getfeature("QC_good") %in% T
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
  file = paste0("results/cell_type/cells_counts_training.csv")
)

infos <- scd$getfeatures
write.csv(
  infos,
  file = paste0("results/cell_type/nfos_training.csv")
)

# single-cells
experiment <- "training"
b_cells <- scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1 &
  scd$getfeature("QC_good") %in% T
DEA_cell_type <- scd$getfeature("phenotype_surface_marker")
levels(DEA_cell_type) <- c(NA, "CM", "CM", "EFF", "EM", "EM", "MEM", "Naive",
                           "Temra", "Temra", "Temra", "TSCM", "TSCM")
DEA_cell_type[b_cells & !( DEA_cell_type %in% c("CM", "EM") )] <- NA
scd$setfeature("DEA_cell_type", DEA_cell_type)

experiment <- "training"
b_cells <- scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1 &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("DEA_cell_type") %in% c("CM", "EM")
system(
  paste0("rm -R results/cell_type/mbatch_", experiment, "_DEA_cell_type_DEA")
)
system(
  paste0("mkdir -p results/cell_type/mbatch_", experiment, "_DEA_cell_type_DEA")
)
mbatch_DEA_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ 1",
  formula_full = "y ~ DEA_cell_type",
  b_cells = b_cells,
  cpus = 10,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_", experiment, "_DEA_cell_type_DEA")
)
save(
  mbatch_DEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", experiment, "_DEA_cell_type_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
write.csv(
  mbatch_DEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_", experiment,
                "_DEA_cell_type_DEA.csv")
)

# bulk
experiment <- "training"
b_cells <- scd$getfeature("experiment") %in% experiment &
  !( scd$getfeature("cell_number") %in% 1  )
DEA_cell_type <- scd$getfeature("phenotype_surface_marker")
levels(DEA_cell_type) <- c(NA, "CM", "CM", "EFF", "EM", "EM", "MEM", "Naive",
                           "Temra", "Temra", "Temra", "TSCM", "TSCM")
DEA_cell_type[b_cells & !( DEA_cell_type %in% c("CM", "EM") )] <- NA
scd$setfeature("DEA_cell_type", DEA_cell_type)

experiment <- "training"
b_cells <- scd$getfeature("experiment") %in% experiment &
  !( scd$getfeature("cell_number") %in% 1  )&
  scd$getfeature("DEA_cell_type") %in% c("CM", "EM")
table(b_cells)

library(DESeq2)
countData <- t(round(
    scd$select(
      b_cells = b_cells,
      genes = ERCC(scd, minus = T))$getcounts))
colData <- scd$select(b_cells = b_cells)$getfeatures

dds <- DESeqDataSetFromMatrix(
  countData = countData + 1,
  colData = colData,
  design = ~ DEA_cell_type)
dds <- dds[ rowSums(counts(dds)) > 1, ]
geoMean <- apply(
  counts(dds),
  1,
  FUN = function(x){
    if(all(x == 0)) {
      0
    } else {
      exp(mean(log(x[x !=0 ])))
    }
  }
)
dds <- estimateSizeFactors(dds, geoMean = geoMean)
dds <- estimateDispersions(dds, fitType = "local")
dds <- nbinomLRT(dds, reduced = ~1)
res <- results(dds)

save(res, file = "results/cell_type/bulk_training_cell_type.Rdata")
write.table(res, file = "results/cell_type/bulk_training_cell_type.csv")

sum(res$padj < 0.05, na.rm=TRUE)

DE_genes <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]

pca_plot(scd$select(b_cells = b_cells, genes = DE_genes), color="DEA_cell_type",
  color_name = "clonality", tmp_file = "results/tmp_bulk_training_pca.Rdata")
ggsave(
  "results/cell_type/pca/pca_bulk_training_cell_type.pdf",
  width = 20, height = 15, units = "cm", dpi = 1200
)

# compaison
load(
  file = paste0("results/cell_type/mbatch_", experiment,
                "_DEA_cell_type_DEA.Rdata")
)
save(file = "results/cell_type/bulk_training_cell_type.Rdata")
DEG_single <- mbatch_DEA_cell_type_DEA$gene[mbatch_DEA_cell_type_DEA$padj < 0.05 &
                             !is.na( mbatch_DEA_cell_type_DEA$padj )]
DEG_bulk <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]
table(DEG_single %in% DEG_bulk)
table(DEG_bulk %in% DEG_single)

b_cells <- scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1 &
  scd$getfeature("QC_good") %in% T
DEA_cell_type <- scd$getfeature("phenotype_surface_marker")
levels(DEA_cell_type) <- c(NA, "CM", "CM", "EFF", "EM", "EM", "MEM", "Naive",
                           "Temra", "Temra", "Temra", "TSCM", "TSCM")
DEA_cell_type[b_cells & !( DEA_cell_type %in% c("CM", "EM") )] <- NA
scd$setfeature("DEA_cell_type", DEA_cell_type)
experiment <- "training"
b_cells <- scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1 &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("DEA_cell_type") %in% c("CM", "EM")

pca_plot(scd$select(b_cells = b_cells, genes = DEG_bulk), color="DEA_cell_type",
  color_name = "clonality", tmp_file = "results/tmp_scd_training_pca.Rdata")
ggsave(
  "results/cell_type/pca/pca_sc_training_cell_type.pdf",
  width = 20, height = 15, units = "cm", dpi = 1200
)
