setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("../scRNAtools/", reset = T)

results_folder <- "results/P8164_bulk/"
system(paste0("mkdir -p ", results_folder))

system("perl -pi -e 's/(P\\d*_\\d*)_\\S*/\\1/g' data/P8164_Bulk/*.tsv")

bcd <- scRNAtools::load_data(
  infos = "data/P8164_Bulk.csv",
  counts = "data/P8164_Bulk/",
  regexp = ".*counts_genes.*",
  infos_sep = ","
)
save(bcd, file = "results/P8164_Bulk_raw_counts.Rdata")
load("results/P8164_Bulk_raw_counts.Rdata")

library(DESeq2)

b_cells <- bcd$getfeature("cell_number") %in% 20 &
  bcd$getfeature("sex") %in% "M"
countData <- t(round(
    bcd$select(
      b_cells = b_cells,
      genes = ERCC(bcd, minus = T))$getcounts))
colData <- bcd$select(b_cells = b_cells)$getfeatures

table(as.factor(as.vector(bcd$select(b_cells = b_cells)$getfeature("clonality"))))

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ clonality)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)

write.table(res, file = paste0(results_folder, "/DE_genes.csv"))

sum(res$padj < 0.05, na.rm=TRUE)

DE_genes <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]

genes_to_excludes <- read.table("data/P8164_Bulk_genes_to_exlude.csv", h = T)
DE_genes <- DE_genes[!(DE_genes %in% genes_to_excludes)]

devtools::load_all("../scRNAtools/", reset = T)
pca_plot(bcd$select(b_cells = b_cells, genes = DE_genes), color="clonality",
  color_name = "clonality")
ggsave(
  paste0(results_folder, "pca_DE_genes.pdf"),
  width = 20, height = 15, units = "cm", dpi = 1200
)
bca_plot(bcd$select(b_cells = b_cells, genes = DE_genes),
  by = bcd$select(b_cells = b_cells, genes = DE_genes)$getfeature("clonality"))
ggsave(
  paste0(results_folder, "bca_DE_genes.pdf"),
  width = 20, height = 15, units = "cm", dpi = 1200
)
