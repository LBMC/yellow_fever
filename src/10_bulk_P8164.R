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

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ clonality)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)

save(res, file = "results/P8164_DEA.Rdata")
load(file = "results/P8164_DEA.Rdata")

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

# 1 vs all comparison
res_1vall <- list()
colData <- bcd$select(b_cells = b_cells)$getfeatures

for (clone in levels(factorize(bcd$select(b_cells = b_cells)$getfeature("clonality")))) {

  print(clone)
  colData <- cbind(
    colData,
    as.factor(
      ifelse(colData$clonality %in% clone, clone, "others")
    )
  )
  colnames(colData) <- c(
    colnames(colData)[1:ncol(colData)-1],
    gsub(" ", "_", clone)
  )
  dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = as.formula(paste("~", gsub(" ", "_", clone))))
  dds <- dds[ rowSums(counts(dds)) > 1,   ]
  dds <- DESeq(dds, test="LRT", reduced=~1)
  res_1vall[[clone]] <- results(dds)

}

length(res_1vall)

save(res_1vall, file = "results/P8164_DEA_1vall.Rdata")
load(file = "results/P8164_DEA_1vall.Rdata")
