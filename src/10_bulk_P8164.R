setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("pkg/", reset = T)
genes_to_excludes <- read.table("data/P8164_Bulk_genes_to_exlude.csv", h = T)
genes_to_excludes <- genes_to_excludes$genes_to_exclude

results_folder <- "results/P8164_bulk/"
system(paste0("mkdir -p ", results_folder))

system("perl -pi -e 's/(P\\d*_\\d*)_\\S*/\\1/g' data/P8164_Bulk/*.tsv")

bcd <- scRNAtools::load_data(
  infos = "data/P8164_Bulk.csv",
  counts = "data/P8164_Bulk/",
  regexp = ".*counts_genes.*",
  infos_sep = ","
)
save(bcd, file = paste0(results_folder, "P8164_Bulk_raw_counts.Rdata"))
load(paste0(results_folder, "P8164_Bulk_raw_counts.Rdata"))

################################################################################
# analysis with ziNB
b_cells <- bcd$getfeature("cell_number") %in% 20 &
  bcd$getfeature("sex") %in% "M"

system("mkdir -p results/P8164_bulk/clonality_DEA")
clonality_DEA <- DEA(
  scd = bcd,
  formula_null = "y ~ 1",
  formula_full = "y ~ clonality",
  b_cells = b_cells,
  cpus = 10,
  v = F,
  folder_name = "results/P8164_bulk/clonality_DEA",
  zi_threshold = 0.99
)
save(
  clonality_DEA,
  file = "results/P8164_bulk/clonality_DEA.Rdata"
)
table(is.na(clonality_DEA$padj))
table(clonality_DEA$padj < 0.05)
write.csv(clonality_DEA, file = "results/P8164_bulk/clonality_DEA.csv")

################################################################################
# analysis with DESeq2

library(DESeq2)

b_cells <- bcd$getfeature("cell_number") %in% 20 &
  bcd$getfeature("sex") %in% "M"
countData <- t(round(
    bcd$select(
      b_cells = b_cells,
      genes = ERCC(bcd, minus = T))$getcounts))
colData <- bcd$select(b_cells = b_cells)$getfeatures

head(colData)

dds <- DESeqDataSetFromMatrix(
  countData = countData + 1,
  colData = colData,
  design = ~ clonality)
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

save(res, file = paste0(results_folder, "P8164_DEA.Rdata"))
load(file = paste0(results_folder, "P8164_DEA.Rdata"))

DE_genes <- res[
    res[["padj"]] < 0.05 &
    !is.na(res[["padj"]]) &
    !(rownames(res) %in% genes_to_excludes)
  , ]

write.table(DE_genes, file = paste0(results_folder, "/DE_genes.csv"))

sum(res$padj < 0.05, na.rm=TRUE)

DE_genes <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]
DE_genes <- DE_genes[!(DE_genes %in% genes_to_excludes)]

devtools::load_all("pkg/", reset = T)
pca_plot(bcd$select(b_cells = b_cells, genes = DE_genes), color="clonality",
  color_name = "clonality", tmp_file = "results/tmp_bulk_pca.Rdata")
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

# check IFNG genes because it's important
norm_counts <- t(counts(dds, normalized = TRUE))
ggplot(
  data = data.frame(
    IFNG = norm_counts[, colnames(norm_counts) %in% "IFNG"],
    clonality = colData$clonality
  ),
  aes(x = clonality, y = IFNG)) +
  geom_point() +
  theme_bw()

# 1 vs all comparison
res_1vall <- list()
colData <- bcd$select(b_cells = b_cells)$getfeatures

for (clone in
  levels(factorize(bcd$select(b_cells = b_cells)$getfeature("clonality")))) {
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

save(res_1vall, file = paste0(results_folder, "P8164_DEA_1vall.Rdata"))
load(file = paste0(results_folder, "P8164_DEA_1vall.Rdata"))


DE_genes <- data.frame()
for (clone in
  levels(factorize(bcd$select(b_cells = b_cells)$getfeature("clonality")))) {
  print(clone)
  DE_genes_clone <- res_1vall[[clone]][
      res_1vall[[clone]][["padj"]] < 0.05 &
      !is.na(res_1vall[[clone]][["padj"]]) &
      !(rownames(res_1vall[[clone]]) %in% genes_to_excludes)
    , ]
  DE_genes_clone$clone <- clone
  if (nrow(DE_genes) != 0) {
    DE_genes <- rbind(
      DE_genes,
      DE_genes_clone
    )
  } else {
    DE_genes <- DE_genes_clone
  }
}
write.table(DE_genes, file = paste0(results_folder, "/DE_genes_1vsall.csv"))


################################################################################
# Analysis with Edger

library(edgeR)
devtools::load_all("pkg/", reset = T)
group <- bcd$getfeature("clonality")
y <- DGEList(
  counts = t(bcd$select(genes = ERCC(bcd, minus = TRUE))$getcounts),
  group = group
)
y <- calcNormFactors(y)
design <- model.matrix(~0+group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)
fit$design
qlf <- glmQLFTest(fit, coef = 1:length(levels(group)))
qlf_padj <- p.adjust(qlf$table$PValue, method = "BH")
sum(qlf_padj < 0.05)
qlf$table[qlf_padj < 0.05, ]


fit <- glmFit(y, design)
lrt <- glmLRT(fit)
lrt_padj <- p.adjust(lrt$table$PValue, method = "BH")
sum(lrt_padj < 0.05)
lrt$table[lrt_padj < 0.05, ]
