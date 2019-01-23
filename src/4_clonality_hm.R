rm(list = ls())
setwd("~/projects/yellow_fever")
devtools::load_all("pkg/", reset = T)
load("results/cycling/cells_counts_QC_cycling.Rdata")

system("mkdir -p results/clonality/")

for (sex in c("M", "F")) {
  days <- c("D15", "D136", "D593")
  if (sex %in% "F") {
    days <- c("D15", "D90")
  }

  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% days[2:length(days)] &
    scd$getfeature("sex") %in% sex

  clone_size <- table(scd$select(b_cells = b_cells)$getfeature("clonality"))
  clone_size <- clone_size[order(clone_size, decreasing = TRUE)]
  clone_size <- clone_size[clone_size > 3][-1]
  clone_list <- names(clone_size)

  b_cells <- b_cells & scd$getfeature("clonality") %in% clone_list

  cells_order <- order_by_groups(
    scd$select(b_cells = b_cells)$getfeature("pDEA_cell_type"),
    by = scd$select(b_cells = b_cells)$getfeature("clonality")
  )

  genes_list <- c("GZMB", "CX3CR1", "CCL4", "GNLY", "GZMH", "KLRD1", "GZMA",
    "PRF1", "HOPX", "CCL5", "GZMK", "SELL", "IL7R", "LEF1", "TCF7", "LTB",
    "NELL2", "CCR7")
  scd_norm <- zinorm(
    scd = scd$select(
      b_cells = b_cells,
      genes = genes_list
    ),
    cpus = 10,
    file = paste0("results/tmp/zi_norm_cells_counts_D100_clone_size_4", sex, ".RData")
  )
  scd_norm$setfeature("clonality",
                      as.factor(as.vector(scd_norm$getfeature("clonality"))))

  hm <- heatmap_genes(
    scd = scd_norm,
    features = c("clonality", "pDEA_cell_type"),
    cells_order = cells_order,
    genes_order = genes_list,
    title = "D100+ clones of size >= 4",
    factor = c(T, F),
    file = paste0( "results/clonality/hm_D100_clone_size_4", sex, ".pdf" ),
    gene_size = 10
  )
  print(hm)
}
