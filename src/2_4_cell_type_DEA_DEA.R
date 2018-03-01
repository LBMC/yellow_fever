setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")

for (day in c("D15", "D136", "D593")) {
  system(
    paste0("mkdir -p results/cell_type/mbatch_", day, "_DEA_cell_type_DEA")
  )
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% day
  mbatch_DEA_cell_type_DEA <- DEA(
    scd = scd,
    formula_null = "y ~ (1|batch)",
    formula_full = "y ~ (1|batch) + DEA_cell_type",
    b_cells = b_cells,
    cpus = 16,
    v = F,
    folder_name = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA")
  )
  save(
    mbatch_DEA_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata")
  )
  system("~/scripts/sms.sh \"DEA done\"")
  print(day)
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
}

system("mkdir -p results/cell_type/heatmap/")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
DEA_cell_type_palette <- cell_type_palette
DEA_genes_M <- c()
for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
  print("GNLY" %in% DEA_genes)
  print("GZMB" %in% DEA_genes)
  print("GZMH" %in% DEA_genes)
  DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
  DEA_gene <- c( DEA_gene, "GNLY", "GZMB", "GZMH")
  DEA_genes_M <- unique(c(DEA_genes, DEA_genes_M))
  DEA_genes_thresholded <- expressed(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day &
        scd$getfeature("DEA_cell_type") %in% "EFF",
      genes = DEA_genes
    ),
    zi_threshold = 0.40
  )
  DEA_genes_thresholded <- unique(c(DEA_genes_thresholded, expressed(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day &
        scd$getfeature("DEA_cell_type") %in% "MEM",
      genes = DEA_genes
    ),
    zi_threshold = 0.30
  )))
  system(paste0("rm results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
      day, ".Rdata"))
  scd_norm <- zinorm(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day &
        scd$getfeature("DEA_cell_type") %in% "MEM",
      genes = DEA_genes_thresholded
    ),
    cpus = 10,
    file = paste0("results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
      day, ".Rdata")
  )
  hm <- heatmap_genes(
    scd = scd_norm,
    features = c("DEA_cell_type", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    genes_order = order_2_groups(
      scd = scd_norm,
      by = scd_norm$getfeature("DEA_cell_type")
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_CB_counts_QC_DEA_cell_type_",
      day, ".pdf"
    ),
    gene_size = 3                      
  )
  print(hm)

  hm_corr <- heatmap_corr_genes(
    scd = scd_norm,
    features = c("DEA_cell_type", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, ".pdf"
    )
  )
  print(hm_corr)
}

system("mkdir -p results/cell_type/pca")
system("mkdir -p results/cell_type/pcmf")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
DEA_cell_type_palette <- cell_type_palette
for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
  system(paste0("rm results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
      day, ".RData"))
  scd_norm <- zinorm(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    cpus = 10,
    file = paste0("results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
      day, ".RData")
  )
  system(paste0("rm results/tmp/pca_zi_norm_counts_QC_DEA_cell_type_", day, ".Rdata"))
  scRNAtools::pca_plot(
    scd_norm,
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_zi_norm_counts_QC_DEA_cell_type_", day, ".Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_zi_norm_counts_QC_DEA_cell_type_pca", day, ".pdf"
  ))
  system(paste0("rm results/tmp/pCMF_CB_counts_QC_DEA_cell_type_", day, ".Rdata"))
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_QC_DEA_cell_type_", day, ".Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pCMF_CB_counts_QC_DEA_cell_type_", day, ".pdf"
  ))
}

################################################################################
## DEA for the F donor

load("results/cell_type/cells_counts_QC_DEA_cell_type_F.Rdata")

for (day in c("D15", "D90")) {
  system(
    paste0("mkdir -p results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F")
  )
  b_cells <- scd$getfeature("QC_good") %in% T &
    scd$getfeature("sex") %in% "F" &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% day
  mbatch_DEA_cell_type_DEA <- DEA(
    scd = scd,
    formula_null = "y ~ (1|batch)",
    formula_full = "y ~ (1|batch) + DEA_cell_type",
    b_cells = b_cells,
    cpus = 16,
    v = F,
    folder_name = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F")
  )
  save(
    mbatch_DEA_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.Rdata")
  )
  system("~/scripts/sms.sh \"DEA done\"")
  print(day)
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
}

for (day in c("D15", "D90")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.Rdata"))
  write.csv(
    mbatch_DEA_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.csv")
  )
}

system("mkdir -p results/cell_type/heatmap/")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "F"
DEA_cell_type_palette <- cell_type_palette
DEA_genes_F <- c()
for (day in c("D15", "D90")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.Rdata"))
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
  DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
  DEA_genes_F <- unique(c(DEA_genes, DEA_genes_F))
  DEA_genes_thresholded <- expressed(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day &
        scd$getfeature("DEA_cell_type") %in% "EFF",
      genes = DEA_genes
    ),
    zi_threshold = 0.40
  )
  DEA_genes_thresholded <- unique(c(DEA_genes_thresholded, expressed(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day &
        scd$getfeature("DEA_cell_type") %in% "MEM",
      genes = DEA_genes
    ),
    zi_threshold = 0.30
  )))
  system(paste0("rm results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
    day, "_thresholded_F.RData"))
  scd_norm <- zinorm(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes_thresholded
    ),
    cpus = 10,
    file = paste0("results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
      day, "_thresholded_F.RData")
  )
  hm <- heatmap_genes(
    scd = scd_norm,
    features = c("DEA_cell_type", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    genes_order = order_2_groups(
      scd = scd_norm,
      by = scd_norm$getfeature("DEA_cell_type")
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_CB_counts_QC_DEA_cell_type_",
      day, "_F.pdf"
    )
  )
  print(hm)

  hm_corr <- heatmap_corr_genes(
    scd = scd_norm,
    features = c("DEA_cell_type", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, "_F.pdf"
    )
  )
  print(hm_corr)
}

system("mkdir -p results/cell_type/pca")
system("mkdir -p results/cell_type/pcmf")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "F"
DEA_cell_type_palette <- cell_type_palette
DEA_genes_F <- c()
for (day in c("D15", "D90")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.Rdata"))
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
  system(paste0("rm results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
      day, "_F.RData"))
  scd_norm <- zinorm(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    cpus = 10,
    file = paste0("results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
      day, "_F.RData")
  )
  system(paste0("rm results/tmp/pca_zi_norm_counts_QC_DEA_cell_type_", day, "_F.Rdata"))
  scRNAtools::pca_plot(
    scd_norm,
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_zi_norm_counts_QC_DEA_cell_type_", day, "_F.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_zi_norm_counts_QC_DEA_cell_type_pca", day, "_F.pdf"
  ))
  system(paste0("rm results/tmp/pCMF_CB_counts_QC_DEA_cell_type_", day, "_F.Rdata"))
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_QC_DEA_cell_type_", day, "_F.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pCMF_CB_counts_QC_DEA_cell_type_", day, "_F.pdf"
  ))
}

################################################################################
## DEA for the F donor

################################################################################

length(DEA_genes_M)

length(DEA_genes_F)

table(DEA_genes_F %in% DEA_genes_M)
table(DEA_genes_M %in% DEA_genes_F)

"GNLY" %in% DEA_genes_M
"GZMB" %in% DEA_genes_M
"GZMH" %in% DEA_genes_M

"GNLY" %in% DEA_genes_F
"GZMB" %in% DEA_genes_F
"GZMH" %in% DEA_genes_F


