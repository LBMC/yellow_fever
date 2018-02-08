setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type"))
mbatch_DEA_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch) + day",
  formula_full = "y ~ (1|batch) + day + DEA_cell_type",
  b_cells = b_cells,
  cpus = 16,
  v = F,
  folder_name = paste0("results/cell_type/mbatch_day_DEA_cell_type_DEA")
)
save(
  mbatch_DEA_cell_type_DEA,
  file = paste0("results/cell_type/mbatch_day_DEA_cell_type_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))

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

load("results/cell_type/CB_counts_QC_DEA_cell_type_weighted.Rdata")
load("results/cell_type/mbatch_day_cell_type_weighted_DEA.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("DEA_cell_type"))
b_genes <- !is.na(mbatch_day_surface_cell_type_weighted_DEA$padj) &
  mbatch_day_surface_cell_type_weighted_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_weighted_DEA$gene[b_genes]

system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells, genes = DEA_genes),
  color = "surface_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_DEA_DEA_cell_type_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_CB_counts_QC_DEA_DEA_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes),
    color = "surface_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day,
      "QC_DEA_DEA_cell_type_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_CB_counts_QC_DEA_DEA_cell_type_", day, ".pdf"
  ))
}

system("mkdir -p results/cell_type/pcmf")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells, genes = DEA_genes),
  color = "surface_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_DEA_DEA_cell_type_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(
  file = "results/cell_type/pcmf/pcmf_CB_counts_QC_DEA_DEA_cell_type.pdf"
)
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes),
    color = "surface_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day,
      "QC_DEA_DEA_cell_type_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pcmf_CB_counts_QC_DEA_DEA_cell_type_", day, ".pdf"
  ))
}

system("mkdir -p results/cell_type/heatmap/")
surface_cell_type_palette <- cell_type_palette
hm <- heatmap_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("surface_cell_type", "day", "psurface_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("psurface_cell_type")
    ))
  ),
  genes_order = order(scd$select(b_cells = b_cells, genes = DEA_genes)$getgenes),
  title = "DE genes between surface_cell_type",
  factor = c(T, T, F),
  file = "results/cell_type/heatmap/hm_CB_counts_QC_DEA_DEA_cell_type.pdf"
)
print(hm)
hm_corr <- heatmap_corr_genes(
  scd = scd$select(b_cells = b_cells, genes = DEA_genes),
  features = c("surface_cell_type", "day", "psurface_cell_type"),
  cells_order = order(
    scd$select(b_cells = b_cells)$getfeature("day"),
    as.numeric(as.vector(
      scd$select(b_cells = b_cells)$getfeature("psurface_cell_type")
    ))
  ),
  title = "corr DE genes between surface_cell_type",
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
    features = c("surface_cell_type", "day", "psurface_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("psurface_cell_type")
      ))
    ),
    genes_order = order(
      scd$select(b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes)$getgenes
    ),
    title = paste0("DE genes between surface_cell_type ", day),
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
    features = c("surface_cell_type", "day", "psurface_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$
          getfeature("psurface_cell_type")
      ))
    ),
    title = paste0("DE genes between surface_cell_type ", day),
    factor = c(T, T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_DEA_cell_type_",
      day, ".pdf"
    )
  )
  print(hm_corr)
}
