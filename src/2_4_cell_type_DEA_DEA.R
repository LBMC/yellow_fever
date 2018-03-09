setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
# options(warn=2)

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
    v = T,
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

for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
  print(day)
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
  write.csv(
    mbatch_DEA_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.csv")
  )
}

system("mkdir -p results/cell_type/heatmap/")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
DEA_cell_type_palette <- cell_type_palette
DEA_genes_inter <- list()
for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  DEA_genes_inter[[day]] <- mbatch_DEA_cell_type_DEA$gene[b_genes]
}

DEA_genes_M <- c()
for (day in c("D15", "D136", "D593")) {
  for (restrict in c("all", "40perc")) {
    load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
    system(paste0("rm results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
        day, "_", restrict, ".Rdata")
)
    b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
      mbatch_DEA_cell_type_DEA$padj < 0.05
    print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
    print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
    DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
    print("GNLY" %in% DEA_genes)
    print("GZMB" %in% DEA_genes)
    print("GZMH" %in% DEA_genes)
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
      zi_threshold = 0.40
    )))
    if (restrict %in% "all") {
      DEA_genes_thresholded <- DEA_genes
    }
    scd_norm <- zinorm(
      scd = scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes_thresholded
      ),
      cpus = 10,
      file = paste0("results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
        day, "_", restrict, ".Rdata")
    )
    hm <- heatmap_genes(
      scd = scd_norm,
      features = c("antigen", "pDEA_cell_type"),
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
        day, "_", restrict, ".pdf"
      )
    )
    print(hm)
    hm_corr <- heatmap_corr_genes(
      scd = scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes_thresholded
      ),
      features = c("antigen", "pDEA_cell_type"),
      cells_order = order(
        as.numeric(as.vector(
          scd_norm$getfeature("pDEA_cell_type")
        ))
      ),
      title = paste0("DE genes between DEA_cell_type ", day),
      factor = c(T, F),
      file = paste0(
        "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
        day, "_", restrict, ".pdf"
      )
    )
    print(hm_corr)
  }
  hm_corr <- heatmap_corr_genes(
    scd = scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes_inter[["D136"]]
      ),
    features = c("antigen", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, "_inter.pdf"
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
  system(paste0("rm results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
      day, ".RData")
  )
  system(paste0("rm results/tmp/pCMF_CB_counts_QC_DEA_cell_type_", day, ".Rdata"))
  system(paste0("rm results/tmp/pca_zi_norm_counts_QC_DEA_cell_type_", day, ".Rdata"))
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
  scd_norm <- zinorm(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    cpus = 11,
    file = paste0("results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
      day, ".RData")
  )
  scRNAtools::pca_plot(
    scd_norm,
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_zi_norm_counts_QC_DEA_cell_type_", day, ".Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_zi_norm_counts_QC_DEA_cell_type_pca", day, ".pdf"
  ))
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

load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
scd_M <- scd
infos_M <- scd$getfeatures
load("results/cell_type/cells_counts_QC_DEA_cell_type_F.Rdata")
b_cells_F <- scd$getfeature("QC_good") %in% T &
  scd$getfeature("sex") %in% "F"
infos_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getfeatures
scd <- scdata$new(
  infos = infos_M,
  counts = scd$getcounts
)

infos_M$pDEA_cell_type[b_cells_F]
infos_M$QC_good[b_cells_F]

write.csv(
  infos_M,
  file = paste0("results/cell_type/cell_type_infos.csv")
)

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

system("mkdir -p results/cell_type/heatmap/")
b_cells <- scd$getfeature("QC_good") %in% T &
infos_M <- scd$getfeatures
load("results/cell_type/cells_counts_QC_DEA_cell_type_F.Rdata")
b_cells_F <- scd$getfeature("QC_good") %in% T &
  scd$getfeature("sex") %in% "F"
infos_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getfeatures
scd <- scdata$new(
  infos = infos_M,
  counts = scd$getcounts
)

write.csv(
  infos_M,
  file = paste0("results/cell_type/cell_type_infos.csv")
)

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
    cpus = 10,
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
  for (restrict in c("all", "40perc")) {
    load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.Rdata"))
    system(paste0("rm results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
        day, "_thresholded_", restrict, "_F.RData"))
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
      zi_threshold = 0.40
    )))
    if (restrict %in% "all") {
      DEA_genes_thresholded <- DEA_genes
    }
    scd_norm <- zinorm(
      scd = scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes_thresholded
      ),
      cpus = 10,
      file = paste0("results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
        day, "_thresholded_", restrict, "_F.RData")
    )
    hm <- heatmap_genes(
      scd = scd_norm,
      features = c("antigen", "pDEA_cell_type"),
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
        day, "_", restrict, "_F.pdf"
      )
    )
    print(hm)
    hm_corr <- heatmap_corr_genes(
      scd = scd_norm,
      features = c("antigen", "pDEA_cell_type"),
      cells_order = order(
        as.numeric(as.vector(
          scd_norm$getfeature("pDEA_cell_type")
        ))
      ),
      title = paste0("DE genes between DEA_cell_type ", day),
      factor = c(T, F),
      file = paste0(
        "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
        day, "_", restrict, "_F.pdf"
      )
    )
    print(hm_corr)
  }
}

DEA_genes_inter <- list()
for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  DEA_genes_inter[[day]] <- mbatch_DEA_cell_type_DEA$gene[b_genes]
}

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "F"
b_cells_M <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
for (day in c("D15", "D136")) {
  for (restrict in c("all", "40perc")) {
    load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
    b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
      mbatch_DEA_cell_type_DEA$padj < 0.05
    print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
    print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
    DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
    DEA_genes_thresholded <- expressed(
      scd = scd$select(
        b_cells = b_cells_M & scd$getfeature("day") %in% day &
          scd$getfeature("DEA_cell_type") %in% "EFF",
        genes = DEA_genes
      ),
      zi_threshold = 0.40
    )
    DEA_genes_thresholded <- unique(c(DEA_genes_thresholded, expressed(
      scd = scd$select(
        b_cells = b_cells_M & scd$getfeature("day") %in% day &
          scd$getfeature("DEA_cell_type") %in% "MEM",
        genes = DEA_genes
      ),
      zi_threshold = 0.40
    )))
    if (restrict %in% "all") {
      DEA_genes_thresholded <- DEA_genes
    }
    table(b_cells & scd$getfeature("day") %in% day)
    system(paste0("rm results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
        day, "_thresholded_", restrict, "_M_genes_F.RData"))
    scd_norm <- zinorm(
      scd = scd$select(
        b_cells = b_cells &
          scd$getfeature("day") %in% ifelse(day %in% "D136", "D90", day),
        genes = DEA_genes_thresholded
      ),
      cpus = 10,
      file = paste0("results/tmp/zi_norm_cells_counts_QC_DEA_cell_type_",
        day, "_thresholded_", restrict, "_M_genes_F.RData")
    )
    hm <- heatmap_genes(
      scd = scd_norm,
      features = c("antigen", "pDEA_cell_type"),
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
        day, "_", restrict, "_M_genes_F.pdf"
      )
    )
    hm_corr <- heatmap_corr_genes(
      scd = scd$select(
        b_cells = b_cells &
          scd$getfeature("day") %in% ifelse(day %in% "D136", "D90", day),
        genes = DEA_genes_thresholded
      ),
      features = c("antigen", "pDEA_cell_type"),
      cells_order = order(
        as.numeric(as.vector(
          scd_norm$getfeature("pDEA_cell_type")
        ))
      ),
      title = paste0("DE genes between DEA_cell_type ", day),
      factor = c(T, F),
      file = paste0(
        "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
        day, "_", restrict, "_M_genes_F.pdf"
      )
    )
    print(hm_corr)
  }
  hm_corr <- heatmap_corr_genes(
    scd = scd$select(
      b_cells = b_cells &
        scd$getfeature("day") %in% ifelse(day %in% "D136", "D90", day),
      genes = DEA_genes_inter[["D136"]]
    ),
    features = c("antigen", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, "_inter_M_genes_F.pdf"
    )
  )
  print(hm_corr)
}

