rm(list=ls())
load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
b_cells <- scd$getfeature("sex") %in% "M" &
  scd$getfeature("day") %in% c("D15", "D136", "D593")
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
options(warn=2)

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

for (day in c("D15", "D136", "D593")) {
  system(
    paste0("mkdir -p results/cell_type/mbatch_", day, "_pDEA_cell_type_DEA")
  )
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% day
  mbatch_pDEA_cell_type_DEA <- DEA(
    scd = scd,
    formula_null = "y ~ (1|batch)",
    formula_full = "y ~ (1|batch) + pDEA_cell_type",
    continuous = "pDEA_cell_type",
    b_cells = b_cells,
    cpus = 11,
    v = T,
    folder_name = paste0("results/cell_type/mbatch_", day, "_pDEA_cell_type_DEA")
  )
  save(
    mbatch_DEA_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day, "_pDEA_cell_type_DEA.Rdata")
  )
  system("~/scripts/sms.sh \"DEA done\"")
  print(day)
  print(table(is.na(mbatch_pDEA_cell_type_DEA$padj)))
  print(table(mbatch_pDEA_cell_type_DEA$padj < 0.05))
}

for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_pDEA_cell_type_DEA.Rdata"))
  print(day)
  print(table(is.na(mbatch_pDEA_cell_type_DEA$padj)))
  print(table(mbatch_pDEA_cell_type_DEA$padj < 0.05))
  write.csv(
    mbatch_pDEA_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day, "_pDEA_cell_type_DEA.csv")
  )
}

for (day in c("D15", "D136", "D593")) {
  system(
    paste0("mkdir -p results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_logit")
  )
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% day
  mbatch_DEA_cell_type_DEA <- DEA(
    scd = scd,
    formula_null = "y ~ (1|batch)",
    formula_full = "y ~ (1|batch) + DEA_cell_type",
    b_cells = b_cells,
    family = "binomial",
    cpus = 16,
    v = T,
    folder_name = paste0("results/cell_type/mbatch_", day,
                         "_DEA_cell_type_DEA_logit")
  )
  save(
    mbatch_DEA_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day,
                  "_DEA_cell_type_DEA_logit.Rdata")
  )
  system("~/scripts/sms.sh \"DEA done\"")
  print(day)
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
}

for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day,
              "_DEA_cell_type_DEA_logit.Rdata"))
  print(day)
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
  write.csv(
    apply(mbatch_DEA_cell_type_DEA,2,as.character),
    file = paste0("results/cell_type/mbatch_", day,
                  "_DEA_cell_type_DEA_logit.csv")
  )
}

load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
system("mkdir -p results/cell_type/heatmap/")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
DEA_cell_type_palette <- cell_type_palette
DEA_genes_inter <- list()
for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.1
  DEA_genes_inter[[day]] <- mbatch_DEA_cell_type_DEA$gene[b_genes]
}

devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
for (day in c("D15", "D136", "D593")) {
  scd_norm <- scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes =c("GZMB", "CX3CR1", "CCL4", "GNLY", "GZMH", "KLRD1", "GZMG",
      "PRF1", "HOPX", "CCL5", "GZMK", "SELL", "IL7R", "LEF1", "TCF7", "LTB",
      "NELL2", "CCR7")
    )
  hm_corr <- heatmap_corr_genes(
    scd = scd_norm,
    features = c("antigen", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd_norm$getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("selected genes between cell-type ", day),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, "_selected_manhattan.pdf"
    ),
    dist_name = "manhattan",
    pca = F,
    pCMF = F,
    ncomp = 5,
    cpus = 10
  )
  print(hm_corr)
}

# same thing with all genes
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
for (day in c("D15", "D136", "D593")) {
    DEA_genes_exp <- expressed(
      scd = scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
      ),
      zi_threshold = 0.90,
    )
    scd_norm <- scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes_exp
      )
  hm_corr <- heatmap_corr_cells(
    scd = scd_norm,
    features = c("antigen", "pDEA_cell_type"),
    title = paste0("DE genes between cell-type ", day),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, "_corr_all_genes.pdf"
    ),
    dist_name = "euclidian",
    gene_size = 3,
    pca = F,
    pCMF = F,
    ncomp = 5,
    cpus = 10
  )
  print(hm_corr)
}


b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
for (day in c("D15", "D136", "D593")) {
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
  b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
    mbatch_DEA_cell_type_DEA$padj < 0.05
  DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
  scd_norm <- scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes
      )
  hm_corr <- heatmap_corr_cells(
    scd = scd_norm,
    features = c("antigen", "pDEA_cell_type"),
    title = paste0("DE genes between cell-type ", day),
    genes_order = 1:length(DEA_genes),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, "_corr_genes.pdf"
    ),
    dist_name = "euclidian",
    gene_size = 3,
    pca = F,
    pCMF = F,
    ncomp = 5,
    cpus = 10
  )
  print(hm_corr)
  for (alt in c("lesser", "greater")) {
    DEA_genes_exp <- expressed(
      scd = scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes
      ),
      zi_threshold = 0.50,
      alt = alt
    )
    scd_norm <- scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes_exp
      )
    hm_title <- paste0("DE genes between cell-type ", day, " (zi <= 0.5)")
    if (alt == "greater") {
      hm_title <- paste0("DE genes between cell-type ", day, " (zi >= 0.5)")
    }
    hm_corr <- heatmap_corr_cells(
      scd = scd_norm,
      features = c("antigen", "pDEA_cell_type"),
      title = paste0("Correlation between genes at ", day),
      genes_order = 1:length(DEA_genes_exp),
      factor = c(T, F),
      file = paste0(
        "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
        day, "_zi0.5_", alt, "_corr_genes.pdf"
      ),
      dist_name = "euclidian",
      gene_size = 3,
      pca = F,
      pCMF = F,
      ncomp = 5,
      cpus = 10
    )
  }
}

DEA_genes_M <- c()
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
for (day in c("D15", "D136", "D593")) {
  for (alt in c("lesser", "greater")) {
    load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata"))
    b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
      mbatch_DEA_cell_type_DEA$padj < 0.05
    DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
    DEA_genes_M <- unique(c(DEA_genes, DEA_genes_M))
    DEA_genes_exp <- expressed(
      scd = scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes
      ),
      zi_threshold = 0.50,
      alt = alt
    )
    scd_norm <- scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes_exp
      )
    hm_title <- paste0("top ", length(DEA_genes_exp),
      "DE genes between cell-type
      ", day, "(zi <= 0.5)")
    if (alt == "greater") {
      hm_title <- paste0("top ", length(DEA_genes_exp),
        "DE genes between cell-type
        ", day, "(zi >= 0.5)")
    }
    system(paste0("rm results/cell_type/gene_cov_",
      day, "_zi0.5_", alt, ".Rdata_norm.Rdata"))
    gene_order <- order_by_factor(
      scd = scd_norm,
      score = scd_norm$getfeature("pDEA_cell_type"),
      tmp_file = paste0("results/cell_type/gene_cov_",
        day, "_zi0.5_", alt, ".Rdata"),
      top = min(100, length(DEA_genes_exp))
    )
    system(paste0("rm results/cell_type/gene_cov_",
        day, "_zi0.5_", alt, "_top100.Rdata_norm.Rdata"))
    gene_order <- order_by_factor(
      scd = scd_norm,
      score = scd_norm$getfeature("pDEA_cell_type"),
      tmp_file = paste0("results/cell_type/gene_cov_",
        day, "_zi0.5_", alt, "_top100.Rdata"),
      top = min(100, length(DEA_genes_exp))
    )

    hm <- heatmap_genes(
      scd = scd_norm,
      features = c("antigen", "pDEA_cell_type"),
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
        day, "_zi0.5_", alt, "_cov.pdf"
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
      title = hm_title,
      factor = c(T, F),
      file = paste0(
        "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
        day, "_zi0.5_", alt, "_cov_manhattan.pdf"
      ),
      dist_name = "manhattan",
      pca = F,
      pCMF = F,
      ncomp = 5,
      cpus = 10
    )
    print(hm_corr)
  }
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

rm(list=ls())
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
infos_M <- scd$getfeatures
load("results/cell_type/cells_counts_QC_DEA_cell_type_F.Rdata")
b_cells_F <- scd$getfeature("sex") %in% "F"
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
    cpus = 16,
    v = T,
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
  day <- "D15"
  load(paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.Rdata"))
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
  write.csv(
    mbatch_DEA_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.csv")
  )
}

for (day in c("D15", "D90")) {
  system(
    paste0("mkdir -p results/cell_type/mbatch_", day,
           "_DEA_cell_type_DEA_F_logit")
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
    family = "binomial",
    cpus = 16,
    v = T,
    folder_name = paste0("results/cell_type/mbatch_", day,
                         "_DEA_cell_type_DEA_F_logit")
  )
  save(
    mbatch_DEA_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day,
                  "_DEA_cell_type_DEA_F_logit.Rdata")
  )
  system("~/scripts/sms.sh \"DEA done\"")
  print(day)
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
}

for (day in c("D15", "D90")) {
  load(paste0("results/cell_type/mbatch_", day,
              "_DEA_cell_type_DEA_F_logit.Rdata"))
  print(table(is.na(mbatch_DEA_cell_type_DEA$padj)))
  print(table(mbatch_DEA_cell_type_DEA$padj < 0.05))
  write.csv(
    apply(mbatch_DEA_cell_type_DEA,2,as.character),
    file = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.csv")
  )
}


system("mkdir -p results/cell_type/heatmap/")

load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
infos_M <- scd$getfeatures
load("results/cell_type/CB_counts_QC_DEA_cell_type_F.Rdata")
b_cells_F <- scd$getfeature("sex") %in% "F"
infos_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getfeatures
scd <- scdata$new(
  infos = infos_M,
  counts = scd$getcounts
)

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
      genes = intersect_3
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

devtools::load_all("../scRNAtools/", reset = T)
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "F"
for (day in c("D15", "D90")) {
  hm_corr <- heatmap_corr_genes(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes =c("GZMB", "CX3CR1", "CCL4", "GNLY", "GZMH", "KLRD1", "GZMG",
      "PRF1", "HOPX", "CCL5", "GZMK", "SELL", "IL7R", "LEF1", "TCF7", "LTB",
      "NELL2", "CCR7")
    ),
    features = c("antigen", "pDEA_cell_type"),
    cells_order = order(
      as.numeric(as.vector(
        scd$select(
          b_cells = b_cells & scd$getfeature("day") %in% day
        )$getfeature("pDEA_cell_type")
      ))
    ),
    title = paste0("DE genes between DEA_cell_type ", day),
    factor = c(T, F),
    file = paste0(
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
      day, "_F_selected_manhattan.pdf"
    ),
    dist_name = "canberra",
    pca = F,
    pCMF = F,
    ncomp = 5,
    cpus = 10
  )
  print(hm_corr)
}

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "F"
for (day in c("D15", "D90")) {
  for (alt in c("lesser", "greater")) {
    load(paste0("results/cell_type/mbatch_",
      ifelse(day %in% "D90", "D136", day), "_DEA_cell_type_DEA.Rdata"))
    b_genes <- !is.na(mbatch_DEA_cell_type_DEA$padj) &
      mbatch_DEA_cell_type_DEA$padj < 0.05
    DEA_genes <- mbatch_DEA_cell_type_DEA$gene[b_genes]
    DEA_genes_exp <- expressed(
      scd = scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes
      ),
      zi_threshold = 0.50,
      alt = alt
    )
    scd_norm <- scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% day,
        genes = DEA_genes_exp
      )
    hm_title <- paste0("top ", length(DEA_genes_exp),
      "DE genes between cell-type
      ", day, "(zi <= 0.5)")
    if (alt == "greater") {
      hm_title <- paste0("top ", length(DEA_genes_exp),
        "DE genes between cell-type
        ", day, "(zi >= 0.5)")
    }
    system(paste0("rm results/cell_type/gene_cov_",
      day, "_zi0.5_", alt, ".Rdata_norm.Rdata"))
    gene_order <- order_by_factor(
      scd = scd_norm,
      score = scd_norm$getfeature("pDEA_cell_type"),
      tmp_file = paste0("results/cell_type/gene_cov_",
        day, "_zi0.5_", alt, "_F.Rdata"),
      top = min(100, length(DEA_genes_exp))
    )
    system(paste0("rm results/cell_type/gene_cov_",
        day, "_zi0.5_", alt, "_top100.Rdata_norm.Rdata"))
    gene_order <- order_by_factor(
      scd = scd_norm,
      score = scd_norm$getfeature("pDEA_cell_type"),
      tmp_file = paste0("results/cell_type/gene_cov_",
        day, "_zi0.5_", alt, "_top100_F.Rdata"),
      top = min(100, length(DEA_genes_exp))
    )

    hm <- heatmap_genes(
      scd = scd_norm,
      features = c("antigen", "pDEA_cell_type"),
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
        day, "_zi0.5_", alt, "_cov_F.pdf"
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
      title = hm_title,
      factor = c(T, F),
      file = paste0(
        "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_cell_type_",
        day, "_zi0.5_", alt, "_cov_manhattan_F.pdf"
      ),
      dist_name = "manhattan",
      pca = F,
      pCMF = F,
      ncomp = 5,
      cpus = 10
    )
    print(hm_corr)
  }
}


load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
b_cells <- scd$getfeature("sex") %in% "M" &
  scd$getfeature("day") %in% c("D15", "D136", "D593")
infos_M <- scd$getfeatures
counts_M <- scd$getcounts
dim(infos_M)
dim(counts_M)
load("results/cell_type/CB_counts_QC_DEA_cell_type_F.Rdata")
b_cells_F <- scd$getfeature("sex") %in% "F" &
  scd$getfeature("day") %in% c("D15", "D90")
infos_M <- rbind(infos_M, scd$select(b_cells = b_cells_F)$getfeatures)
counts_M <- rbind(counts_M, scd$select(b_cells = b_cells_F)$getcounts)
infos_M <- t(infos_M)
counts_M <- t(counts_M)
dim(infos_M)
dim(counts_M)

normalized_counts <- rbind(infos_M, counts_M)
dim(normalized_counts)
write.csv(
  normalized_counts,
  file = paste0("results/cell_type/cell_type_CB_counts.csv")
)



load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
b_cells <- scd$getfeature("sex") %in% "M" &
  scd$getfeature("day") %in% c("D15", "D136", "D593") &
  scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type"))

scd_data <- scd$select(b_cells = b_cells)$getfeatures
data <- data.frame(
  ccr7 = as.numeric(as.vector(scd_data$ccr7)),
  il7ra = as.numeric(as.vector(scd_data$il7ra)),
  pMEM = ifelse(
    scd_data$phenotype_surface_cell_type %in% "MEM",
    1,
    ifelse(
      scd_data$phenotype_surface_cell_type %in% "EFF",
      0,
      NA
    )
  ),
  day = factor(scd_data$day, levels = c("D15", "D136", "D593")),
  classification = "manual"
)
data <- rbind(
  data,
  data.frame(
    ccr7 = as.numeric(as.vector(scd_data$ccr7)),
    il7ra = as.numeric(as.vector(scd_data$il7ra)),
    pMEM = scd_data$psurface_cell_type,
    day = factor(scd_data$day, levels = c("D15", "D136", "D593")),
    classification = "PLS from manual"
  )
)
data <- rbind(
  data,
  data.frame(
    ccr7 = as.numeric(as.vector(scd_data$ccr7)),
    il7ra = as.numeric(as.vector(scd_data$il7ra)),
    pMEM = scd_data$pDEA_cell_type,
    day = factor(scd_data$day, levels = c("D15", "D136", "D593")),
    classification = "PLS from DEA"
  )
)

ggplot(data = data,
       aes(x = classification, y = ccr7, color = pMEM)) +
  scale_y_log10() +
  geom_jitter(height = 0) +
  facet_wrap(~day) +
  theme_bw()
ggsave("results/cell_type/M_in_vivo_ccr7_pMEM.pdf",
       width = 11,
       height = 10)

ggplot(data = data,
       aes(x = classification, y = il7ra, color = pMEM)) +
  scale_y_log10() +
  geom_jitter(height = 0) +
  facet_wrap(~day) +
  theme_bw()
ggsave("results/cell_type/M_in_vivo_il7ra_pMEM.pdf",
       width = 11,
       height = 10)

load("results/cell_type/CB_counts_QC_DEA_cell_type_F.Rdata")
b_cells <- scd$getfeature("sex") %in% "F" &
  scd$getfeature("day") %in% c("D15", "D90") &
  scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type"))

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
ggsave("results/cell_type/F_in_vivo_ccr7_pMEM.pdf",
       width = 11,
       height = 10)

ggplot(data = data,
       aes(x = classification, y = il7ra, color = pMEM)) +
  scale_y_log10() +
  geom_jitter(height = 0) +
  facet_wrap(~day) +
  theme_bw()
ggsave("results/cell_type/F_in_vivo_il7ra_pMEM.pdf",
       width = 11,
       height = 10)

load("results/cell_type/CB_counts_QC_DEA_cell_type_F.Rdata")
b_cells <- scd$getfeature("sex") %in% "F" &
  scd$getfeature("day") %in% c("D15", "D90") &
  scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type"))

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
ggsave("results/cell_type/F_in_vivo_ccr7_pMEM.pdf",
       width = 11,
       height = 10)

ggplot(data = data,
       aes(x = classification, y = il7ra, color = pMEM)) +
  scale_y_log10() +
  geom_jitter(height = 0) +
  facet_wrap(~day) +
  theme_bw()
ggsave("results/cell_type/F_in_vivo_il7ra_pMEM.pdf",
       width = 11,
       height = 10)

