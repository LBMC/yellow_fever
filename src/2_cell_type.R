setwd("~/projects/yellow_fever/")
devtools::load_all("../scRNAtools/", reset = T)
load("results/QC/CB_counts_QC.Rdata")
system("mkdir -p results/cell_type")

# csv table to select feature to classify the cell-type one
feature_to_select = c("tcr_found", "norm_fsc", "all_events_fsc_h_mean",
  "norm_ssc", "all_events_ssc_h_mean", "cd57", "fas", "ptprc_cd45ra", "cd4",
  "il7ra", "dextramer", "cd3e", "cd8a", "ccr7", "lineage_neg", "itga6.cd49f",
  "pdcd1", "cd27", "cfse", "quality", "phenotype_surface_marker")
to_select <- data.frame(features = c(feature_to_select, rep(NA, scd$getngenes - length(feature_to_select))),
  genes = scd$getgenes)
write.csv(to_select, file = "results/cell_type/feature_to_select.csv")

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

# build cell_type factor
phenotype_surface_marker <- scd$getfeature("phenotype_surface_marker")
phenotype_surface_marker[phenotype_surface_marker == ""] <- NA
levels(phenotype_surface_marker) <- c("", "MEM", "MEM", "EFF", "EFF", "EFF",
  "MEM", "Naive", "EFF", "EFF", "EFF", "MEM", "MEM")
phenotype_surface_marker <- as.factor(as.vector(phenotype_surface_marker))
scd$setfeature("phenotype_surface_cell_type", phenotype_surface_marker)
b_cells <- scd$getfeature("QC_good") %in% T


################################################################################
# classification on surface_marker

surface_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "phenotype_surface_cell_type",
  features = surface_marker,
  genes = genes_marker,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/surface_cell_types"
)
save(
  surface_cell_type_classification,
  file = "results/cell_type/surface_cell_types_all_smplscv.Rdata"
)
system("~/scripts/sms.sh \"DEA done\"")
load("results/cell_type/surface_cell_types_all_smplscv.Rdata")

length(surface_cell_type_classification$groups)
surface_cell_type_classification$classification$fit_spls$fit$selected
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- surface_cell_type_classification$groups
scd$setfeature("surface_cell_type", cell_type_groups)
cell_type_pgroups <- rep(NA, scd$getncells)
cell_type_pgroups[b_cells] <- surface_cell_type_classification$pgroups
scd$setfeature("psurface_cell_type", cell_type_pgroups)

save(scd, file = "results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
scd_norm <- scd
load("results/QC/cells_counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/cells_counts_QC_surface_cell_type.Rdata")

load("results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
devtools::load_all("../scRNAtools/", reset = T)
b_cells <- scd$getfeature("QC_good") %in% T
data_gplot <- data.frame(
  ccr7 = scd$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd$select(b_cells = b_cells)$getfeature("phenotype_surface_cell_type")
)
cell_type_color <- scRNAtools::cell_type_palette(levels(data_gplot$cell_type))
data_gplot$cell_type <- as.vector(data_gplot$cell_type)
data_gplot$cell_type[is.na(data_gplot$cell_type)] <- "unknown"
data_gplot$cell_type <- as.factor(data_gplot$cell_type)
data_gplot$ccr7 <- as.numeric(as.vector(data_gplot$ccr7))
data_gplot$il7ra <- as.numeric(as.vector(data_gplot$il7ra))
summary(data_gplot)
cell_type_color <- c(cell_type_color, "gray")
names(cell_type_color) <- c(names(cell_type_color)[-3], "unknown")
ggplot() +
  scale_fill_manual( values = cell_type_color ) +
  scale_color_manual( values = cell_type_color ) +
  geom_point(data = data_gplot[data_gplot$cell_type == "unknown", ],
    aes(x = ccr7, y = il7ra, color = cell_type)) +
  geom_point(data = data_gplot[data_gplot$cell_type != "unknown", ],
    aes(x = ccr7, y = il7ra, color = cell_type)) +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_phenotype_surface_cell_type.pdf")
data_gplot <- data.frame(
  ccr7 = scd$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd$select(b_cells = b_cells)$getfeature("surface_cell_type")
)
data_gplot$ccr7 <- as.numeric(as.vector(data_gplot$ccr7))
data_gplot$il7ra <- as.numeric(as.vector(data_gplot$il7ra))
ggplot(data = data_gplot, aes(x = ccr7, y = il7ra, color = cell_type)) +
  geom_point() +
  scale_fill_manual(
      values = scRNAtools::cell_type_palette(levels(data_gplot$cell_type))
  ) +
  scale_color_manual(
    values = scRNAtools::cell_type_palette(levels(data_gplot$cell_type))
  ) +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_surface_cell_type.pdf")


system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells),
  color = "surface_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_counts_QC_surface_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "surface_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_counts_QC_surface_cell_type_", day, ".pdf"
  ))
}

system("mkdir -p results/cell_type/pcmf")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells),
  color = "surface_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/cell_type/pcmf/pcmf_counts_QC_surface_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "surface_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pcmf_counts_QC_surface_cell_type_", day, ".pdf"
  ))
}

################################################################################
# DEA on surface_cell_type

setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_surface_cell_type.Rdata")

b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
DEA_paraload_parameters(
  paraload_file = "results/cell_type/paraload_mbatch_surface_cell_type_DEA.txt",
  scd = scd,
  job_DEA_number = 5,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + surface_cell_type",
  b_cells = b_cells,
  cpus = 1,
  folder_name = "results/cell_type/mbatch_surface_cell_type_DEA"
)
table(is.na(batch_surface_cell_type_DEA$padj))
table(batch_surface_cell_type_DEA$padj < 0.05)

system("mkdir -p results/cell_type/batch_surface_cell_type_DEA")
b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
devtools::load_all("../scRNAtools/", reset = T)
batch_surface_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ batch",
  formula_full = "y ~ batch + surface_cell_type",
  b_cells = b_cells,
  cpus = 10,
  v = F,
  folder_name = "results/cell_type/batch_surface_cell_type_DEA"
)
save(
  batch_surface_cell_type_DEA,
  file = "results/cell_type/batch_surface_cell_type_DEA.Rdata"
)
system("~/scripts/sms.sh \"DEA done\"")
table(is.na(mbatch_surface_cell_type_DEA$padj))
table(mbatch_surface_cell_type_DEA$padj < 0.05)

system("mkdir -p results/cell_type/batch_day_surface_cell_type_DEA")
b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
devtools::load_all("../scRNAtools/", reset = T)
batch_day_surface_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ batch + day",
  formula_full = "y ~ batch + day + surface_cell_type",
  b_cells = b_cells,
  cpus = 10,
  v = F,
  folder_name = "results/cell_type/batch_day_surface_cell_type_DEA"
)
save(
  batch_day_surface_cell_type_DEA,
  file = "results/cell_type/batch_day_surface_cell_type_DEA.Rdata"
)
system("~/scripts/sms.sh \"DEA done\"")
table(is.na(batch_day_surface_cell_type_DEA$padj))
table(batch_day_surface_cell_type_DEA$padj < 0.05)

system("mkdir -p results/cell_type/mbatch_day_surface_cell_type_DEA")
b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
devtools::load_all("../scRNAtools/", reset = T)
mbatch_day_surface_cell_type_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch) + day",
  formula_full = "y ~ (1|batch) + day + surface_cell_type",
  b_cells = b_cells,
  cpus = 10,
  v = F,
  folder_name = "results/cell_type/mbatch_day_surface_cell_type_DEA"
)
save(
  mbatch_day_surface_cell_type_DEA,
  file = "results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata"
)
system("~/scripts/sms.sh \"DEA done\"")
table(is.na(mbatch_day_surface_cell_type_DEA$padj))
table(mbatch_day_surface_cell_type_DEA$padj < 0.05)

system("mkdir -p results/cell_type/mbatch_day_surface_cell_type_DEA")
b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
devtools::load_all("../scRNAtools/", reset = T)
DEA_paraload_parameters(
  paraload_file = "results/cell_type/mbatch_day_surface_cell_type_DEA/paraload.csv",
  scd = scd,
  job_DEA_number = 5,
  formula_null = "y ~ (1|batch) + day",
  formula_full = "y ~ (1|batch) + day + surface_cell_type",
  b_cells = b_cells,
  cpus = 1,
  folder_name = "results/cell_type/mbatch_day_surface_cell_type_DEA"
)

# launch paraload server
system("
bin/paraload --server \
--port 13469 \
--input results/cell_type/mbatch_day_surface_cell_type_DEA/paraload.csv \
--output results/cell_type/mbatch_day_surface_cell_type_DEA/paraload_run.txt \
--log results/cell_type/mbatch_day_surface_cell_type_DEA/paraload.log \
--report results/cell_type/mbatch_day_surface_cell_type_DEA/paraload_report.txt \
--conf src/pbs/DEA/DEA.conf
")

# launch paralod clients
system("
bin/paraload --client --port 13469 --host pbil-deb
")
system("
while [ $(ps -u modolo | grep paraload | wc -l) -gt 0 ]
do
stat | wc -l
iter=$(echo 200 - $(qstat -u modolo | grep -e \"[RQ]\" | wc -l) | bc)
for ((i = 1;i <= $iter;i += 1))
do
qsub src/pbs/DEA/DEA_cell_type.pbs &
/bin/sleep 0.5
done
/bin/sleep 3600
done
")

load("results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
load("results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
table(scd$getgenes %in% expressed(scd$select(b_cells = b_cells)))
table(is.na(mbatch_day_surface_cell_type_DEA$padj))
table(mbatch_day_surface_cell_type_DEA$padj < 0.05)

b_genes <- !is.na(mbatch_day_surface_cell_type_DEA$padj) &
  mbatch_day_surface_cell_type_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_DEA$gene[b_genes]

system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells, genes = DEA_genes),
  color = "surface_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_DEA_surface_cell_type_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_CB_counts_QC_DEA_surface_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes),
    color = "surface_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day,
      "QC_DEA_surface_cell_type_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_CB_counts_QC_DEA_surface_cell_type_", day, ".pdf"
  ))
}

system("mkdir -p results/cell_type/pcmf")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells, genes = DEA_genes),
  color = "surface_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_DEA_surface_cell_type_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(
  file = "results/cell_type/pcmf/pcmf_CB_counts_QC_DEA_surface_cell_type.pdf"
)
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes),
    color = "surface_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day,
      "QC_DEA_surface_cell_type_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pcmf_CB_counts_QC_DEA_surface_cell_type_", day, ".pdf"
  ))
}


load("results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
load("results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
b_genes <- !is.na(mbatch_day_surface_cell_type_DEA$padj) &
  mbatch_day_surface_cell_type_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_DEA$gene[b_genes]

system("mkdir -p results/cell_type/heatmap/")
surface_cell_type_palette <- cell_type_palette
devtools::load_all("../scRNAtools/", reset = T)
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
  file = "results/cell_type/heatmap/hm_CB_counts_QC_DEA_surface_cell_type.pdf"
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
  file = "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_surface_cell_type.pdf"
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
      "results/cell_type/heatmap/hm_CB_counts_QC_DEA_surface_cell_type_",
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
      "results/cell_type/heatmap/hm_corr_CB_counts_QC_DEA_surface_cell_type_",
      day, ".pdf"
    )
  )
  print(hm_corr)
}


################################################################################
# classification on DEA genes for surface_cell_type

devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
load("results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata")

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
load(
  file = "results/cell_type/mbatch_day_surface_cell_type_DEA.Rdata",
  v = T
)
b_genes <- !is.na(mbatch_day_surface_cell_type_DEA$padj) &
  mbatch_day_surface_cell_type_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_DEA$gene[b_genes]

DEA_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "surface_cell_type",
  features = surface_marker,
  genes = DEA_genes,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/DEA_cell_types"
)
save(
  DEA_cell_type_classification,
  file = "results/cell_type/DEA_cell_types_all_splsstab.Rdata"
)

cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- cell_type_classification$groups
scd$setfeature("DEA_cell_type", cell_type_groups)

save(scd, file = "results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
scd_norm <- scd
load("results/QC/cells_counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")


load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
ggplot(data = data.frame(
  ccr7 = scd_norm$select(b_cells = b_cells)$getfeature("ccr7"),
  il7ra = scd_norm$select(b_cells = b_cells)$getfeature("il7ra"),
  cell_type = scd_norm$select(b_cells = b_cells)$getfeature("DEA_cell_type")
), aes(x = ccr7, y = il7ra, color = cell_type)) +
  geom_point() +
  theme_bw()
ggsave(file = "results/cell_type/counts_QC_DEA_cell_type.pdf")


system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_counts_QC_DEA_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0(
    "results/cell_type/pca/pca_counts_QC_DEA_cell_type_", day, ".pdf"
  ))
}

system("mkdir -p results/cell_type/pcmf")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells),
  color = "DEA_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "DEA_cell_type", color_name = "cell_type",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0(
    "results/cell_type/pcmf/pcmf_counts_QC_DEA_cell_type_", day, ".pdf"
  ))
}
