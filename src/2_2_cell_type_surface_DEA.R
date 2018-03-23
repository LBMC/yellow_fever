################################################################################
# DEA on surface_cell_type

setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)

load("results/cell_type/cells_counts_QC_surface_cell_type_weighted.Rdata")
system("mkdir -p results/cell_type/mbatch_day_surface_cell_type_weighted_DEA")
b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
devtools::load_all("../scRNAtools/", reset = T)
mbatch_day_surface_cell_type_weighted_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch) + day",
  formula_full = "y ~ (1|batch) + day + surface_cell_type",
  b_cells = b_cells,
  cpus = 10,
  v = F,
  folder_name = "results/cell_type/mbatch_day_surface_cell_type_weighted_DEA"
)
save(
  mbatch_day_surface_cell_type_weighted_DEA,
  file = "results/cell_type/mbatch_day_surface_cell_type_weighted_DEA.Rdata"
)
system("~/scripts/sms.sh \"DEA done\"")
table(is.na(mbatch_day_surface_cell_type_weighted_DEA$padj))
table(mbatch_day_surface_cell_type_weighted_DEA$padj < 0.05)

b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
devtools::load_all("../scRNAtools/", reset = T)
DEA_paraload_parameters(
  paraload_file = "results/cell_type/mbatch_day_surface_cell_type_weighted_DEA/paraload.csv",
  scd = scd,
  job_DEA_number = 5,
  formula_null = "y ~ (1|batch) + day",
  formula_full = "y ~ (1|batch) + day + surface_cell_type",
  b_cells = b_cells,
  cpus = 1,
  folder_name = "results/cell_type/mbatch_day_surface_cell_type_weighted_DEA"
)

# launch paraload server
system("
bin/paraload --server \
--port 13469 \
--input results/cell_type/mbatch_day_surface_cell_type_weighted_DEA/paraload.csv \
--output results/cell_type/mbatch_day_surface_cell_type_weighted_DEA/paraload_run.txt \
--log results/cell_type/mbatch_day_surface_cell_type_weighted_DEA/paraload.log \
--report results/cell_type/mbatch_day_surface_cell_type_weighted_DEA/paraload_report.txt \
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

load("results/cell_type/CB_counts_QC_surface_cell_type_weighted.Rdata")
load("results/cell_type/mbatch_day_surface_cell_type_weighted_DEA.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T & !is.na(scd$getfeature("surface_cell_type"))
b_genes <- !is.na(mbatch_day_surface_cell_type_weighted_DEA$padj) &
  mbatch_day_surface_cell_type_weighted_DEA$padj < 0.05
DEA_genes <- mbatch_day_surface_cell_type_weighted_DEA$gene[b_genes]

system("mkdir -p results/cell_type/pca")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells, genes = DEA_genes),
  color = "surface_cell_type", color_name = "cell_type",
  tmp_file = "results/tmp/pca_CB_counts_QC_DEA_surface_cell_type_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/cell_type/pca/pca_CB_counts_QC_DEA_surface_cell_type.pdf")
for (day in c("D15", "D136", "D593")) {
  scd_norm <- zinorm(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    cpus = 10,
    file = paste0("results/tmp/zi_norm_cells_counts_surface_DEA_cell_type_",
      day, ".RData")
  )
  scRNAtools::pca_plot(
    scd_norm,
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
  scd_norm <- zinorm(
    scd = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% day,
      genes = DEA_genes
    ),
    cpus = 10,
    file = paste0("results/tmp/zi_norm_cells_counts_surface_DEA_cell_type_",
      day, ".RData")
  )
  scRNAtools::pCMF_plot(
    scd_norm,
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
