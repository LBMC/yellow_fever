source("src/00_functions.R")
load(file = "results/sce_sctrf.Rdata", verbose = T)

genes_PLS <- read_csv("data/2017_11_28_List_Laurent_Genes_PLS.csv") %>% 
  pivot_longer(cols = c("Genes_EFF", "Genes_MEM"),
               names_to = "type",
               values_to = "genes") %>% 
  dplyr::rename(proteins = `Protein_Markers`)

# build cell_type factor
sce$manual_cell_type <- NA
sce$manual_cell_type <- factor(sce$manual_cell_type, levels = c("EFF", "MEM"))
sce$manual_cell_type[sce$male_invivo] <- sce$phenotype_surface_marker[sce$male_invivo] %>%
  as_factor() %>%
  lvls_revalue(c("MEM", "MEM", "EFF", "EFF", "EFF", "MEM"))
no_surface <- colData(sce)[, genes_PLS %>%
                 pull(proteins) %>%
                 na.omit() %>%
                 unique()
               ] %>%
  as_tibble() %>% rowSums() %>% is.na()
sce$manual_cell_type[no_surface] <- NA
genes_PLS %>% pull(proteins) %>% na.omit() %>% unique()
genes_PLS

fit <- PLS_fit(
  sce = sce,
  group_by = sce$manual_cell_type,
  genes = rownames(sce)[
    match(genes_PLS %>% pull(genes) %>% na.omit(),rowData(sce)$gene)
  ],
  features = genes_PLS %>% pull(proteins) %>% na.omit() %>% unique(),
  assay_name = "counts_vst",
  altExp_name = "PLS_surface_cell_type",
  cpus = 10
)
save(fit, file = "results/03_PLS_fit.Rdata")
load(file = "results/03_PLS_fit.Rdata", verbose = T)

sce <- PLS_predict(
  sce = fit$sce,
  fit = fit,
  group_by = sce$manual_cell_type,
  genes = rownames(sce)[
    match(genes_PLS %>% pull(genes) %>% na.omit(),rowData(sce)$gene)
  ],
  features = genes_PLS %>% pull(proteins) %>% na.omit() %>% unique(),
  assay_name = "counts_vst",
  altExp_name = "PLS_surface_cell_type",
  cpus = 10
)
save(sce, file = "results/sce_PLS_surface_cell_type.Rdata")
load(file = "results/sce_PLS_surface_cell_type.Rdata", verbose = T)

altExp(sce, "PLS_surface_cell_type") %>% 
  rowData() %>% 
  as_tibble(rownames = "id") %>% 
  select(id, predictor) %>% 
  filter(predictor) %>% 
  left_join(rowData(sce) %>% as_tibble())

colData(sce) %>% 
  as_tibble() %>% 
  filter(male_invivo) %>% 
  ggplot() +
  geom_point(aes(x = p_PLS_surface_cell_type, y = sum, color = PLS_surface_cell_type)) +
  facet_wrap(~day) +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes detected",
       x = " sum")

DEA_surface_cell_type <- DEA(
  sce[, sce$male_invivo],
  test = "~ p_PLS_surface_cell_type",
  formula = "count ~ p_PLS_surface_cell_type + day + (1|batch)",
  assay_name = "counts_vst",
  cpus = 10
)
save(DEA_surface_cell_type, file = "results/03_DEA.Rdata")
load(file = "results/03_DEA.Rdata", verbose = T)

assay(sce[
    match(genes_PLS %>% pull(genes) %>% na.omit(),rowData(sce)$gene),
    sce$male_invivo
  ],
  "counts_vst") %>% 
  t() %>% 
  as.matrix() %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything()) %>% 
  left_join(rowData(sce) %>% 
              as_tibble() %>% 
              dplyr::rename(name = id)) %>% 
  ggplot() +
  geom_histogram(aes(x = value)) +
  facet_wrap(~gene)

rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo <- NA
rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo <- 
  get_genes_pval(DEA_surface_cell_type, sce)

rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo %>% 
  is.na() %>% 
  table()

DEA_surface_cell_type_fixed <- DEA(
  sce[
    rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo %>% is.na(),
    sce$male_invivo
  ],
  test = "~ p_PLS_surface_cell_type",
  formula = "count ~ p_PLS_surface_cell_type + day + batch",
  assay_name = "counts_vst",
  cpus = 10
)

save(DEA_surface_cell_type_fixed, file = "results/03_DEA_fixed.Rdata")

get_genes_pval(DEA_surface_cell_type_fixed, sce[
  rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo %>% is.na(),
]) %>%
  is.na() %>%
  table()

rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo[
  rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo %>% is.na()
] <- get_genes_pval(DEA_surface_cell_type_fixed, sce[
  rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo %>% is.na(),
])

rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo %>% 
  is.na() %>% 
  table()

rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo_adj <- p.adjust(
  rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo,
  method = "BH"
)

rowData(sce) %>% 
  as_tibble() %>% 
  ggplot(aes(log_mean, log_var, color = is.na(pval_DEA_p_PLS_surface_cell_type_male_invivo_adj))) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_density_2d(size = 0.3) +
  geom_abline(intercept = 0,
              slope = 1,
              color = 'red')

save(sce, file = "results/sce_DEA_surface_cell_type.Rdata")
load(file = "results/sce_DEA_surface_cell_type.Rdata")

colData(sce) %>% 
  as_tibble() %>% 
  write_csv(path = "results/2020_12_15_sce_DEA_cell_type_cellData.csv")
rowData(sce) %>% 
  as_tibble(rownames = "id") %>% 
  write_csv(path = "results/2020_12_15_sce_DEA_cell_type_geneData.csv")

sce <- logNormCounts(sce, exprs_values = "counts_vst")
sce_DEA <- runPCA(sce[
    rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo_adj <= 0.05 &
      !is.na(rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo_adj), 
    colData(sce)$male_invivo & !is.na(colData(sce)$p_PLS_surface_cell_type)
  ])
sce_DEA <- runTSNE(sce_DEA)

plotPCA(sce_DEA,
  colour_by = "p_PLS_surface_cell_type")
plotTSNE(sce_DEA,
  colour_by = "p_PLS_surface_cell_type")
plotPCA(sce_DEA,
  colour_by = "manual_cell_type")
plotTSNE(sce_DEA,
  colour_by = "manual_cell_type")
plotPCA(sce_DEA,
  colour_by = "day")
plotTSNE(sce_DEA,
  colour_by = "day")

rownames(sce_DEA) <- rowData(sce_DEA)$gene
plotHeatmap(
  sce_DEA,
  features = rownames(sce_DEA) %in% (
    rownames(sce_DEA)[order(rowData(sce_DEA)$pval_DEA_p_PLS_surface_cell_type_male_invivo_adj)] %>% head(30)
  ),
  order_columns_by = "p_PLS_surface_cell_type",
  order_rows_by = "pval_DEA_p_PLS_surface_cell_type_male_invivo_adj",
  colour_columns_by = c("day", "antigen", "manual_cell_type"),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-5, 5)
) 

plotHeatmap(
  sce_DEA,
  features = rownames(sce_DEA) %in% (
    genes_PLS %>% pull(genes) %>% na.omit()   
  ),
  order_columns_by = "p_PLS_surface_cell_type",
  order_rows_by = "pval_DEA_p_PLS_surface_cell_type_male_invivo_adj",
  colour_columns_by = c("day", "antigen", "manual_cell_type"),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-5, 5)
) 

rowData(sce[
  rowData(sce)$gene %in% (genes_PLS %>% pull(genes) %>% na.omit()),
  ]
)

rowData(sce[
  rowData(sce)$gene %in% (genes_PLS %>% pull(genes) %>% na.omit()),
  ]
)
