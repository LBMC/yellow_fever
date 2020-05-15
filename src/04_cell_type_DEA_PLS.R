rm(list = ls())
source("src/00_functions.R")

load(file = "results/sce_DEA_surface_cell_type.Rdata", verbose = T)

DEA_genes <- rownames(sce)[
  rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo_adj <= 0.05 &
  !is.na(rowData(sce)$pval_DEA_p_PLS_surface_cell_type_male_invivo_adj)
]
length(DEA_genes)

genes_PLS <- read_csv("data/2017_11_28_List_Laurent_Genes_PLS.csv") %>% 
  pivot_longer(cols = c("Genes_EFF", "Genes_MEM"),
               names_to = "type",
               values_to = "genes") %>% 
  dplyr::rename(proteins = `Protein_Markers`)

fit <- PLS_fit(
  sce = sce,
  group_by = sce$manual_cell_type,
  genes = c(
    rownames(sce)[
      match(
        genes_PLS %>%
          pull(genes) %>%
          na.omit(),
        rowData(sce)$gene
      )
    ], DEA_genes) %>%
    unique(),
  assay_name = "counts_vst",
  altExp_name = "PLS_DEA_cell_type",
  force = rownames(sce)[
      match(
        genes_PLS %>%
          pull(genes) %>%
          na.omit(),
        rowData(sce)$gene
      )
    ],
  cpus = 10
)
save(fit, file = "results/04_PLS_fit.Rdata")
load(file = "results/04_PLS_fit.Rdata", verbose = T)

sce <- PLS_predict(
  sce = fit$sce,
  fit = fit,
  group_by = fit$sce$manual_cell_type,
  genes = rownames(altExp(fit$sce, "PLS_DEA_cell_type")),
  assay_name = "counts_vst",
  altExp_name = "PLS_DEA_cell_type",
  cpus = 10
)
colData(sce) %>%
  as_tibble() %>%
  filter(male_invivo) %>%
  summary()

colData(sce) %>%
  as_tibble() %>%
  filter(YVF2003_D1401) %>%
  summary()
sce[rownames(sce) %in% (
    altExp(fit$sce, "PLS_DEA_cell_type") %>%
      assay("counts") %>% 
      as_tibble(rownames = "gene_id") %>% 
      dplyr::filter_all(any_vars(is.nan(.))) %>% 
      pull(gene_id)
  ), ] %>%
  scale_zi_nb_sce(sce = ., assay_name = "counts_vst", cpus = 4)
  

  
  pivot_longer(cols = -c("gene_id")) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~gene_id)
  

save(sce, file = "results/sce_PLS_DEA_cell_type.Rdata")
load(file = "results/sce_PLS_DEA_cell_type.Rdata", verbose = T)

altExp(sce, "PLS_DEA_cell_type") %>%
  rowData() %>%
  as_tibble(rownames = "id") %>%
  dplyr::select(id, predictor, predictor_force) %>%
  left_join(
    rowData(sce) %>%
      as_tibble()
    ) %>%
  dplyr::filter(predictor_force)

colData(sce) %>%
  as_tibble() %>%
  filter(male_invivo) %>%
  ggplot() +
  geom_point(aes(x = p_PLS_DEA_cell_type, y = sum, color = PLS_DEA_cell_type)) +
  facet_wrap(~day) +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes detected",
       x = " sum")

DEA_DEA_cell_type <- DEA(
  sce[, sce$male_invivo],
  test = "~ p_PLS_DEA_cell_type",
  formula = "count ~ p_PLS_DEA_cell_type + day + (1|batch)",
  assay_name = "counts_vst",
  cpus = 10
)

save(DEA_DEA_cell_type, file = "results/04_DEA.Rdata")
load(file = "results/04_DEA.Rdata", verbose = T)

rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo <- NA
rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo[sce$male_invivo] <- 
  get_genes_pval(DEA_DEA_cell_type)

rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo %>% 
  is.na() %>% 
  table()

DEA_DEA_cell_type_fixed <- DEA(
  sce[
    rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo %>% is.na(),
    sce$male_invivo
  ],
  test = "~ p_PLS_DEA_cell_type",
  formula = "count ~ p_PLS_DEA_cell_type + day + batch",
  assay_name = "counts_vst",
  cpus = 10
)

save(DEA_DEA_cell_type_fixed, file = "results/04_DEA_fixed.Rdata")

get_genes_pval(DEA_DEA_cell_type_fixed) %>%
  is.na() %>%
  table()

rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo[
  rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo %>% is.na()
] <- get_genes_pval(DEA_DEA_cell_type_fixed)

rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo %>% 
  is.na() %>% 
  table()

rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo_adj <- p.adjust(
  rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo,
  method = "BH"
)


save(sce, file = "results/sce_DEA_DEA_cell_type.Rdata")
load(file = "results/sce_DEA_DEA_cell_type.Rdata")

sce <- logNormCounts(sce, exprs_values = "counts_vst")
sce_DEA <- runPCA(sce[
    rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo_adj <= 0.05 &
      !is.na(rowData(sce)$pval_DEA_p_PLS_DEA_cell_type_male_invivo_adj), 
    sce$male_invivo
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

colData(sce) %>% 
  as_tibble() %>% 
  mutate(day = fct_relevel(day, "D15", "D90", "D136", "D593", "D1401")) %>%  
  filter(!is.na(ccr7), !is.na(il7ra)) %>% 
  ggplot(aes(
      x = ccr7,
      y = il7ra,
  )) +
  scale_color_manual(values=c("#0571b0", "#c90000")) +
  geom_point(aes(color = manual_cell_type), size = 2) +
  ggnewscale::new_scale_color() +
  scale_colour_gradient2(low = "blue", mid = "gray", high = "red",
            midpoint = 0.5) +
  geom_point(aes(color = p_PLS_DEA_cell_type), size = 1) +
  facet_wrap(~sex + day, ncol = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "il7ra",
       x = "ccr7")

colData(sce_DEA) %>% 
  as_tibble() %>% 
  left_join(
    reducedDim(sce_DEA, "TSNE") %>%
      as_tibble(rownames = "id")
    ) %>% 
  mutate(day = fct_relevel(day, "D15", "D136", "D593")) %>% 
  rename(c("V1" = "TSNE1", "V2" = "TSNE2")) %>% 
  pivot_longer(cols = c("p_PLS_DEA_cell_type", "p_PLS_surface_cell_type"),
               names_to = "PLS_type",
               values_to = "PLS_proba") %>%  
  mutate(PLS_type = fct_relevel(PLS_type, "p_PLS_surface_cell_type", "p_PLS_DEA_cell_type")) %>% 
  ggplot(aes(
      x = TSNE1,
      y = TSNE2,
  )) +
  scale_color_manual(values=c("#0571b0", "#c90000")) +
  labs(color = "manual cell-type") +
  geom_point(aes(color = manual_cell_type), size = 2) +
  ggnewscale::new_scale_color() +
  scale_colour_gradient2(low = "blue", mid = "gray", high = "red",
            midpoint = 0.5) +
  labs(color = "PLS cell-type probability") +
  geom_point(aes(color = PLS_proba), size = 1) +
  facet_wrap(~PLS_type + sex + day, ncol = 3) +
  theme_bw() +
  labs(y = "t-SNE1",
       x = "t-SNE2")

sce_hm <- sce

sce_DEA_hm <- sce[
  rowData(sce_DEA_hm)$gene %in% c(
    "GZMB",
    
    
    "CX3CR1",
    "CCL4",
    "GNLY",
    "GZMH",
    "KLRD1",
    "GZMG",
    "PRF1",
    "HOPX",
    "CCL5",
    "GZMK",
    "SELL",
    "IL7R",
    "LEF1",
    "TCF7",
    "LTB",
    "NELL2",
    "CCR7"
  ),
  sce$male_invivo]
rownames(sce_DEA_hm) <- rowData(sce_DEA_hm)$gene

plotHeatmap(
  sce_DEA_hm,
  features = rownames(sce_DEA_hm),
  order_columns_by = "p_PLS_DEA_cell_type",
  colour_columns_by = c("p_PLS_DEA_cell_type", "manual_cell_type", "day"),
  center = TRUE,
  symmetric = TRUE,
  zlim = c(-5, 5)
) 

colData(sce) %>% 
  as_tibble() %>% 
  mutate(pMEM = ifelse(PLS_DEA_cell_type %in% "EFF", 0, 1),
         pMEM = ifelse(is.na(PLS_DEA_cell_type), NA, pMEM)
         ) %>% 
  write_csv(path = "results/sce_DEA_DEA_cell_type_cellData.csv")
rowData(sce) %>% 
  as_tibble(rownames = "id") %>% 
  write_csv(path = "results/sce_DEA_DEA_cell_type_geneData.csv")
