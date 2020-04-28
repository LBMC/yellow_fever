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

fit <- PLS_fit(
  sce = sce,
  group_by = sce$manual_cell_type,
  genes = genes_PLS %>% pull(genes) %>% na.omit(),
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
  genes = genes_PLS %>% pull(genes) %>% na.omit(),
  features = genes_PLS %>% pull(proteins) %>% na.omit() %>% unique(),
  assay_name = "counts_vst",
  altExp_name = "PLS_surface_cell_type",
  cpus = 10
)
save(sce, file = "results/sce_PLS_surface_cell_type.Rdata")
load(file = "results/sce_PLS_surface_cell_type.Rdata", verbose = T)

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

save(fit, file = "results/03_DEA.Rdata")
load(file = "results/03_DEA.Rdata", verbose = T)

rowData(sce)$p_PLS_surface_cell_type_male_invivo <- NA
rowData(sce)$p_PLS_surface_cell_type_male_invivo[sce$male_invivo] <- 
  get_genes_pval(DEA_surface_cell_type)

save(sce, file = "results/sce_DEA_surface_cell_type.Rdata")