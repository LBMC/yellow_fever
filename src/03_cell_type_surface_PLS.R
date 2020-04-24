library(tidyverse)
library(SingleCellExperiment)
library(SummarizedExperiment)
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

source("src/00_functions.R")
fit <- PLS_fit(
  sce = sce,
  group_by = sce$manual_cell_type,
  genes = genes_PLS %>% pull(genes) %>% na.omit(),
  features = genes_PLS %>% pull(proteins) %>% na.omit() %>% unique(),
  assay_name = "counts_vst",
  cpus = 10
)
save(fit, file = "results/03_PLS_fit.Rdata")

load(file = "results/03_PLS_fit.Rdata", verbose = T)



