source("src/00_functions.R")

load(file = "results/sce_DEA_DEA_cell_type.Rdata")
large_clone <-
  readxl::read_xlsx("data/2020_05_12_In_Vivo_Clones_DEA_Test_LM_May2020.xlsx") %>% 
    rename(c(
      Donor = "sex",
      Day = "day",
      `Sample:ID` = "id",
      Clone_ID = "clonality"
    ))

big_M_clone <- colData(sce) %>% 
  as_tibble() %>% 
  mutate(sex = ifelse(day == "D1401", "F", sex)) %>% 
  dplyr::filter(!is.na(clonality), !is.na(day), !is.na(sex)) %>% 
  dplyr::group_by(clonality, day, sex)  %>% 
  dplyr::mutate(clone_size = n()) %>% 
  dplyr::arrange(desc(clone_size)) %>% 
  dplyr::filter(male_invivo, clone_size > 3) %>% 
  pull(id)

DEA_large_clone_M <- DEA(
  sce = sce[
      sce$id %in% big_M_clone & sce$male_invivo
    ],
  test = "(1|clonality)",
  formula = "count ~ (1|clonality) + p_PLS_DEA_cell_type + day + (1|batch)",
  zi_formula = "count ~ (1|batch)",
  assay_name = "counts_vst",
  cpus = 6
)

save(DEA_large_clone_M, file = "results/05_DEA_large_clone_M.Rdata")
load(file = "results/05_DEA_large_clone_M.Rdata", verbose = T)

rowData(sce)$pval_DEA_large_clone_M <- NA
rowData(sce)$pval_DEA_large_clone_M <- 
  get_genes_pval(DEA_large_clone_M, sce[, sce$id %in% big_M_clone & sce$male_invivo])

rowData(sce)$pval_DEA_large_clone_M_adj <- p.adjust(
  rowData(sce)$pval_DEA_large_clone_M,
  method = "BH"
)

rowData(sce)$pval_DEA_large_clone_M %>% 
  is.na() %>% 
  table()

rowData(sce) %>% 
  as_tibble() %>% 
  mutate(test = pval_DEA_large_clone_M_adj < 0.05) %>% 
  pull(test) %>% 
  table()

sce <- logNormCounts(sce, exprs_values = "counts_vst")
sce_DEA <- runPCA(sce[
    rowData(sce)$pval_DEA_large_clone_M_adj <= 0.05 &
      !is.na(rowData(sce)$pval_DEA_large_clone_M_adj), 
    sce$id %in% big_M_clone & sce$male_invivo
  ])
sce_DEA <- runTSNE(sce_DEA)

plotPCA(sce_DEA,
  colour_by = "clonality")
plotTSNE(sce_DEA,
  colour_by = "clonality")
plotPCA(sce_DEA,
  colour_by = "day")


# plot BCA
library(ade4)
for (day in (sce_DEA$day %>% as.factor() %>% levels())){
  g <- assay(sce_DEA[, sce_DEA$day %in% day], "logcounts") %>%
    as.matrix() %>% 
    t() %>% 
    as_tibble() %>% 
    dudi.pca(scannf = F, nf = 4) %>%
    bca(sce_DEA[, sce_DEA$day %in% day]$clonality %>% as.factor(), scannf = F) %>% 
    .$ls %>% 
    as_tibble() %>% 
    mutate(clonality = sce_DEA[, sce_DEA$day %in% day]$clonality) %>% 
    ggplot() + 
      stat_ellipse(
          aes(
            x = CS1,
            y = CS2,
            color = clonality,
            fill = clonality
          ),
          geom = "polygon", alpha = 0.25, level = 0.8) +
      geom_point(
          aes(
            x = CS1,
            y = CS2,
            color = clonality,
            fill = clonality
          )) +
      labs(title = day) +
      theme_bw()
  print(g)
}

# In vitro data
DEA_large_clone_P1902 <- DEA(
  sce = sce[
      sce$experiment %in% "P1902" & sce$male_in_vitro
    ],
  test = "(1|clonality)",
  formula = "count ~ (1|clonality) + (1|batch)",
  assay_name = "counts_vst",
  cpus = 6
)

save(DEA_large_clone_P1902, file = "results/05_DEA_large_clone_P1902.Rdata")
load(file = "results/05_DEA_large_clone_P1902.Rdata", verbose = T)

rowData(sce)$pval_DEA_large_clone_P1902 <- NA
rowData(sce)$pval_DEA_large_clone_P1902 <- 
  get_genes_pval(DEA_large_clone_P1902, sce[, sce$experiment %in% "P1902" & sce$male_in_vitro])

rowData(sce)$pval_DEA_large_clone_M_adj <- p.adjust(
  rowData(sce)$pval_DEA_large_clone_M,
  method = "BH"
)

rowData(sce)$pval_DEA_large_clone_M %>% 
  is.na() %>% 
  table()

rowData(sce) %>% 
  as_tibble() %>% 
  mutate(test = pval_DEA_large_clone_M_adj < 0.05) %>% 
  pull(test) %>% 
  table()