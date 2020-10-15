source("src/00_functions.R")
load(file = "results/2020_01_02_clonality_paper_sce.Rdata", v = T)

genes_PLS <- read_csv("data/2017_11_28_List_Laurent_Genes_PLS.csv") %>% 
  pivot_longer(cols = c("Genes_EFF", "Genes_MEM"),
               names_to = "type",
               values_to = "genes") %>% 
  dplyr::rename(proteins = `Protein_Markers`) %>% 
  pull(genes)

random_id <- runif(1) * 1e8
DEA_clone_PCA_cell_type_size <- list()
min_clone_size <- 10
for (day in names(sce_day)[-1]) {
  colData(sce_day[[day]])$clone_id <- colData(sce_day[[day]]) %>% 
    as_tibble() %>% 
    mutate(clone_id = ifelse(
      clone_id == 0,
      (max(clone_id) + 1):(max(clone_id) + length(which(clone_id == 0))),
      clone_id)) %>% 
    pull(clone_id)
  colData(sce_day[[day]])$clone_size <- colData(sce_day[[day]]) %>% 
    as_tibble() %>% 
    left_join(
      colData(sce_day[[day]]) %>%
        as_tibble() %>%
        group_by(clone_id) %>%
        dplyr::summarise(clone_size = n())
    ) %>%
    pull(clone_size)
  colData(sce_day[[day]])$cell_type_pca_a <- prcomp(
    assay(sce_day[[day]], "counts_vst")[
      rowData(sce_day[[day]])$gene_name %in% genes_PLS, 
    ]
  )$rotation[, 1]
  colData(sce_day[[day]])$cell_type_pca_b <- prcomp(
    assay(sce_day[[day]], "counts_vst")[
      rowData(sce_day[[day]])$gene_name %in% genes_PLS, 
    ]
  )$rotation[, 2]
  sce_tmp <- sce_day[[day]][,
    colData(sce_day[[day]])$clone_size >= min_clone_size
  ]
  colData(sce_tmp)$clone_id <- as.factor(colData(sce_tmp)$clone_id)
  colData(sce_tmp)$clone_suffle <- sample(
    colData(sce_tmp)$clone_id,
    size = length(colData(sce_tmp)$clone_id)
  )
  DEA_clone_PCA_cell_type_size[[day]] <- DEA(
    sce_tmp,
    test = "~ (1|clone_id)",
    formula = "count ~ cell_type_pca_a + cell_type_pca_b + (1|clone_suffle)",
    assay_name = "counts_vst",
    cpus = 32
  )
  save(
    DEA_clone_PCA_cell_type_size,
    file = str_c(
      "results/2020_01_01_DEA_clone_shuffled_PCA_cell_type_size_",
      min_clone_size, 
      "_",
      random_id,
      ".Rdata"
    )
  )
}

