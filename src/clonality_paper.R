# Analysis for the clonality paper
source("src/00_functions.R")

# load data
count_data <- readr::read_csv("data/2020_09_15_SmartSeq3/JH_dex_umi.csv") %>% 
  dplyr::select(-X1) %>% 
  dplyr::left_join(
    readr::read_csv("data/2020_09_15_SmartSeq3/YFV2003_DEX_UMI.csv") %>% 
    dplyr::select(-X1),
    by = "Geneid"
  ) %>% 
  dplyr::left_join(
    readr::read_tsv("data/2020_09_15_SmartSeq3/P3128_merged_gene_counts.csv"),
    by = "Geneid"
  ) %>% 
  dplyr::rename(id = `Geneid`)
# %>% 
#   dplyr::filter(!(Geneid %in% (readr::read_delim(
#       "data/2020_09_15_SmartSeq3/Genes_Exclude_Sept2020_LM.csv",
#       delim = ";"
#     ) %>% pull(Geneid)))
#   )

row_data <- count_data %>%
  dplyr::select(id, gene_name)
count_data <- count_data %>% 
  dplyr::select(-c(id, gene_name))
col_data <- readr::read_csv("data/2020_09_15_SmartSeq3/Laurent_Meta_Dex_Late_Only.csv") %>% 
  dplyr::rename(id = `X1`)

col_data <- dplyr::tibble(
    id = count_data %>% colnames()
  ) %>% 
  dplyr::left_join(col_data) %>% 
  tidyr::drop_na() %>% 
  dplyr::rename(day = Day)

count_data <- count_data %>% 
  dplyr::select(all_of(col_data %>% pull(id)))

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(
    counts_raw = count_data %>% as.matrix() %>% Matrix::Matrix(sparse = T)
  ),
  colData = col_data,
  rowData = row_data
)
rm(count_data, col_data, row_data)

colData(sce)$donor_id <- 
  str_replace(colData(sce)$id, "([^_]*)_.*", "\\1")
save(sce, file = "results/2020_09_16_clonality_paper_sce.Rdata")
load(file = "results/2020_09_16_clonality_paper_sce.Rdata")

# QC
colData(sce)$gene_mean <- colMeans(assays(sce)$counts_raw)
colData(sce)$gene_var <- apply(assays(sce)$counts_raw, 2, var)
colData(sce)$detection_rate <- colMeans(assays(sce)$counts_raw > 0)
colData(sce)$detected <- colSums(assays(sce)$counts_raw) > 0

rowData(sce)$gene_mean <- rowMeans(assays(sce)$counts_raw)
rowData(sce)$gene_var <- apply(assays(sce)$counts_raw, 1, var)
rowData(sce)$detection_rate <- rowMeans(assays(sce)$counts_raw > 0)
rowData(sce)$detected <- rowSums(assays(sce)$counts_raw) > 0

colData(sce) %>%
  as_tibble() %>%
  ggplot(aes(x = gene_mean, y = detection_rate, color = donor_id)) +
  facet_wrap(~day) +
  geom_point(size = 0.3) +
  geom_density_2d(size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "detection rate",
       x = "genes mean expression")

colData(sce) %>%
  as_tibble() %>%
  ggplot(aes(x = gene_mean, y = gene_var, color = donor_id)) +
  facet_wrap(~day, scales = "free") +
  geom_point(size = 0.3) +
  geom_smooth() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes var",
       x = "genes mean")

# we compute the variance stabilizing transformation model
vst_out <- list()
rowData(sce) <- tibble(
  mean = rowMeans(assays(sce)$counts_raw),
  var = apply(assays(sce)$counts_raw, 1, var),
  detection_rate = rowMeans(assays(sce)$counts_raw > 2),
  log_mean = log10(mean),
  log_var = log10(var),
  n_counts =  rowSums(assays(sce)$counts_raw),
  log_counts = log10(n_counts)) %>%
  as_tibble() %>%
  cbind(rowData(sce), .)
colData(sce) <- tibble(
  mean = colMeans(assays(sce)$counts_raw),
  var = apply(assays(sce)$counts_raw, 2, var),
  detection_rate = colMeans(assays(sce)$counts_raw > 2),
  log_mean = log10(mean),
  log_var = log10(var),
  n_counts =  colSums(assays(sce)$counts_raw),
  log_counts = log10(n_counts)) %>% 
  as_tibble() %>% 
  cbind(colData(sce), .)

colnames(sce) <- colData(sce)$id
rownames(sce) <- rowData(sce)$id

save(sce, file = "results/2020_09_21_clonality_paper_sce.Rdata")
load(file = "results/2020_09_21_clonality_paper_sce.Rdata")

for(day in colData(sce)$day %>% as.factor() %>% levels()){
  vst_out[[day]] <- sctransform::vst(
    assay(sce, "counts_raw")[, colData(sce)$day %in% day],
    method = "nb",
    cell_attr = colData(sce)[colData(sce)$day %in% day, ],
    latent_var = c("log_counts"),
    return_gene_attr = T,
    return_cell_attr = T,
    n_genes = assays(sce)$counts_raw[, colData(sce)$day %in% day] %>% nrow(),
    n_cells = assays(sce)$counts_raw[, colData(sce)$day %in% day] %>% ncol(),
    show_progress = T)
}
save(vst_out, file = "results/vst_out.Rdata")
load(file = "results/vst_out.Rdata")

for(day in colData(sce)$day %>% as.factor() %>% levels()){
  sctransform::plot_model_pars(vst_out[[day]]) %>% print()
}

sce_day <- list()
for(day in colData(sce)$day %>% as.factor() %>% levels()){
  print(day)
  vst_norm <- correct(vst_out[[day]]) %>% Matrix::Matrix(sparse = T)
  vst_norm %>% dim() %>% print()
  sce_day[[day]] <- sce[rownames(vst_norm), colnames(vst_norm)]
  assays(sce_day[[day]])$counts_vst <- vst_norm
}
rm(vst_norm, vst_out)

for (day in names(sce_day)){
  rowData(sce_day[[day]]) <- tibble(
    mean = rowMeans(assays(sce_day[[day]])$counts_vst),
    var = apply(assays(sce_day[[day]])$counts_vst, 1, var),
    detection_rate = rowMeans(assays(sce_day[[day]])$counts_vst > 2),
    log_mean = log10(mean),
    log_var = log10(var),
    n_counts =  rowSums(assays(sce_day[[day]])$counts_vst),
    log_counts = log10(n_counts)) %>%
    as_tibble() %>%
    cbind(rowData(sce_day[[day]]), .)
  colData(sce_day[[day]]) <- tibble(
    mean = colMeans(assays(sce_day[[day]])$counts_vst),
    var = apply(assays(sce_day[[day]])$counts_vst, 2, var),
    detection_rate = colMeans(assays(sce_day[[day]])$counts_vst > 2),
    log_mean = log10(mean),
    log_var = log10(var),
    n_counts =  colSums(assays(sce_day[[day]])$counts_vst),
    log_counts = log10(n_counts)) %>% 
    as_tibble() %>% 
    cbind(colData(sce_day[[day]]), .)
}

col_data <- tibble() 
for (day in names(sce_day)){
  col_data <- col_data %>% 
    bind_rows(
      colData(sce_day[[day]]) %>%
      as_tibble()
    )
}
col_data %>% 
  ggplot(aes(x = mean, y = detection_rate, color = donor_id)) +
  facet_wrap(~day) +
  geom_point(size = 0.3) +
  geom_density_2d(size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "detection rate",
       x = "genes mean expression")

col_data %>%
  as_tibble() %>%
  ggplot(aes(x = mean, y = var, color = donor_id)) +
  facet_wrap(~day, scales = "free") +
  geom_point(size = 0.3) +
  geom_smooth() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes var",
       x = "genes mean")
rm(col_data)

save(sce_day, file = "results/2020_09_22_clonality_paper_sce.Rdata")
load(file = "results/2020_09_22_clonality_paper_sce.Rdata")

# PLS classification

fit_day <- list()
predict_day <- list()
load(file = "results/04_PLS_fit.Rdata", verbose = T)
load(file = "results/sce_DEA_surface_cell_type.Rdata", verbose = T)
load(file = "results/2020_09_30_fit_day.Rdata")
load(file = "results/2020_09_30_predict_day.Rdata")

for (day in names(sce_day)) {
  colData(sce_day[[day]])$manual_cell_type <- NA
  PLS_genes <- tibble(
    gene_id = rownames(altExp(fit$sce, "PLS_DEA_cell_type"))
    ) %>% 
    filter(gene_id %in% rownames(sce_day[[day]])) %>% 
    pull(gene_id)
  
  PLS_counts <- assay(sce_day[[day]], "counts_vst") %>%
    as.matrix() %>%
    as_tibble() %>% 
    dplyr::mutate(id = rownames(sce_day[[day]])) %>% 
    dplyr::filter(id %in% PLS_genes) %>% 
    left_join(
      assay(sce,
        "counts_vst")[
        rownames(sce) %in% PLS_genes, 
        !is.na(colData(sce)$manual_cell_type) &
          !(colData(sce)$id %in% colData(sce_day[[day]])$id)] %>% 
        as.matrix() %>%
        as_tibble() %>% 
        dplyr::mutate(id = rownames(sce[rownames(sce) %in% PLS_genes, ]))
    )
  
  PLS_coldata <- colData(sce_day[[day]]) %>% 
    as_tibble() %>% 
    mutate(day = as.character(day)) %>% 
    bind_rows(
      colData(sce)[!is.na(colData(sce)$manual_cell_type) &
          !(colData(sce)$id %in% colData(sce_day[[day]])$id), ] %>%
        as_tibble()
    )
  
  PLS_rowdata <- rowData(sce_day[[day]]) %>% 
    as_tibble() %>% 
    dplyr::filter(id %in% PLS_counts$id)
  
  PLS_counts <- PLS_counts %>%
      select(-c(id)) %>%
      as.matrix()
    
  rownames(PLS_counts) <- PLS_rowdata$id
  colnames(PLS_counts) <- PLS_coldata$id
  
  if (!(day %in% names(predict_day))){
    predict_day[[day]] <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
          counts_vst = PLS_counts
        ),
        colData = PLS_coldata,
        rowData = PLS_rowdata
      )
  }
  if (!(day %in% names(fit_day))){
    fit_day[[day]] <- PLS_fit(
      sce = sce,
      group_by = sce$manual_cell_type,
      genes = PLS_genes,
      assay_name = "counts_vst",
      altExp_name = "PLS_DEA_cell_type",
      force = PLS_genes,
      cpus = 8
    )
    save(fit_day, file = "results/2020_09_30_fit_day.Rdata")
  }
  if (!(day %in% names(predict_day))){
    predict_day[[day]] <- PLS_predict(
      sce = predict_day[[day]],
      fit = fit_day[[day]],
      group_by = colData(predict_day[[day]])$manual_cell_type,
      genes = PLS_genes,
      assay_name = "counts_vst",
      altExp_name = "PLS_DEA_cell_type",
      cpus = 8
    )
    save(predict_day, file = "results/2020_09_30_predict_day.Rdata")
  }
}

for (day in names(sce_day)) {
  colData(sce_day[[day]]) <- colData(sce_day[[day]]) %>% 
    as_tibble() %>% 
    left_join(
      colData(predict_day[[day]]) %>% 
        as_tibble() %>% 
        select(id, ends_with("PLS_DEA_cell_type"))
    ) %>% 
    select(ends_with("PLS_DEA_cell_type")) %>% 
    cbind(colData(sce_day[[day]]), .)
}

rm(fit, fit_day, PLS_coldata, PLS_counts, PLS_rowdata, predict_day, sce)
  
save(sce_day, file = "results/2020_01_01_clonality_paper_sce.Rdata")
load(file = "results/2020_01_01_clonality_paper_sce.Rdata")

for (day in names(sce_day)) {
  p <- colData(sce_day[[day]]) %>% 
    as_tibble() %>% 
    ggplot(aes(x = mean, y = p_PLS_DEA_cell_type, color = PLS_DEA_cell_type)) +
    facet_wrap(~day) +
    geom_point(size = 0.3) +
    geom_density_2d(size = 0.3) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    annotation_logticks() +
    labs(y = "detection rate",
         x = "genes mean expression") 
  print(p)
}

# heatmap


for (day in names(sce_day)) {
  assays(sce_day[[day]])$logcounts <- scater::logNormCounts(
      sce_day[[day]],
      exprs_values = "counts_raw",
      log = T
    ) %>% 
    assay(., "logcounts") %>% 
    Matrix::Matrix(sparse = T)
  sce_DEA_hm <- sce_day[[day]][
    rowData(sce_day[[day]])$gene_name %in% c(
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
    ]
  rownames(sce_DEA_hm) <- rowData(sce_DEA_hm)$gene_name
  plotHeatmap(
    sce_DEA_hm,
    features = rownames(sce_DEA_hm),
    order_columns_by = "p_PLS_DEA_cell_type",
    colour_columns_by = c("p_PLS_DEA_cell_type"),
    center = TRUE,
    symmetric = TRUE,
    zlim = c(-5, 5),
    main = day
  ) 
}

# DEA clonality ################################################################
load(file = "results/2020_01_01_clonality_paper_sce.Rdata")

DEA_DEA_clone_cell_type <- list()
for (day in names(sce_day)) {
  colData(sce_day[[day]])$clone_id <- colData(sce_day[[day]]) %>% 
    as_tibble() %>% 
    mutate(clone_id = ifelse(
      clone_id == 0,
      (max(clone_id) + 1):(max(clone_id) + length(which(clone_id == 0))),
      clone_id)) %>% 
    pull(clone_id)
  DEA_DEA_clone_cell_type[[day]] <- DEA(
    sce_day[[day]],
    test = "~ (1|clone_id)",
    formula = "count ~ p_PLS_DEA_cell_type + (1|clone_id)",
    assay_name = "counts_vst",
    cpus = 10
  )
  save(DEA_DEA_clone_cell_type, file = "results/2020_01_01_DEA_DEA_clone_cell_type.Rdata")
}

load("results/2020_01_01_DEA_DEA_clone_cell_type.Rdata", v = T)

## adj pvalue

for (day in names(sce_day)) {
  print(day)
  rowData(sce_day[[day]])$pval_DEA_clone_cell_type <- NA
  rowData(sce_day[[day]])$pval_DEA_clone_cell_type <- 
    get_genes_pval(DEA_DEA_clone_cell_type[[day]], sce_day[[day]])
  
  rowData(sce_day[[day]])$pval_DEA_clone_cell_type %>% 
    is.na() %>% 
    table() %>%
    print()
  
  rowData(sce_day[[day]])$pval_DEA_clone_cell_type_adj <- p.adjust(
    rowData(sce_day[[day]])$pval_DEA_clone_cell_type,
    method = "BH"
  )
  table(rowData(sce_day[[day]])$pval_DEA_clone_cell_type_adj < 0.05) %>% print()
}

save(sce_day, file = "results/2020_01_02_clonality_paper_sce.Rdata")
load(file = "results/2020_01_02_clonality_paper_sce.Rdata")

genes_PLS <- read_csv("data/2017_11_28_List_Laurent_Genes_PLS.csv") %>% 
  pivot_longer(cols = c("Genes_EFF", "Genes_MEM"),
               names_to = "type",
               values_to = "genes") %>% 
  dplyr::rename(proteins = `Protein_Markers`) %>% 
  pull(genes)

gc()
DEA_clone_PCA_cell_type_size <- list()
min_clone_size <- 3
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
  DEA_clone_PCA_cell_type_size[[day]] <- DEA(
    sce_tmp,
    test = "~ (1|clone_id)",
    formula = "count ~ cell_type_pca_a + cell_type_pca_b + (1|clone_id)",
    assay_name = "counts_vst",
    cpus = 10
  )
  save(
    DEA_clone_PCA_cell_type_size,
    file = "results/2020_01_01_DEA_clone_PCA_cell_type_size_3.Rdata"
  )
}
load(file = "results/2020_01_01_DEA_clone_PCA_cell_type_size_3.Rdata", v =T)

################################################################################
################## shuffling exp ###############################################

load(file = "results/2020_01_02_clonality_paper_sce.Rdata", v = T)

genes_PLS <- read_csv("data/2017_11_28_List_Laurent_Genes_PLS.csv") %>% 
  pivot_longer(cols = c("Genes_EFF", "Genes_MEM"),
               names_to = "type",
               values_to = "genes") %>% 
  dplyr::rename(proteins = `Protein_Markers`) %>% 
  pull(genes)

DEA_clone_PCA_cell_type_size <- list()
for (min_clone_size in c(10, 3)) {
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
    cbind(
      colData(sce_tmp)$clone_id,
      colData(sce_tmp)$clone_suffle
    )
    table(colData(sce_tmp)$clone_id)
    table(colData(sce_tmp)$clone_suffle)
    DEA_clone_PCA_cell_type_size[[day]] <- DEA(
      sce_tmp[1:3],
      test = "~ (1|clone_suffle)",
      formula = "count ~ cell_type_pca_a + cell_type_pca_b + (1|clone_suffle)",
      assay_name = "counts_vst",
      cpus = 2
    )
    save(
      DEA_clone_PCA_cell_type_size,
      file = str_c(
        "results/2020_01_01_DEA_clone_shuffled_PCA_cell_type_size_",
        min_clone_size, 
        ".Rdata"
      )
    )
  }
}

# load results
load(file = "results/2020_01_02_clonality_paper_sce.Rdata", v = T)
files_results <- list.files(path = "results/") %>% 
  .[str_detect(., "2020_01_01_DEA_clone_shuffled_PCA_cell_type_size_*")]
gene_to_keep <- read_csv("data/2020_10_19_JH_dex_Normalized_MW.csv") %>%
  select(-X1) %>%
  colnames()
test_infos <- tibble()
for (file_results in files_results) {
  for (day in names(sce_day)[-1]) {
    print(file_results)
    load(str_c("results/", file_results))
    print(day)
    test_infos <- tibble(
      p_value = get_genes_pval(DEA_clone_PCA_cell_type_size[[day]], sce_day[[day]])[
        rownames(sce_day[[day]]) %in% gene_to_keep
      ]
    ) %>% 
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        converge = !is.na(p_value),
        day = day,
        min_clone_size = file_results %>% str_replace(".*size_(.*)_.*", "\\1"),
        id =  file_results %>% str_replace(".*size_.*_(.*)\\..*", "\\1")
      ) %>% 
      bind_rows(test_infos)
  }
}

test <- test_infos %>%
  filter(converge) %>% 
  select(-c(converge, p_value)) %>% 
  group_by(day, min_clone_size, id) %>% 
  nest() %>% 
  ungroup() %>%
  bind_cols(
    rep(seq(from = 0.01, to = 0.5, by = 0.01), time = nrow(.)) %>%
      enframe() %>%
      pivot_wider(id_col = value)
  ) %>% 
  pivot_longer(cols = !c(day:data), values_to = "fdr") %>% 
  select(-name) %>% 
  mutate(
    error = lapply(X = as.list(1:nrow(.)), FUN = function(x, data, fdr){
      mean(data[[x]]$p_adj < fdr[x], na.rm = T)
    }, data = data, fdr = fdr)
  ) %>% 
  unnest(error)

test %>% 
  ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_line(aes(x = fdr, y = error, color = id)) +
  facet_wrap(~day + min_clone_size) +
  theme_bw()




################################################################################



genes_PLS <- read_csv("data/2017_11_28_List_Laurent_Genes_PLS.csv") %>% 
  pivot_longer(cols = c("Genes_EFF", "Genes_MEM"),
               names_to = "type",
               values_to = "genes") %>% 
  dplyr::rename(proteins = `Protein_Markers`) %>% 
  pull(genes)
rm_genes <- readr::read_delim(
    "data/2020_09_15_SmartSeq3/Genes_Exclude_Sept2020_LM.csv",
    delim = ";"
  ) %>% 
  janitor::clean_names()
load(file = "results/2020_01_02_clonality_paper_sce.Rdata", v = T)
logit_DEA_clone_PCA_cell_type_size <- list()
min_clone_size <- 10
options(future.globals.maxSize= 32000*1024^2)
# future::plan("multicore", workers = 4)
future::plan("sequential")

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
  print(day)
  print(table(colData(sce_tmp)$clone_id))
  logit_DEA_clone_PCA_cell_type_size[[day]] <- assay(sce_tmp, "counts_vst") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "gene_name") %>% 
    filter(!(gene_name %in% rm_genes$geneid)) %>% 
    tidyr::nest(counts = !c(gene_name)) %>% 
    mutate(
      counts = furrr::future_map(.x = counts, .f = function(.x){
        tibble(
          id = colnames(.x),
          count = t(.x)[, 1],
          cell_type_pca_a = colData(sce_tmp)$cell_type_pca_a,
          cell_type_pca_b = colData(sce_tmp)$cell_type_pca_b,
          clone_id = colData(sce_tmp)$clone_id
          ) %>% 
          mutate(expressed = as.numeric(count > 0))
        }
      , .progress = T)
    ) %>%
    mutate(
      count_var = purrr::map(.x = counts, .f = function(.x){
          var(.x$count)
      }),
    ) %>%  
    unnest(count_var) %>% 
    filter(count_var > 0) %>% 
    mutate(
      models = furrr::future_map(.x = counts, .f = function(.x){
        tryCatch({
          list(model0 = glm(
              expressed ~ cell_type_pca_a + cell_type_pca_b,
              data = .x,
              family = binomial
            ),
            model = lme4::glmer(
              expressed ~ cell_type_pca_a + cell_type_pca_b + (1|clone_id),
              data = .x,
              family = binomial,
              nAGQ = 0 
            )
          )
        },
        error = function(e){ list(model0 = NULL, model = NULL) }
        )
      }, .progress = T)
    ) %>%
    select(-c(counts, count_var)) %>%
    mutate(
      test = furrr::future_map(.x = models, .f = function(.x){
        if (is.null(.x$model0)) {
          return(NULL)
        }
        anova(.x$model0, .x$model1, test = "Chisq") %>% 
          as_tibble() %>% 
          janitor::clean_names()
      }, .progress = T)
    ) %>%
    select(c(gene_name, test)) %>%
    unnest(test) %>%
    filter(!is.na(pr_chi)) %>%
    mutate(
      pval_logit_DEA_clonality_PCA_cell_type_size_10 = pr_chi,
      pval_logit_DEA_clonality_PCA_cell_type_size_10_adj = p.adjust(
        pr_chi,
        method = "BH"
      )
    ) %>%
    select(c(gene_name, starts_with("pval")))
  save(
    logit_DEA_clone_PCA_cell_type_size,
    file = "results/2020_10_13_logit_DEA_clone_PCA_cell_type_size_10.Rdata"
  )
}

load(
  file = "results/2020_10_13_logit_DEA_clone_PCA_cell_type_size_10.Rdata",
  v = T
)

for (day in names(sce_day)[-1]) {
  logit_DEA_clone_PCA_cell_type_size[[day]] %>%
    mutate(signif = pval_logit_DEA_clonality_PCA_cell_type_size_10_adj < 0.05) %>%
    summary() %>%
    print()
}



## heatmap

min_clone_size <- 10
load(file = "results/2020_01_02_clonality_paper_sce.Rdata")
load(file = "results/2020_01_01_DEA_DEA_clone_size.Rdata")
load(file = "results/2020_01_01_DEA_clone_PCA_cell_type_size_3.Rdata", v = T)
load(file = "results/2020_01_01_DEA_clone_PCA_cell_type_size_10.Rdata", v = T)
## adj pvalue

for (min_clone_size in c(3, 10)) {
  load(
    file = str_c(
      "results/2020_01_01_DEA_clone_PCA_cell_type_size_", 
      min_clone_size,
      ".Rdata"
      ), v = T
  )
  
for (day in names(sce_day)[-1]) {
  print(day)
  rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size <- NA
  rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size <- 
    get_genes_pval(DEA_clone_PCA_cell_type_size[[day]], sce_day[[day]])
  
  rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size %>% 
    is.na() %>% 
    table() %>%
    print()
  
  rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size_adj <- p.adjust(
    rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size,
    method = "BH"
  )
  table(rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size_adj[!is.na(rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size_adj)] < 0.05) %>% print()
}

for (day in names(sce_day)[-1]) {
  assays(sce_day[[day]])$logcounts <- scater::logNormCounts(
      sce_day[[day]],
      exprs_values = "counts_raw",
      log = T
    ) %>% 
    assay(., "logcounts") %>% 
    Matrix::Matrix(sparse = T)
  sce_DEA_hm <- sce_day[[day]][
      !is.na(rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size_adj) &
      rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size_adj < 0.05
    ]
  sce_DEA_hm %>% dim()
  colData(sce_DEA_hm) <- colData(sce_DEA_hm) %>% 
    as_tibble() %>% 
    left_join(
      colData(sce_DEA_hm) %>%
        as_tibble() %>%
        group_by(clone_id) %>%
        dplyr::summarise(clone_size = n())
    ) %>%
    select(clone_size) %>% 
    cbind(colData(sce_DEA_hm), .)
  sce_DEA_hm <- sce_DEA_hm[, colData(sce_DEA_hm)$clone_size >= min_clone_size]
  colData(sce_DEA_hm)$clone_id <-  colData(sce_DEA_hm)$clone_id %>% as.factor()
  
  library("tidymodels")
  library("tidyverse")
  library("glmnet")
  library("furrr")
  
  rm_genes <- readr::read_delim(
      "data/2020_09_15_SmartSeq3/Genes_Exclude_Sept2020_LM.csv",
      delim = ";"
    ) %>% 
    janitor::clean_names()

  cluster_row <- assay(sce_DEA_hm, "logcounts") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "gene_name") %>% 
    filter(!(gene_name %in% rm_genes$geneid)) %>% 
    tidyr::nest(counts = !c(gene_name)) %>% 
    mutate(
      counts = purrr::map(.x = counts, .f = function(.x){
        tibble(
          id = colnames(.x),
          count = t(.x)[, 1],
          clone_id = colData(sce_DEA_hm)$clone_id)
        }
      )
    ) %>%
    mutate(
      count_var = purrr::map(.x = counts, .f = function(.x){
          var(.x$count)
      })
    ) %>%  
    unnest(count_var) %>% 
    mutate(
      model = furrr::future_map(.x = counts, .f = function(.x){
        lm(count ~ clone_id, data = .x)
      },
      .progress = TRUE)
    ) %>%
    mutate(coefs = map(model, tidy)) %>%
    select(-c(counts, model)) %>%
    tidyr::unnest(coefs) %>% 
    janitor::clean_names() %>% 
    dplyr::mutate(term = ifelse(
      term == "(Intercept)",
      colData(sce_DEA_hm)$clone_id %>% as.factor() %>% levels() %>% .[1] %>% 
        str_c("clone_id", .),
      term)
    ) %>% 
    dplyr::mutate(
      term = str_replace(term, "clone_id(.*)", "\\1")
    ) %>% 
    dplyr::select(gene_name, term, estimate) %>% 
    tidyr::pivot_wider(
      id_cols = gene_name,
      names_from = term,
      values_from = estimate,
      values_fill = 0,
      values_fn = sum
    ) %>% 
    as.data.frame()
  rownames(cluster_row) <- cluster_row$gene_name
  
  sce_DEA_hm_plot <- sce_DEA_hm[rownames(sce_DEA_hm) %in% cluster_row$gene_name, ]
  rowData(sce_DEA_hm_plot)$gene_order <- cluster_row[, -1] %>%
    dist(method = "manhattan") %>% 
    hclust() %>% .$order
  rownames(sce_DEA_hm_plot) <- rowData(sce_DEA_hm_plot)$gene_name
  sce_DEA_hm_plot <- sce_DEA_hm_plot[rowData(sce_DEA_hm_plot)$gene_order, ]
  plotHeatmap(
    sce_DEA_hm_plot[!(rowData(sce_DEA_hm_plot)$gene_name %in% rm_genes),],
    features = rownames(sce_DEA_hm_plot),
    order_columns_by = c("clone_id", "p_PLS_DEA_cell_type"),
    colour_columns_by = c("clone_id", "p_PLS_DEA_cell_type"),
    center = TRUE,
    symmetric = TRUE,
    zlim = c(-5, 5),
    main = str_c("DEG for clonality (size >= ",
                 min_clone_size,
                 ") accounting for cell_type",
                 day),
    cluster_rows = F,
  ) 
  
  bet <- assay(sce_DEA_hm, "logcounts") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "gene_name") %>% 
    filter(!(gene_name %in% rm_genes$geneid)) %>% 
    select(-gene_name) %>% 
    as.matrix() %>% 
    t() %>% 
    ade4::dudi.pca(scannf = F, nf = 100) %>% 
    ade4::bca(dudi, fac = colData(sce_DEA_hm)$clone_id, scannf = F, nf = 2)
  bca <- bet$ls %>% 
    janitor::clean_names() %>% 
    tibble(rownames = colData(sce_DEA_hm)$id) %>% 
    mutate(clone_id = colData(sce_DEA_hm)$clone_id) %>% 
    group_by(clone_id) %>% 
    nest() %>% 
    mutate(
      clone_x = map(data, function(x){
        sum(x$cs1) / nrow(x)
      }) %>% unlist(),
      clone_y = map(data, function(x){
        sum(x$cs2) / nrow(x)
      }) %>% unlist())
  bca <- bca %>%
    unnest(data) %>% 
    ggplot(aes(x = cs1, y = cs2, color = clone_id, group = clone_id)) +
    geom_point() +
    stat_ellipse(type = "t") +
    geom_segment(
      aes(x = cs1, y = cs2, xend = clone_x, yend = clone_y), color = "gray"
    ) +
    ggrepel::geom_label_repel(
      data = bca,
      aes(x = clone_x, y = clone_y, color = clone_id, group = clone_id, label = clone_id)
    ) +
    labs(
      title = str_c("Between Class Analysis ", day, " (clone size >= ",
                    min_clone_size, ")"),
      x = "PC1",
      y = "PC2"
    ) +
    theme_bw()
  bca %>% print()
  
  rowData(sce_day[[day]]) %>% 
    as_tibble(rownames = "id") %>% 
    write_csv(
      path = str_c(
        "results/2020_10_08_DEA_",
        day,
        "_clone_PCA_cell_type_size_",
        min_clone_size,
        ".csv")
    )
}
}

load(file = "results/2020_01_02_clonality_paper_sce.Rdata", v = T)
files_results <- list.files(path = "results/") %>% 
  .[str_detect(., "2020_01_01_DEA_clone.*_PCA_cell_type_size_.*")]
gene_to_keep <- read_csv("data/2020_10_19_JH_dex_Normalized_MW.csv") %>%
  select(-X1) %>%
  colnames()
for (file_results in files_results) {
  for (day in names(sce_day)[-1]) {
    print(file_results)
    load(str_c("results/", file_results))
    print(day)
    pval_name <- file_results %>%
      str_replace(".*(DEA_.*)\\.Rdata", "\\1")
    rowData(sce_day[[day]])[[str_c(pval_name, "_pval")]] <- get_genes_pval(
      DEA_clone_PCA_cell_type_size[[day]], sce_day[[day]]
    )
    rowData(sce_day[[day]])[[str_c(pval_name, "_padj")]] <- NA
    rowData(sce_day[[day]])[[str_c(pval_name, "_padj")]][
        rownames(sce_day[[day]]) %in% gene_to_keep
      ] <- p.adjust(
          rowData(sce_day[[day]])[[str_c(pval_name, "_pval")]][
            rownames(sce_day[[day]]) %in% gene_to_keep
          ],
          method = "BH"
        )
  }
}

for (day in names(sce_day)[-1]) {
  rowData(sce_day[[day]]) %>% 
    as_tibble(rownames = "id") %>% 
    write_csv(
      path = str_c(
        "results/2020_10_08_DEA_",
        day,
        "_clone_PCA_cell_type.csv")
    )
}
