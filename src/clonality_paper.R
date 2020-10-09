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
  DEA_clone_PCA_cell_type_size[[day]] <- DEA(
    sce_tmp,
    test = "~ (1|clone_id)",
    formula = "count ~ cell_type_pca_a + cell_type_pca_b + (1|clone_id)",
    assay_name = "counts_vst",
    cpus = 10
  )
  save(
    DEA_clone_PCA_cell_type_size,
    file = "results/2020_01_01_DEA_clone_PCA_cell_type_size_10.Rdata"
  )
}


rm_genes <- readr::read_delim(
    "data/2020_09_15_SmartSeq3/Genes_Exclude_Sept2020_LM.csv",
    delim = ";"
  ) %>% 
  janitor::clean_names()
load(file = "results/2020_01_02_clonality_paper_sce.Rdata", v = T)
logit_DEA_clone_PCA_cell_type_size <- list()
min_clone_size <- 10
future::plan("multiprocess", workers = 10)
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
  logi_DEA_clone_PCA_cell_type_size[[day]] <- assay(sce_tmp, "counts_vst") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "gene_name") %>% 
    filter(!(gene_name %in% rm_genes$geneid)) %>% 
    tidyr::nest(counts = !c(gene_name)) %>% 
    mutate(
      counts = purrr::map(.x = counts, .f = function(.x){
        tibble(
          id = colnames(.x),
          count = t(.x)[, 1],
          cell_type_pca_a = colData(sce_tmp)$cell_type_pca_a,
          cell_type_pca_b = colData(sce_tmp)$cell_type_pca_b,
          clone_id = colData(sce_tmp)$clone_id
          ) %>% 
          mutate(expressed = count > 0)
        }
      )
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
      .progress = TRUE)
    ) %>%
    mutate(
      test = furrr::future_map(.x = models, .f = function(.x){
        anova(.x$model, .x$model0, test = "Chisq") %>% 
          as_tibble() %>% 
          janitor::clean_names()
      },
      .progress = TRUE)
    ) %>% 
    tidyr::unnest(test) %>% 
    filter(!is.na(chisq))
  save(
    logit_DEA_clone_PCA_cell_type_size,
    file = "results/2020_01_01_logit_DEA_clone_PCA_cell_type_size_10.Rdata"
  )
}

## heatmap

min_clone_size <- 10
load(file = "results/2020_01_02_clonality_paper_sce.Rdata")
load(file = "results/2020_01_01_DEA_DEA_clone_size.Rdata")
load(file = "results/2020_01_01_DEA_clone_PCA_cell_type_size.Rdata", v=T)
d()
 e = "results/2020_01_01_DEA_clone_PCA_cell_type_size_10.Rdata"
##) adj pvalue

for (day in names(sce_day)) {[-1]
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
  table(rowData(sce_day[[day]])$pval_DEA_clone_PCA_cell_type_size_adj < 0.05) %>% print()
}
  day <- "593"
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

  future::plan("multiprocess", workers = 10)
  test <- assay(sce_DEA_hm, "logcounts") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "gene_name") %>% 
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
      }),
      count_mean = purrr::map(.x = counts, .f = function(.x){
          (.x$count)
      })
    ) %>%  
    unnest(count_var) %>% 
    unnest(count_mean) %>% 
    select(count_var, count_mean) %>% summary()
    filter(!(gene_name %in% rm_genes$geneid)) %>% 
    filter(count_var > 0.1) %>% 
    mutate(
      model = furrr::future_map(.x = counts, .f = function(.x){
        model.matrix(~-1 + clone_id, data = .x) %>% 
        glmnet::glmnet(
          x = .,
          y = (.x %>% pull(count)),
          lambda = glmnet::cv.glmnet(
            .,
            (.x %>% pull(count))
          )$lambda.1se
        )
      },
      .progress = TRUE)
    ) %>%
    mutate(coefs = map(model, tidy)) %>%
    select(-c(counts, model)) %>%
    tidyr::unnest(coefs)
    
    
  
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
  
cluster_row <- 
  assay(sce_DEA_hm, "logcounts") %>% 
  as.matrix() %>% 
  as_tibble(rownames = "gene_name") %>% 
  filter(!(gene_name %in% rm_genes$geneid)) %>% 
  tidyr::nest(counts = !c(gene_name)) %>% 
  mutate(
    counts = purrr::map(.x = counts, .f = function(.x){
      tibble(
        id = colnames(.x),
        count = t(.x)[, 1],
        clone_id = as.factor(colData(sce_DEA_hm)$clone_id)
        ) %>% 
        group_by(clone_id) %>% 
        dplyr::summarise(
          id = id,
          count = count,
          clone_id = clone_id,
          count_mean = max(count)
        ) %>% 
        ungroup()
      }
    )
  ) %>% 
  tidyr::unnest(counts) %>% 
  janitor::clean_names() %>% 
  dplyr::select(gene_name, clone_id, count_mean) %>% 
  tidyr::pivot_wider(
    id_cols = gene_name,
    names_from = clone_id,
    values_from = count_mean,
    values_fill = 0,
    values_fn = sum
  ) %>% 
  as.data.frame()
rownames(cluster_row) <- cluster_row$gene_name
  
sce_DEA_hm_plot <- sce_DEA_hm[rownames(sce_DEA_hm) %in% cluster_row$gene_name, ]
rowData(sce_DEA_hm_plot)$gene_order <- cluster_row[, -1] %>%
  dist(method = "canberra") %>% 
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
  main = day,
  cluster_rows = F,
) 
