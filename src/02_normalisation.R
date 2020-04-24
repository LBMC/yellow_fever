library(tidyverse)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(sctransform)

load(file = "results/sce_QC.Rdata", verbose = T)
colnames(sce) <- sce$id
sce$batch[sce$experiment %in% "P9997"] <- 40
sce$batch[sce$experiment %in% "P9998"] <- 41
sce$day[sce$experiment %in% c("P9997", "P9998")] <- "D1401"
sce$donor[sce$experiment %in% c("P9997", "P9998")] <- "YVF2003"
sce$antigen[sce$experiment %in% c("P9997", "P9998")] <- "A2"

sce$male_invivo <- sce$sex %in% "M" &
  sce$cell_number <= 1 &
  sce$day %in% c("D15", "D136", "D593") &
  sce$sequencing %in% "paired" &
  !(sce$batch %in% c(6:8)) &
  !(sce$id %in% c("P1299_1797", "P1299_1896"))
sce$female_invivo <- sce$sex %in% "F" &
  sce$cell_number <= 1 &
  !(sce$id %in% str_c("P1292_", 1097:1192)) &
  sce$day %in% c("D15", "D90")
sce$YVF2003_D1401 <- sce$day %in% c("D1401")
sce$male_in_vitro <- sce$day %in% c("InVitro") &
  sce$cell_number <= 1
sce$male_in_vitro_restim <- sce$day %in% c("In_Vitro_Restim") &
  sce$cell_number <= 1 &
  sce$sex %in% "M"
sce$female_in_vitro_restim <- sce$day %in% c("In_Vitro_Restim") &
  sce$cell_number <= 1 &
  sce$sex %in% "F"

colData(sce) %>% 
  as_tibble() %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(color = sex)) +
  facet_wrap(~day) +
  geom_density_2d(size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes detected",
       x = "counts sum")

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

rowData(sce) %>% 
  as_tibble() %>% 
  ggplot(aes(log_mean, log_var)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_density_2d(size = 0.3) +
  geom_abline(intercept = 0,
              slope = 1,
              color = 'red')

vst_out <- sctransform::vst(
  assays(sce)$counts_raw,
  cell_attr = colData(sce),
  latent_var = c("log_counts"),
  return_gene_attr = T,
  return_cell_attr = T,
  n_genes = NULL,
  show_progress = T)
sctransform::plot_model_pars(vst_out)
counts_norm <- correct(vst_out)

rowData(sce)$gene_mean <- rowMeans(assays(sce)$counts_raw[, colData(sce)$to_QC])
rowData(sce)$gene_var <- apply(assays(sce)$counts_raw[, colData(sce)$to_QC], 1, var)
rowData(sce)$detection_rate <- rowMeans(assays(sce)$counts_raw[, colData(sce)$to_QC] > 0)
rowData(sce)$detected <- rowSums(assays(sce)$counts_raw[, colData(sce)$to_QC]) > 0

rowData(sce) %>% 
  as_tibble() %>% 
  mutate(log_mean = log10(gene_mean),
         log_var = log10(gene_var),
         keep = rownames(sce) %in% rownames(counts_norm)) %>% 
  ggplot(aes(x = gene_mean, y = detection_rate, color = keep)) +
  geom_point(alpha = 0.3) +
  geom_density_2d(size = 0.3) +
  geom_line(data = poisson_model,
            aes(x = mean, y = detection_rate), color = "red") +
  scale_x_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes detected",
       x = "counts sum",
       color = "gene keeped by SCtransform")

sce <- sce[rownames(counts_norm), ]
assays(sce)$counts_vst <- counts_norm %>% Matrix::Matrix(sparse = T)

assays(sce)$counts_scaled <- scater::logNormCounts(
    sce,
    exprs_values = "counts_raw",
    log = F
  ) %>% 
  assay(., "normcounts") %>% 
  Matrix::Matrix(sparse = T)

rowData(sce)$gene_mean_vst <- rowMeans(assays(sce)$counts_vst)
rowData(sce)$gene_var_vst <- apply(assays(sce)$counts_vst, 1, var)
rowData(sce)$detection_rate_vst <- rowMeans(assays(sce)$counts_vst > 0)
rowData(sce)$detected_vst <- rowSums(assays(sce)$counts_vst) > 0
rowData(sce)$gene_mean_scaled <- rowMeans(assays(sce)$counts_scaled %>% as.matrix())
rowData(sce)$gene_var_scaled <- apply(assays(sce)$counts_scaled, 1, var)
rowData(sce)$detection_rate_scaled <- rowMeans(assays(sce)$counts_scaled > 0)
rowData(sce)$detected_scaled <- rowSums(assays(sce)$counts_scaled) > 0

x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- tibble(
  log_mean = x,
  mean = 10 ^ x,
  detection_rate = 1 - dpois(0, lambda = 10 ^ x)
)

rowData(sce) %>% 
  as_tibble() %>% 
  filter(detected > 0) %>% 
  ggplot(aes(x = gene_mean, y = detection_rate)) +
  geom_point(alpha = 0.3) +
  geom_density_2d(size = 0.3) +
  geom_point(aes(x = gene_mean_scaled, y = detection_rate_scaled), alpha = 0.3, color = "green") +
  geom_density_2d(aes(x = gene_mean_scaled, y = detection_rate_scaled), size = 0.3, color = "green") +
  geom_point(aes(x = gene_mean_vst, y = detection_rate_vst), alpha = 0.3, color = "red") +
  geom_density_2d(aes(x = gene_mean_vst, y = detection_rate_vst), size = 0.3, color = "red") +
  geom_line(data = poisson_model,
            aes(x = mean, y = detection_rate), color = "red") +
  scale_x_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes detected",
       x = "counts sum")

rowData(sce) %>% 
  as_tibble() %>% 
  mutate(log_mean = log10(gene_mean),
         log_var = log10(gene_var)) %>% 
  filter(detected > 0) %>% 
  ggplot(aes(x = gene_mean, y = gene_var, color = gene_var/gene_mean >= 1)) +
  geom_point(alpha = 0.3) +
  geom_density_2d(size = 0.3) +
  geom_point(aes(x = gene_mean_scaled, y = gene_var_scaled), alpha = 0.3, color = "green") +
  geom_density_2d(aes(x = gene_mean_scaled, y = gene_var_scaled), size = 0.3, color = "green") +
  geom_point(aes(x = gene_mean_vst, y = gene_var_vst), alpha = 0.3, color = "red") +
  geom_density_2d(aes(x = gene_mean_vst, y = gene_var_vst), size = 0.3, color = "red") +
  geom_abline(intercept = 0,
              slope = 1,
              color = 'red') +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "genes var",
       x = "gene mean")

save(sce, file = "results/sce_sctrf.Rdata", verbose = T)
