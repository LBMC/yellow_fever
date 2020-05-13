source("src/00_functions.R")

load(file = "results/sce_QC.Rdata", verbose = T)

# annotation of the P9997 and P9998 experiment
sce$batch[sce$experiment %in% "P9997"] <- 40
sce$batch[sce$experiment %in% "P9998"] <- 41
sce$day[sce$experiment %in% c("P9997", "P9998")] <- "D1401"
sce$donor[sce$experiment %in% c("P9997", "P9998")] <- "YVF2003"
sce$antigen[sce$experiment %in% c("P9997", "P9998")] <- "A2"

# sortcut annotation for the different subsets of data
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


# We add localisation information to the genes

hub_infos <- AnnotationHub::AnnotationHub()
hub_ids <- mcols(hub_infos) %>%
  data.frame() %>% 
  rownames_to_column(var = "id") %>% # we keep the rownames
  as_tibble() %>% 
  dplyr::filter(
    dataprovider %in% "Ensembl" & # Ensembl annotation
      species %in% c("Homo sapiens"), # for the two species we want
    genome %in% c("GRCh38"), # on the right genome
    str_detect(title, "99"), # on the right revision
    rdataclass %in% "EnsDb",
  ) %>% 
  dplyr::select(id, species)

pull_loc <- function(id, sce, hub_infos){
  AnnotationDbi::mapIds(
    hub_infos[[id]],
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SEQNAME")
}

merge_loc <- function(sce, id, hub_infos){
  sapply(id %>% pull(id), pull_loc, sce = sce, hub_infos = hub_infos) %>% 
    as_tibble() %>% 
    unite(col = "chr_pos") %>% 
    mutate(chr_pos = ifelse(chr_pos == "NA", NA, chr_pos)) %>% 
    pull(chr_pos)
}

rowData(sce)$chr_pos = merge_loc(sce, hub_ids, hub_infos)
rowData(sce)$is_genomic <- rowData(sce)$chr_pos %in% c(as.character(1:22), "X", "Y")

# MT genes
rowData(sce)$gene[!rowData(sce)$is_genomic & !is.na(rowData(sce)$chr_pos)]

altExp(sce, "MT") <- sce[!rowData(sce)$is_genomic & !is.na(rowData(sce)$chr_pos), ]
sce <- sce[!(!rowData(sce)$is_genomic & !is.na(rowData(sce)$chr_pos)), ]

# we look at the unannotated genes
table(is.na(rowData(sce)$chr_pos))
rowData(sce)$gene[is.na(rowData(sce)$chr_pos)]

# among those which ones are in another ensembl id ?
duplicated_annotation <- rowData(sce)$gene[is.na(rowData(sce)$chr_pos)][
  rowData(sce)$gene[is.na(rowData(sce)$chr_pos)] %in%
    rowData(sce)$gene[rowData(sce)$is_genomic]
]

assay(sce, "logcounts") <- NULL
for (x in duplicated_annotation) {
  print(x)
  assay(sce, "counts_raw")[
    rowData(sce)$gene %in% x,
  ] <- colSums(
    assay(sce, "counts_raw")[
      rowData(sce)$gene %in% x,
    ]
  )
}
sce <- sce[!is.na(rowData(sce)$chr_pos), ]

# QC plot
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

# update QC metrix on the new selection of cell / genes after 01_QC.R
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
colnames(sce) <- sce$id

rowData(sce) %>% 
  as_tibble() %>% 
  ggplot(aes(log_mean, log_var)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_density_2d(size = 0.3) +
  geom_abline(intercept = 0,
              slope = 1,
              color = 'red')

# we compute the variance stabilizing transformation model
vst_out <- sctransform::vst(
  assay(sce, "counts_raw"),
  method = "nb",
  cell_attr = colData(sce),
  latent_var = c("log_counts"),
  return_gene_attr = T,
  return_cell_attr = T,
  n_genes = assays(sce)$counts_raw %>% nrow(),
  n_cells = assays(sce)$counts_raw %>% ncol(),
  show_progress = T)
save(vst_out, file = "results/vst_out.Rdata")
load(file = "results/vst_out.Rdata")
sctransform::plot_model_pars(vst_out)

# we can look at the most variable genes
vst_out$gene_attr %>% 
  as_tibble(rownames = "id") %>% 
  mutate(gene = rowData(sce)$gene[match(id, rownames(sce))]) %>% 
  arrange(-residual_variance) %>% 
  head(20)

counts_norm <- correct(vst_out)

rowData(sce)$gene_mean <- rowMeans(assays(sce)$counts_raw[, colData(sce)$to_QC])
rowData(sce)$gene_var <- apply(assays(sce)$counts_raw[, colData(sce)$to_QC], 1, var)
rowData(sce)$detection_rate <- rowMeans(assays(sce)$counts_raw[, colData(sce)$to_QC] > 0)
rowData(sce)$detected <- rowSums(assays(sce)$counts_raw[, colData(sce)$to_QC]) > 0

x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- tibble(
  log_mean = x,
  mean = 10 ^ x,
  detection_rate = 1 - dpois(0, lambda = 10 ^ x)
)
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
rowData(sce)$id <- rownames(sce)

rowData(sce) <- rowData(sce) %>%
  as_tibble() %>% 
  right_join(
    vst_out$gene_attr %>% 
    as_tibble(rownames = "id") %>% 
    rename_all(funs(str_replace(., "^", "vst_"))),
    by = c("id" = "vst_id"))

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
  dplyr::filter(detected > 0) %>% 
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
  dplyr::filter(detected > 0) %>% 
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
load(file = "results/sce_sctrf.Rdata", verbose = T)
