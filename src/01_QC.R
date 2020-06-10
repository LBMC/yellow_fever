source("src/00_functions.R")

load(file = "results/sce_annot.Rdata", verbose = T)

# We move the ERCC data to an altExp in the sce object
altExp(sce, "ERCC") <- sce[stringr::str_detect(rownames(sce), "ERCC-\\d+"), ]
sce <- sce[!stringr::str_detect(rownames(sce), "ERCC-\\d+"), ]
# 92 ERCC

# add quality metric
sce <- scater::addPerFeatureQC(
  sce,
  exprs_values = "counts_raw"
)
rowData(sce)$detection_rate = rowMeans(
  assays(sce)$counts_raw > 0)
sce <- scater::addPerCellQC(
  sce,
  exprs_values = "counts_raw"
)

# add quality metric for the ERCC
altExp(sce, "ERCC") <- scater::addPerFeatureQC(
  altExp(sce, "ERCC"),
  exprs_values = "counts_raw"
)
rowData(altExp(sce, "ERCC"))$detection_rate = rowMeans(
  assays(altExp(sce, "ERCC"))$counts_raw > 0)
altExp(sce, "ERCC") <- scater::addPerCellQC(
  altExp(sce, "ERCC"),
  exprs_values = "counts_raw"
)

# QC plots

colData(sce) %>% 
  as_tibble() %>% 
  ggplot() +
  geom_histogram(aes(x = sum, fill = day)) +
  theme_bw() +
  labs(x = "count sum")

colData(sce) %>% 
  as_tibble() %>% 
  bind_rows(colData(altExp(sce, "ERCC")) %>% 
              as_tibble()) %>% 
  ggplot() +
  geom_point(aes(x = sum, y = altexps_ERCC_sum, color = day)) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(x = "count sum",
       y = "ERCC count sum")

colData(sce) %>% 
  as_tibble() %>% 
  bind_rows(colData(altExp(sce, "ERCC")) %>% 
              as_tibble()) %>% 
  ggplot() +
  geom_point(aes(x = detected, y = altexps_ERCC_detected, color = day)) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(x = "genes detected",
       y = "ERCC detected")

colData(sce) %>% 
  as_tibble() %>% 
  bind_rows(colData(altExp(sce, "ERCC")) %>% 
              as_tibble()) %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(color = experiment)) +
  geom_density_2d(size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes detected",
       x = "counts sum")


# knee plot
bcrank <- DropletUtils::barcodeRanks(assays(sce)$counts_raw[rowData(sce)$detected, ])
colData(sce)
tibble(
  rank = rank(-colSums(assays(sce)$counts_raw[rowData(sce)$detected, ] > 0)),
  expressed = colSums(assays(sce)$counts_raw[rowData(sce)$detected, ] > 0),
  rank_cutoff = max(bcrank$rank[bcrank$total > metadata(bcrank)[["inflection"]]]),
  cell_number = colData(sce)$cell_number
  ) %>%
  distinct() %>%
  ggplot(aes(y = rank, x = expressed)) +
  geom_vline(aes(xintercept = rank_cutoff), col = "red") +
  geom_point(aes(color = cell_number == 0)) +
  geom_line() +
  theme_bw() +
  labs(x = "expressed genes",
       y = "rank")

# arbitraty cell filter
weird_D15 <- c("P1299_1105", "P1299_1106", "P1299_1111", "P1299_1112", "P1299_1117", "P1299_1129", "P1299_1133", "P1299_1150", "P1299_1151", "P1299_1185", "P1299_1222", "P1299_1263", "P1299_1284", "P1299_1297", "P1299_1299", "P1299_1313", "P1299_1328", "P1299_1336", "P1299_1345", "P1299_1356", "P1299_1364", "P1299_1371", "P1299_1390", "P1299_1397", "P1299_1404", "P1299_1416", "P1299_1429", "P1299_1432", "P1299_1437", "P1299_1445", "P1299_1457", "P1299_1465", "P1299_1466", "P1299_1473", "P1299_1478", "P1299_1770", "P1299_1772", "P1299_1781", "P1299_1795", "P1299_1802", "P1299_1803", "P1299_1810", "P1299_1818", "P1299_1819", "P1299_1826", "P1299_1838", "P1299_1843", "P1299_1847", "P1299_1850", "P1299_1861", "P1299_1881", "P1299_1882", "P1299_1884", "P1299_1908", "P1299_1913", "P1299_1921", "P1299_1922", "P1299_1928", "P1299_1949", "P1299_2012", "P1299_2017", "P1299_2035", "P1299_2052", "P1299_2054", "P1299_2056")
sce <- sce[, !(colData(sce)$id %in% weird_D15)]
bad_F_cells <- str_c("P1292_", 1097:1192)
sce <- sce[, !(colData(sce)$id %in% bad_F_cells)]

# we add a logcounts matrix
logcounts(sce) <- log1p(assays(sce)$counts_raw)

# svm bagging QC analysis for Male invivo data

# we work on the male in vivo data
colData(sce)$to_QC <- colData(sce) %>% 
  as_tibble() %>% 
  dplyr::mutate(
    to_QC = sex %in% "M" &
    cell_number <= 1 &
    day %in% c("D15", "D136", "D593") &
    sequencing %in% "paired" &
    !(batch %in% c(6:8)) &
    !(id %in% c("P1299_1797", "P1299_1896"))
  ) %>% pull(to_QC)

colData(sce)$QC_score <- -1
colData(sce)$QC_score[colData(sce)$to_QC] <- QC_score(
  sce = sce[, colData(sce)$to_QC],
  assay_name = "logcounts",
  ncell_name = "cell_number",
) - 1

save(sce, file = "results/sce_QC_M.Rdata")
load(file = "results/sce_QC_M.Rdata", verbose = T)

colData(sce)[colData(sce)$to_QC, ] %>% 
  as_tibble() %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(color = QC_score)) +
  geom_density_2d(size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes detected",
       x = "counts sum")

colData(sce)[colData(sce)$to_QC, ] %>%
  as_tibble() %>%
  ggplot(aes(x = QC_score)) +
  geom_histogram() +
  geom_vline(aes(xintercept = quantile(QC_score, 0.8))) +
  theme_bw()

# From the male in vivo data we extend the QC classification to all cells
colData(sce)$QC_good <- QC_predict(
  fit = QC_fit(
    sce = sce[, colData(sce)$QC_score == 0 | colData(sce)$QC_score >= 1]
  ),
  sce = sce) != 1

save(sce, file = "results/sce_QC_all.Rdata")
load(file = "results/sce_QC_all.Rdata")

# file name for Joana

data_dir <- "data/salmon_output/"
cell_id <- tibble(
    file_name = list.files(
      path = data_dir,
      pattern = ".*quant\\.sf",
      recursive = T
    ),
    id = list.files(
      path = data_dir,
      pattern = ".*quant\\.sf",
      recursive = T
    ) %>%
    stringr::str_replace(
      pattern = ".*[_/]P(\\d+)_(\\w{0,1}\\d+)[_/].*$",
      replacement = "P\\1_\\2"
    ) %>%
    stringr::str_replace(
      pattern = "output/(.).*_(\\d+)_S(\\d+).*$",
      replacement = "P999\\1_\\2\\3"
    )
  )
colData(sce) <- colData(sce) %>% 
  as_tibble() %>% 
  left_join(cell_id, by = "id") %>%
  filter(!(id %>% duplicated())) %>% 
  DataFrame()
colData(sce) %>% 
  as_tibble() %>% 
  write_csv(path = "results/sce_QC_cellData.csv")
# end file name for Joana

colData(sce)[colData(sce)$to_QC, ] %>% 
  as_tibble() %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(color = QC_good)) +
  geom_density_2d(size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  annotation_logticks() +
  labs(y = "genes detected",
       x = "counts sum")

sce$to_QC[sce$day %in% "D1401"] %>% table()

colData(sce) %>% 
  as_tibble() %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(color = QC_good)) +
  geom_density_2d(size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~experiment) +
  annotation_logticks() +
  labs(y = "genes detected",
       x = "counts sum")

colData(sce) %>% 
  as_tibble() %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(color = QC_good)) +
  geom_density_2d(size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~experiment) +
  annotation_logticks() +
  labs(y = "genes detected",
       x = "counts sum")

colData(sce) %>% 
  as_tibble() %>% 
  ggplot(aes(x = detected)) +
  geom_histogram() +
  facet_wrap(~experiment) +
  theme_bw() +
  labs(x = "genes detected")
  

bcrank <- DropletUtils::barcodeRanks(assays(sce)$counts_raw[rowData(sce)$detected, ])
tibble(
  rank = rank(-colSums(assays(sce)$counts_raw[rowData(sce)$detected, ] > 0)),
  expressed = colSums(assays(sce)$counts_raw[rowData(sce)$detected, ] > 0),
  rank_cutoff = max(bcrank$rank[bcrank$total > metadata(bcrank)[["inflection"]]]),
  cell_number = colData(sce)$cell_number,
  QC_good = colData(sce)$QC_good,
  experiment = sce$experiment,
  ) %>%
  distinct() %>%
  ggplot(aes(x = rank, y = expressed)) +
  geom_vline(aes(xintercept = rank_cutoff), col = "red") +
  geom_point(aes(color = QC_good)) +
  facet_wrap(~experiment) +
  geom_line() +
  theme_bw() +
  labs(y = "expressed genes",
       x = "rank")

sce <- sce[, colData(sce)$QC_good]

# we filter the genes
x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- tibble(
  log_mean = x,
  mean = 10 ^ x,
  detection_rate = 1 - dpois(0, lambda = 10 ^ x)
)

rowData(sce)$gene_mean <- rowMeans(assays(sce)$counts_raw[, colData(sce)$to_QC])
rowData(sce)$gene_var <- apply(assays(sce)$counts_raw[, colData(sce)$to_QC], 1, var)
rowData(sce)$detection_rate <- rowMeans(assays(sce)$counts_raw[, colData(sce)$to_QC] > 0)
rowData(sce)$detected <- rowSums(assays(sce)$counts_raw[, colData(sce)$to_QC]) > 0
rowData(sce) %>% 
  as_tibble() %>% 
  filter(detected > 0) %>% 
  ggplot(aes(x = gene_mean, y = detection_rate, color = gene_var/gene_mean >= 1)) +
  geom_point(alpha = 0.3) +
  geom_density_2d(size = 0.3) +
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
  geom_abline(intercept = 0,
              slope = 1,
              color = 'red') +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "genes var",
       x = "gene mean")

rowData(sce)$QC_good <- rowData(sce)$gene_mean > 0 &
  (rowData(sce)$gene_var / rowData(sce)$gene_mean) >= 1
table(rowData(sce)$QC_good)
sce <- sce[rowData(sce)$QC_good, ]

save(sce, file = "results/sce_QC.Rdata")
load(file = "results/sce_QC.Rdata")

