library(tidyverse)
library(SingleCellExperiment)
library(SummarizedExperiment)

load(file = "results/sce_QC.Rdata", verbose = T)
colnames(sce) <- sce$id
sce$batch[sce$experiment %in% "P9997"] <- 40
sce$batch[sce$experiment %in% "P9998"] <- 41
sce$day[sce$experiment %in% c("P9997", "P9998")] <- "D1401"
sce$donor[sce$experiment %in% c("P9997", "P9998")] <- "YVF2003"
sce$antigen[sce$experiment %in% c("P9997", "P9998")] <- "A2"

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


norm_counts <- list()

sce$male_invivo <- sce$sex %in% "M" &
  sce$cell_number <= 1 &
  sce$day %in% c("D15", "D136", "D593") &
  sce$sequencing %in% "paired" &
  !(sce$batch %in% c(6:8)) &
  !(sce$id %in% c("P1299_1797", "P1299_1896"))

norm_counts[["male_invivo"]] <- SCnorm::SCnorm(
  Data = assays(sce[, sce$male_invivo])$counts_raw,
  Conditions = rep(c(1), ncol(sce[, sce$male_invivo])),
  NCores = 1
)

save(norm_counts, file = "results/sce_norm_list.Rdata")

sce$female_invivo <- sce$sex %in% "F" &
  sce$cell_number <= 1 &
  sce$day %in% c("D15", "D90")

norm_counts[["female_invivo"]] <- SCnorm::SCnorm(
  Data = assays(sce[, sce$male_invivo])$counts_raw,
  Conditions = rep(c(1), ncol(sce[, sce$male_invivo])),
  NCores = 1
)

save(norm_counts, file = "results/sce_norm_list.Rdata")

sce$YVF2003_D1401 <- sce$day %in% c("D1401")

norm_counts[["YVF2003_D1401"]] <- SCnorm::SCnorm(
  Data = assays(sce[, sce$male_invivo])$counts_raw,
  Conditions = rep(c(1), ncol(sce[, sce$male_invivo])),
  NCores = 1
)

save(norm_counts, file = "results/sce_norm_list.Rdata")

sce$male_in_vitro <- sce$day %in% c("InVitro") &
  sce$cell_number <= 1

norm_counts[["male_in_vitro"]] <- SCnorm::SCnorm(
  Data = assays(sce[, sce$male_in_vitro])$counts_raw,
  Conditions = rep(c(1), ncol(sce[, sce$male_in_vitro])),
  NCores = 1
)

save(norm_counts, file = "results/sce_norm_list.Rdata")

sce$male_in_vitro_restim <- sce$day %in% c("In_Vitro_Restim") &
  sce$cell_number <= 1 &
  sce$sex %in% "M"

norm_counts[["male_in_vitro"]] <- SCnorm::SCnorm(
  Data = assays(sce[, sce$male_in_vitro_restim])$counts_raw,
  Conditions = rep(c(1), ncol(sce[, sce$male_in_vitro_restim])),
  NCores = 1
)

save(norm_counts, file = "results/sce_norm_list.Rdata")

sce$male_in_vitro_restim <- sce$day %in% c("In_Vitro_Restim") &
  sce$cell_number <= 1 &
  sce$sex %in% "F"

norm_counts[["female_in_vitro_restim"]] <- SCnorm::SCnorm(
  Data = assays(sce[, sce$female_in_vitro_restim])$counts_raw,
  Conditions = rep(c(1), ncol(sce[, sce$female_in_vitro_restim])),
  NCores = 1
)

save(norm_counts, file = "results/sce_norm_list.Rdata")

assays(sce)$counts <- norm_counts$NormalizedData %>%
  Matrix::Matrix(sparse = T)

save(sce, file = "results/sce_scnorm.Rdata")
load(file = "results/sce_scnorm.Rdata", verbose = T)
