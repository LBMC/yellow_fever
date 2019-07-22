rm(list=ls())
setwd("~/projects/mold/yellow_fever/")
devtools::load_all("pkg/", reset = T)
library(tidyr)
library(tidyverse)
Sys.setenv("DISPLAY"=":0")

system("mkdir -p results/prot_classification")

data <- read_delim("data/190719_YF_flow.csv", delim=";") %>%
  mutate(Day = as.factor(Day),
         Donor = as.factor(Donor),
         HLA = as.factor(HLA)
  )

markers = c("CCR7", "GZMB", "TCF7", "LEF1", "CD45RA")
prots <- data %>% colnames() %>% .[5:19]

annotate_marker <- function(data, marker, nfirst = 500) {
  new_col <- paste0(marker, "p")
  qmarker <- sym(marker)
  data <- data %>% arrange(!! qmarker) %>% as.data.frame()
  data[[new_col]] <- NA
  data[[new_col]][1:nfirst] <- T
  data[[new_col]][(nrow(data) - nfirst + 1):nrow(data)] <- F
  data %>% as_tibble()
}

for(marker in markers) {
  data <- data %>%
    annotate_marker(marker)
}

data %>%
  gather(prots, key = prot, value = level) %>%
  mutate(prot = as.factor(prot)) %>%
  group_by(prot) %>%
  mutate(level = scale(level)) %>%
  ungroup() %>%
  ggplot(aes(x = level)) +
  theme_bw() +
  geom_density() +
  facet_wrap(~prot, scale = "free")

data %>%
  gather(prots, key = prot, value = level) %>%
  mutate(prot = as.factor(prot)) %>%
  group_by(prot) %>%
  mutate(level = level) %>%
  ungroup() %>%
  ggplot(aes(x = level+abs(min(level)))) +
  theme_bw() +
  geom_density() +
  facet_wrap(~prot, scale = "free")

data <- data %>%
  gather(prots, key = prot, value = raw) %>%
  mutate(prot = as.factor(prot)) %>%
  group_by(prot) %>%
  mutate(scaled = scale(raw)) %>%
  ungroup() %>%
  pivot_wider(names_from = prot, values_from = c(raw, scaled))

# PLS training
for ( marker in markers) {
  print(marker)
  pmarker <- paste0(marker, "p")
  training_file <- paste0("results/prot_classification/",
                          pmarker,
                          "_training.Rdata")
  if (file.exists(training_file)){
    load(training_file)
  } else {
    qpmarker <- sym(pmarker)
    by_var <- list(by = (data %>%
                        filter(!is.na(!! qpmarker)) %>%
                        pull(!! qpmarker) %>%
                        as.factor() %>%
                        as.numeric() %>%
                        as.vector())-1)
    training <- logistic_spls_stab_training(
      by = by_var,
      data = data %>%
        filter(!is.na(!! qpmarker)) %>%
               select(paste0("scaled_", prots)) %>%
               as.matrix(),
      ncores = 8,
      file = paste0("results/prot_classification/", pmarker, "_training"),
      force = F
    )
    save(training, file = training_file)
    system(paste0("~/scripts/sms.sh \"PLS done", marker,"\""))
  }
}

# PLS classification
for ( marker in markers) {
  print(marker)
  pmarker <- paste0(marker, "p")
  training_file <- paste0("results/prot_classification/",
                          pmarker,
                          "_training.Rdata")
  classification_file <- paste0("results/prot_classification/",
                          pmarker,
                          "_classification.Rdata")
  if (file.exists(classification_file)){
    load(classification_file)
  } else {
    load(training_file)
    qpmarker <- sym(pmarker)
    classification <- logistic_spls_stab_classification(
      fit = training,
      data = data %>%
        filter(is.na(!! qpmarker)) %>%
        select(paste0("scaled_", prots)),
      data_train = data %>%
        filter(!is.na(!! qpmarker)) %>%
        select(paste0("scaled_", prots)),
      ncores = 8,
      file = paste0("results/prot_classification/",
                    pmarker, "_classification"),
      force = training$classification$fit_spls$fit$selected
    )
    save(classification, file = classification_file)
    system(paste0("~/scripts/sms.sh \"PLS done", marker,"\""))
  }
}

for ( marker in markers) {
  print(marker)
  marker <- "CCR7"
  pmarker <- paste0(marker, "p")
  classification_file <- paste0("results/prot_classification/",
                          pmarker,
                          "_classification.Rdata")
  load(classification_file)
  print(classification$fit_spls$fit$selected)
}

for ( marker in markers) {
  print(marker)
  marker <- "CCR7"
  pmarker <- paste0(marker, "p")
  training_file <- paste0("results/prot_classification/",
                          pmarker,
                          "_training.Rdata")
  load(training_file)
  print(plsgenomics::stability.selection(training$fit)$selected.predictor)
}

for ( marker in markers) {
  print(marker)
  pmarker <- paste0("p", marker)
  pmarkerg <- paste0("p", marker, "g")
  markerp <- paste0(marker, "p")
  classification_file <- paste0("results/prot_classification/",
                          markerp,
                          "_classification.Rdata")
  if (file.exists(classification_file)){
    rm(classification)
    load(classification_file)
    data[[pmarker]] <- NA
    data[[pmarkerg]] <- NA
    data[[pmarkerg]][!is.na(data[[markerp]])] <- classification$model$hatY
    data[[pmarker]][!is.na(data[[markerp]])] <- classification$model$proba
    data[[pmarkerg]][is.na(data[[markerp]])] <- classification$model$hatYtest
    data[[pmarker]][is.na(data[[markerp]])] <- classification$model$proba.test
  }
}

save(data, file = "results/prot_classification/data.Rdata")
load("results/prot_classification/data.Rdata")
data %>%
  ggplot(aes(x = UMAPX, y = UMAPY, color = pCCR7)) +
  theme_bw() +
  geom_point()
data %>%
  ggplot(aes(x = raw_CCR7, y = scaled_CCR7, color = pCCR7)) +
  theme_bw() +
  geom_point()
data %>%
  ggplot(aes(x = UMAPX, y = UMAPY, color = pGZMB)) +
  theme_bw() +
  geom_point()
data %>%
  ggplot(aes(x = UMAPX, y = UMAPY, color = pTCF7)) +
  theme_bw() +
  geom_point()
data %>%
  ggplot(aes(x = UMAPX, y = UMAPY, color = pLEF1)) +
  theme_bw() +
  geom_point()
data %>%
  ggplot(aes(x = UMAPX, y = UMAPY, color = pCD45RA)) +
  theme_bw() +
  geom_point()

write_delim(data,
            path = "results/prot_classification/prot_classification.csv",
            delim = ";"
)

library(ComplexHeatmap)

pdf(file = "results/prot_classification/CCR7.pdf", width=16, height=8)
col <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256))
col_fun <- colorRamp2(seq(from = 0.2, to = 0.8, length = 256), col)
ha = HeatmapAnnotation(Donor = data %>% arrange(pCCR7) %>% pull(Donor),
                       D15 = data %>% arrange(pCCR7) %>% mutate(D15 = ifelse(Day == 'd15', "d15", "d90")) %>% pull(D15),
                       D90 = data %>% arrange(pCCR7) %>% mutate(D90 = ifelse(Day == 'd90', "d90", "d15")) %>% pull(D90),
                       HLA = data %>% arrange(pCCR7) %>% pull(HLA),
                       pCCR7 = data %>% arrange(pCCR7) %>% pull(pCCR7),
                       pCCR7_group = data %>% arrange(pCCR7) %>%
                         mutate(pCCR7g = ifelse(pCCR7g, "CCR7p", "CCR7m")) %>%
                         pull(pCCR7g) %>% as.factor(),
                       col = list(pCCR7 = col_fun,
                                  HLA = c("B7" = "#9AC5DE", "A2" = "#C93737"),
                                  pCCR7_group = c("CCR7p" = "#C5C659", "CCR7m" = "#2D2950"),
                                  D15 = c("d15" = "black", "d90" = "white"),
                                  D90 = c("d90" = "black", "d15" = "white"),
                                  Donor = c("donor1" = "pink", "donor2" = "orange", "donor3" = "black")
                       ),
                       show_annotation_name = T
)
ascb <- function(x){
  scale(apply(x, 2, FUN=function(x){sqrt(x + abs(min(x)) + 3/8)}))
}
data2 <- data %>%
  arrange(pCCR7) %>%
  select(starts_with("raw_")) %>%
  rename_all(funs(str_replace(., "raw_", ""))) %>%
  t() %>%
  ascb()
col_fun_2 <- colorRamp2(seq(from = -2, to = 2, length = 256), col)
Heatmap(data2,
        top_annotation = ha,
        cluster_columns = F,
        col = col_fun_2,
        color_space = "LAB",
)
dev.off()


pdf(file = "results/prot_classification/GZMB.pdf", width=16, height=8)
col <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256))
col_fun <- colorRamp2(seq(from = 0.2, to = 0.8, length = 256), col)
ha = HeatmapAnnotation(Donor = data %>% arrange(pGZMB) %>% pull(Donor),
                       D15 = data %>% arrange(pGZMB) %>% mutate(D15 = ifelse(Day == 'd15', "d15", "d90")) %>% pull(D15),
                       D90 = data %>% arrange(pGZMB) %>% mutate(D90 = ifelse(Day == 'd90', "d90", "d15")) %>% pull(D90),
                       HLA = data %>% arrange(pGZMB) %>% pull(HLA),
                       pGZMB = data %>% arrange(pGZMB) %>% pull(pGZMB),
                       pGZMB_group = data %>% arrange(pGZMB) %>%
                         mutate(pGZMBg = ifelse(pGZMBg, "GZMBp", "GZMBm")) %>%
                         pull(pGZMBg) %>% as.factor(),
                       col = list(pGZMB = col_fun,
                                  HLA = c("B7" = "#9AC5DE", "A2" = "#C93737"),
                                  pGZMB_group = c("GZMBp" = "#C5C659", "GZMBm" = "#2D2950"),
                                  D15 = c("d15" = "black", "d90" = "white"),
                                  D90 = c("d90" = "black", "d15" = "white"),
                                  Donor = c("donor1" = "pink", "donor2" = "orange", "donor3" = "black")
                       ),
                       show_annotation_name = T
)
ascb <- function(x){
  scale(apply(x, 2, FUN=function(x){sqrt(x + abs(min(x)) + 3/8)}))
}
data2 <- data %>%
  arrange(pGZMB) %>%
  select(starts_with("raw_")) %>%
  rename_all(funs(str_replace(., "raw_", ""))) %>%
  t() %>%
  ascb()
col_fun_2 <- colorRamp2(seq(from = -2, to = 2, length = 256), col)
Heatmap(data2,
        top_annotation = ha,
        cluster_columns = F,
        col = col_fun_2,
        color_space = "LAB",
)
dev.off()


pdf(file = "results/prot_classification/TCF7.pdf", width=16, height=8)
col <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256))
col_fun <- colorRamp2(seq(from = 0.2, to = 0.8, length = 256), col)
ha = HeatmapAnnotation(Donor = data %>% arrange(pTCF7) %>% pull(Donor),
                       D15 = data %>% arrange(pTCF7) %>% mutate(D15 = ifelse(Day == 'd15', "d15", "d90")) %>% pull(D15),
                       D90 = data %>% arrange(pTCF7) %>% mutate(D90 = ifelse(Day == 'd90', "d90", "d15")) %>% pull(D90),
                       HLA = data %>% arrange(pTCF7) %>% pull(HLA),
                       pTCF7 = data %>% arrange(pTCF7) %>% pull(pTCF7),
                       pTCF7_group = data %>% arrange(pTCF7) %>%
                         mutate(pTCF7g = ifelse(pTCF7g, "TCF7p", "TCF7m")) %>%
                         pull(pTCF7g) %>% as.factor(),
                       col = list(pTCF7 = col_fun,
                                  HLA = c("B7" = "#9AC5DE", "A2" = "#C93737"),
                                  pTCF7_group = c("TCF7p" = "#C5C659", "TCF7m" = "#2D2950"),
                                  D15 = c("d15" = "black", "d90" = "white"),
                                  D90 = c("d90" = "black", "d15" = "white"),
                                  Donor = c("donor1" = "pink", "donor2" = "orange", "donor3" = "black")
                       ),
                       show_annotation_name = T
)
ascb <- function(x){
  scale(apply(x, 2, FUN=function(x){sqrt(x + abs(min(x)) + 3/8)}))
}
data2 <- data %>%
  arrange(pTCF7) %>%
  select(starts_with("raw_")) %>%
  rename_all(funs(str_replace(., "raw_", ""))) %>%
  t() %>%
  ascb()
col_fun_2 <- colorRamp2(seq(from = -2, to = 2, length = 256), col)
Heatmap(data2,
        top_annotation = ha,
        cluster_columns = F,
        col = col_fun_2,
        color_space = "LAB",
)
dev.off()


pdf(file = "results/prot_classification/LEF1.pdf")
col <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256))
col_fun <- colorRamp2(seq(from = 0.2, to = 0.8, length = 256), col)
ha = HeatmapAnnotation(Donor = data %>% arrange(pLEF1) %>% pull(Donor),
                       D15 = data %>% arrange(pLEF1) %>% mutate(D15 = ifelse(Day == 'd15', "d15", "d90")) %>% pull(D15),
                       D90 = data %>% arrange(pLEF1) %>% mutate(D90 = ifelse(Day == 'd90', "d90", "d15")) %>% pull(D90),
                       HLA = data %>% arrange(pLEF1) %>% pull(HLA),
                       pLEF1 = data %>% arrange(pLEF1) %>% pull(pLEF1),
                       pLEF1_group = data %>% arrange(pLEF1) %>%
                         mutate(pLEF1g = ifelse(pLEF1g, "LEF1p", "LEF1m")) %>%
                         pull(pLEF1g) %>% as.factor(),
                       col = list(pLEF1 = col_fun,
                                  HLA = c("B7" = "#9AC5DE", "A2" = "#C93737"),
                                  pLEF1_group = c("LEF1p" = "#C5C659", "LEF1m" = "#2D2950"),
                                  D15 = c("d15" = "black", "d90" = "white"),
                                  D90 = c("d90" = "black", "d15" = "white"),
                                  Donor = c("donor1" = "pink", "donor2" = "orange", "donor3" = "black")
                       ),
                       show_annotation_name = T
)
ascb <- function(x){
  scale(apply(x, 2, FUN=function(x){sqrt(x + abs(min(x)) + 3/8)}))
}
data2 <- data %>%
  arrange(pLEF1) %>%
  select(starts_with("raw_")) %>%
  rename_all(funs(str_replace(., "raw_", ""))) %>%
  t() %>%
  ascb()
col_fun_2 <- colorRamp2(seq(from = -2, to = 2, length = 256), col)
Heatmap(data2,
        top_annotation = ha,
        cluster_columns = F,
        col = col_fun_2,
        color_space = "LAB",
)
dev.off()

pdf(file = "results/prot_classification/CD45RA.pdf", width=16, height=8)
col <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256))
col_fun <- colorRamp2(seq(from = 0.2, to = 0.8, length = 256), col)
ha = HeatmapAnnotation(Donor = data %>% arrange(pCD45RA) %>% pull(Donor),
                       D15 = data %>% arrange(pCD45RA) %>% mutate(D15 = ifelse(Day == 'd15', "d15", "d90")) %>% pull(D15),
                       D90 = data %>% arrange(pCD45RA) %>% mutate(D90 = ifelse(Day == 'd90', "d90", "d15")) %>% pull(D90),
                       HLA = data %>% arrange(pCD45RA) %>% pull(HLA),
                       pCD45RA = data %>% arrange(pCD45RA) %>% pull(pCD45RA),
                       pCD45RA_group = data %>% arrange(pCD45RA) %>%
                         mutate(pCD45RAg = ifelse(pCD45RAg, "CD45RAp", "CD45RAm")) %>%
                         pull(pCD45RAg) %>% as.factor(),
                       col = list(pCD45RA = col_fun,
                                  HLA = c("B7" = "#9AC5DE", "A2" = "#C93737"),
                                  pCD45RA_group = c("CD45RAp" = "#C5C659", "CD45RAm" = "#2D2950"),
                                  D15 = c("d15" = "black", "d90" = "white"),
                                  D90 = c("d90" = "black", "d15" = "white"),
                                  Donor = c("donor1" = "pink", "donor2" = "orange", "donor3" = "black")
                       ),
                       show_annotation_name = T
)
ascb <- function(x){
  scale(apply(x, 2, FUN=function(x){sqrt(x + abs(min(x)) + 3/8)}))
}
data2 <- data %>%
  arrange(pCD45RA) %>%
  select(starts_with("raw_")) %>%
  rename_all(funs(str_replace(., "raw_", ""))) %>%
  t() %>%
  ascb()
col_fun_2 <- colorRamp2(seq(from = -2, to = 2, length = 256), col)
Heatmap(data2,
        top_annotation = ha,
        cluster_columns = F,
        col = col_fun_2,
        color_space = "LAB",
)
dev.off()
