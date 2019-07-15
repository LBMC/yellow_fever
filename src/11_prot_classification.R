rm(list=ls())
setwd("~/projects/mold/yellow_fever/")
devtools::load_all("pkg/", reset = T)
library(tidyverse)

system("mkdir -p results/prot_classification")

data <- read_csv2("data/190710_YF_flow.csv") %>%
  mutate(Day = as.factor(Day),
         Donor = as.factor(Donor),
         HLA = as.factor(HLA),
         CD16 = as.numeric(CD16),
         CD45RA = as.numeric(CD45RA),
         CCR7 = as.numeric(CCR7),
         CD57 = as.numeric(CD57),
         CD127 = as.numeric(CD127),
         TCF7 = as.numeric(TCF7),
         CD94 = as.numeric(CD94),
         UMAPX = as.numeric(UMAPX),
         UMAPY = as.numeric(UMAPY)
  )

markers = c("CCR7", "GZMB", "TCF7", "LEF1", "IL7R")
prots <- data %>% colnames() %>% .[5:19]

annotate_marker <- function(data, marker, nfirst = 500) {
  new_col <- paste0(marker, "p")
  data <- data %>% as.data.frame()
  data[[new_col]] <- NA
  data[[new_col]][1:nfirst] <- T
  data[[new_col]][(nrow(data) - nfirst + 1):nrow(data)] <- F
  data %>% as.tibble()
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

data <- data %>%
  gather(prots, key = prot, value = level) %>%
  mutate(prot = as.factor(prot)) %>%
  group_by(prot) %>%
  mutate(level = scale(level)) %>%
  ungroup() %>%
  spread(prot, level)

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
               select(prots) %>%
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
      data = data %>% filter(is.na(!! qpmarker)) %>% select(prots),
      data_train = data %>% filter(!is.na(!! qpmarker)) %>% select(prots),
      ncores = 8,
      file = paste0("results/prot_classification/",
                    pmarker, "_classification"),
      force = training$classification$fit_spls$fit$selected
    )
    save(classification, file = classification_file)
    system(paste0("~/scripts/sms.sh \"PLS done", marker,"\""))
  }
}

