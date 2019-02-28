rm(list = ls())
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
system("mkdir -p results/cycling/")

# we try to refine the regev cell-cycle genes list
regev_genes <- c("MCM5", "HMGB2", "PCNA", "CDK1", "TYMS", "NUSAP1", "FEN1",
  "UBE2C", "MCM2", "BIRC5", "MCM4", "TPX2", "RRM1", "TOP2A", "UNG", "NDC80",
  "GINS2", "CKS2", "MCM6", "NUF2", "CDCA7", "CKS1B", "DTL", "MKI67", "PRIM1",
  "TMPO", "UHRF1", "CENPF", "MLF1IP", "TACC3", "HELLS", "FAM64A", "RFC2",
  "SMC4", "RPA2", "CCNB2", "NASP", "CKAP2L", "RAD51AP1", "CKAP2", "GMNN",
  "AURKB", "WDR76", "BUB1", "SLBP", "KIF11", "CCNE2", "ANP32E", "UBR7",
  "TUBB4B", "POLD3", "GTSE1", "MSH2", "KIF20B", "ATAD2", "HJURP", "RAD51",
  "HJURP", "RRM2", "CDCA3", "CDC45", "HN1", "CDC6", "CDC20", "EXO1", "TTK",
  "TIPIN", "CDC25C", "DSCC1", "KIF2C", "BLM", "RANGAP1", "CASP8AP2", "NCAPD2",
  "USP1", "DLGAP5", "CLSPN", "CDCA2", "POLA1", "CDCA8", "CHAF1B", "ECT2",
  "BRIP1", "KIF23", "E2F8", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR",
  "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")
save(regev_genes, file="results/cycling/regev_genes.RData")

cycling_score <- rep(0, scd$getncells)
cycling <- rep("SLC", scd$getncells)
pcycling <- rep(0, scd$getncells)
genes_cycling <- regev_genes
require("mixtools")
for (sex in c("M", "F")) {
  time_range <- c("D15", "D136", "D593")
  if (sex %in% "F") {
    time_range <- c("D15", "D90")
  }
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% time_range &
    scd$getfeature("sex") %in% sex
  cycling_score[b_cells] <- pca_loading(
    scd = scd$select(
      b_cells = b_cells,
      genes = genes_cycling
    ),
    cells = TRUE
  )[, 1]
  cycling_score[b_cells] <- log(abs(
    cycling_score[b_cells] - max(cycling_score[b_cells])
  ) + 1)
  print(summary(cycling_score[b_cells & scd$getfeature("day") %in% "D136"]))
  model <- normalmixEM(
    cycling_score[b_cells],
    lambda = .5,
    mu = c(0, 2), sigma = c(1,2)
  )
  cycling[b_cells] <- ifelse(model$posterior[,2]>0.5, "cycling", "SLC")
  pcycling[b_cells] <- model$posterior[,2]
}
scd$setfeature("cycling", cycling)
scd$setfeature("pcycling", as.vector(pcycling))
scd$setfeature("cycling_score", as.vector(cycling_score))

save(scd, file = "results/cycling/cells_counts_QC_cycling.Rdata")
infos <- scd$getfeatures
write.csv(
  infos,
  file = paste0("results/cycling/cell_type_infos.csv")
)

b_cells <- scd$getfeature("day") %in% c("D15", "D136", "D593", "D90") &
  scd$getfeature("QC_good") %in% T &
  scd$getfeature("cell_number") %in% 1
infos_M <- scd$select(b_cells = b_cells)$getfeatures
infos_M <- apply(infos_M, 2, as.vector)
counts_M <- scd$select(b_cells = b_cells)$getcounts
infos_M <- t(infos_M)
counts_M <- t(counts_M)
dim(infos_M)
dim(counts_M)
normalized_counts <- rbind(infos_M, counts_M)
dim(normalized_counts)
write.csv(
  normalized_counts,
  file = paste0("results/cell_type/cell_type_cells_counts.csv")
)

# plots

for (sex in c("M", "F")) {
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("sex") %in% sex
  time_range <- c("D15", "D136", "D593")
  if (sex == "F") {
    time_range <- c("D15", "D90")
  }
  g <- ggplot(
    data = scd$select(
      b_cells = b_cells & scd$getfeature("day") %in% time_range
      )$getfeatures,
    aes(x = cycling_score,
        y = pDEA_cell_type,
        color = day)
    ) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = day_palette(levels(as.factor(as.vector(
      scd$select(
        b_cells = b_cells & scd$getfeature("day") %in% time_range
      )$getfeature("day")
    ))))) +
    labs(title = paste0("cycling score vs Memory probability ", sex),
      x = "cycling score", y = "Memory probability", color = "day")
  print(g)
  ggsave(
    file = paste0("results/cycling/cells_counts_QC_DEA_cell_type_vs_cycling_score_", sex, ".pdf")
  )
}

for (sex in c("M", "F")) {
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("sex") %in% sex
  for (time_range in list(c("D15"), c("D15", "D136", "D593"))) {
    if (sex == "F" & length(time_range) == 3) {
      time_range <- c("D15", "D90")
    }
    system(ifelse(time_range[1] == "D15",
      paste0("rm results/tmp/pca_cells_counts_QC_D15_", sex, ".Rdata"),
      paste0("rm results/tmp/pca_cells_counts_QC_all_day_", sex, ".Rdata")
    ))
    for (score_type in c("cycling", "pcycling")) {
      scRNAtools::pca_plot(
        scd$select(b_cells = b_cells &
          scd$getfeature("day") %in% time_range,
          genes = genes_cycling
        ),
        color = score_type,
        color_name = ifelse(score_type == "cycling",
          "cycling", "cycling_score"),
        tmp_file = ifelse(time_range[1] == "D15",
          paste0("results/tmp/pca_cells_counts_QC_D15_", sex, ".Rdata"),
          paste0("results/tmp/pca_cells_counts_QC_all_day_", sex, ".Rdata")
        ),
        main = score_type
      )
      ggsave(file = ifelse(time_range[1] == "D15",
          paste0("results/cycling/pca_cells_counts_QC_", score_type, "_D15_", sex, ".pdf"),
          paste0("results/cycling/pca_cells_counts_QC_", score_type, "_all_day_", sex, ".pdf")
        )
      )
    }
  }
}

# we find all the genes that covariate with the regev genes
load("results/cycling/regev_genes.RData")
regev_genes_cov <- list()
for(gene in regev_genes[-c(38,50,60,66,68,85,98)]){
  print(gene)
  regev_genes_cov[[gene]] <- gene_cov(data_genes_norm[r_select,],
                                      gene)
}
save(regev_genes_cov, file=paste0(outdir, "regev_genes_cov.RData"))

# from this list we remove which genes are expressed at D136 and D908
load(paste0(outdir, "regev_genes_cov.RData"))
regev_genes_cov <- unique(c(regev_genes, unlist(regev_genes_cov)))
r_select <- scell_QC_M(data_infos, T) & data_infos$day %in% c("D136", "D908") & !(data_infos$batch %in% c(1,6:8))
c_select <- which(colnames(data_genes) %in% unique(regev_genes_cov))
regev_genes_cov_expressed <- colnames(del_zero(data_genes[r_select, c_select],
                                               n_cell = 0.05,
                                               min_count = 10,
                                               percentage = TRUE,
                                               conditions = data_infos$day[r_select]))
cell_cycle_genes <- regev_genes_cov[!(regev_genes_cov %in% regev_genes_cov_expressed)]
length(regev_genes_cov)
length(regev_genes_cov_expressed)
length(cell_cycle_genes)
summary(cell_cycle %in% cell_cycle_genes)

cell_cycle[cell_cycle %in% cell_cycle_genes]
cell_cycle[!(cell_cycle %in% cell_cycle_genes)]

c_select <- which(colnames(data_genes) %in% "MKI67")
data_genes[r_select, c_select]

save(cell_cycle_genes, file=paste0(outdir, "cell_cycle_genes.RData"))
load(paste0(outdir, "cell_cycle_genes.RData"))

library("tidyverse")
library("readxl")

load(file = "results/cycling/cells_counts_QC_cycling.Rdata")
infos <- read_xlsx("data/Column_Headers.xlsx")
infos <- infos %>% colnames()
infos[1] <- infos[1] %>% tolower()
infos[3] <- infos[3] %>% tolower()

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type"))

infos <- cbind(
      scd$select(b_cells = b_cells)$getfeatures[, ( infos %>% .[1:10]  )] ,
      scd$select(b_cells = b_cells)$getcounts[, ( infos %>% .[11:60] )],
      scd$select(b_cells = b_cells)$getfeatures[, ( infos %>% .[63:73]  )]
      )
write.csv(infos, file = "results/cycling/infos_cell_type_formated.csv")
