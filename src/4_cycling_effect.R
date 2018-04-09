setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
infos_M <- scd$getfeatures
load("results/cell_type/CB_counts_QC_DEA_cell_type_F.Rdata")
b_cells_F <- scd$getfeature("sex") %in% "F"
infos_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getfeatures
scd <- scdata$new(
  infos = infos_M,
  counts = scd$getcounts
)

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
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"
genes_cycling <- expressed(
  scd = scd$select(
    b_cells = b_cells &
      scd$getfeature("day") %in% c("D136", "D593"),
    genes = regev_genes[-c(38,50,60,66,68,85,98)]
  ),
  zi_threshold = 0.90,
)
table(genes_cycling %in% regev_genes)
genes_cycling <- regev_genes
require("mixtools")
for (sex in c("M", "F")) {
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) & 
    scd$getfeature("day") %in% ifelse(sex == "M",
     c("D15", "D136", "D593"), c("D15", "D90")) &
    scd$getfeature("sex") %in% sex
  cycling_score[b_cells] <- pca_loading(
    scd = scd$select(
      b_cells = b_cells,
      genes = genes_cycling
    ),
    cells = TRUE
  )[, 1]
  model <- normalmixEM(
    log(abs(cycling_score[b_cells])+1),
    lambda = .5,
    mu = c(0, 2), sigma = c(1,2)
  )
  cycling[b_cells] <- ifelse(model$posterior[,2]>0.5, "cycling", "SLC")
  pcycling[b_cells] <- model$posterior[,2]
}
scd$setfeature("cycling", cycling)
scd$setfeature("pcycling", as.vector(pcycling))
scd$setfeature("cycling_score", as.vector(cycling_score))

# cycling phase with cyclone
cycling_phase <- rep(NA, scd$getncells)
cycling_G1_score <- rep(NA, scd$getncells)
cycling_G2M_score <- rep(NA, scd$getncells)
cycling_S_score <- rep(NA, scd$getncells)
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
for (sex in c("M", "F")) {
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) & 
    scd$getfeature("day") %in% "D15" &
    scd$getfeature("sex") %in% sex
  hs_pairs <- list()
  for (phase in names(hs.pairs)) {
    hs_pairs[[phase]] <- hs.pairs[[phase]]
    for (pairs in names(hs.pairs[[phase]])) {
      conversion <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "hgnc_symbol"),
        values = hs.pairs[[phase]][[pairs]],
        mart = mart
      )
      rosette <- list()
      for (ensembl_gene_id in conversion$ensembl_gene_id) {
        rosette[[ensembl_gene_id]] <- conversion$hgnc_symbol[
          conversion$ensembl_gene_id %in% ensembl_gene_id
        ]
      } 
      ensembl_gene_ids <- hs_pairs[[phase]][[pairs]] 
      hs_pairs[[phase]][[pairs]] <- unlist(rosette[ ensembl_gene_ids ])
    }
  }
  cycling_infos <- scran::cyclone(
    x = t(scd$select(
      b_cells = b_cells)$getcounts),
    pairs = hs_pairs,
    verbose = T
  )
  cycling_phase[b_cells] <- cycling_infos$phases
  cycling_G1_score[b_cells] <- cycling_infos$normalized.scores$G1
  cycling_G2M_score[b_cells] <- cycling_infos$normalized.scores$G2M
  cycling_S_score[b_cells] <- cycling_infos$normalized.scores$S
}
scd$setfeature("cycling_phase", cycling_phase)
scd$setfeature("cycling_G1_score", cycling_G1_score)
scd$setfeature("cycling_G2M_score", cycling_G2M_score)
scd$setfeature("cycling_S_score", cycling_S_score)

save(scd, file = "results/cycling/CB_counts_QC_cycling.Rdata")
infos <- scd$getfeatures
write.csv(
  infos,
  file = paste0("results/cycling/cell_type_infos.csv")
)


load(file = "results/cycling/CB_counts_QC_cycling.Rdata")

# plots
for (sex in c("M", "F")) {
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("sex") %in% sex
  for (time_range in list(c("D15"), c("D15", "D136", "D593"))) {
    if (sex == "F" & length(time_range) == 3) {
      time_range <- c("D15", "D90")
    }
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
          paste0("results/tmp/pca_CB_counts_QC_D15_", sex, ".Rdata"),
          paste0("results/tmp/pca_CB_counts_QC_all_day_", sex, ".Rdata")
        ),
        main = score_type
      )
      ggsave(file = ifelse(time_range[1] == "D15",
          paste0("results/tmp/pca_CB_counts_QC_", score_type, "_D15_", sex, ".pdf"),
          paste0("results/tmp/pca_CB_counts_QC_", score_type, "_all_day_", sex, ".pdf")
        )
      )
    }
  }
}

for (sex in c("M", "F")) {
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("sex") %in% sex
  for (time_range in list(c("D15"), c("D15", "D136", "D593"))) {
    if (sex == "F" & length(time_range) == 3) {
      time_range <- c("D15", "D90")
    }
    for (score_type in c("cycling", "pcycling")) {
      g <- ggplot(
        data = scd$select(
          b_cells = b_cells & scd$getfeature("day") %in% time_range
          )$getfeatures,
        aes(x = cycling_score,
            y = pDEA_cell_type,
            color = day)
        ) +
        geom_point() +
        theme_bw()
      print(g)
    }
  }
}

# we find all the genes that covariate with the regev genes
load("results/cycling/regev_genes.RData")
regev_genes_cov <- list()
for(gene in regev_genes[-c(38,50,60,66,68,85,98)]){
  print(gene)
  regev_genes_cov[[gene]] <- pls_on_one_genes(data_genes_norm[r_select,], gene)
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
