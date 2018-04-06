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
cycling <- rep(0, scd$getncells)
cycling[b_cells = b_cells &
      scd$getfeature("day") %in% "D15"] <- pca_loading(
  scd = scd$select(
    b_cells = b_cells &
      scd$getfeature("day") %in% "D15",
    genes = genes_cycling
  ),
  cells = TRUE
)[, 1]

hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
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

cycling <- scran::cyclone(
  x = t(scd$select(
    b_cells = b_cells & scd$getfeature("day") %in% "D15")$getcounts),
  pairs = hs_pairs,
  verbose = T
)
print(cycling$phases)


b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "F"
cycling[b_cells = b_cells &
      scd$getfeature("day") %in% "D15"] <- pca_loading(
  scd = scd$select(
    b_cells = b_cells &
      scd$getfeature("day") %in% "D15",
    genes = genes_cycling
  ),
  cells = TRUE
)[, 1]

scd$setfeature("cycling", as.vector(cycling))

save(scd, file = "results/cycling/CB_counts_QC_cycling.Rdata")
infos <- scd$getfeatures

write.csv(
  infos,
  file = paste0("results/cycling/cell_type_infos.csv")
)


b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M"

system("rm results/tmp/pca_zi_norm_counts_QC_all_day.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells &
    scd$getfeature("day") %in% c("D15", "D136", "D593")),
  color = "cycling", 
  color_name = "cycling_score",
  shape = "day",
  tmp_file = "results/tmp/pca_zi_norm_counts_QC_all_day.Rdata",
  main = "all day cycling"
)
ggsave(file = 
  "results/cycling/pca_zi_norm_counts_QC_cycling_pca_all_day.pdf"
)

system("rm results/tmp/pca_zi_norm_counts_QC_D15.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells &
    scd$getfeature("day") %in% "D15"),
  color = "cycling", 
  color_name = "cycling_score",
  tmp_file = "results/tmp/pca_zi_norm_counts_QC_D15.Rdata",
  main = "D15 cycling"
)
ggsave(file = 
  "results/cycling/pca_zi_norm_counts_QC_cycling_pca_D15.pdf"
)

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "F"

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells &
    scd$getfeature("day") %in% "D15"),
  color = "cycling", 
  color_name = "cycling_score",
  tmp_file = "results/tmp/pca_zi_norm_counts_QC_D15_F.Rdata",
  main = "D15 cycling"
)
ggsave(file = 
  "results/cycling/pca_zi_norm_counts_QC_cycling_pca_D15_F.pdf"
)

load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
scd <- scdata$new(
  infos = infos,
  counts = scd$getcounts
)
save(scd, file = "results/cycling/cells_counts_QC_cycling.Rdata")


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
