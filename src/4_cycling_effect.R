setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")

system("mkdir -p results/cycling/")

b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% "D15" &
  scd$getfeature("sex") %in% "M"


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
save(regev_genes, file=paste0(outdir, "cycling/regev_genes.RData"))

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
