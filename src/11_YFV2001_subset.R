rm(list=ls())
load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")
b_cells <- scd$getfeature("sex") %in% "M" &
  scd$getfeature("day") %in% c("D15", "D136", "D593")
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")
options(warn=1)

infos <- read.csv("data/sampleinfo_P1_P2.csv", sep = ";", h = T)
b_cells <- scd$getfeature("id") %in% infos$id

s_scd <- scdata$new(
  infos = infos,
  counts = scd$select(b_cells = b_cells)$getcounts
)

mbatch_P1_P2_DEA <- DEA(
  scd = s_scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + GATE",
  b_cells = rep(T, s_scd$getncells),
  cpus = 8,
  v = T,
  folder_name = paste0("results/cell_type/mbatch_sampleinfos_P1_P2_DEA")
)
save(
  mbatch_P1_P2_DEA,
  file = paste0("results/cell_type/mbatch_sampleinfos_P1_P2_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_P1_P2_DEA$padj)))
print(table(mbatch_P1_P2_DEA$padj < 0.05))

load(paste0("results/cell_type/mbatch_sampleinfos_P1_P2_DEA.Rdata"))
print(day)
print(table(is.na(mbatch_P1_P2_DEA$padj)))
print(table(mbatch_P1_P2_DEA$padj < 0.05))
write.csv(
  mbatch_P1_P2_DEA,
  file = paste0("results/cell_type/mbatch_sampleinfos_P1_P2_DEA.csv")
)
system("sh src/dump_dropbox.sh")

DEA_genes <- mbatch_P1_P2_DEA$gene[mbatch_P1_P2_DEA$padj < 0.05]

pca_plot(s_scd$select(gene = DEA_genes),
         color = "GATE",
         color_name = "antigen",
         tmp_file = "results/tmp/pca_P1_VS_P2_DEG.Rdata")
ggsave(file = "results/cell_type/pca/pca_P1_vs_P2_DEG.pdf")

pca_plot(s_scd,
         color = "GATE",
         color_name = "antigen",
         tmp_file = "results/tmp/pca_P1_VS_P2.Rdata")
ggsave(file = "results/cell_type/pca/pca_P1_vs_P2_all_genes.pdf")

pCMF_plot(s_scd$select(gene = DEA_genes),
          color = "GATE",
          color_name = "antigen",
          tmp_file = "results/tmp/pcmf_P1_VS_P2_DEG.Rdata",
          file = "results/cell_type/pcmf/pcmf_P1_vs_P2_all_genes_DEG.pdf")

system("sh src/dump_dropbox.sh")
