rm(list = ls())
setwd("~/projects/mold/yellow_fever")
devtools::load_all("pkg/", reset = T)
load("results/cycling/cells_counts_QC_cycling.Rdata")

days <- c("D15", "D136", "D593")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "M"

is_LLC <- function(scd) {
  clone_names <- table(
    factorize(scd$getfeature("clonality")),
    factorize(scd$getfeature("day"))
  )
  clone_names <- clone_names > 1
  clone_names <- clone_names[, c(2,1,3)]
  presence_matrix <- c(111,
                       101)
  clone_names_LLC <- rownames(clone_names)[apply(clone_names, 1, function(x){
    x <- as.numeric(paste0(as.numeric(x), collapse = ""))
    x %in% presence_matrix
    })]
  presence_matrix <- c(100)
  clone_names_SLC <- rownames(clone_names)[apply(clone_names, 1, function(x){
    x <- as.numeric(paste0(as.numeric(x), collapse = ""))
    x %in% presence_matrix
    })]
  return(list(LLC = clone_names_LLC[!(clone_names_LLC %in% "")],
              SLC = clone_names_SLC[!(clone_names_SLC %in% "")]))
}
clone <- is_LLC(scd$select(b_cells = b_cells))
LLC <- ifelse(scd$getfeature("clonality") %in% clone$LLC, "LLC",
              ifelse(scd$getfeature("clonality") %in% clone$SLC, "SLC", NA))
scd$setfeature("LLC", LLC)

day <- "D15"
system(
  paste0("rm -R results/LLC/mbatch_", day, "_LLC_DEA")
)
system(
  paste0("mkdir -p results/LLC/mbatch_", day, "_LLC_DEA")
)
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("sex") %in% "M" &
  scd$getfeature("day") %in% day &
  !is.na(scd$getfeature("LLC")) 


table(scd$select(b_cells = b_cells)$getfeature("LLC"), scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"))

mbatch_LLC_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + LLC",
  b_cells = b_cells,
  cpus = 11,
  v = T,
  folder_name = paste0("results/LLC/mbatch_", day, "_LLC_DEA")
)
save(
  mbatch_LLC_DEA,
  file = paste0("results/LLC/mbatch_", day, "_LLC_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(day)
print(table(is.na(mbatch_LLC$padj)))
print(table(mbatch_LLC_DEA$padj < 0.05))
write.csv(
  mbatch_LLC_DEApadj,
  file = paste0("results/LLC/mbatch_", day, "_LLC_DEA.csv")
)

################################################################################

day <- "D15"
b_genes <- !is.na(mbatch_LLC_DEA$padj) &
  mbatch_LLC_DEA$padj < 0.05
DEA_genes <- mbatch_LLC_DEA$gene[b_genes]
scd_norm <- zinorm(
  scd = scd$select(
    b_cells = b_cells & scd$getfeature("day") %in% day,
    genes = DEA_genes
  ),
  cpus = 11,
  file = paste0("results/tmp/zi_norm_cells_counts_QC_DEA_LLC_",
    day, ".RData")
)
system(paste0("rm results/tmp/pca_zi_norm_counts_QC_DEA_LLC_", day, ".Rdata"))

length(DEA_genes)
dim(scd_norm$getcounts)

scRNAtools::pca_plot(
  scd_norm,
  color = "LLC", color_name = "cell_type_color",
  tmp_file = paste0("results/tmp/pca_zi_norm_counts_QC_DEA_LLC_", day, ".Rdata"),
  main = "LLC vs SLC"
)
ggsave(file = paste0(
  "results/cell_type/pca/pca_zi_norm_counts_QC_DEA_cell_type_pca", day, ".pdf"
))
scRNAtools::pCMF_plot(
  scd_norm,
  color = "LLC", color_name = "cell_type_color",
  tmp_file = paste0("results/tmp/pCMF_zi_norm_counts_QC_DEA_LLC_", day, ".Rdata"),
  main = "LLC vs SLC",
  ncores = 11
)
ggsave(file = paste0(
  "results/cell_type/pca/pcmf_zi_norm_counts_QC_DEA_cell_type_pca", day, ".pdf"
))

