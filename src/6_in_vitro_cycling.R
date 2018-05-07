setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)

load("results/cycling/CB_counts_QC_cycling.Rdata")
b_cells <- scd$getfeature("sex") %in% "M" &
  scd$getfeature("day") %in% c("D15", "D136", "D593")
infos <- scd$getfeatures
CB_counts <- scd$getcounts
load("results/QC/CB_counts_QC.Rdata")
cells_counts <- scd$getcounts

load("results/cycling/CB_counts_QC_cycling.Rdata")
b_cells <- scd$getfeature("sex") %in% "F" &
  scd$getfeature("day") %in% c("D15", "D90")
infos[b_cells, ] <- scd$select(b_cells = b_cells)$getfeatures
CB_counts[b_cells, ] <- scd$select(b_cells = b_cells)$getcounts
load("results/QC/CB_counts_QC_F.Rdata")
cells_counts[b_cells, ] <- scd$select(b_cells = b_cells)$getcounts

load("results/QC/CB_counts_QC_in_vitro_P1902_P3128.Rdata")
experiment <- c("P1902", "P3128")
day <- "InVitro"
b_cells <- scd$getfeature("day") %in% day &
  scd$getfeature("experiment") %in% experiment &
  scd$getfeature("cell_number") %in% 1
c_cells <- which(colnames(infos) %in% colnames(scd$getfeatures))
c_cells_bis <- which(colnames(scd$getfeatures) %in% colnames(infos))
infos[b_cells, c_cells] <- scd$select(b_cells = b_cells)$getfeatures[, c_cells_bis]
c_cells_bis <- which(colnames(scd$getfeatures) %in%
  setdiff(colnames(scd$getfeatures), colnames(infos)))
infos <- cbind(infos, scd$getfeatures[, c_cells_bis])
CB_counts[b_cells, ] <- scd$select(b_cells = b_cells)$getcounts
load("results/QC/cells_counts_QC_in_vitro_P1902_P3128.Rdata")
cells_counts[b_cells, ] <- scd$select(b_cells = b_cells)$getcounts

founder_phenotype <- rep(NA, scd$getncells)
founder_infos <- read.csv("data/FounderCellP1902_P3128.csv")
for (i in 1:nrow(founder_infos)) {
  b_select <- b_cells <- scd$getfeature("day") %in% day &
    scd$getfeature("experiment") %in% founder_infos$Project[i] &
    scd$getfeature("clonality") %in% founder_infos$Clone.ID[i]
  founder_phenotype[b_select] <- as.vector(founder_infos$Founder_Type_Simple[i])
}

scd <- scdata$new(
  infos = infos,
  counts = CB_counts
)

scd$setfeature(
  "founder_phenotype",
  founder_phenotype
)
scd$setfeature(
  "experiment",
  gsub("(P\\d+)_\\d+", "\\1", scd$getfeature("id"), perl = T)
)

save(scd, file = "results/cell_type/CB_counts_QC_all.Rdata")

scd <- scdata$new(
  infos = infos,
  counts = cells_counts
)

scd$setfeature(
  "founder_phenotype",
  founder_phenotype
)
scd$setfeature(
  "experiment",
  gsub("(P\\d+)_\\d+", "\\1", scd$getfeature("id"), perl = T)
)

save(scd, file = "results/cell_type/cells_counts_QC_all.Rdata")



system("mkdir -p results/cycling/")

# we try to refine the regev cell-cycle genes list
load(file="results/cycling/regev_genes.RData", v = T)
load("results/cell_type/CB_counts_QC_all.Rdata")

cycling_score <- scd$getfeature("cycling_score")
cycling <- scd$getfeature("cycling")
pcycling <- scd$getfeature("pcycling")
genes_cycling <- regev_genes
require("mixtools")
for (experiment in c("P1902", "P3128")) {
  day <- "InVitro"
  b_cells <- scd$getfeature("day") %in% day &
    scd$getfeature("experiment") %in% experiment &
    scd$getfeature("cell_number") %in% 1
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

save(scd,
  file = "results/cycling/CB_counts_QC_cycling_invitro_P1902_P3128.Rdata"
)
infos <- scd$getfeatures
write.csv(
  infos,
  file = paste0("results/cycling/cell_type_infos_invitro_P1902_P3128.csv")
)


load("results/cycling/CB_counts_QC_cycling_invitro_P1902_P3128.Rdata")
day <- "InVitro"
for (experiment in c("P1902", "P3128")) {
  b_cells <- scd$getfeature("day") %in% day &
    scd$getfeature("experiment") %in% experiment &
    scd$getfeature("cell_number") %in% 1

  tmp_infos <- scd$select(b_cells = b_cells)$getfeatures
  tmp_infos$clonality <- as.factor( as.vector(tmp_infos$clonality) )
  clone_cycling <- by(
    tmp_infos$cycling_score, tmp_infos$clonality, mean
  )
  tmp_infos$clonality <- factor(
    tmp_infos$clonality,
    levels = names(clone_cycling)[order(clone_cycling)])
  g <- ggplot(tmp_infos,
    aes(x = clonality,
        y = cycling_score,
        color = founder_phenotype)) +
    geom_jitter() +
    geom_violin(alpha = 0.5) +
    theme_bw() +
    scale_color_manual(
      values = cell_type_palette(
        levels(factorize(tmp_infos$founder_phenotype)))) +
    labs(x = "clones",
        y =  "cell-cycle score",
        color = "founder phenotype",
        title = day)
  print(g)
  ggsave(file = paste0(
    "results/cycling/violing_invitro_",
    experiment,
    "_cycling_vs_founder_phenotype.pdf"
  ))
}

