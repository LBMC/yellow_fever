rm(list = ls())
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cycling/cells_counts_QC_cycling.Rdata")

system("mkdir -p results/clonality/")

for (sex in c("M", "F")) {
  days <- c("D15", "D136", "D593")
  if (sex %in% "F") {
    days <- c("D15", "D90")
  }

  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% days[1] &
    scd$getfeature("sex") %in% sex
  tmp_infos <- scd$select(b_cells = b_cells)$getfeatures
  tmp_infos$clonality <- factorize(tmp_infos$clonality)
  clone_size <- table(tmp_infos$clonality)
  clone_size <- clone_size[order(clone_size, decreasing = TRUE)]
  clone_size <- clone_size[clone_size > 3][-1]
  clone_list_D15 <- names(clone_size)
  r_select_clone <- tmp_infos$clonality %in% clone_list_D15
  tmp_infos$clonality[!r_select_clone] <- NA
  tmp_infos$clonality <- factor(tmp_infos$clonality,
    levels = clone_list_D15)
  clones_names <- levels(tmp_infos$clonality)
  clones_av_p_MEM <- by(tmp_infos$pDEA_cell_type, tmp_infos$clonality, median)
  clones_av_p_MEM <- as.vector(clones_av_p_MEM)
  clones_names <- clones_names[order(clones_av_p_MEM,
    decreasing = TRUE)]
  clones_color <- clonality_MEM_palette(
    levels(tmp_infos$clonality),
    clones_av_p_MEM)
  tmp_infos$clonality <- factor(tmp_infos$clonality,
    levels = clones_names)
  tmp_infos$wilcox_pval <- NULL
  for (clone in levels(tmp_infos$clonality)){
    r_select <- tmp_infos$clonality %in% clone
    tmp_infos$wilcox_pval[r_select] <-
      wilcox.test(tmp_infos$pDEA_cell_type[!r_select], tmp_infos$pDEA_cell_type[r_select])$p.value
  }
  tmp_infos$wilcox_pval <- paste("wilcox p-value :",
    signif(tmp_infos$wilcox_pval,
      digits = 3))
  g <- ggplot(tmp_infos[!is.na(tmp_infos$clonality), ],
    aes(x = pDEA_cell_type,
      color = clonality,
      fill = clonality)) +
    geom_density(fill="NA") +
    geom_histogram() +
    scale_color_manual(values = clones_color) +
    scale_fill_manual(values = clones_color) +
    geom_text(aes(x = 0.5, y = 6, label=wilcox_pval, color = NULL), size = 2) +
    facet_wrap(~clonality,
      ncol = 1,
      switch = "y", scales = "free_y") +
    theme_bw() +
    labs(title = "D15",
      x = "Memory score",
      color = "clone",
      fill = "clone") +
    theme(
      axis.text.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      panel.margin = unit(0, "lines"),
      panel.margin.x = unit(0, "lines"),
      panel.margin.y = unit(0, "lines"),
      axis.ticks = element_blank(),
      strip.text.y = element_blank(),
      strip.background = element_blank())
  print(g)
  ggsave(file = paste0(
    "results/clonality/pMEM_density_by_clone_D15_", sex, ".pdf"
  ))

  g <- ggplot(tmp_infos,
    aes(x = pDEA_cell_type)) +
    geom_histogram(color="NA") +
    theme_bw() +
    labs(title = "D15",
      x = "Memory score") +
    theme(
      axis.text.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      panel.margin = unit(0, "lines"),
      panel.margin.x = unit(0, "lines"),
      panel.margin.y = unit(0, "lines"),
      axis.ticks = element_blank(),
      strip.text.y = element_blank(),
      strip.background = element_blank())
  print(g)
  ggsave(file = paste0("results/clonality/pMEM_density_D15_", sex, ".pdf"))

  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% days[2:length(days)] &
    scd$getfeature("sex") %in% sex
  tmp_infos <- scd$select(b_cells = b_cells)$getfeatures
  tmp_infos$clonality <- factorize(tmp_infos$clonality)
  clone_size <- table(tmp_infos$clonality)
  clone_size <- clone_size[order(clone_size, decreasing = TRUE)]
  clone_size <- clone_size[clone_size > 3][-1]
  clone_list_D100 <- names(clone_size)
  r_select_clone <-  tmp_infos$clonality %in% clone_list_D100
  tmp_infos$clonality[!r_select_clone] <- NA
  tmp_infos$clonality <- factor(tmp_infos$clonality,
    levels = clone_list_D100)
  clones_names <- levels(tmp_infos$clonality)
  clones_av_p_MEM <- by(tmp_infos$pDEA_cell_type, tmp_infos$clonality, median)
  clones_av_p_MEM <- as.vector(clones_av_p_MEM)
  clones_names <- clones_names[order(clones_av_p_MEM,
    decreasing = TRUE)]
  clones_color <- clonality_MEM_palette(
    levels(tmp_infos$clonality),
    clones_av_p_MEM)
  tmp_infos$clonality <- factor(tmp_infos$clonality,
    levels = clones_names)
  tmp_infos$wilcox_pval <- NULL
  for (clone in levels(tmp_infos$clonality)){
    r_select <- tmp_infos$clonality %in% clone
    tmp_infos$wilcox_pval[r_select] <-
      wilcox.test(tmp_infos$pDEA_cell_type[!r_select],
                  tmp_infos$pDEA_cell_type[r_select])$p.value
  }
  tmp_infos$wilcox_pval <- paste("wilcox p-value :",
    signif(tmp_infos$wilcox_pval,
      digits = 3))
  g <- ggplot(tmp_infos[!is.na(tmp_infos$clonality), ],
    aes(x = pDEA_cell_type,
      color = clonality,
      fill = clonality)) +
    geom_density(fill="NA") +
    geom_histogram() +
    scale_color_manual(values = clones_color) +
    scale_fill_manual(values = clones_color) +
    geom_text(aes(x = 0.5, y = 6, label=wilcox_pval, color = NULL), size = 2) +
    facet_wrap(~clonality,
      ncol = 1,
      switch = "y", scales = "free_y") +
    theme_bw() +
    labs(title = "D100+",
      x = "Memory score",
      color = "clone",
      fill = "clone") +
    theme(
      axis.text.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      panel.margin = unit(0, "lines"),
      panel.margin.x = unit(0, "lines"),
      panel.margin.y = unit(0, "lines"),
      axis.ticks = element_blank(),
      strip.text.y = element_blank(),
      strip.background = element_blank())
  print(g)
  ggsave(file = paste0(
    "results/clonality/pMEM_density_by_clone_D100+_", sex, ".pdf"))

  g <- ggplot(tmp_infos,
    aes(x = pDEA_cell_type)) +
    geom_histogram(color="NA") +
    theme_bw() +
    labs(title = "D100+",
      x = "Memory score") +
    theme(
      axis.text.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      panel.margin = unit(0, "lines"),
      panel.margin.x = unit(0, "lines"),
      panel.margin.y = unit(0, "lines"),
      axis.ticks = element_blank(),
      strip.text.y = element_blank(),
      strip.background = element_blank())
  print(g)
  ggsave(file = paste0("results/clonality/pMEM_density_D100+", sex, ".pdf"))
}
