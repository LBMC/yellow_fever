rm(list = ls())
setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cycling/cells_counts_QC_cycling.Rdata")

system("mkdir -p results/clonality/")

for (sex in c("M", "F")) {
  for (days in list(c("D15"), c("D136", "D593"))) {
    if (sex %in% "F" & days[1] %in% "D136") {
      days <- c("D90")
    }
    b_cells <- scd$getfeature("QC_good") %in% T &
      !is.na(scd$getfeature("DEA_cell_type")) &
      scd$getfeature("day") %in% days &
      scd$getfeature("sex") %in% sex
    # we get the list of genes to display
    infos_file <- "data/DotPlotOrder_Male_d136_d593.csv"
    if (days[1] %in% "D15") {
      infos_file <- "data/DotPlotOrder_Male_d15.csv"
    }
    genes_lists <- read.table(infos_file,
      fill = T,
      h = TRUE,
      sep = ",")
    # the first column of the table correspond to the clone to display
    clone_list <- as.vector(genes_lists[!is.na(genes_lists[, 1]), 1])
    # for each gene category in the previous genes table
    for (gene_type in colnames(genes_lists)[-1]){
      # we select the right column in the genes_list table
      c_select <- which(colnames(genes_lists) %in% gene_type)
      genes_list <- genes_lists[, c_select]
      genes_list <- genes_list[!(genes_list %in% "")]
      scRNAtools::dotplot(
        scd$select(b_cells = b_cells),
        title = paste("D15", gene_type),
        clones_order = clone_list,
        genes_order = genes_list,
        file = paste0(
          "results/clonality/cells_counts_QC_dotplot_",
          gene_type, "_",
          ifelse(length(days) == 2, "D100+", days[1]),
          "_", sex, ".pdf"
        )
      )
    }
  }
}

