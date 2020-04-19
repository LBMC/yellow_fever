load(file = "results/sce_scnorm.Rdata", verbose = T)


genes_PLS <- read.csv("~/data/yellow_fever/2017_11_28_List_Laurent_Genes_PLS.csv")
surface_marker <- c()
genes_marker <- c()
for (marker_type in colnames(genes_PLS)) {
  for (marker in genes_PLS[[marker_type]]) {
    if (marker %in% scd$getgenes) {
      genes_marker <- c(genes_marker, marker)
    }
    if (marker %in% colnames(scd$getfeatures)) {
      surface_marker <- c(surface_marker, marker)
    }
  }
}
