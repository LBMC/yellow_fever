#' normalize counts in a scdata object
#'
#' @param scd scdata object
#' @param method (default = "SCnorm") to use for the normalization
#' @param cpus (default = 4) number of cpus to use (if possible)
#' @param v (default = FALSE) verbose mode
#' @return return scdata object
#' @examples
#' \dontrun{
#' scd_norm <- normalize(scd)
#' }
#' @import SCnorm
#' @export SCnorm
normalize <- function(
    scd,
    method = "SCnorm",
    cpus = 4,
    tmp_file,
    v = F
  ) {
  if (method == "SCnorm"){
    if (!missing(tmp_file) & file.exists(tmp_file)) {
      print("tmp file found skipping pca...")
      load(tmp_file)
    } else {
      countDeptEst <- plotCountDepth(
        Data = t(scd$select(b_cells = b_cells)$getcounts),
        Conditions = rep(1, scd$select(b_cells = b_cells)$getncells),
        FilterCellProportion = .1,
        NCores=cpus)
      DataNorm <- SCnorm(
        Data = t(scd$select(b_cells = b_cells)$getcounts),
        Conditions = rep(1, scd$select(b_cells = b_cells)$getncells),
        PrintProgressPlots = TRUE,
        FilterCellNum = 10,
        NCore=cpus)
      if (v) {
        GenesNotNormalized <- results(DataNorm, type="GenesFilteredOut")
        print("genes not normalized:")
        print(str(GenesNotNormalized))
      }
      if (!missing(tmp_file)) {
        save(countDeptEst, DataNorm, file = file)
      }
    }
    b_cells = scd$getfeature('to_QC') & scd$getfeature("QC_good") %in% T
    counts <- scd$getcounts
    counts[b_cells] <- t(results(DataNorm))
    return(scdata$new(
      infos = scd$getfeatures,
      counts = counts,
      v = v
    )
    )
  }
}
