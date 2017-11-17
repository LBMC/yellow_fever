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
#' @export normalize
normalize <- function(
    scd,
    b_cells = scd$getfeature("QC_good") %in% T,
    method = "SCnorm",
    cpus = 4,
    tmp_file,
    v = F
  ) {
  algo_norm <- get(paste0(method, "_normalize"))
  return(algo_norm(
    scd = scd,
    b_cells = b_cells,
    cpus = cpus,
    tmp_file = tmp_file,
    v = v
  ))
}

#' @import SCnorm
SCnorm_normalize <- function(
    scd,
    b_cells = scd$getfeature("QC_good") %in% T,
    cpus = 4,
    tmp_file,
    v = F
  ) {
  if (!missing(tmp_file) & file.exists(tmp_file)) {
    print("tmp file found skipping SCnorm...")
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
  counts <- scd$getcounts
  counts[b_cells] <- t(results(DataNorm))
  rownames(counts) <- rownames(scd$getcounts)
  colnames(counts) <- colnames(scd$getcounts)
  return(
    scdata$new(
      infos = scd$getfeatures,
      counts = counts,
      v = v
    )
  )
}

#' @importFrom sva ComBat
combat_normalize <- function(
    scd,
    b_cells = scd$getfeature("QC_good") %in% T,
    cpus = 4,
    tmp_file,
    v = F
  ) {
  if (!missing(tmp_file) & file.exists(tmp_file)) {
    print("tmp file found skipping SCnorm...")
    load(tmp_file)
  } else {
    DataNorm <- ComBat(
      dat = t(ascb(scd$select(b_cells = b_cells)$getcounts)),
      batch =  scd$select(b_cells = b_cells)$getfeature("batch"),
      BPPARAM = bpparam("SerialParam")
    )
    if (!missing(tmp_file)) {
      save(DataNorm, file = file)
    }
  }
  counts <- scd$getcounts
  counts[b_cells] <- round(ascb_inv(t(DataNorm)))
  rownames(counts) <- rownames(scd$getcounts)
  colnames(counts) <- colnames(scd$getcounts)
  return(
    scdata$new(
      infos = scd$getfeatures,
      counts = counts,
      v = v
    )
  )
}
