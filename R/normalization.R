#' normalize counts in a scdata object
#'
#' @param scd scdata object
#' @param method (default = "SCnorm") to use for the normalization
#' @param cpus (default = 4) number of cpus to use (if possible)
#' @param tmp_file temporary file to save intermediate results
#' @param v (default = FALSE) verbose mode
#' @param ... other arguments for the method function
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
    v = F,
    ...
  ) {
  algo_norm <- get(paste0(method, "_normalize"))
  return(algo_norm(
    scd = scd,
    b_cells = b_cells,
    cpus = cpus,
    tmp_file = tmp_file,
    v = v,
    ...
  ))
}

#' @importFrom SCnorm SCnorm
SCnorm_normalize <- function(
    scd,
    b_cells = scd$getfeature("QC_good") %in% T,
    cpus = 4,
    tmp_file,
    v = F,
    ...
  ) {
  if (!missing(tmp_file) & file.exists(tmp_file)) {
    print("tmp file found skipping SCnorm...")
    load(tmp_file)
  } else {
    scnorm_arg = list( PrintProgressPlots = TRUE,
                      FilterCellNum = 10)
    new_args <- list(...)
    for (new_arg in names(new_args)) {
      scnorm_arg[[new_arg]] <- new_args[[new_arg]]
    }
    DataNorm <- SCnorm(
      Data = t(scd$select(b_cells = b_cells,
                          genes = ERCC(scd, minus = T))$getcounts),
      Conditions = rep(1, scd$select(b_cells = b_cells)$getncells),
      NCore=cpus,
      new_args)
    if (v) {
      GenesNotNormalized <- results(DataNorm, type="GenesFilteredOut")
      print("genes not normalized:")
      print(str(GenesNotNormalized))
    }
    if (!missing(tmp_file)) {
      save(DataNorm, file = tmp_file)
    }
  }
  counts <- scd$getcounts
  counts[b_cells, ERCC(scd, minus = T)] <- t(results(DataNorm))
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
ComBat_normalize <- function(
    scd,
    b_cells = scd$getfeature("QC_good") %in% T,
    cpus = 4,
    tmp_file,
    v = F,
    ...
  ) {
  if (!missing(tmp_file) & file.exists(tmp_file)) {
    print("tmp file found skipping ComBat...")
    load(tmp_file)
  } else {
    expressed <- scd$getgenes[
      colSums(scd$select(b_cells = b_cells,
                         genes = ERCC(scd, minus = T))$getcounts) > 0
    ]
    combat_arg = list( par.prior = F,
                      BPPARAM = bpparam("SerialParam"))
    new_args <- list(...)
    for (new_arg in names(new_args)) {
      combat_arg[[new_arg]] <- new_args[[new_arg]]
    }
    DataNorm <- ComBat(
      dat = t(ascb(scd$select(b_cells = b_cells, genes = expressed)$getcounts)),
      batch =  scd$select(b_cells = b_cells)$getfeature("batch"),
      combat_arg
    )
    if (!missing(tmp_file)) {
      save(DataNorm, expressed, file = tmp_file)
    }
  }
  b_genes <- colnames(scd$getcounts) %in% expressed
  counts <- scd$getcounts
  norm_counts <- round(ascb_inv(t(DataNorm)))
  print( norm_counts[1:10, 1:10] )
  norm_expressed <- colSums( norm_counts ) > 0
  print(summary(norm_expressed))
  norm_counts[, !norm_expressed] <- counts[b_cells, b_genes][, !norm_expressed]
  counts[b_cells, b_genes] <- norm_counts
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

#' @importFrom scran mnnCorrect
mnnCorrect_normalize <- function(
    scd,
    b_cells = scd$getfeature("QC_good") %in% T,
    cpus = 4,
    tmp_file,
    v = F,
    ...
  ) {
  if (!missing(tmp_file) & file.exists(tmp_file)) {
    print("tmp file found skipping mnnCorrect...")
    load(tmp_file)
  } else {
    expressed <- scd$getgenes[
      colSums(scd$select(b_cells = b_cells,
                         genes = ERCC(scd, minus = T))$getcounts) > 0
    ]
    batchs <- levels(as.factor(scd$select(b_cells = b_cells)$
                               getfeature("batch")))
    arg_matrix <- c()
    for (i in batchs) {
      b_batch <- b_cells & scd$getfeature("batch") %in% i
      assign(paste0("batch_", i), t(round(scd$select(b_cells = b_batch,
                                                     genes = expressed)
                                             $getcounts)))
      arg_matrix <- c(arg_matrix, paste0("batch_", i))
    }

    arg_matrix <- paste(arg_matrix, collapse = ", ")
    mnnCorrect_cmd <- paste0("DataNorm <- mnnCorrect(",
                             arg_matrix,
                             ", svd.dim = 0 ",
                             ")")
    eval(parse(text = mnnCorrect_cmd))
    if (!missing(tmp_file)) {
      save(DataNorm, expressed, file = tmp_file)
    }
  }
  b_genes <- colnames(scd$getcounts) %in% expressed
  tDataNorm <- lapply(DataNorm$corrected, t)
  norm_counts <- as.data.frame(do.call(rbind, tDataNorm))
  counts <- as.data.frame(scd$getcounts)
  counts[b_cells, b_genes] <- norm_counts
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
