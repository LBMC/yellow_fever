#' QC for scRNASeq data
#'
#' @param scd is a scRNASeq data object
#' @param iter is the number of bootstrap sample to perform
#' @param sample_size of the cells sample to draw (the number of blanks for examples)
#' @param is_blank boolean vector set to TRUE if a cells is a blank
#' @param output_file if non empty file to save the classification results
#' @return a list() with the number of time each cells is classied as 'good' (
#' non is_blank looking) or 'bad' (blank looking)
#' @examples
#' \dontrun{
#' classification <- QC_boot(scd, 500000)
#' }
#' @importFrom stats predict
#' @export QC_boot
QC_boot <- function(
  scd,
  iter,
  sample_size = length(which(scd$getfeature("cell_number") == 0)),
  is_blank = scd$getfeature("cell_number") == 0,
  output_file = ""
) {
  cell_id <- 1:scd$getncells
  classification <- list(
    good = rep(0, scd$getncells),
    bad = rep(0, scd$getncells))
  for (i in 1:iter) {
    print(paste0("iteration ", i))
    # index of non-blank cells to select for the QC_fit
    cell_sample <- sample(cell_id[!is_blank], replace = T, sample_size)
    selected <- ifelse(cell_id %in% cell_sample, TRUE, FALSE)
    # we also select the blanks for the positive sample
    selected[is_blank] <- TRUE
    # we define the label for the set of selected cells
    class_labels <- rep(TRUE, sum(selected))
    class_labels[is_blank[selected]] <- FALSE
    # fit
    fit <- QC_fit(class_labels, scd$getcounts[selected, ])
    prediction <- as.vector(stats::predict(fit, scd$getcounts[!selected, ]))
    # counts the number of time a given cell is classified as good or bad
    classification$bad[!selected] <- classification$bad[!selected] +
      ifelse(prediction, 0, 1)
    classification$good[!selected] <- classification$good[!selected] +
      ifelse(prediction, 1, 0)
  }
  if (output_file != "") {
    print(paste0("saving output in: ", output_file))
    save(classification, file = output_file)
  }
  return(classification)
}

#' @importFrom e1071 svm
QC_fit <- function(class_labels, counts) {
  class_weights <- 100 / table(as.factor(class_labels))
  return(
    suppressWarnings(
      e1071::svm(
        x = counts,
        y = as.factor(class_labels),
        kernel = "linear",
        class.weights = class_weights,
        type = "C-classification",
        scale = TRUE
  )))
}

#' helper function for QC for scRNASeq data
#'
#' @param paraload_file paraload table file
#' @param bootstraps total number of bootstrap to perform
#' @param job_boot_number  of bootstrap per paraload jobs
#' @return write a paraload ready table with the parameters (as line) on which
#' to run every QC_boot() commands
#' @examples
#' \dontrun{
#' QC_paraload_parameters('results/QC/paraload_file.txt')
#' }
#' @export QC_paraload_parameters
QC_paraload_parameters <- function(
  paraload_file,
  bootstraps = 500000,
  job_boot_number = 50
) {
  boot_number <- bootstraps / job_boot_number
  print(boot_number)
  parameters <- c()
  parameters <- rbind(parameters,
                     cbind(rep("QC", boot_number),
                           rep(job_boot_number, boot_number),
                           1:boot_number)
                     )
  utils::write.table(parameters,
    file = paraload_file,
    append = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE)
}

#' helper function to run QC in a script
#'
#' @param scd_file RData file with an scTools data object 'data' saved in it
#' @param QC_folder folder where to save QC_boot results
#' @param Args agument passsed when used within a script
#' @return write a paraload ready table with the parameters (as line) on which
#' to run every QC_boot() commands
#' @examples
#' \dontrun{
#' QC_paraload_parameters('results/QC/paraload_file.txt')
#' }
#' @export QC_pbs
QC_pbs <- function(scd_file, QC_folder, Args = commandArgs()) {
  print(getwd())
  if (length(Args) > 4) {
      args        <- utils::read.table(Args[6])
      boot_number <- as.numeric(args[2])
      file         <- args[3]
      print(date())
      load(scd_file)
      genes <- ERCC(scd, minus = T)
      b_cells <- scd$getfeature("cell_number") == 1 |
        scd$getfeature("cell_number") == 0
      print(paste0(
        "QC_boot with ",
        length(genes),
        " genes and ",
        length(which(b_cells)),
        " cells"))
      if (!("to_QC") %in% colnames(scd$getfeatures)){
        scd$setfeature("to_QC", rep(T, scd$getncells))
      }
      classification <- QC_boot(
        scd = scd$select(
          genes = genes,
          b_cells = b_cells & scd$getfeature("to_QC")),
        iter = boot_number,
        output_file = paste0(QC_folder, "QC_", file, ".Rdata"))
  }
  utils::write.table(
    c(file, boot_number, classification$good),
    file = Args[7],
    append = F, row.names = F, col.names = F, quote = F
  )
  print("QC_pbs done.")
}

#' helper function for QC for scRNASeq data
#'
#' @param scd scdata object
#' @param paraload_folder dir where the QC_pbs output are saved
#' @return return scdata object with a "QC_score" feature
#' @examples
#' \dontrun{
#' QC_paraload_parameters('results/QC/paraload_file.txt')
#' }
#' @export QC_load_bootstraps
QC_load_bootstraps <- function(scd, paraload_folder, rt_result = F) {
  scd$addfeature("QC_score")
  files_list <- get_files(path = paraload_folder, regexp = ".*Rdata")
  print(paste0("loading QC results..."))
  print(head(files_list))
  b_cells <- scd$getfeature("cell_number") == 1 |
    scd$getfeature("cell_number") == 0
  if (!("to_QC") %in% colnames(scd$getfeatures)){
    scd$setfeature("to_QC", rep(T, scd$getncells))
  }
  b_cells <- b_cells & scd$getfeature("to_QC")
  classification_summary <- list(
    bad = rep(0, scd$getncells),
    good = rep(0, scd$getncells)
  )
  for (classif_file in files_list) {
    load(classif_file)
    classification_summary$bad[b_cells] =
      classification_summary$bad[b_cells] + classification$bad
    classification_summary$good[b_cells] =
      classification_summary$good[b_cells] + classification$good
  }
  QC_score <- rep(NA, scd$getncells)
  QC_score[b_cells] <- classification_summary$good[b_cells] /
    (classification_summary$good[b_cells] + classification_summary$bad[b_cells])
  QC_score[scd$getfeature("cell_number") == 0] <- 0
  scd$setfeature("QC_score", QC_score)
  print("done.")
  if (rt_result){
    return(scd)
  }
}

#' helper function for QC for scRNASeq data
#'
#' @param scd scdata object
#' @param is_blank boolean vector indicating lines where there is not cells
#' @return return scdata object with a "QC_good" feature
#' @examples
#' \dontrun{
#' QC_classification(scd)
#' }
#' @importFrom stats predict
#' @export QC_classification
QC_classification <- function(
  scd,
  is_blank = scd$getfeature("cell_number") == 0,
  rt_result = F,
  quant = ""
) {
  # select cells to classify
  b_cells <- (scd$getfeature("cell_number") == 1 |
              scd$getfeature("cell_number") == 0) &
                             !is.na(scd$getfeature("QC_score"))
  if (!("to_QC") %in% colnames(scd$getfeatures)) {
    scd$setfeature("to_QC", rep(T, scd$getncells))
  }
  b_cells <- b_cells & scd$getfeature("to_QC")
  data_fit <- data.frame(pos = which(b_cells),
                         label = ifelse((b_cells & is_blank)[b_cells], F, NA),
                         score = scd$select(b_cells = b_cells)$getfeature("QC_score")
  )
  data_fit <- data_fit[order(data_fit$score, decreasing = T), ]
  true_sample_size <- length(which(!data_fit$label))
  if (!missing(quant)) {
    true_sample_size <- round(length(which(is.na(data_fit$label))) * quant)
  }
  data_fit$label[1:true_sample_size] <- TRUE
  fit <- QC_fit(data_fit$label[!is.na(data_fit$label)],
                scd$getcounts[data_fit$pos[!is.na(data_fit$label)], ])
  to_predict <- data_fit$pos[is.na(data_fit$label)]
  prediction <- stats::predict(
    fit,
    scd$getcounts[to_predict, ],
    na.action = na.fail
  )
  QC_good <- rep(NA, scd$getncells)
  QC_good[data_fit$pos[data_fit$label %in% F]] <- F
  QC_good[data_fit$pos[data_fit$label %in% T]] <- T
  QC_good[to_predict] <- as.vector(prediction)
  scd$setfeature("QC_good", as.factor(QC_good))
  if (rt_result){
    return(scd)
  }
}
