#' QC for scRNASeq data
#'
#' @param data is a scRNASeq data object
#' @param iter is the number of bootstrap sample to perform
#' @param size of the cells sample to draw (the number of blanks for examples)
#' @param is_blank boolean vector set to TRUE if a cells is a blank
#' @param output_file if non empty file to save the classification results
#' @return a list() with the number of time each cells is classied as 'good' (
#' non is_blank looking) or 'bad' (blank looking)
#' @examples
#' \dontrun{
#' classification <- QC_boot(data, 500000)
#' }
#' @export QC_boot
QC_boot <- function(
  data,
  iter,
  sample_size = length(which(data$getfeature["cell_number"] == 0)),
  is_blank = data$getfeature["cell_number"] == 0,
  output_file = ""
) {
  cell_id <- 1:data$getncells
  classification <- list(
    good = rep(0, data$getncells),
    bad = rep(0, data$getncells))
  for (i in 1:iter) {
    print(paste0("iteration ", i))
    # index of non-blank cells to select for the QC_fit
    selected <- ifelse(
      cell_id %in% sample(cell_id[!is_blank], replace = T, sample_size),
      TRUE, FALSE)
    # we also select the blanks for the positive sample
    selected[is_blank] <- TRUE
    # we define the label for the set of selected cells
    class_labels <- rep(TRUE, sum(selected))
    class_labels[is_blank[selected]] <- FALSE
    # fit
    fit <- QC_fit(class_labels, data$getcounts[selected, ])
    prediction <- as.vector(predict(fit, data$getcounts[!selected, ]))
    # counts the number of time a given cell is classified as good or bad
    classification$bad[!selected] <- classification$bad[!selected] +
      ifelse(prediction, 0, 1)
    classification$good[!selected] <- classification$good[!selected] +
      ifelse(prediction, 1, 0)
  }
  if (output_file != "") {
    save(classification, output_file)
  }
  return(classification)
}

#' @importFrom e1071 svm
QC_fit <- function(class_labels, counts) {
  class_weights <- 100 / table(as.factor(class_labels))
  return(
    suppressWarnings(
      svm(
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
#' @param paraload table file
#' @param total number of bootstrap to perform
#' @param number of bootstrap per paraload jobs
#' @return write a paraload ready table with the parameters (as line) on which
#' to run every QC_boot() commands
#' @examples
#' \dontrun{
#' QC_paraload_parameters('results/QC/paraload_file.txt')
#' }
#' @export QC_boot
QC_paraload_parameters <- function(
  paraload_file,
  bootstraps=500000,
  job_boot_number=50
) {
  boot_number <- bootstraps / job_boot_number
  print(boot_number)
  parameters <- c()
  parameters <- rbind(parameters,
                     cbind(rep("SVM", boot_number),
                           rep(job_boot_number, boot_number),
                           1:boot_number)
                     )
  write.table(parameters,
    file = paraload_file,
    append = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE)
}

#' helper function to run QC in a script
#'
#' @param data_file RData file with an scTools data object 'data' saved in it
#' @param args
#' @param number of bootstrap per paraload jobs
#' @return write a paraload ready table with the parameters (as line) on which
#' to run every QC_boot() commands
#' @examples
#' \dontrun{
#' QC_paraload_parameters('results/QC/paraload_file.txt')
#' }
#' @export QC_boot
QC_pbs <- function(data_file, Args = commandArgs()) {
  if (length(Args) > 4) {
      args        <- read.table(Args[6])
      boot_number <- as.numeric(args[2])
      col         <- as.numeric(args[3])
      file        <- Args[7]
      print(date())
      load(data_file, v = T)
      genes <- ERCC(data, minus = T)
      b_cells <- data$getfeature("cell_number") == 1 | \
        data$getfeature("cell_number") == 0
      QC_boot(
        data = data$select(genes = genes, b_cells = b_cells),
        iter = boot_number,
        sample_size = length(which(data$getfeature("cell_number") == 0)),
        is_blank = data$getfeature("cell_number") == 0,
        output_file = file)
  }
}

#' helper function for QC for scRNASeq data
#'
#' @param data scdata object
#' @param paraload_folder dir where the QC_pbs output are saved
#' @return return scdata object with a "QC_score" feature
#' @examples
#' \dontrun{
#' QC_paraload_parameters('results/QC/paraload_file.txt')
#' }
#' @export QC_load_bootstraps
QC_load_bootstraps <- function(data, paraload_folder) {
  data$addfeature("QC_score")
  files_list <- get_files(path = paraload_folder, regexp = "*.Rdata")
  print(paste0("loading QC results...")
  classification_summary <- list(
    bad = rep(0, data$getncells),
    good = rep(0, data$getncells)
  )
  for (classif_file in files_list) {
    load(classif_file)
    classification_summary$bad = \
      classification_summary$bad + calssification$bad
    classification_summary$good = \
      classification_summary$good + calssification$good
  }
  b_cells <- data$getfeature("cell_number") == 1 | \
        data$getfeature("cell_number") == 0
  QC_score <- rep(NA, data$getncells)
  QC_good <- rep(NA, data$getncells)
  QC_score[b_cells] <- classification_summary$good / \
    (classification_summary$good + classification_summary$bad)
  data$setfeature("QC_score", QC_score)
  return(data)
}

#' helper function for QC for scRNASeq data
#'
#' @param data scdata object
#' @param is_blank boolean vector indicating lines where there is not cells
#' @return return scdata object with a "QC_good" feature
#' @examples
#' \dontrun{
#' QC_classification(data)
#' }
#' @export QC_classification
QC_classification <- function(
  data,
  is_blank = data$getfeature["cell_number"] == 0
) {
  data$addfeature("QC_good")
  cell_id <- 1:data$getncells
  classification <- list(
    good = rep(0, data$getncells),
    bad = rep(0, data$getncells)
  )
  b_cells <- data$getfeature("cell_number") == 1 | \
        data$getfeature("cell_number") == 0
  sample_size <- length(which(data$getfeature("cell_number") == 0))
  is_good <- order(data$getfeature("QC_score"), decreasing = T)[1:sample_size]
  is_bad <- which(is_blank)
  class_labels <- rep(TRUE, sample_size * 2)
  class_labels[1:sample_size] <- FALSE
  fit <- QC_fit(class_labels, data$getcounts[c(is_bad, is_good), ])
  prediction <- as.vector(predict(
    fit,
    data$getcounts[which(!(1:data$getncells %in% c(is_bad, is_good))), ]))
  QC_good <- rep(NA, data$getncells)
  QC_good[is_bad] <- F
  QC_good[is_good] <- T
  QC_good[which(!(1:data$getncells %in% c(is_bad, is_good))), ])] <- prediction
  data$setfeature("QC_good", QC_good)
  return(data)
}
