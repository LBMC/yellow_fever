#' scRNASeq data object
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export scdata
#' @format An \code{\link{R6Class}} generator object
#' @keywords infos counts
scdata <- R6::R6Class("scdata",
  private = list(
    genes = NULL,
    cells = NULL,
    counts = NULL,
    features = NULL,
    get_rc_counts = function(cells = NULL, genes = NULL) {
      c_num <- 1:length(private$genes)
      if (!is.null(genes)) {
        c_num <- which(private$genes %in% genes)
      }
      r_num <- 1:length(private$cells)
      if (!is.null(cells)) {
        r_num <- which(private$cells %in% cells)
      }
      cat(paste0("ncol: ", length(c_num), ".\n"))
      cat(paste0("nrow: ", length(r_num), ".\n"))
      return(list(c_num = c_num, r_num = r_num))
    },
    get_rc_features = function(cells = NULL, features = NULL) {
      c_num <- 1:length(private$features)
      if (!is.null(features)) {
        c_num <- which(colnames(private$features) %in% features)
      }
      r_num <- 1:length(private$cells)
      if (!is.null(cells)) {
        r_num <- which(private$cells %in% cells)
      }
      cat(paste0("ncol: ", length(c_num), ".\n"))
      cat(paste0("nrow: ", length(r_num), ".\n"))
      return(list(c_num = c_num, r_num = r_num))
    }
    ),
  public = list(
    initialize = function(infos = NA, counts = NA) {
      private$genes <- colnames(counts)
      private$cells <- rownames(counts)
      private$features <- as.data.frame(infos)
      private$counts <- as.matrix(counts)
      self$summary()
    },
    summary = function() {
      cat(paste0("ncol: ", ncol(private$counts), ".\n"))
      cat(paste0("nrow: ", nrow(private$counts), ".\n"))
      cat(paste0("features: ", ncol(private$features), ".\n"))
    },
    getfeature = function(feature) {
      return(private$features[[feature]])
    },
    getgene = function(gene) {
      c_num <- which(private$genes %in% gene)
      return(private$counts[, c_num])
    },
    getcountsw = function(cells = NULL, genes = NULL) {
      rc_num <- private$get_rc_counts(cells = cells, genes = genes)
      return(private$counts[rc_num$r_num, rc_num$c_num])
    },
    getfeaturesw = function(cells = NULL, features = NULL) {
      rc_num <- private$get_rc_features(cells = cells, features = features)
      return(private$features[rc_num$r_num, rc_num$c_num])
    },
    getcountso = function(cells = NULL, genes = NULL) {
      rc_num <- private$get_rc_counts(cells = cells, genes = genes)
      return(private$counts[-rc_num$r_num, -rc_num$c_num])
    },
    getfeatureso= function(cells = NULL, features = NULL) {
      rc_num <- private$get_rc_features(cells = cells, features = features)
      return(private$features[-rc_num$r_num, -rc_num$c_num])
    }
  ),
  active = list(
    getcounts = function() {
      return(private$counts)
    },
    getfeatures = function() {
      return(private$features)
    },
    getcells = function() {
      return(private$cells)
    },
    getgenes = function() {
      return(private$genes)
    }
  )
)

#' load scRNASeq data
#'
#' @param infos a text tabular information file on cells with cells in rows and
#' features in columns
#' @param counts a text tabular count file for cells in rows and genes in
#' columns
#' @param ... additional argument to be passed to read.table function
#' @return a list with infos and counts associated and formated for the others
#' scRNASeq functions
#' @examples
#' \dontrun{
#' data <- load_data('data/infos.csv', 'data/counts.csv')
#' }
#' @export load_data
load_data <- function(infos, counts, ...) {
  data <- scdata$new(
    infos = utils::read.table(infos, ...),
    counts = utils::read.table(counts, ...)
    )
  return(data)
}
