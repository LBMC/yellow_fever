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
    features = NULL
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
    }
  ),
  active = list(
    getcounts = function() {
      return(private$counts)
    }
  )
)

#' load scRNASeq data
#'
#' @param infos a text tabular information file on cells with cells in rows and
#' features in columns
#' @param counts a text tabular count file for cells in rows and genes in
#' columns
#' @return a list with infos and counts associated and formated for the others
#' scRNASeq functions
#' @examples
#' \dontrun{
#' data <- load_data('data/infos.csv', 'data/counts.csv')
#' }
#' @export load_data
load_data <- function(infos, counts) {
  data <- scdata$new(
    infos = read.table(infos),
    counts = read.table(counts)
    )
  return(data)
}
