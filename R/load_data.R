#' scRNASeq data object
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export scdata
#' @format An \code{\link{R6Class}} generator object
#' @keywords infos counts
scdata <- R6::R6Class("scdata",
  private = list(
    genes = NULL, # a vector of genes names
    cells = NULL, # a vector of cells id
    counts = NULL, # a matrix of counts
    features = NULL, # a data.frame of cells features
    # private methods go get position of cells and genes in rows and columns
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
    # private methods go get position of cells and features in rows and columns
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
    },
    order_by_cells = function() {
      private$counts <- private$counts[order(private$cells), ]
      private$features <- private$features[order(private$cells), ]
      private$cells <- private$cells[order(private$cells)]
    }
    ),
  public = list(
    initialize = function(infos = NA, counts = NA) {
      private$genes <- colnames(counts)
      private$cells <- rownames(counts)
      private$features <- as.data.frame(infos)
      private$counts <- as.matrix(counts)
      if (nrow(private$features) != nrow(private$counts)) {
        stop("error: number of cells differ between infos and counts")
      }
      if ("id" %in% colnames(private$features)) {
        if (length(intersect(private$cells, as.vector(private$features$id))) !=
          length(private$cells)) {
          stop("error: id's in infos don't match cells name in counts")
        }
        private$order_by_cells()
      }
      self$summary()
    },
    add = function(infos = NA, counts = NA) {
      genes <- colnames(counts)
      cells <- rownames(counts)
      features <- as.data.frame(infos)
      counts <- as.matrix(counts)
      if (nrow(features) != nrow(counts)) {
        stop("error: number of cells differ between infos and counts")
      }
      if (length(intersect(private$cells, cells)) != 0) {
        stop("error: trying to add cells already present")
      }
      if (length(intersect(private$genes, genes)) != length(private$genes)) {
        stop("error: genes set don't match existing genes set")
      }

      if ("id" %in% colnames(features)) {
        if (length(intersect(cells, as.vector(features$id))) !=
          length(cells)) {
          stop("error: id's in infos don't match cells name in counts")
        }
      }
      private$cells <- c(private$cells, cells)
      private$counts <- rbind(private$counts, counts)
      private$features <- rbind(private$features, features)
      private$order_by_cells()
    },
    summary = function() {
      cat(paste0("ncol: ", ncol(private$counts), ".\n"))
      cat(paste0("nrow: ", nrow(private$counts), ".\n"))
      cat(paste0("features: ", ncol(private$features), ".\n"))
    },
    # accessors methods
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
    getfeatureso = function(cells = NULL, features = NULL) {
      rc_num <- private$get_rc_features(cells = cells, features = features)
      return(private$features[-rc_num$r_num, -rc_num$c_num])
    }
  ),
  active = list(
    # accessors methods without arguments
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
load_data <- function(infos, counts, regexp = ".*", ...) {
  if (base::file.info(counts)$isdir) {
    return(load_multiple_file(
      infos = infos,
      counts = counts,
      regexp = regexp,
      ...
      )
    )
  } else {
    return(scdata$new(
        infos = utils::read.table(infos, ...),
        counts = utils::read.table(counts, ...)
      )
    )
  }
}

get_files <- function(counts, regexp) {
  file_list <- list.files(path = counts, full.names = TRUE, recursive = TRUE)
  file_list <- file_list[grepl(regexp, file_list, perl = T)]
  return(file_list)
}

load_multiple_file <- function(infos, counts, regexp, ...) {
  features <- utils::read.table(infos, ...)
  files_list <- get_files(counts = counts, regexp = regexp)
  data <- scdata$new(
      infos = features,
      counts = utils::read.table(files_list[1], ...),
    )
  for (counts_file in files_list[-1]) {
    data$add(
        infos = features,
        counts = utils::read.table(counts_file, ...),
      )
  }
  return(data)
}
