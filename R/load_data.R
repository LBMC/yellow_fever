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
    # private method to get position of cells and genes in rows and columns
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
    # private method to get position of cells and features in rows and columns
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
    # private method to get position of cells in common between features and
    # counts
    get_common_cells = function(features, counts) {
      common_cells <- intersect(as.vector(features$id), rownames(counts))
      print(
        paste0(length(common_cells),
          " cells in common between infos and counts"))
      r_features <- which(as.vector(features$id) %in% common_cells)
      r_counts <- which(rownames(counts) %in% common_cells)
      return(list(features = r_features, counts = r_counts))
    },
    order_by_cells = function() {
      private$counts <- private$counts[order(private$cells), ]
      private$features <- private$features[order(private$cells), ]
      private$cells <- private$cells[order(private$cells)]
    },
    transpose_count = function(counts) {
      counts <- as.matrix(counts)
      if (nrow(private$features) < nrow(counts)) {
        print("transposing counts...")
        return(t(counts))
      }
      return(counts)
    },
    display_dim = function(infos, counts) {
      print(
        paste0("dim of infos :",
          nrow(infos),
          " x ",
          ncol(infos)
        ))
      print(
        paste0("dim of counts :",
          nrow(counts),
          " x ",
          ncol(counts)
        ))
    }
    ),
  public = list(
    initialize = function(infos = NA, counts = NA) {
      private$features <- as.data.frame(infos)
      private$counts <- private$transpose_count(counts)
      private$genes <- colnames(private$counts)
      private$cells <- rownames(private$counts)
      private$display_dim(private$features, private$counts)
      if ("id" %in% colnames(private$features)) {
        print("id column found in infos")
        if (length(intersect(private$cells, as.vector(private$features$id)))
          == 0) {
            print(private$cells)
            print(as.vector(private$features$id))
          stop("error: id's in infos don't match cells name in counts")
        }
        rr_num <- private$get_common_cells(private$features, private$counts)
        private$features <- private$features[rr_num$features, ]
        private$counts <- private$counts[rr_num$counts, ]
        private$order_by_cells()
      } else {
        if (nrow(private$features) != nrow(private$counts)) {
          stop(
            paste0("error: number of cells (",
              nrow(private$features),
              ") differ between infos and counts (",
              nrow(private$counts)))
        }
      }
      self$summary()
    },
    add = function(infos = NA, counts = NA) {
      features <- as.data.frame(infos)
      counts <- private$transpose_count(counts)
      genes <- colnames(counts)
      cells <- rownames(counts)
      private$display_dim(features, counts)
      if (length(intersect(private$cells, cells)) != 0) {
        stop("error: trying to add cells already present")
      }
      if (length(intersect(private$genes, genes)) != length(private$genes)) {
        stop("error: genes set don't match existing genes set")
      }
      if ("id" %in% colnames(features)) {
        print("id column found in infos")
        if (length(intersect(cells, as.vector(features$id))) == 0) {
          stop("error: id's in infos don't match cells name in counts")
        }
      } else {
        if (nrow(features) != nrow(counts)) {
          stop(
            paste0("error: number of cells (",
              nrow(private$features),
              ") differ between infos and counts (",
              nrow(private$counts)))
        }
      }
      rr_num <- private$get_common_cells(features, counts)
      features <- features[rr_num$features, ]
      counts <- counts[rr_num$counts, ]
      genes <- colnames(counts)
      cells <- rownames(counts)
      not_here <- !(cells %in% private$cells)
      private$cells <- c(private$cells, cells[not_here])
      private$counts <- rbind(private$counts, counts[not_here, ])
      private$features <- rbind(private$features, features[not_here, ])
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
    addfeature = function(feature) {
      private$features[[feature]] <- NA
    },
    setfeature = function(feature, value) {
      private$features[[feature]] <- value
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
    },
    copy = function(
      cells = NULL, genes = NULL, features = NULL, b_cells = NULL){
      if (is.null(b_cells)){
        b_select <- T
      }
      return(
        scdata$new(
          infos = self$getfeaturesw(
            cells = cells, features = features)[b_cells, ],
          counts = self$getcountsw(
            cells = cells, genes = genesm)[b_cells, ]
        )
      )
    },
    select = function(
      cells = NULL, genes = NULL, features = NULL, b_cells = NULL){
      if (is.null(b_cells)){
        b_select <- T
      }
      return(
        scdata$new(
          infos = self$getfeaturesw(
            cells = cells, features = features)[b_cells, ],
          counts = self$getcountsw(
            cells = cells, genes = genesm)[b_cells, ]
        )
      )
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
    },
    getncells = function() {
      return(length(private$cells))
    },
    getngenes = function() {
      return(length(private$genes))
    }
  )
)

#' load scRNASeq data
#'
#' @param infos a text tabular information file on cells with cells in rows and
#' features in columns
#' @param counts a text tabular count file for cells in rows and genes in
#' columns. Can also be a folder in which case every files in the subtree will
#' be loaded'
#' @param regexp a regular expression marching the file we want to load in case
#' of counts is a folder
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
    print(paste0("loading ", infos))
    print(paste0("loading ", counts))
    return(scdata$new(
        infos = utils::read.table(infos, fill = T, h = T, ...),
        counts = utils::read.table(counts, h = T, ...)
      )
    )
  }
}

load_multiple_file <- function(infos, counts, regexp, ...) {
  print(paste0("loading ", infos))
  features <- utils::read.table(infos, fill = T, h = T, ...)
  files_list <- get_files(path = counts, regexp = regexp)
  print(paste0("loading ", files_list[1]))
  scd <- scdata$new(
      infos = features,
      counts = utils::read.table(files_list[1], h = T, ...)
    )
  for (counts_file in files_list[-1]) {
    print(paste0("loading ", counts_file))
    scd$add(
        infos = features,
        counts = utils::read.table(counts_file, h = T, ...)
      )
  }
  return(scd)
}

#' random sample of scRNASeq data
#'
#' @param data a scdata object
#' @param number number of cells to sample
#' @param replace should sampling be with replacement?
#' @return a scdata object with number of cells
#' @examples
#' \dontrun{
#' rdata <- random_sample_data(data, 150)
#' }
#' @export load_data
random_sample_data <- function(data, number, replace = FALSE) {
  number <- min(number, data$getncells)
  cells <- sample(data$getcells, number, replace = FALSE)
  rdata <- scdata$new(
      infos = data$getfeaturesw(cells = cells),
      counts = data$getcountsw(cells = cells)
    )
  return(rdata)
}
