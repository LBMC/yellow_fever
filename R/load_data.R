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
    get_rc_counts = function(cells = NULL, genes = NULL, v = F) {
      c_num <- 1:self$getngenes
      if (!is.null(genes)) {
        c_num <- which(private$genes %in% genes)
      }
      r_num <- 1:self$getncells
      if (!is.null(cells)) {
        r_num <- which(private$cells %in% cells)
      }
      return(list(c_num = c_num, r_num = r_num))
    },
    # private method to get position of cells and features in rows and columns
    get_rc_features = function(cells = NULL, features = NULL, v = F) {
      c_num <- 1:self$getnfeatures
      if (!is.null(features)) {
        c_num <- which(colnames(private$features) %in% features)
      }
      r_num <- 1:self$getncells
      if (!is.null(cells)) {
        r_num <- which(private$cells %in% cells)
      }
      return(list(c_num = c_num, r_num = r_num))
    },
    # private method to get position of cells in common between features and
    # counts
    get_common_cells = function(features, counts, v = F) {
      id_features <- as.vector(features$id)
      id_counts <- rownames(counts)
      common_cells <- intersect(id_features, id_counts)
      if (v) {
        print(
          paste0(length(common_cells),
            " cells in common between infos and counts"))
        print("cells present in infos, but not counts")
        print(setdiff(id_features, id_counts))
        print("cells present in counts, but not infos")
        print(setdiff(id_counts, id_features))
      }
      r_features <- which(id_features %in% common_cells & !duplicated(id_features))
      r_counts <- which(id_counts %in% common_cells & !duplicated(id_counts))
      return(list(features = r_features, counts = r_counts))
    },
    order_by_cells = function(v = F) {
      if (nrow(private$counts) == length(private$cells)) {
        private$counts <- private$counts[order(private$cells), ]
      } else {
          stop(paste0(
            "error: number of cells (",
            length(private$cells),
            ") don't match number of counts rows (",
            nrow(private$counts),
          ")"))
      }
      if (nrow(private$features) == length(private$cells)) {
        private$features <- private$features[order(self$getfeature("id")), ]
      } else {
          stop(paste0(
            "error: number of cells (",
            length(private$cells),
            ") don't match number of feature rows (",
            nrow(private$features),
          ")"))
      }
      private$cells <- rownames(private$counts)
    },
    transpose_count = function(counts, v = F) {
      counts <- as.matrix(counts)
      if (nrow(private$features) < nrow(counts)) {
        if (v) {
          print("transposing counts...")
        }
        return(t(counts))
      }
      return(counts)
    },
    set_na_to_zero = function(v = F) {
      if (v) {
        print("cells with NA's counts:")
        print(self$getcells[rowSums(is.na(self$getcounts)) != 0])
        print("genes with NA's counts:")
        print(self$getgenes[colSums(is.na(self$getcounts)) != 0])
      }
      private$counts[is.na(private$counts)] <- 0
    },
    display_dim = function(infos, counts, v = F) {
      if (v) {
        print(
          paste0("dim of infos :",
            nrow(infos),
            " x ",
            ncol(infos)
          )
        )
        print(
          paste0("dim of counts :",
            nrow(counts),
            " x ",
            ncol(counts)
          )
        )
      }
    }
    ),
  public = list(
    initialize = function(infos = NA, counts = NA, v = F) {
      private$features <- as.data.frame(infos)
      private$counts <- private$transpose_count(counts, v = v)
      private$genes <- colnames(private$counts)
      private$cells <- rownames(private$counts)
      private$display_dim(private$features, private$counts, v = v)
      if ("id" %in% colnames(private$features)) {
        if (v) {
          print("id column found in infos")
        }
        if (length(intersect(private$cells, as.vector(private$features$id)))
          == 0) {
          print(setdiff(private$cells, as.vector(private$features$id)))
          stop("error: id's in infos don't match cells name in counts")
        }
        rr_num <- private$get_common_cells(
          private$features,
          private$counts,
          v = v
        )
        private$features <- private$features[rr_num$features, ]
        private$counts <- private$counts[rr_num$counts, ]
        private$cells <- rownames(private$counts)
        private$display_dim(private$features, private$counts, v = v)
        private$order_by_cells(v = v)
        if (any(private$getfeature['id'] != private$cells) ){
          stop("error : features order don't match counts order")
        }
      } else {
        if (nrow(private$features) != nrow(private$counts)) {
          stop(
            paste0("error: number of cells (",
              nrow(private$features),
              ") differ between infos and counts (",
              nrow(private$counts)))
        }
      }
      private$set_na_to_zero(v = v)
      self$summary(v = v)
    },
    add = function(infos = NA, counts = NA, v = F) {
      features <- as.data.frame(infos)
      counts <- private$transpose_count(counts, v = v)
      genes <- colnames(counts)
      cells <- rownames(counts)
      private$display_dim(features, counts, v = v)
      if (length(intersect(private$cells, cells)) != 0) {
        stop("error: trying to add cells already present")
      }
      if (length(intersect(private$genes, genes)) != length(private$genes)) {
        print(private$cells)
        print(as.vector(private$features$id))
        stop("error: genes set don't match existing genes set")
      }
      if ("id" %in% colnames(features)) {
        if (v) {
          print("id column found in infos")
        }
        if (length(intersect(cells, as.vector(features$id))) == 0) {
          print(setdiff(private$cells, as.vector(private$features$id)))
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
      rr_num <- private$get_common_cells(features, counts, v = v)
      features <- features[rr_num$features, ]
      counts <- counts[rr_num$counts, ]
      genes <- colnames(counts)
      cells <- rownames(counts)
      not_here <- !(cells %in% private$cells)
      private$counts <- rbind(private$counts, counts[not_here, ])
      private$features <- rbind(private$features, features[not_here, ])
      private$cells <- rownames(private$counts)
      private$order_by_cells(v = v)
      private$set_na_to_zero(v = v)
      if (any(private$getfeature['id'] != private$cells) ){
        stop("error : features order don't match counts order")
      }
    },
    summary = function(v = F) {
      if (v) {
        cat(paste0("ncol: ", ncol(private$counts), ".\n"))
        cat(paste0("nrow: ", nrow(private$counts), ".\n"))
        cat(paste0("features: ", ncol(private$features), ".\n"))
      }
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
        b_cells <- T
      }
      if(is.null(cells) | length(cells) > 1){
        return(
          scdata$new(
            infos = self$getfeaturesw(
              cells = cells, features = features)[b_cells, ],
            counts = self$getcountsw(
              cells = cells, genes = genes)[b_cells, ]
          )
        )
      } else {
        return(
          scdata$new(
            infos = self$getfeaturesw(
              cells = cells, features = features),
            counts = self$getcountsw(
              cells = cells, genes = genes)
          )
        )
      }
    },
    select = function(
      cells = NULL, genes = NULL, features = NULL, b_cells = NULL) {
      return(
        self$copy(
          cells = cells,
          genes =  genes,
          features = features,
          b_cells = b_cells
      ))
    },
    order = function(cells = NULL, genes = NULL) {
      if (!is.null(cells)) {
        private$features <- private$features[cells, ]
        private$counts <- private$counts[cells, ]

      }
      if (!is.null(genes)) {
        private$counts <- private$counts[, genes]
      }
      private$genes <- colnames(private$counts)
      private$cells <- rownames(private$counts)
    },
    transform = function(FUN = function(x){x}) {
      private$counts <- FUN(private$counts)
    },
    update = function(counts = NULL, features = NULL) {
      if (!is.null(counts)) {
        private$counts <- counts
        private$genes <- colnames(private$counts)
        private$cells <- rownames(private$counts)
      }
      if (!is.null(features)) {
        private$features <- features
      }
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
    },
    getnfeatures = function() {
      return(ncol(private$features))
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
load_data <- function(
  infos, counts, regexp = ".*", infos_sep = "", counts_sep = "", ...) {
  if (base::file.info(counts)$isdir) {
    return(load_multiple_file(
      infos = infos,
      counts = counts,
      regexp = regexp,
      infos_sep = infos_sep,
      counts_sep = counts_sep,
      ...
      )
    )
  } else {
    print(paste0("loading ", infos))

    print(paste0("loading ", counts))
    return(scdata$new(
        infos = utils::read.table(infos, fill = T, h = T, sep = infos_sep, ...),
        counts = utils::read.table(counts, h = T, sep = counts_sep, ...)
      )
    )
  }
}

load_multiple_file <- function(
  infos, counts, regexp, infos_sep = "", counts_sep = "", ...) {
  print(paste0("loading ", infos))
  features <- utils::read.table(infos, fill = T, h = T, sep = infos_sep, ...)
  files_list <- get_files(path = counts, regexp = regexp)
  print(paste0("loading ", files_list[1]))
  scd <- scdata$new(
      infos = features,
      counts = utils::read.table(files_list[1], h = T, sep = counts_sep, ...)
    )
  for (counts_file in files_list[-1]) {
    print(paste0("loading ", counts_file))
    scd$add(
        infos = features,
        counts = utils::read.table(counts_file, h = T, sep = counts_sep, ...)
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

#' load_data_salmon scRNASeq data
#'
#' @param infos a tabular csv file that contains information on the cells
#' @param counts path to kallisto quant.sf output
#' @param id_regexp a regular expression caputing the id of the cells in the
#' salmon folder output paths
#' @param id_regexp_b a regular expression match the path of the quant.sf files
#' @param ERCC_regexp a regular expression caputing genes names of the ERCC
#' @param feature the feature to extract from salmon output
#' @param ... additional argument to be passed to read.table function
#' @return scdata object
#' @examples
#' \dontrun{
#' data <- load_data_salmon('data/infos.csv', 'data/salmon_output')
#' }
#' @import tximport readr
#' @export load_data_salmon
load_data_salmon <- function(
  infos, counts,
  id_regexp = ".*_(P[0-9]{4}_[0-9]{1,4})_.*",
  id_regexp_b = "P[0-9]{4}_[0-9]{1,4}",
  ERCC_regexp = "^ERCC.*",
  feature = "counts",
  infos_sep = ",",
  tximport_obj = "",
  ...) {

  if(missing(tximport_obj) |
    (!missing(tximport_obj) &
    !file.exists(paste0(tximport_obj, ".Rdata")))) {

    print("loading quant.sf files...")
    dir_list <- list.dirs(counts)
    dir_list <- paste0(dir_list, "/quant.sf")
    dir_list <- dir_list[file.exists(dir_list)]
    names(dir_list) <- gsub(id_regexp, "\\1", dir_list, perl = T)
    dir_list <- dir_list[grepl(id_regexp_b, names(dir_list))]
    transcript_to_genes <- read.table(dir_list[1], h = T)[1]
    scd_paired_name <- data.frame(
      transcript = as.vector(transcript_to_genes$Name),
      genes = unlist(
        sapply(
          as.vector(transcript_to_genes$Name),
          FUN = function(x){
            if (grepl(ERCC_regexp, x, perl = TRUE)){
              return(x)
            } else {
              return(unlist(strsplit(x, "|", fixed = T))[6])
            }
          }
        )
      )
    )
    scd <- tximport(
      dir_list,
      type = "salmon",
      tx2gene = scd_paired_name,
      countsFromAbundance = "lengthScaledTPM"
    )
    print(names(scd))
    if (!missing(tximport_obj)) {
      save(scd,
        file = paste0(tximport_obj, ".Rdata"))
    }
  } else {
    load(file = paste0(tximport_obj, ".Rdata"))
  }
  counts <- scd[[feature]]
  print(dim(counts))
  if (any(duplicated(rownames(counts)))) {
    print("Warning: duplicated cell id's, merging them...")
    count_list_b <- by(
      data = counts,
      INDICES = as.factor(rownames(counts)),
      FUN = colMeans
    )
    counts <- data.frame(matrix(
      unlist(count_list_b),
      nrow = length(count_list_b),
      byrow = T)
    )
    rownames(counts) <- names(count_list_b)
    colnames(counts) <- names(count_list_b[[1]])
  }
  print(dim(counts))
  if (any(duplicated(rownames(counts)))) {
    stop("error: duplicated genes id's")
  }
  infos_table <- utils::read.table(
    infos, fill = T, h = T, sep = infos_sep, ...)
  return(scdata$new(
    infos = infos_table,
    counts = counts
  ))
}
