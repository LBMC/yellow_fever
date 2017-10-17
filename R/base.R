#' Select ERCC from gene list
#'
#' @param data scdata object
#' @param minus (default: false) if true return non ERCC
#' @return return a vector of ERCC gene names
#' @examples
#' \dontrun{
#' genes = ERCC(data, minus = true)
#' data = data$select(genes = genes)
#' }
#' @export ERCC
ERCC <- function(scd, minus = FALSE) {
  if (minus){
    return(scd$getgenes[!grepl("ERCC\\.", scd$getgenes)])
  }
  return(scd$getgenes[grepl("ERCC\\.", scd$getgenes)])
}

get_files <- function(path, regexp) {
  file_list <- base::list.files(
    path = path, full.names = TRUE, recursive = TRUE
  )
  file_list <- file_list[grepl(regexp, file_list, perl = T)]
  return(file_list)
}

#' anscombe transform of counts
#'
#' @param x vector of counrs
#' @return return anscomb transform of x
#' @examples
#' \dontrun{
#' x = anscombe(scd$getgene("gene_a"))
#' }
#' @export anscombe
anscombe <- function(x){
  return(2 * sqrt(x + 3 / 8))
}

#' anscombe transform of counts
#'
#' @param x vector of counrs
#' @param to_zero (bool default: TRUE) set ascb(x) to 0 when x = 0
#' @param na.rm remove NA values
#' @return return anscomb transform of x
#' @examples
#' \dontrun{
#' x = ascb(scd$getgene("gene_a"))
#' }
#' @export ascb
ascb <- function(x, to_zero=TRUE, na.rm=FALSE){
  if (na.rm){
    if (!is.null(dim(x))){
      x <- x[, apply(x, 2, function(x){
        !any(is.na(x))
      })]
    }
  }
  if (to_zero){
    if (min(x, na.rm = T) < 0){
      x[!is.na(x)] <- x[!is.na(x)] - min(x, na.rm = TRUE)
    }
    return(anscombe(x) - anscombe(0))
  }
  return(anscombe(x))
}

#' factorize a vector or data.frame
#'
#' @param x vector data.frame to vectorize
#' @param columns number to vectorize in case of a data.frame
#' @return return x with as.factor(as.vector(x)) applied
#' @examples
#' \dontrun{
#' x = factorize(x)
#' }
#' @export factorize
factorize <- function(x, columns) {
  if (is.null(ncol(x))) {
    return(as.factor(as.vector(x)))
  } else {
    if (missing(columns)) {
      for (i in seq_len(ncol(x))) {
        x[, i] <- as.factor(as.vector(x[,i]))
      }
    } else {
      for (i in columns) {
        x[, i] <- as.factor(as.vector(x[,i]))
      }
    }
    return(x)
  }
}

#' vectorize a vector or data.frame
#'
#' @param x vector data.frame to vectorize
#' @param columns number to vectorize in case of a data.frame
#' @return return x with as.numeric(as.vector(x)) applied
#' @examples
#' \dontrun{
#' x = vectorize(x)
#' }
#' @export vectorize
vectorize <- function(x, columns) {
  if (is.null(ncol(x))){
    return(as.numeric(as.vector(x)))
  } else {
    if( missing(columns)) {
      for (i in seq_len(ncol(x))) {
        x[, i] <- as.numeric(as.vector(x[,i]))
      }
    } else {
      for (i in columns) {
        x[, i] <- as.numeric(as.vector(x[,i]))
      }
    }
    return(x)
  }
}

#' compute bca loading
#'
#' @param x vector data.frame to vectorize
#' @param by a factor defining the groups
#' @param norm_by a factor defining the groups with want to normalize
#' @param ncomp (default=5) number of component to compute
#' @param cells_l (default=FALSE) shall the cells loading be also returned
#' @param loading (default=FALSE) return only the genes contribution to the axes
#' @return return a matrix of cells coordinates or a list including this matrix
#' @examples
#' \dontrun{
#' x = bca_loading(scd, scd$getfeature("sex"))
#' }
#' @import ade4
#' @export bca_loading
bca_loading <- function(scd, by, norm_by, ncomp=5, cells_l=FALSE, loading = F){
  by <- scRNAtools::factorize(by)
  pca_out <- dudi.pca(scd$getcounts,
                      scan = F,
                      nf = ncomp)
  if (!missing(norm_by)){
    norm_by <- factorize(norm_by)
    wca_out <- wca(pca_out,
                   as.factor(as.vector(norm_by)),
                   scan = F,
                   nf = ncomp)
    pca_out <- dudi.pca(wca_out$tab,
                        scan = F,
                        nf = ncomp)
  }
  bca_out <- bca(pca_out,
                 by,
                 scan = F,
                 nf = ncomp)
  if (loading){
    return(bca_out$c1)
  }
  gene_order <- order(abs(bca_out$c1[, 1]), decreasing = TRUE)
  results <- as.data.frame(bca_out$c1[gene_order, ])
  rownames(results) <- scd$getgenes
  names(results) <- names(bca_out$c1)[gene_order]
  if (cells_l){
    return(list(res = results,
                cells_l = bca_out$ls,
                c1 = bca_out$c1))
  }
  return(results)
}
