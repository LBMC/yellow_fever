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
    return(scd$getgenes[!grepl("ERCC.*", scd$getgenes, perl = TRUE)])
  }
  return(scd$getgenes[grepl("ERCC.*", scd$getgenes, perl = TRUE)])
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

#' inverse anscombe transform of counts
#'
#' @param x vector of anscomb(counrs)
#' @return return anscomb transform of x
#' @examples
#' \dontrun{
#' x = anscombe_inv(anscombe(scd$getgene("gene_a")))
#' }
#' @export anscombe_inv
anscombe_inv <- function(x){
  return((x / 2)^2 - 3 / 8)
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
  } else {
    if (is.null(dim(x))) {
      x <- ifelse(is.na(x), 0, x)
    } else {
      x <- apply(x, 2, FUN = function(y) {
        ifelse(is.na(y), 0, y)
      })
    }
  }
  if (is.null(dim(x))) {
    x <- ifelse(x < 0, 0, x)
  } else {
    x <- apply(x, 2, FUN = function(y) {
      ifelse(y < 0, 0, y)
    })
  }
  if (to_zero){
    if (min(x, na.rm = T) < 0){
      x[!is.na(x)] <- x[!is.na(x)] - min(x, na.rm = TRUE)
    }
    return(anscombe(x) - anscombe(0))
  }
  return(anscombe(x))
}

#' inverse anscombe transform of counts
#'
#' @param x vector of ascb(counrs)
#' @param to_zero (bool default: TRUE) set ascb(x) to 0 when x = 0
#' @param na.rm remove NA values
#' @return return anscomb transform of x
#' @examples
#' \dontrun{
#' x = ascb_inv(scd$getgene("gene_a"))
#' }
#' @export ascb_inv
ascb_inv <- function(x, to_zero=TRUE, na.rm=FALSE){
  if (to_zero){
    return(anscombe_inv(x + anscombe(0)))
  }
  return(anscombe_inv(x))
}
anscombe_inv(anscombe(3))
ascb_inv(ascb(3))
ascb_inv(ascb(3, to_zero = F), to_zero = F)

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
#' @importFrom ade4 dudi.pca wca bca
#' @export bca_loading
bca_loading <- function(scd, by, norm_by, ncomp=5, cells_l=FALSE, loading = F){
  by <- scRNAtools::factorize(by)
  pca_out <- ade4::dudi.pca(scd$getcounts,
                      scan = F,
                      nf = ncomp)
  if (!missing(norm_by)){
    norm_by <- factorize(norm_by)
    wca_out <- ade4::wca(pca_out,
                   as.factor(as.vector(norm_by)),
                   scan = F,
                   nf = ncomp)
    pca_out <- ade4::dudi.pca(wca_out$tab,
                        scan = F,
                        nf = ncomp)
  }
  bca_out <- ade4::bca(pca_out,
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

#' compute wca normalized pca
#' @param x vector data.frame to vectorize
#' @param by a factor defining the groups
#' @param keep_zero (default=FALSE) should the zero values be keept
#' @return return cells coordinates normalized for 'by' variance
#' @examples
#' \dontrun{
#' x = wca_loading(scd, scd$getfeature("sex"))
#' }
#' @importFrom ade4 dudi.pca wca
#' @export wca_norm
wca_norm <- function(scd, by, ncomp=2, keep_zero=FALSE){
  pca_out <- ade4::dudi.pca(scd$getcounts,
                      scan = F,
                      nf = scd$getngenes - 1)
  if (is.null(ncol(by))){
    wca_out <- ade4::wca(pca_out,
                   factorize(by),
                   scan = F,
                   nf = ncomp)
  } else {
    for (i in seq_len(ncol(by))){
      wca_out <- ade4::wca(pca_out,
                   factorize(by[, i]),
                   scan = F,
                   nf = ncomp)
    }
  }
  if (keep_zero){
    wca_out$tab[scd$getcounts == 0] <- 0
  }
  return(wca_out$tab)
}

#' compute TPM
#' @param counts counts vector
#' @param len genes length vector
#' @return return TPM
#' @examples
#' \dontrun{
#' x = TPM_vector(scd$getcounts[1, ], len)
#' }
#' @export TPM_vector
TPM_vector <- function(counts, len) {
  if (!is.null(ncol(counts))){
    stop("error: TPM_vector(counts, len), counts is not a vector")
  }
  if (!is.null(ncol(len))){
    stop("error: TPM_vector(counts, len), len is not a vector")
  }
  denum <- sum( counts / len )
  return(
    apply(
      data.frame(
        counts = counts,
        len = len
      ),
      MARGIN = 1,
      FUN = function(x, denum){
        if (x[1] == 0 | x[2] == 0) {
          return(0)
        } else {
          return( 10^6 * ( ( x[1] / x[2] ) / denum ) )
        }
      },
      denum = denum
    )
  )
}

#' compute TPM
#' @param counts matrix
#' @param len genes length vector
#' @return return TPM
#' @examples
#' \dontrun{
#' x = TPM_matrix(scd$getcounts, len)
#' }
#' @export TPM_matrix
TPM_matrix <- function(counts, len) {
  if (is.null(ncol(counts))){
    stop("error: TPM_matrix(counts, len), counts is not a matrix")
  }
  if (!is.null(ncol(len))){
    len = len[1 ,]
  }
  return(
    apply(
      counts,
      MARGIN = 1,
      FUN = function(x, len){
        TPM_vector(x, len)
      },
      len = len
    )
  )
}

#' compute TPM
#' @param counts matrix or vector
#' @param len genes length vector
#' @return return TPM
#' @examples
#' \dontrun{
#' x = TPM(scd$getcounts, len)
#' }
#' @export TPM_matrix
TPM <- function(counts, len) {
  if (is.null(ncol(counts))) {
    return(TPM_vector(counts, len))
  } else {
    return(TPM_matrix(counts, len))
  }
}

#' return genes expressed with the following criteria
#' @param scd scdata object
#' @param zi_threshold persentage of cells == 0
#' @return return list of genes names
#' @examples
#' \dontrun{
#' expressed_genes = expressed(scd, 0.90)
#' }
#' @export expressed
expressed <- function(scd, zi_threshold = 0.90){
  zi_rate <- colSums(round(scd$getcounts) == 0) / scd$getncells
  return(scd$getgenes[zi_rate <= zi_threshold])
}

#' return order base on the rank FUN in groups
#' @param score score to order on
#' @param by factor to group on
#' @param FUN (default: mean) function on apply on factor in each groups
#' @return return order
#' @examples
#' \dontrun{
#' cells_order = order_by_groups(scd$get_feature("pDEA_cell_type"), scd$get_feature("DEA_cell_type"))
#' }
#' @export order_by_groups
order_by_groups <- function(score, by, FUN = mean){
  score_av <- aggregate(score, by = list(as.factor(by)), FUN = FUN)
  names(score_av) <- c("by", "mean")
  score_av_list <- as.list(score_av$mean)
  names(score_av_list) <- score_av$by
  score_by <- unlist(score_av_list[as.vector(by)])
  return(order(score_by))
}

#' return scalled counts for zidata
#' @param scd scdata object to scale
#' @return return scdata object with scaled counts
#' @examples
#' \dontrun{
#' scd_norm = zinorm(scd)
#' }
#' @export zinorm
zinorm <- function(scd, cpu = 4, v = F){
  weight <- scRNAtools::get_weights(
    scd = scd,
    genes = scd$getgenes,
    cpus = cpus,
    v = v
  )
  weighted_counts <- apply(
    scRNAtools::get_genes(scd, genes),
    1,
    FUN = function(x, weight) {
      x / as.numeric(weight$gene_scale)
    },
    weight
  )
  weighted_counts <- apply(
    t(weighted_counts),
    1,
    FUN = function(x, weight) {
      x * as.numeric(weight$gene_weight)
    },
    weight
  )
  return(scdata$new(
    infos = scd$getfeatures,
    counts = t(weighted_counts)
  ))
  return(order(score_by))
}

#' return order base on the rank FUN in groups
#' @param scd an scdata object
#' @param by factor to group on
#' @param FUN (default: mean) function on apply on factor in each groups
#' @return return order
#' @examples
#' \dontrun{
#' genes_order = order_by_groups(scd$getcounts, scd$get_feature("DEA_cell_type"))
#' }
#' @export order_2_groups
order_2_groups <- function(
  scd, b_cells = NULL, cells = NULL, genes = NULL,
  by, FUN = function(x){mean(x)}, top = ncol(data), min_cell = 5) {
  data <- scd$select(
    genes = genes,
    cells = cells,
    b_cells = b_cells)$getcounts
  if(top > ncol(data)){
    top <- ncol(data)
  }
  by <- factorize(by)
  by_dist <- list()
  for(i in levels(by)){
    r_select <- by %in% i
    by_dist[[i]] <- apply(data[r_select,], 2, FUN=function(x, FUNx){
        FUNx(x)
      }, FUNx=FUN)
  }
  by_order <- as.vector(by_dist[[1]])
  by_order_abs <- as.vector(by_dist[[1]])
  for(i in 2:length(by_dist)){
    by_order <- by_order - as.vector(by_dist[[i]])
    by_order_abs <- abs(by_order_abs - as.vector(by_dist[[i]]))
  }
  results <- which(by_order_abs %in% by_order_abs[order(by_order_abs,
                                                        decreasing=TRUE)][1:top])
  results <- results[order(by_order[results], decreasing=TRUE)]
  if(top < ncol(data)){
    results <- results[colSums(data[, results] >= 5) >= 5]
    return(results[1:top])
  }else{
    return(results)
  }
}
