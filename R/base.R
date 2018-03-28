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

#' compute pca loading
#'
#' @param x vector data.frame to vectorize
#' @param cells (default=FALSE) shall the cells loading be returned
#' @return return a matrix of cells (or genes) coordinates or a list including
#' this matrix
#' @examples
#' \dontrun{
#' x = pca_loading(scd)
#' }
#' @importFrom ade4 dudi.pca
#' @export pca_loading
pca_loading <- function(scd, cells = FALSE, ncomp = 5){
  pca_out <- ade4::dudi.pca(
    ascb(scd$getcounts, to_zero = TRUE),
    scan = F,
    nf = ncomp
  )
  if (cells) {
    return(pca_out$l1)
  }
  return(pca_out$c1)
}

#' compute pCMF loading
#'
#' @param x vector data.frame to vectorize
#' @param cells (default=FALSE) shall the cells loading be returned
#' @return return a matrix of cells (or genes) coordinates or a list including
#' this matrix
#' @examples
#' \dontrun{
#' x = pCMF_loading(scd)
#' }
#' @importFrom ade4 dudi.pca
#' @export pca_loading
pCMF_loading <- function(scd, cells = FALSE, ncomp = 5, cpus = 4, tmp_file){
  if (!missing(tmp_file) & file.exists(tmp_file)) {
    print("tmp file found skipping pCMF...")
    load(tmp_file)
  } else {
    pCMF_out <- pCMF(
      X = scd$getcounts,
      K = ncomp,
      iterMax = 500,
      iterMin = 100,
      epsilon = 1e-3,
      verbose = TRUE,
      sparse = TRUE,
      ZI = TRUE,
      ncores = cpus
    )
    if (!missing(tmp_file)){
      save(pCMF_out, file = tmp_file)
    }
  }
  if (cells) {
    return(pCMF_out$stats$ElogU)
  }
  return(pCMF_out$stats$ElogV)
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
  return(order(score_by, score))
}

#' return order of genes by a covariate
#' @param score covariate to order on
#' @param scd scdata with genes to order
#' @return return order
#' @examples
#' \dontrun{
#' genes_order = order_by_factor(scd$get_feature("pDEA_cell_type"), scd)
#' }
#' @export order_by_groups
order_by_factor <- function(score, scd, tmp_file, top = FALSE){
  score_cov <- gene_cov(
    scd = scd, 
    score = score,
    sparse = F,
    ncomp = 1,
    tmp_file = tmp_file)
  la_score <- log(abs(score_cov$B))
  if (top) {
    la_score_100 <- seq(
      from = min(la_score),
      to = max(la_score),
      length.out = 100
    )
    la_score_cdf <- ecdf(la_score)
    la_score_cdf <- la_score_cdf(la_score_100)
    dd_score_cdf <- diff(diff(la_score_cdf))
    score_min <- la_score_100[which(dd_score_cdf == min(dd_score_cdf))]
    return(order(la_score)[order(la_score) %in% which(la_score > score_min)])
  } else {
    return(order(score_cov$B))
  }
}

#' return scalled counts for zidata
#' @param scd scdata object to scale
#' @return return scdata object with scaled counts
#' @examples
#' \dontrun{
#' scd_norm = zinorm(scd)
#' }
#' @export zinorm
zinorm <- function(scd, cpus = 4, v = F, file, zi_scale = TRUE, sd_scale = TRUE){
  if(!missing(file) & file.exists(file)){
    print("cache found...")
    load(file)
  } else {
    print("computing weight...")
    genes <- scd$getgenes
    weight <- scRNAtools::get_weights(
      scd = scd,
      genes = genes,
      cpus = cpus,
      v = v
    )
    weighted_counts <- apply(
      scRNAtools::get_genes(scd, genes),
      1,
      FUN = function(x, weight, sd_scale) {
        if (sd_scale) {
          return( x / as.numeric(weight$gene_scale) )
        } else {
          return( x )
        }
      },
      weight,
      sd_scale
    )
    weighted_counts <- apply(
      t(weighted_counts),
      1,
      FUN = function(x, weight, zi_scale) {
        if (zi_scale) {
          return( x * as.numeric(weight$gene_weight) )
        } else {
          return( x )
        }
      },
      weight,
      zi_scale
    )
    scd_norm <- scdata$new(
      infos = scd$getfeatures,
      counts = t(weighted_counts)
    )
  }
  if(!missing(file) & !file.exists(file)){
    save(scd_norm, file = file)
  }
  return(scd_norm)
}


#' return order base on the rank FUN in groups
#' @param scd an scdata object
#' @param by factor to group on
#' @param FUN (default: mean) function on apply on factor in each groups
#' @return return order
#' @examples
#' \dontrun{
#' genes_order = order_TMP(scd$getcounts, scd$get_feature("DEA_cell_type"))
#' }
#' @export order_TMP
order_TMP <- function(scd, by, FUN=mean, top=scd$getngenes){
  data <- scd$getcounts
  by <- factorize(by)
  by_list <- list()
  by_keep <- c()
  top <- top / length(levels(by))
  for(i in levels(by)){
    tmp <- data[by==i,]
    by_list[[i]] <- apply(tmp, 2, FUN=function(x){
      FUN(x)
    })
    by_keep <- unique(c(by_keep,
                        order(by_list[[i]], decreasing=TRUE)[1:top]))
  }
  by_dist <- by_list[[1]]
  for(i in 2:length(by_list)){
    by_dist <- by_dist - by_list[[i]]
  }
  return(order(by_dist)[by_keep])
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
  scd, b_cells = NULL, cells = NULL, genes = NULL, by,
  FUN = function(x){mean(x) / sd(x)},
  top = ncol(data), min_cell = 5) {
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

#' return top and bottom of list for order base on the rank FUN in groups
#' @param scd an scdata object
#' @param by factor to group on
#' @param top number of item to return
#' @param FUN (default: mean) function on apply on factor in each groups
#' @return return indice of top and bottom item
#' @examples
#' \dontrun{
#' genes_top = top_by_groups(scd$getcounts, scd$get_feature("DEA_cell_type"), 100)
#' }
#' @export order_2_groups
top_2_groups <- function(
  scd, b_cells = NULL, top = 100, cells = NULL, genes = NULL, by,
  FUN = function(x){mean(x[x != 0]) * length(x[x != 0])/length(x) / sd(x[x!=0])}
) {
  indices <- order_2_groups(
    scd = scd, b_cells = b_cells, cells = cells, genes = genes,
    by = by, FUN = FUN, top = ncol(scd$getcounts), min_cell = 5
  )
  top_indice <- 1:round(top/2)
  bottom_indice <- (length(indices) - round(top/2)):length(indices)
  gene_indices <- indices[c(top_indice, bottom_indice)]
  return(scd$getgenes[gene_indices])
}
