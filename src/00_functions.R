needed_packages <- c(
  "optparse",
  "tidyverse",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "sctransform",
  "broom",
  "broom.mixed",
  "DHARMa",
  "glmmTMB",
  "lmtest",
  "parallel",
  "pbmcapply",
  "plsgenomics",
  "scater",
  "plyr"
)

install_if_needed <- function(packages){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
  sapply(
    packages,
    function(x) {
      if (!requireNamespace(x, quietly = TRUE))
        tryCatch({
          BiocManager::install(x, update = F, ask = F)
        }, error = function(e) {
          devtools::install_github(x, upgrade = "never")
        })
      if (!stringr::str_detect(x, "/")){
        library(x, character.only = TRUE)
      } else {
        library(stringr::str_replace(x, ".*/(.*)", "\\1"),
                character.only = TRUE)
      }
    }
  )
  return(str_c(str_c(packages, collapse = ", "), " installed and loaded."))
}
install_if_needed(needed_packages)

#' compute anscombe transform
#' compute the anscombe transform of x (sqrt(x + 3/8))
#' @param x a vector of numeric
#' @return x a vector of numeric
#' @examples
#' \dontrun {
#' anscombe(1:10)
#' }
#' @export anscombe
anscombe <- function(x){
  2.0 * sqrt(x + 3.0 / 8.0)
}

#' compute anscombe transform shifted to zero
#' compute the anscombe transform of x (anscombe(x) + anscombe(0)
#' @param x a vector of numeric
#' @return x a vector of numeric
#' @examples
#' \dontrun {
#' anscombe(1:10)
#' }
#' @export ascb
ascb <- function(x){
  anscombe(x) - anscombe(0.0)
}

logspace <- function(d1, d2, n) {
  return(exp(log(10) * seq(d1, d2, length.out = n)))
}

cell_type_palette <- function(cell_type){
  cell_type_color <- list(EM = "#c90000",
                          TEMRA = "#f4a582",
                          TSCM = "#92c5de",
                          CM = "#0571b0",
                          Naive = "gray",
                          NAIVE = "gray",
                          Effector = "#c90000",
                          EFF = "#c90000",
                          Memory = "#0571b0",
                          MEM = "#0571b0",
                          UNK = "gray")
  cell_types_color <- lapply(as.list(cell_type),
  FUN = function(x, cell_type_color){
      cell_type_color[[x]]
    }
    , cell_type_color)
  cell_types_color <- unlist(cell_types_color)
  names(cell_types_color) <- cell_type
  return(cell_types_color)
}

#' compute anscombe transform shifted to zero
#' compute the anscombe transform of x (anscombe(x) + anscombe(0)
#' @param x a vector of numeric
#' @return x a vector of numeric
#' @examples
#' \dontrun {
#' anscombe(1:10)
#' }
#' @export ascb
QC_sample <- function(is_good) {
  c(
    base::sample(
      which(is_good),
      size = sum(!is_good),
      replace = T
    ),
    which(!is_good)
  )
}

#' compute anscombe transform shifted to zero
#' compute the anscombe transform of x (anscombe(x) + anscombe(0)
#' @param x a vector of numeric
#' @return x a vector of numeric
#' @examples
#' \dontrun {
#' anscombe(1:10)
#' }
#' @export ascb
QC_fit <- function(sce, assay_name = "logcounts", ncell_name = "cell_number") {
  is_good <- colData(sce)[[ncell_name]] > 0
  select_sample <- QC_sample(is_good)
  suppressWarnings(
    e1071::svm(
      x = t(assays(sce)[[assay_name]][, select_sample]),
      y = as.factor(is_good[select_sample]),
      kernel = "linear",
      class.weights = 100 / table(is_good),
      type = "C-classification",
      scale = TRUE
  ))
}
#' compute anscombe transform shifted to zero
#' compute the anscombe transform of x (anscombe(x) + anscombe(0)
#' @param x a vector of numeric
#' @return x a vector of numeric
#' @examples
#' \dontrun {
#' anscombe(1:10)
#' }
#' @export ascb
QC_predict <- function(fit, sce, assay_name = "logcounts") {
  stats::predict(
    fit,
    t(assays(sce)[[assay_name]]),
    na.action = na.fail
  ) %>%
    as.numeric()
}
QC_score <- function(sce, assay_name = "logcounts", ncell_name = "cell_number",
                     run = ncol(sce), ncpus = 6){
  pbmcapply::pbmclapply(
    as.list(1:run),
    FUN = function(x, sce, assay_name, ncell_name){
      message(x)
      QC_predict(
        fit = QC_fit(
          sce = sce,
          assay_name = assay_name,
          ncell_name = ncell_name
        ),
        sce = sce,
        assay_name = assay_name)
    }, sce = sce, assay_name = assay_name, ncell_name = ncell_name,
    mc.cores = ncpus) %>% 
    do.call(cbind, .) %>% 
    rowMeans(. - 1)
}

test_zi <- function(model_nozi, model_zi, threshold = 0.05){
  test <- stats::anova(model_nozi, model_zi)
  if (!is.na(test["model_zi", "Pr(>Chisq)"])) {
    if (test["model_zi", "Pr(>Chisq)"] <= threshold) {
      return(model_zi)
    }
  }
  return(model_nozi)
}

fit_zi_nb <- function(data, formula, zi_formula = formula, threshold = 0.05){
  test_zi( 
    model_nozi = glmmTMB::glmmTMB(
      as.formula(formula),
      family = glmmTMB::nbinom2,
      data = data
    ),
    model_zi = glmmTMB::glmmTMB(
      as.formula(formula),
      family = glmmTMB::nbinom2,
      ziformula = as.formula(zi_formula),
      data = data
    ),
    threshold = threshold
  )
}

residuals_zi_nb <- function(model) {
  tryCatch({
      sim_resid <- DHARMa::simulateResiduals(model, n = 100)
      tibble(uniformity = DHARMa::testUniformity(
               simulationOutput = sim_resid, plot = F)$p.value,
             dispersion = DHARMa::testDispersion(
               simulationOutput = sim_resid, plot = F)$p.value,
             outlier = DHARMa::testOutliers(
               simulationOutput = sim_resid, plot = F)$p.value,
             zi = !is.null(model$modelInfo$terms$zi)
      )
    },
    error = function(e){ NULL }
  )
}

drop_term <- function(formula, test) {
  formula %>%  
    str_remove(fixed(str_remove(test, "\\s*~\\s*"))) %>% 
    str_replace("(~)\\s*[+-:*]", "\\1") %>% 
    str_replace("[:]\\s*([+-])", "\\1") %>% 
    str_replace("([+-])\\s*[:]", "\\1") %>% 
    str_replace("([+-:])\\s*[+-:]\\s*", "\\1") %>% 
    str_remove("\\s*[+-:*]\\s*$")
}

LRT_zi_nb <- function(data, test, formula, zi_formula = formula, threshold = 0.05){
  model <- fit_zi_nb(
    data = data,
    formula = formula,
    zi_formula = zi_formula,
    threshold = threshold
  )
  model0 <- fit_zi_nb(
    data = data,
    formula = drop_term(formula, test),
    zi_formula = drop_term(zi_formula, test),
    threshold = threshold
  )
  list(
    parameters = model %>%
      broom.mixed::tidy(),
    LRT = anova(model0, model) %>%
      broom.mixed::tidy(),
    residuals = model %>%
      residuals_zi_nb()
  )
}

DEA_data <- function(x, sce, assay_name, formula, zi_formula){
  colData(sce) %>%
    as_tibble() %>%
    mutate(
      count = as.vector(assay(sce, assay_name)[x, ]),
      gene_name = x
    ) %>% 
    dplyr::select(
      any_of(
        str_split(
          str_c(formula, zi_formula),
          "[ 1~+()|]"
        )[[1]]
      ),
      gene_name
    ) %>%
    drop_na()
}

DEA <- function(sce, 
                test, formula, zi_formula = formula, zi_threshold = 0.05,
                assay_name = "counts",
                cpus = 4){
  results <- pbmcapply::pbmclapply(
    rownames(sce),
    FUN = function(x, sce, assay_name, test, formula, zi_formula, zi_threshold){
      LRT_zi_nb(
        data = DEA_data(
          x = x,
          sce = sce,
          assay_name = assay_name,
          formula = formula,
          zi_formula = zi_formula
          ),
        test = test,
        formula = formula,
        zi_formula = zi_formula,
        threshold = zi_threshold
      )
    },
    sce = sce,
    assay_name = assay_name,
    test = test,
    formula = formula,
    zi_formula = zi_formula,
    zi_threshold = zi_threshold,
    mc.preschedule = F,
    mc.cores = cpus,
    ignore.interactive = T
  )
  names(results) <- rownames(sce)
  return(results)
}

get_genes_pval <- function(results, sce) {
  lapply(results,
       FUN = function(x){
         if ("LRT" %in% names(x)) {
           return(x$LRT$p.value[2])
         } else {
           return(NA)
         }
       }) %>%
    do.call(c, .) %>% 
    tibble(
      id = names(.),
      pval = .
    ) %>% 
    right_join(tibble(id = rownames(sce))) %>%
    pull(pval)
}

#' scaling for for scRNASeq data
#' Helper function that scale the genes expression matrix for the PLS function
#' scale a scRNA genes counts vector by dividing by the estimate of a glmmTMB 
#' nbinom2 variance estimate and if a zero inflation is detected by multipling
#' by the proportion of counts outside the dirac in zero
#' @param data a data.frame with a *count* column
#' @return a data.frame with a *count* column
#' @examples
#' \dontrun {
#'scaled_counts <- scale_zi_nb(
#'  data
#')
#' }
#' @import tidyverse boot
scale_zi_nb <- function(data) {
  model <- fit_zi_nb(
      data = data,
      formula = "count ~ 1"
    )
  zi_prop <- ifelse(
    "betazi" %in% names(model$fit$par),
    boot::inv.logit(model$fit$par["betazi"]),
    0.0)
  scaling <-
    exp(model$fit$par["beta"]) + exp(model$fit$par["beta"])^2 / sigma(model)
  scaling <- ifelse(
    scaling != 0.0,
    scaling,
    1.0)
  data$count <- (data$count / scaling) * (1.0 - zi_prop)
  return(data)
}

#' scaling for for scRNASeq data
#' Helper function that scale the genes expression matrix for the PLS function
#' apply the function scale_zi_nb on a list of genes
#' @param sce SingleCellExperiment object
#' @param genes vector of gene name to scale (rownames of sce)
#' @param assay_name (default: "counts") name of the assay of genes counts to use
#' @param cpus (default: 4) number of cpus to use
#' @return a matrix with the result of the scaling 
#' @examples
#' \dontrun {
#'scaled_counts <- scale_zi_nb_sce(
#'  sce = sce,
#'  genes = rownames(sce),
#'  assay_name = "counts",
#'  cpus = 10
#')
#' }
#' @import plsgenomics pbmcappl SingleCellExperiment Matrix tidyverse
scale_zi_nb_sce <- function(sce,
                            genes = rownames(sce),
                            assay_name = "counts",
                            cpus = 4) {
  pbmcapply::pbmclapply(
    as.list((names(genes) <- genes)),
    function(x, sce, assay_name){
      tibble(count = assay(sce, assay_name)[x, ]) %>%
        scale_zi_nb() %>%
        plyr::rename(x = ., replace = c("count" = x))
    },
    sce = sce,
    assay_name = assay_name,
    mc.cores = cpus
  ) %>%
    do.call(bind_cols, .)
}

# classification function

#' classification for scRNASeq data
#' Helper function that scale the genes expression matrix for the PLS function
#' @examples
#' @import plsgenomics pbmcappl SingleCellExperiment Matrix tidyverse
PLS_scaling <- function(sce, genes, features = NULL, assay_name = "counts", cpus = 4) {
  if (length(colnames(colData(sce)[features])) == 0) {
    scale_zi_nb_sce(
      sce = sce,
      genes = genes,
      assay_name = assay_name,
      cpus = cpus
    ) %>% 
      ascb() %>%
      as.matrix() %>%
      t() %>%
      SingleCellExperiment::SingleCellExperiment(assays = list(counts = .)) %>%
      return()
  }
  cbind(
    scale_zi_nb_sce(
      sce = sce,
      genes = genes,
      assay_name = assay_name,
      cpus = cpus
    ),
    scale(colData(sce)[features])
  ) %>%
    ascb() %>%
    as.matrix() %>%
    t() %>%
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = .))
}

#' classification for scRNASeq data
#' Helper function that compute the X and Xtest matrix for the PLS function
#' @examples
#' @import plsgenomics pbmcappl SingleCellExperiment Matrix tidyverse
PLS_filter <- function(sce,  group_by, genes , features = NULL,
                       assay_name = "counts",
                       altExp_name = "PLS_scaled",
                       cpus = 4){
  if (altExp_name %in% altExpNames(sce)) {
    print(str_c(altExp_name, " altExp, found, skipping PLS filter & scalling"))
    colData(altExp(sce, altExp_name))$group_predict <- NA
    colData(altExp(sce, altExp_name))$group_predict <- 
      !apply(assay(altExp(sce, altExp_name), "counts"), 2, anyNA)
    return(sce)
  }
  altExp(sce, altExp_name) <- PLS_scaling(
      sce = sce,
      genes = genes,
      features = features,
      assay_name = assay_name,
      cpus = cpus
    )
  colData(altExp(sce, altExp_name))$group_name <- group_by
  colData(altExp(sce, altExp_name))$group_train <-
    !is.na(group_by) &
    !apply(assay(altExp(sce, altExp_name), "counts"), 2, anyNA)
  colData(altExp(sce, altExp_name))$group_by <-
    as.numeric(as.factor(group_by)) - 1
  return(sce)
}

#' classification for scRNASeq data
#' copy of the plsgenomics:logit.spls.stab function using pbmcapply to display
#' a progress bar
#' @examples
#' @import plsgenomics pbmcappl 
#' @export plsgenomics_logit_spls_stab
plsgenomics_logit_spls_stab <- function (X, Y, lambda.ridge.range, lambda.l1.range, ncomp.range, 
          adapt = TRUE, maxIter = 100, svd.decompose = TRUE, ncores = 1, 
          nresamp = 100, center.X = TRUE, scale.X = FALSE, weighted.center = TRUE, 
          seed = NULL, verbose = TRUE) 
{
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  index.p <- c(1:p)
  if (is.factor(Y)) {
    Y <- as.numeric(levels(Y))[Y]
  }
  Y <- as.integer(Y)
  Y <- as.matrix(Y)
  q <- ncol(Y)
  one <- matrix(1, nrow = 1, ncol = n)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  cnames <- NULL
  if (!is.null(colnames(X))) {
    cnames <- colnames(X)
  }
  else {
    cnames <- paste0(1:p)
  }
  if (length(table(Y)) > 2) {
    warning("message from logit.spls.stab: multicategorical response")
    results <- multinom.spls.stab(X = X, Y = Y, lambda.ridge.range = lambda.ridge.range, 
                                  lambda.l1.range = lambda.l1.range, ncomp.range = ncomp.range, 
                                  adapt = adapt, maxIter = maxIter, svd.decompose = svd.decompose, 
                                  ncores = ncores, nresamp = nresamp, center.X = center.X, 
                                  scale.X = scale.X, weighted.center = weighted.center, 
                                  seed = seed, verbose = verbose)
    return(results)
  }
  if ((!is.matrix(X)) || (!is.numeric(X))) {
    stop("message from logit.spls.stab: X is not of valid type")
  }
  if (p == 1) {
    stop("message from logit.spls.stab: p=1 is not valid")
  }
  if ((!is.matrix(Y)) || (!is.numeric(Y))) {
    stop("message from logit.spls.stab: Y is not of valid type")
  }
  if (q != 1) {
    stop("message from logit.spls.stab: Y must be univariate")
  }
  if (nrow(Y) != n) {
    stop("message from logit.spls.stab: the number of observations in Y is not equal to the number of row in X")
  }
  if (sum(is.na(Y)) != 0) {
    stop("message from logit.spls.stab: NA values in Ytrain")
  }
  if (sum(!(Y %in% c(0, 1))) != 0) {
    stop("message from logit.spls.stab: Y is not of valid type")
  }
  if (sum(as.numeric(table(Y)) == 0) != 0) {
    stop("message from logit.spls.stab: there are empty classes")
  }
  if (any(!is.numeric(lambda.ridge.range)) || any(lambda.ridge.range < 
                                                  0) || any(!is.numeric(lambda.l1.range)) || any(lambda.l1.range < 
                                                                                                 0) || any(lambda.l1.range > 1)) {
    stop("Message from logit.spls.stab: lambda is not of valid type")
  }
  if (any(!is.numeric(ncomp.range)) || any(round(ncomp.range) - 
                                           ncomp.range != 0) || any(ncomp.range < 0) || any(ncomp.range > 
                                                                                            p)) {
    stop("Message from logit.spls.stab: ncomp is not of valid type")
  }
  if ((!is.numeric(maxIter)) || (round(maxIter) - maxIter != 
                                 0) || (maxIter < 1)) {
    stop("message from logit.spls.stab: maxIter is not of valid type")
  }
  if ((!is.numeric(ncores)) || (round(ncores) - ncores != 0) || 
      (ncores < 1)) {
    stop("message from logit.spls.stab: ncores is not of valid type")
  }
  if ((!is.numeric(nresamp)) || (round(nresamp) - nresamp != 
                                 0) || (nresamp < 1)) {
    stop("message from logit.spls.stab: nresamp is not of valid type")
  }
  grid.resampling <- as.matrix(Reduce("rbind", pbmcapply::pbmclapply(1:nresamp, 
                                                        function(id.samp) {
                                                          ntrain = floor(0.5 * n)
                                                          index.train = sort(sample(1:n, size = ntrain))
                                                          Xtrain = X[index.train, ]
                                                          Ytrain = Y[index.train]
                                                          condition = length(table(Ytrain)) < 2
                                                          test = 0
                                                          while (condition & test < 100) {
                                                            index.train = sort(sample(1:n, size = ntrain))
                                                            Xtrain = X[index.train, ]
                                                            Ytrain = Y[index.train]
                                                            condition = length(table(Ytrain)) < 2
                                                            test = test + 1
                                                          }
                                                          if (test == 100) {
                                                            ind0 = which(Y == 0)
                                                            ind1 = which(Y == 1)
                                                            index.train = sample(ind0, size = 1)
                                                            index.train = c(index.train, sample(ind1, size = 1))
                                                            index.train = c(index.train, sample((1:n)[which(!(1:n) %in% 
                                                                                                              index.train)], size = ntrain - 2))
                                                            Xtrain = X[index.train, ]
                                                            Ytrain = Y[index.train]
                                                          }
                                                          paramGrid <- expand.grid(lambdaL1 = lambda.l1.range, 
                                                                                   lambdaL2 = lambda.ridge.range, ncomp = ncomp.range, 
                                                                                   KEEP.OUT.ATTRS = FALSE)
                                                          grid_out <- as.matrix(Reduce("rbind", lapply(1:nrow(paramGrid), 
                                                                                                       function(gridRow) {
                                                                                                         lambdaL1 <- paramGrid$lambdaL1[gridRow]
                                                                                                         lambdaL2 <- paramGrid$lambdaL2[gridRow]
                                                                                                         ncomp <- paramGrid$ncomp[gridRow]
                                                                                                         fit_out <- logit.spls(Xtrain = Xtrain, Ytrain = Ytrain, 
                                                                                                                               lambda.ridge = lambdaL2, lambda.l1 = lambdaL1, 
                                                                                                                               ncomp = ncomp, Xtest = NULL, adapt = adapt, 
                                                                                                                               maxIter = maxIter, svd.decompose = svd.decompose, 
                                                                                                                               center.X = center.X, scale.X = scale.X, weighted.center = weighted.center)
                                                                                                         sel_var <- fit_out$Anames
                                                                                                         status_var <- rep(0, length(cnames))
                                                                                                         status_var[which(cnames %in% sel_var)] <- rep(1, 
                                                                                                                                                       length(which(cnames %in% sel_var)))
                                                                                                         tmp <- c(lambdaL1, lambdaL2, ncomp, id.samp, 
                                                                                                                  sum(status_var), status_var)
                                                                                                         return(tmp)
                                                                                                       })))
                                                          rownames(grid_out) <- NULL
                                                          return(grid_out)
                                                        }, mc.cores = ncores, mc.silent = !verbose)))
  grid.resampling <- data.frame(grid.resampling)
  colnames(grid.resampling) <- c("lambdaL1", "lambdaL2", "ncomp", 
                                 "id", "nbVar", cnames)
  grid.resampling$point <- paste0(grid.resampling$lambdaL1, 
                                  "_", grid.resampling$lambdaL2, "_", grid.resampling$ncomp)
  o.grid <- order(grid.resampling$nbVar)
  grid.resampling <- grid.resampling[o.grid, ]
  if (any(table(grid.resampling$point) < nresamp)) {
    warning("message from logit.spls.stab: empty classe in a resampling")
    print(table(grid.resampling$point))
  }
  tmp_qLambda <- as.matrix(Reduce("rbind", pbmcapply::pbmclapply(1:nresamp, 
                                                    function(id.samp) {
                                                      tmp1 <- subset(grid.resampling, grid.resampling$id == 
                                                                       id.samp)
                                                      tmp2 <- apply(t(tmp1)[-c(1:5, tail(1:ncol(grid.resampling), 
                                                                                         1)), ], 1, cumsum)
                                                      tmp3 <- apply(t(tmp2), 2, function(x) return(sum(x != 
                                                                                                         0)))
                                                      tmp4 <- cbind(tmp1[, c(1:4)], unname(tmp3))
                                                      return(tmp4)
                                                    }, mc.cores = ncores, mc.silent = !verbose)))
  tmp_qLambda <- data.frame(tmp_qLambda)
  colnames(tmp_qLambda) <- c("lambdaL1", "lambdaL2", "ncomp", 
                             "id.samp", "qLambda")
  qLambda <- plyr::ddply(tmp_qLambda, c("lambdaL1", "lambdaL2", "ncomp"), 
                   function(x) colMeans(x[c("qLambda")], na.rm = TRUE))
  o.qLambda <- order(qLambda$qLambda)
  qLambda <- qLambda[o.qLambda, ]
  probs_lambda <- plyr::ddply(grid.resampling, c("lambdaL1", "lambdaL2", 
                                           "ncomp"), function(x) colMeans(x[cnames], na.rm = TRUE))
  probs_lambda <- probs_lambda[o.qLambda, ]
  return(list(q.Lambda = qLambda, probs.lambda = probs_lambda, 
              p = p))
}

#' classification for scRNASeq data
#'
#' @param sce SingleCellExperiment object
#' @param group_by factor definifing the training group, with NA for the cells to classify
#' @param genes genes names to use
#' @param features colData() of sce to use
#' @param assay_name (default: "counts") name of the assay of genes counts to use
#' @param altExp_name (default: "PLS_scaled") name of the altExp to store the scaled data for the PLS
#' @param force vector of genes names to force the usage of for the classification
#' @param cpus (default: 4) number of cpus to use
#' @return a list() with the result of the classification
#' @examples
#' \dontrun {
#'model <- PLS_fit(
#'  sce = sce,
#'  group_by = sce$manual_cell_type,
#'  genes = genes_PLS %>% pull(genes) %>% na.omit(),
#'  features = genes_PLS %>% pull(proteins) %>% na.omit() %>% unique(),
#'  assay_name = "counts_vst",
#'  altExp_name = "PLS_surface_cell_type",
#'  cpus = 10
#')
#' }
#' @import plsgenomics tidyverse SingleCellExperiment Matrix
#' @export PLS_fit
PLS_fit <- function(sce,
                    group_by,
                    genes,
                    features = NULL,
                    assay_name = "counts",
                    altExp_name = "PLS_scaled",
                    force = c(),
                    cpus = 4) {
  print("scaling data...")
  sce <- PLS_filter(
    sce = sce,
    group_by = group_by,
    genes = c(genes, force) %>% unique(),
    features = features,
    assay_name = assay_name,
    altExp_name = altExp_name,
    cpus = cpus
  )
  print("training sparse PLS...")
  fit_sparse <- assay(
          altExp(sce, altExp_name)[ ,
            colData(altExp(sce, altExp_name))$group_train
          ],
          "counts") %>%
        as.matrix() %>%
        t() %>%
    plsgenomics_logit_spls_stab(
      X = .,
      Y = colData(
        altExp(sce, altExp_name))$group_by[
          colData(altExp(sce, altExp_name))$group_train
      ],
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(
        nrow(assay(altExp(sce, altExp_name), "counts")),
        3)
      ),
      svd.decompose = FALSE,
      center.X = FALSE,
      scale.X = FALSE,
      weighted.center = FALSE,
      ncores = cpus,
      verbose = TRUE
    )
  print("training PLS with selected genes...")
  rowData(altExp(sce, altExp_name))$predictor <-
    rownames(altExp(sce, altExp_name)) %in%
    stability.selection(fit_sparse)$selected.predictors
  rowData(altExp(sce, altExp_name))$predictor_force <-
    rownames(altExp(sce, altExp_name)) %in%
    c(
      stability.selection(fit_sparse)$selected.predictors,
      force
    )
  print(rownames(altExp(sce, altExp_name))[
      rowData(altExp(sce, altExp_name))$predictor
    ])
  print(rownames(altExp(sce, altExp_name))[
      rowData(altExp(sce, altExp_name))$predictor_force
    ])
  fit <- assay(
      altExp(sce, altExp_name)[
          rowData(altExp(sce, altExp_name))$predictor_force,
          colData(altExp(sce, altExp_name))$group_train
        ],
        "counts") %>%
      as.matrix() %>%
      t() %>%
    plsgenomics::rpls.cv(
      Xtrain = .,
      Ytrain = colData(
        altExp(sce, altExp_name))$group_by[
          colData(altExp(sce, altExp_name))$group_train
        ],
      LambdaRange = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      ncompMax = min(3, nrow(assay(altExp(sce, altExp_name), "counts"))),
      ncores = cpus
    )
  return(list(
    fit = fit,
    fit_sparse = fit_sparse,
    sce = sce
  ))
}

#' classification for scRNASeq data
#'
#' @param fit output of PLS_fit
#' @param sce SingleCellExperiment object
#' @param group_by factor definifing the training group, with NA for the cells to classify
#' @param genes genes names to use
#' @param features colData() of sce to use
#' @param assay_name (default: "counts") name of the assay of genes counts to use
#' @param altExp_name (default: "PLS_scaled") name of the altExp to store the scaled data for the PLS
#' @param cpus (default: 4) number of cpus to use
#' @return a list() with the result of the classification
#' @examples
#' \dontrun {
#'model <- PLS_predict(
#'  fit = fit,
#'  sce = fit$sce,
#'  group_by = sce$manual_cell_type,
#'  genes = genes_PLS %>% pull(genes) %>% na.omit(),
#'  features = genes_PLS %>% pull(proteins) %>% na.omit() %>% unique(),
#'  assay_name = "counts_vst",
#'  altExp_name = "PLS_surface_cell_type",
#'  cpus = 10
#')
#' }
#' @import plsgenomics tidyverse SingleCellExperiment Matrix
#' @export PLS_predict
PLS_predict <- function(fit, sce, group_by, genes, features,
                        assay_name = "counts",
                        altExp_name = "PLS_scaled",
                        cpus = 4) {
  print("scaling data...")
  sce <- PLS_filter(
    sce = sce,
    group_by = group_by,
    genes = genes,
    features = features,
    assay_name = assay_name,
    altExp_name = altExp_name,
    cpus = cpus
  )
  model <- assay(altExp(sce, altExp_name)[
          rowData(altExp(sce, altExp_name))$predictor_force,
          colData(altExp(sce, altExp_name))$group_train
        ],
        "counts") %>%
      t() %>%
      as.matrix() %>%
    plsgenomics::rpls(
      Ytrain = colData(altExp(sce, altExp_name))$group_by[
        colData(altExp(sce, altExp_name))$group_train],
      Xtrain = .,
      Lambda = fit$fit$Lambda,
      ncomp = fit$fit$ncomp,
      Xtest = assay(altExp(sce, altExp_name)[
            rowData(altExp(sce, altExp_name))$predictor_force,
            colData(altExp(sce, altExp_name))$group_predict
          ],
          "counts") %>%
        as.matrix() %>%
        t()
    )
  group_names <- colData(altExp(sce, altExp_name)) %>%
    as_tibble() %>%
    drop_na() %>%
    select(group_name, group_by) %>%
    unique()
  colData(sce) <- colData(sce) %>%
    cbind(.,
      tibble(rownames = rownames(colData(sce))) %>%
      mutate(
        rownames = as.vector(rownames)
      ) %>%
      left_join(
        tibble(
          rownames = rownames(model$proba.test),
          !!altExp_name := group_names$group_name[
            match(model$hatYtest, group_names$group_by)
            ],
          !!(str_c("p_", altExp_name)) := model$proba.test[, 1]
        ),
        by = "rownames"
      ) %>%
      select(-rownames)
    )
  return(sce)
}
