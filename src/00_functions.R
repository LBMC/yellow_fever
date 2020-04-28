library(tidyverse)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(broom)
library(broom.mixed)
library(DHARMa)
library(glmmTMB)
library(parallel)
library(plsgenomics)

anscombe <- function(x){
  2 * sqrt(x + 3 / 8)
}

ascb <- function(x){
  anscombe(x) - 1.224745
}

logspace <- function(d1, d2, n) {
  return(exp(log(10) * seq(d1, d2, length.out = n)))
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

LRT_zi_nb <- function(data, test, formula, zi_formula = formula, threshold = 0.05){
  model <- fit_zi_nb(
    data = data,
    formula = formula,
    zi_formula = zi_formula,
    threshold = threshold
  )
  list(
    parameters = model %>%
      broom.mixed::tidy(),
    LRT = model %>%
      stats::drop1(as.formula(test), test = "Chisq") %>%
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
    select(
      any_of(
        str_split(
          str_c(formula, zi_formula),
          "[ ~+()|]"
        )[[1]]
      ), gene_name
    ) %>%
    drop_na()
}

DEA <- function(sce, 
                test, formula, zi_formula = formula, zi_threshold = 0.05,
                assay_name = "counts",
                cpus = 4){
  results <- parallel::mclapply(
    rownames(sce),
    FUN = function(x, sce, assay_name, test, formula, zi_formula, zi_threshold){
      print(x)
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
    mc.cores = cpus
  )
  names(results) <- rownames(sce)
  return(results)
}

get_genes_pval <- function(results) {
  lapply(results,
       FUN = function(x){
         x$LRT$p.value[2]
       }) %>% do.call(c, .)
}

scale_zi_nb <- function(data) {
  model <- fit_zi_nb(
      data = data,
      formula = "count ~ 1"
    )
  zi_prop <- ifelse(
    "betazi" %in% names(model$fit$par),
    exp(model$fit$par["betazi"]),
    0)
  data$count <- (data$count / sigma(model)) * (1 - zi_prop)
  return(data)
}

# classification function

PLS_scaling <- function(sce, genes, features, assay_name = "counts", cpus = 4) {
  cbind(
    scale(colData(sce)[features]),
    parallel::mclapply(
      as.list(genes),
      function(x, sce, assay_name){
        tibble(count = assay(sce, assay_name)[x, ]) %>% 
          scale_zi_nb() %>% 
          plyr::rename(x = ., replace = c("count" = x))
      },
      sce = sce,
      assay_name = assay_name,
      mc.cores = cpus
    ) %>% 
      do.call(cbind, .)
  ) %>%
    ascb() %>% 
    as.matrix() %>% 
    t() %>% 
    SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = .)
    )
}

PLS_filter <- function(sce,  group_by, genes , features,
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

PLS_fit <- function(sce,
                    group_by,
                    genes,
                    features,
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
  print("training sparse PLS...")
  fit_sparse <- assay(
          altExp(sce, altExp_name)[ ,
            colData(altExp(sce, altExp_name))$group_train
          ],
          "counts") %>%
        as.matrix() %>%
        t() %>% 
    plsgenomics::logit.spls.stab(
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
  print(stability.selection(fit_sparse)$selected.predictors)
  rowData(altExp(sce, altExp_name))$predictor <-
    rownames(altExp(sce, altExp_name)) %in% 
    stability.selection(fit_sparse)$selected.predictors
  fit <- assay(
      altExp(sce, altExp_name)[
          rowData(altExp(sce, altExp_name))$predictor,
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
#'  sce = fit$sce,
#'  fit = fit,
#'  group_by = sce$manual_cell_type,
#'  genes = genes_PLS %>% pull(genes) %>% na.omit(),
#'  features = genes_PLS %>% pull(proteins) %>% na.omit() %>% unique(),
#'  assay_name = "counts_vst",
#'  altExp_name = "PLS_surface_cell_type",
#'  cpus = 10
#')
#' }
#' @import plsgenomics
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
          rowData(altExp(sce, altExp_name))$predictor,
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
            rowData(altExp(sce, altExp_name))$predictor,
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