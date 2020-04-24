library(tidyverse)
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

logspace <- function( d1, d2, n) {
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
      family = glmmTMB::nbinom1,
      data = data
    ),
    model_zi = glmmTMB::glmmTMB(
      as.formula(formula),
      family = glmmTMB::nbinom1,
      ziformula = as.formula(zi_formula),
      data = data
    ),
    threshold = threshold
  )
}

residuals_zi_nb <- function(model) {
  sim_resid <- DHARMa::simulateResiduals(model)
  tibble(uniformity = DHARMa::testUniformity(
           simulationOutput = sim_resid, plot = F)$p.value,
         dispersion = DHARMa::testDispersion(
           simulationOutput = sim_resid, plot = F)$p.value,
         outlier = DHARMa::testOutliers(
           simulationOutput = sim_resid, plot = F)$p.value,
         zi = !is.null(model$modelInfo$terms$zi)
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

PLS_scaling <- function(sce, genes, features, assay_name = "counts", cpus = 4) {
  cbind(
    scale(colData(sce)[features]),
    parallel::mclapply(
      as.list(genes),
      function(x, sce, assay_name){
        tibble(count = assay(sce, assay_name)[x, ]) %>% 
          scale_zi_nb() %>% 
          plyr::rename(replace = c("count" = x))
      },
      sce = sce,
      assay_name = assay_name,
      mc.cores = cpus
    ) %>% 
      do.call(cbind, .)
  ) %>%
    ascb() %>% 
    as.matrix()
}

PLS_filter <- function(data, group_by){
  names(group_by) <- rownames(data) 
  group_by <- group_by[!is.na(group_by) & !apply(data, 1, anyNA)]
  data <- data[names(group_by), ]
  names(group_by) <- group_by
  group_by <- as.numeric(as.factor(group_by)) - 1
  return(list(data = data, group_by = group_by))
}

PLS_fit <- function(sce, group_by, genes, features, assay_name = "counts", cpus = 4) {
  print("scaling data...")
  data <- PLS_scaling(
    sce = sce,
    genes = genes,
    features = features,
    assay_name = assay_name,
    cpus = cpus
  ) %>%
    PLS_filter(group_by = group_by)
  group_by = data$group_by
  data = data$data
  print("training sparse PLS...")
  fit_sparse <- plsgenomics::logit.spls.stab(
      X = data,
      Y = group_by,
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(ncol(data), 3)),
      svd.decompose = FALSE,
      center.X = FALSE,
      scale.X = FALSE,
      weighted.center = FALSE,
      ncores = cpus,
      verbose = TRUE
    )
  print("training PLS with selected genes...")
  print(stability.selection(fit_sparse)$selected.predictors)
  fit <- plsgenomics::rpls.cv(
    Xtrain = data[, stability.selection(fit_sparse)$selected.predictors],
    Ytrain = group_by,
    LambdaRange = signif(
      logspace(d1 = -2, d2 = 3, n = 21), digits = 3
    ),
    ncompMax = min(3, max(ncol(data))),
    ncores = cpus
  )
  return(list(
    fit = fit,
    fit_sparse = fit_sparse,
    selected = stability.selection(fit_sparse)$selected.predictors,
    data = data,
    group_by = group_by
  ))
}

PLS_predict <- function(fit, sce, genes, features, assay_name = "counts", cpus = 4) {
  print("scaling data...")
  data <- PLS_scaling(
    sce = sce,
    genes = genes,
    features = features,
    assay_name = assay_name,
    cpus = cpus
  )
  data_train <- PLS_filter(data = data, group_by = group_by)
}