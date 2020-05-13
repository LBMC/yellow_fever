library(tidyverse)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(sctransform)
library(broom)
library(broom.mixed)
library(DHARMa)
library(glmmTMB)
library(parallel)
library(pbmcapply)
library(plsgenomics)

anscombe <- function(x){
  2.0 * sqrt(x + 3.0 / 8.0)
}

ascb <- function(x){
  anscombe(x) - anscombe(0.0)
}

logspace <- function(d1, d2, n) {
  return(exp(log(10) * seq(d1, d2, length.out = n)))
}

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
    mc.cores = cpus
  )
  names(results) <- rownames(sce)
  return(results)
}

get_genes_pval <- function(results) {
  lapply(results,
       FUN = function(x){
         if ("LRT" %in% names(x)) {
           return(x$LRT$p.value[2])
         } else {
           return(NA)
         }
       }) %>% do.call(c, .)
}

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

PLS_filter <- function(sce,  group_by, genes , features = NULL,
                       assay_name = "counts",
                       altExp_name = "PLS_scaled",
                       cpus = 4){
  if (altExp_name %in% altExpNames(sce)) {
    print(str_c(altExp_name, " altExp, found, skipping PLS filter & scalling"))
    colData(altExp(sce, altExp_name))$group_predict <- NA
    colData(altExp(sce, altExp_name))$group_predict <- 
      !apply(assay(altExp(sce, altExp_name), "counts"), 2, anyNA)
    assay(altExp(sce, altExp_name), "counts") %>% t() %>% summary() %>% print()
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
  rowData(altExp(sce, altExp_name))$predictor <-
    rownames(altExp(sce, altExp_name)) %in%
    stability.selection(fit_sparse)$selected.predictors
  rowData(altExp(sce, altExp_name))$predictor_force <-
    rownames(altExp(sce, altExp_name)) %in%
    c(
      stability.selection(fit_sparse)$selected.predictors,
      force
    )
  rowData(altExp(sce, altExp_name))$predictor %>% table() %>% print()
  rowData(altExp(sce, altExp_name))$predictor_force %>% table() %>% print()
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
