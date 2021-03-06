#' classification for scRNASeq data
#'
#' @param scd is a scRNASeq data object
#' @param feature the feature defining the groups to learn on
#' @param ncores (default:4) number of cpus to use
#' @param features a vector of the names of the feature to classifiy on
#' @param genes a vector of the names of the genes to classifiy on
#' @param training is a vector of cell names to train on
#' @param algo is the algorithm to use
#' @param output_file if non empty file to save the classification results
#' @param force a list of features and genes
#' @return a list() with the result of the classification
#' @examples
#' \dontrun {
#' classif <- classification(scd, 500000)
#' }
#' @import plsgenomics
#' @export classification

classification <- function(
  scd,
  feature,
  ncores = 4,
  features,
  genes,
  algo = "spls_pls",
  selection = "stab",
  output_file = "",
  force = c(),
  weight = TRUE,
  v = F) {

  if (weight) {
    scd <- weight_genes_features(
      scd = scd,
      genes = genes,
      features = features,
      file = output_file,
      cpus = min(6, ncores),
      v = v
    )
  }
  to_train_on <- !is.na(scd$getfeature(feature))
  print("training on:")
  group_by <- group_by_transform(
    scd$select(b_cells = to_train_on)$getfeature(feature)
  )
  model_type <- scRNAtools::pls_type(as.factor(group_by$names))

  print("building training set...")
  rm_data_train <- scRNAtools::get_data(
    scd = scd,
    b_cells = to_train_on,
    features = features,
    genes = genes,
    group_by = group_by,
    v = v
  )
  if ( length(force) > 0) {
    if (any(!(force %in% colnames(rm_data_train$data)))) {
      stop("error: genes or features to force not in the features and genes to
      learn on")
    }
  }
  print("training PLS...")
  algo_training <- get(paste0(model_type, "_", algo, "_training"))
  training <- algo_training(
    by = rm_data_train$group_by,
    data = rm_data_train$data,
    ncores = ncores,
    file = paste0(output_file, "_training"),
    force = force
  )
  print("building set to predict...")
  rm_data <- scRNAtools::get_data(
    scd = scd,
    b_cell = !to_train_on,
    features = features,
    genes = genes,
    group_by = group_by,
    v = v
  )
  print("predicting from PLS model...")
  algo_classification <- get(paste0(model_type, "_", algo , "_classification"))
  if (length(training$group$by) != nrow(rm_data_train$data)) {
    print("reusing old model for prediction, adjusting grouping factor...")
    training$group_by <- rm_data_train$group_by
  }
  classification <- algo_classification(
    fit = training,
    data = rm_data$data,
    data_train = rm_data_train$data,
    ncores = ncores,
    file = paste0(output_file, "_classification"),
    force = force
  )
  groups <- rep(NA, scd$getncells)
  groups[rm_data_train$b_cells] <- as.vector(
    scd$select(b_cells = rm_data_train$b_cells)$getfeature(feature)
  )
  groups[rm_data$b_cells] <- as.vector(
    classification$model$groups
  )
  pgroups <- rep(NA, scd$getncells)
  pgroups[rm_data_train$b_cells] <- classification$model$proba
  pgroups[rm_data$b_cells] <- classification$model$proba.test
  return(list(
    classification = classification,
    scd = scd,
    feature = feature,
    groups = groups,
    pgroups = pgroups
  ))
}

################################################################################
# helper functions

logspace <- function( d1, d2, n) {
  return(exp(log(10) * seq(d1, d2, length.out = n)))
}

weight_regression <- function( gene, scd, v) {
  data <- data.frame(y = round(scd$getgene(gene)))
  is_zi <- zi_test(
    data = data,
    formula = "y ~ 1",
    gene_name = gene,
    model_family = "nbinom1",
    link = "log",
    threshold = 0.05,
    v = v)
  models_result <- scRNAtools::ziNB_fit(
    data = data,
    formula = "y ~ 1",
    gene_name = gene,
    model_family = "nbinom1",
    link = "log",
    zi = is_zi,
    v = v
  )
  zi_weight <- 1
  if (is_zi & !is.na(models_result$residuals[1])) {
    zi_weight <- 1 - models_result$pz
  }
  gene_scale <- 1
  if (!is.na(models_result$residuals[1])) {
    gene_scale <- exp(models_result$b) * models_result$alpha
  }
  return(list(
    gene = gene,
    gene_weight = zi_weight,
    gene_scale = gene_scale
  ))
}

get_weights <- function(scd, genes, cpus = 1, v = TRUE) {
  genes_list <- as.list(unique(genes))
  results <- list()
  if (cpus > 1) {
    results <- parallel::mclapply(
      X = genes_list,
      FUN  = function(x, scd, v){
        scRNAtools::weight_regression(x, scd, v)
      },
      mc.cores = cpus,
      scd = scd,
      v = v
    )
  } else {
    results <- lapply(
      X = genes_list,
      FUN  = function(x, scd, v){
        scRNAtools::weight_regression(x, scd, v)
      },
      scd = scd,
      v = v
    )
  }
  results_unlisted <- as.data.frame(do.call(rbind, results))
  return(results_unlisted)
}

weight_genes_features <- function(scd, genes, features, file, cpus = 1, v = T) {
  if (!missing(file) & file.exists(paste0(file, "_weighting.Rdata"))) {
    print(paste0(file, "_weighting.Rdata found. skipping training step..."))
    load(paste0(file, "_weighting.Rdata"))
  } else {
    print("scaling genes...")
    weight <- scRNAtools::get_weights(scd, genes, cpus, v)
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
    infos <- scd$getfeatures
    if (length(features) > 0) {
      print("scaling features...")
      infos[, features] <- apply(
        infos[, features],
        2,
        FUN = function(x) {
          scale(as.numeric(as.vector(x)))
        })
    }
    genes_na <-apply(t(weighted_counts), 2, FUN = function(x) {
                      any(is.na(x))
        })
    if (any(!genes_na)) {
      print(paste0("warning: ", genes[genes_na],
                  " contain NA after weighting, droping them"))
      print(weight[genes_na,])
    }
    scd_new <- scdata$new(
      infos = infos,
      counts = t(weighted_counts)[ ,!genes_na]
    )
    if (!missing(file)) {
      save(scd_new, file = paste0(file, "_weighting.Rdata"))
    }
  }
  return(scd_new)
}

get_data <- function(scd, b_cells, features, genes, group_by, v = TRUE) {
  if (missing(features) & missing(genes)) {
    stop("error: missing features and / or genes to classify on")
  }
  data <- scRNAtools::get_features(
    scd = scd, features = features
  )
  data <- cbind(
    data,
    scRNAtools::get_genes(scd = scd, genes = genes)
  )
  data <- apply(as.matrix(data), c(1, 2), as.numeric)
  data <- scRNAtools::ascb(data)
  if (!any(b_cells)) {
    if (v) {
      print("warning: predicted set empty. Using predicting set.")
    }
    b_cells <- rep(TRUE, nrow(data))
    group_by$by <- group_by$by[b_cells]
  } else {
    group_by$by <- group_by$by
  }
  print(dim(data[b_cells, ]))
  print(length(group_by$by))
  return(
    list(
      data = data[b_cells, ],
      b_cells = b_cells,
      group_by = group_by
    )
  )
}

get_features <- function(scd, features, v = TRUE) {
  data <- c()
  if (!missing(features)) {
    if (length(features) == 1) {
      data <- scd$getfeature(features)
    }
    if (length(features) > 1) {
      data <- scd$select(features = unique(features))$getfeatures
    }
  }
  return(data)
}

get_genes <- function(scd, genes, v = TRUE) {
  data <- c()
  if (!missing(genes)) {
    if (length(genes) == 1) {
      data <- scd$getgene(genes)
    } else {
      data <- scd$getcounts[,
          which(scd$getgenes %in% unique(genes))
        ]
    }
  }
  return(data)
}

group_by_transform <- function(by) {
  by <- as.factor(as.vector(by))
  print(summary(by))
  return(list(by = as.numeric(by) - 1, names = levels(by)))
}

pls_type <- function(by){
  if (length(levels(by)) == 2) {
    print(paste0(length(levels(by)), " levels detected using logistic model"))
    return("logistic")
  }
  print(paste0(length(levels(by)), " levels detected using multinomial model"))
  return("multinomial")
}

################################################################################
# multinomial functions

#' @importFrom plsgenomics multinom.spls.stab
multinomial_spls_stab_training <- function(by, data, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_msplsstab.Rdata"))) {
    print(paste0(file, "_msplsstab.Rdata found. skipping training step..."))
    load(paste0(file, "_msplsstab.Rdata"))
  } else {
    print(paste0(file, "_msplsstab.Rdata not found. training step..."))
    fit <- plsgenomics::multinom.spls.stab(
      X = data,
      Y = by$by,
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(ncol(data), 3)),
      svd.decompose = FALSE,
      center.X = FALSE,
      scale.X = FALSE,
      weighted.center = FALSE,
      ncores = ncores,
      verbose = TRUE
    )
    fit$group_by <- by
    if (!missing(file)) {
      save(fit, file = paste0(file, "_msplsstab.Rdata"))
    }
  }
  return(list(fit = fit))
}

#' @importFrom plsgenomics multinom.spls.stab
multinomial_spls_cv_training <- function(by, data, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_msplscv.Rdata"))) {
    print(paste0(file, "_msplscv.Rdata found. skipping training step..."))
    load(paste0(file, "_msplscv.Rdata"))
  } else {
    print(paste0(file, "_msplscv.Rdata not found. training step..."))
    fit <- plsgenomics::multinom.spls.cv(
      X = data,
      Y = by$by,
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(ncol(data), 3)),
      svd.decompose = FALSE,
      center.X = FALSE,
      scale.X = FALSE,
      weighted.center = FALSE,
      ncores = ncores,
      verbose = TRUE
    )
    fit$group_by <- by
    if (!missing(file)) {
      save(fit, file = paste0(file, "_msplscv.Rdata"))
    }
  }
  return(list(fit = fit))
}


#' @importFrom plsgenomics stability.selection
multinomial_spls_stab_classification <- function(
    fit, data, data_train, ncores, file, force) {
  fit$fit$selected <- unique(c(
    plsgenomics::stability.selection(fit$fit)$selected.predictor,
    force
  ))
  print("training PLS with selected genes...")
  fit_pls <- multinomial_pls_cv_training(
    by = fit$fit$group_by,
    data = data_train[, fit$fit$selected],
    ncores = ncores,
    file = file,
    force = force,
  )
  model <- multinomial_pls_cv_classification(
    fit = fit_pls$fit,
    data = data[, fit$fit$selected],
    data_train = data_train[, fit$fit$selected],
    file = file,
    force = force
  )
  return(list(
    fit_spls = fit,
    fit_pls = fit_pls,
    model = model$model,
  ))
}

#' @importFrom plsgenomics multinom.spls
multinomial_spls_classification <- function(
    fit, data, data_train, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_mspls.Rdata"))) {
    print(paste0(file, "_mspls.Rdata found. skipping classification step..."))
    load(paste0(file, "_mspls.Rdata"))
  } else {
    model <- plsgenomics::multinom.spls(
      Xtrain = data_train,
      Ytrain = fit$by,
      lambda.ridge = fit$lambda.ridge.opt,
      lambda.l1 = fit$lambda.l1.opt,
      ncomp = fit$ncomp.opt,
      Xtest = data,
      svd.decompose = FALSE,
      X.center = FALSE,
      X.scale = TRUE,
      weighted.center = FALSE
    )
    model$groups_names <- fit$group_by$names
    model$groups <- model$groups_names[model$hatYtest + 1]
    if (!missing(file)) {
      save(model, file = paste0(file, "_mspls.Rdata"))
    }
  }
  return(list(
    fit = fit,
    model = model
  ))
}


multinomial_spls_cv_classification <- function(
    fit, data, data_train, ncores, file, force) {
  model_spls <- multinomial_spls_classification(
    fit = fit$fit,
    data = data,
    data_train = data_train,
    ncores = ncores,
    file = file,
    force = force
  )
  model_spls$selected <- unique(c(
    colnames(data)[model_spls$model$A.full],
    force
  ))
  fit_pls <- multinomial_pls_cv_training(
    by = model_spls$fit$group_by,
    data = data_train[, model_spls$selected],
    ncores = ncores,
    file = file,
    force = force
  )
  model <- multinomial_pls_cv_classification(
    fit = fit_pls$fit,
    data = data[, model_spls$selected],
    data_train = data[, model_spls$selected],
    file = file,
    force = force
  )
  return(list(
    fit_spls = model_spls$fit,
    model_spls = model_spls$model,
    fit_pls = fit_pls,
    model = model$model
  ))
}

multinomial_pls_stab_training <- function(by, data, ncores, file, force) {
  print("no stability selection for multinomial pls switching to CV")
  return(multinomial_pls_cv_training(data, by, ncores, file, force))
}

#' @importFrom plsgenomics mrpls.cv
multinomial_pls_cv_training <- function(by, data, ncores, file, force){
  if (!missing(file) & file.exists(paste0(file, "_mplscv.Rdata"))) {
    print(paste0(file, "_mplscv.Rdata found. skipping training step..."))
    load(paste0(file, "_mplscv.Rdata"))
  } else {
    print(paste0(file, "_mplscv.Rdata not found. training step..."))
    fit <- plsgenomics::mrpls.cv(
      Xtrain = data,
      Ytrain = by$by,
      LambdaRange = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      ncompMax = min(3, max(ncol(data))),
      ncores = ncores
    )
    fit$group_by <- by
    if (!missing(file)) {
      save(fit, file = paste0(file, "_mplscv.Rdata"))
    }
  }
  return(list(fit = fit))
}

multinomial_pls_stab_classification <- function(
    fit, data, data_train, ncores, file, force) {
  multinomial_pls_cv_classification(fit, data, data_train, ncores, file, force)
}

#' @importFrom plsgenomics mrpls
multinomial_pls_cv_classification <- function(
    fit, data, data_train, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_mpls.Rdata"))) {
    print(paste0(file, "_mpls.Rdata found. skipping classification step..."))
    load(paste0(file, "_mpls.Rdata"))
  } else {
    print(paste0(file, "_mpls.Rdata not found. training step..."))
    model <- plsgenomics::mrpls(
      Xtrain = data_train,
      Ytrain = fit$group_by$by,
      Lambda = fit$Lambda,
      ncomp = fit$ncomp,
      Xtest = data
    )
    model$groups_names <- fit$group_by$names
    model$groups <- model$groups_names[model$hatYtest + 1]
    if (!missing(file)) {
      save(model, file = paste0(file, "_mpls.Rdata"))
    }
  }
  return(list(fit = fit, model = model))
}

################################################################################
# logistic functions

#' @importFrom plsgenomics logit.spls.stab
logistic_spls_stab_training <- function(by, data, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_lsplsstab.Rdata"))) {
    print(paste0(file, "_lsplsstab.Rdata found. skipping training step..."))
    load(paste0(file, "_lsplsstab.Rdata"))
  } else {
    print(paste0(file, "_lsplsstab.Rdata not found. training step..."))
    fit <- plsgenomics::logit.spls.stab(
      X = data,
      Y = by$by,
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(ncol(data), 3)),
      svd.decompose = FALSE,
      center.X = FALSE,
      scale.X = FALSE,
      weighted.center = FALSE,
      ncores = ncores,
      verbose = TRUE
    )
    fit$group_by <- by
    if (!missing(file)) {
      save(fit, file = paste0(file, "_lsplsstab.Rdata"))
    }
  }
  return(list(fit = fit))
}

#' @importFrom plsgenomics logit.spls.cv
logistic_spls_cv_training <- function(by, data, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_lsplscv.Rdata"))) {
    print(paste0(file, "_lsplscv.Rdata found. skipping training step..."))
    load(paste0(file, "_lsplscv.Rdata"))
  } else {
    print(paste0(file, "_lsplscv.Rdata not found. training step..."))
    fit <- plsgenomics::logit.spls.cv(
      X = data,
      Y = by$by,
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(ncol(data), 3)),
      svd.decompose = FALSE,
      center.X = FALSE,
      scale.X = FALSE,
      weighted.center = FALSE,
      ncores = ncores,
      verbose = TRUE
    )
    fit$group_by <- by
    if (!missing(file)) {
      save(fit, file = paste0(file, "_lsplscv.Rdata"))
    }
  }
  return(list(fit = fit))
}


#' @importFrom plsgenomics stability.selection
logistic_spls_stab_classification <- function(
    fit, data, data_train, ncores, file, force) {
  fit$fit$selected <- unique(c(
    plsgenomics::stability.selection(fit$fit)$selected.predictor,
    force
  ))
  print("training PLS with selected genes...")
  print(fit$fit$selected)
  print(str(fit$fit$group_by))
  fit_pls <- logistic_pls_cv_training(
    by = fit$fit$group_by,
    data = data_train[, fit$fit$selected],
    ncores = ncores,
    file = file,
    force = force
  )
  if (length(fit_pls$group_by$by) != length(fit$fit$group_by$by)) {
    fit_pls$group_by$by <- fit$fit$group_by$by
    fit_pls$fit$group_by$by <- fit$fit$group_by$by
  }
  model <- logistic_pls_cv_classification(
    fit = fit_pls$fit,
    data = data[, fit$fit$selected],
    data_train = data_train[, fit$fit$selected],
    file = file,
    force = force
  )
  return(list(
    fit_spls = fit,
    fit_pls = fit_pls,
    model = model$model,
    force = force
  ))
}

#' @importFrom plsgenomics logit.spls
logistic_spls_classification <- function(
    fit, data, data_train, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_lspls.Rdata"))) {
    print(paste0(file, "_lspls.Rdata found. skipping classification step..."))
    load(paste0(file, "_lspls.Rdata"))
  } else {
    print(paste0(file, "_lspls.Rdata not found. training step..."))
    model <- plsgenomics::logit.spls(
      Xtrain = data_train,
      Ytrain = fit$by,
      lambda.ridge = fit$lambda.ridge.opt,
      lambda.l1 = fit$lambda.l1.opt,
      ncomp = fit$ncomp.opt,
      Xtest = data,
      svd.decompose = FALSE,
      center.X = FALSE,
      scale.X = FALSE,
      weighted.center = FALSE
    )
    model$groups_names <- fit$group_by$names
    model$groups <- model$groups_names[model$hatYtest + 1]
    if (!missing(file)) {
      save(model, file = paste0(file, "_lspls.Rdata"))
    }
  }
  return(list(
    fit = fit,
    model = model
  ))
}


logistic_spls_cv_classification <- function(
    fit, data, data_train, ncores, file, force) {
  model_spls <- logistic_spls_classification(
    fit = fit$fit,
    data = data,
    data_train = data_train,
    ncores = ncores,
    file = file,
    force = force
  )
  model_spls$selected <- unique(c(
    colnames(data)[model_spls$model$A.full],
    force
  ))
  fit_pls <- logistic_pls_cv_training(
    by = model_spls$fit$group_by,
    data = data_train[, model_spls$selected],
    ncores = ncores,
    file = file,
    force = force
  )
  model <- logistic_pls_cv_classification(
    fit = fit_pls$fit,
    data = data[, model_spls$selected],
    data_train = data[, model_spls$selected],
    file = file,
    force = force
  )
  return(list(
    fit_spls = model_spls$fit,
    model_spls = model_spls$model,
    fit_pls = fit_pls,
    model = model$model,
    force = force
  ))
}

logistic_pls_stab_training <- function(by, data, ncores, file, force) {
  print("no stability selection for logistic pls switching to CV")
  return(logistic_pls_cv_training(data, by, ncores, file, force))
}

#' @importFrom plsgenomics rpls.cv
logistic_pls_cv_training <- function(by, data, ncores, file, force){
  if (!missing(file) & file.exists(paste0(file, "_lplscv.Rdata"))) {
    print(paste0(file, "_lplscv.Rdata found. skipping training step..."))
    load(paste0(file, "_lplscv.Rdata"))
  } else {
    print(paste0(file, "_lplscv.Rdata not found. training step..."))
    print(dim(data))
    print(length(by$by))
    print(str(by))
    fit <- plsgenomics::rpls.cv(
      Xtrain = data,
      Ytrain = by$by,
      LambdaRange = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      ncompMax = min(3, max(ncol(data))),
      ncores = ncores
    )
    fit$group_by <- by
    if (!missing(file)) {
      save(fit, file = paste0(file, "_lplscv.Rdata"))
    }
  }
  return(list(fit = fit))
}

logistic_pls_stab_classification <- function(
    fit, data, data_train, ncores, file, force) {
  logistic_pls_cv_classification(fit, data, data_train, ncores, file, force)
}

#' @importFrom plsgenomics mrpls
logistic_pls_cv_classification <- function(
    fit, data, data_train, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_lpls.Rdata"))) {
    print(paste0(file, "_lpls.Rdata found. skipping classification step..."))
    load(paste0(file, "_lpls.Rdata"))
  } else {
    print(paste0(file, "_lpls.Rdata not found. training step..."))
    print(dim(data))
    print(dim(data_train))
    model <- plsgenomics::rpls(
      Xtrain = data_train,
      Ytrain = fit$group_by$by,
      Lambda = fit$Lambda,
      ncomp = fit$ncomp,
      Xtest = data
    )
    model$groups_names <- fit$group_by$names
    model$groups <- model$groups_names[model$hatYtest + 1]
    if (!missing(file)) {
      save(model, file = paste0(file, "_lpls.Rdata"))
    }
  }
  return(list(fit = fit, model = model))
}


#' classification for scRNASeq data
#'
#' @param scd is a scRNASeq data object
#' @param gene gene_name to find other correlated genes
#' @param score factor to find other correlated genes
#' @param cpus (default:4) number of cpus to use
#' @return a vector with the name of the genes selected
#' @examples
#' \dontrun {
#' gene_list <- gene_cov(scd, "LTS")
#' }
#' @import plsgenomics
#' @export gene_cov
gene_cov <- function(
    scd, gene, score, cpus = 4, tmp_file, sparse = T, ncomp = 3,
    FUN = function(x){
      x <- log( x + 1 )
      x <- (x - mean(x) ) / sd(x)
    }){
  if (missing(gene) & missing(score)) {
    return(c())
  }
  scd_norm <- zinorm(
    scd = scd,
    cpus = cpus,
    file = paste0(tmp_file, "_norm.Rdata")
  )
  if (missing(score)) {
    b_genes <- scd_norm$getgenes %in% gene
    gene_exp <- scd_norm$getgene(gene)
  } else {
    gene_exp <- score
    b_genes <- rep(F, scd_norm$getngenes)
  }
  data_exp <- as.matrix(scd_norm$getcounts[, !b_genes])
  data_exp <- data_exp[, !apply( data_exp, 2, function(x){ any(is.na(x) | is.infinite(x)) } )]
  data_exp <- FUN(data_exp)
  if (sparse) {
    if (!missing(tmp_file) & file.exists(paste0(tmp_file, "fit.Rdata"))) {
      print(paste0(tmp_file, "fit.Rdata, found, loading..."))
      load(paste0(tmp_file, "fit.Rdata"))
    } else {
      fit <- spls.cv(
        X = data_exp,
        Y = gene_exp,
        lambda.l1.range = seq(0.05, 0.95, by=0.3),
        ncomp.range = 1:2,
        weight.mat = NULL, adapt = TRUE, center.X = FALSE,
        center.Y = FALSE, scale.X = FALSE, scale.Y = FALSE,
        weighted.center = FALSE,
        return.grid = TRUE, ncores = cpus, nfolds = 10)
    }
    if (!missing(tmp_file)) {
      save(fit, file = paste0(tmp_file, "fit.Rdata"))
    }
    if (!missing(tmp_file) & file.exists(paste0(tmp_file, "predict.Rdata"))) {
      print(paste0(tmp_file, "predict.Rdata, found, loading..."))
      load(paste0(tmp_file, "predict.Rdata"))
    } else {
      predict <- spls(
        Xtrain = data_exp,
        Ytrain = gene_exp,
        lambda.l1 = fit$lambda.l1.opt, ncomp=fit$ncomp.opt,
        weight.mat = NULL,
        Xtest = data_exp,
        adapt = TRUE, center.X = FALSE, center.Y = FALSE, scale.X = FALSE,
        scale.Y = FALSE, weighted.center = FALSE
      )
    }
    if (!missing(tmp_file)) {
      save(predict, file = paste0(tmp_file, "predict.Rdata"))
    }
    genes_names <- colnames(scd$getcounts[, !b_genes])[predict$A]
    return(genes_names)
  } else {
    fit <- pls.regression(
      Xtrain = data_exp,
      Ytrain = gene_exp,
      ncomp = ncomp
    )
    return(fit)
  }
}
