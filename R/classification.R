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
  v = F) {

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
  data_train <- rm_data_train$data
  if ( length(force) > 0) {
    if (any(!(force %in% colnames(data_train)))) {
      stop("error: genes or features to force not in the features and genes to
      learn on")
    }
  }
  group_by <- rm_data_train$group_by
  print("training PLS...")
  algo_training <- get(paste0(model_type, "_", algo , "_training"))
  training <- algo_training(
    by = group_by,
    data = data_train,
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
  data <- rm_data$data
  to_predict_on <- rm_data$b_cells
  print("predicting from PLS model...")
  algo_classification <- get(paste0(model_type, "_", algo , "_classification"))
  classification <- algo_classification(
    fit = training,
    data = data,
    data_train = data_train,
    ncores = ncores,
    file = paste0(output_file, "_classification"),
    force = force
  )
  classification$rm_NA_training <- rm_data_train$row_rm
  classification$rm_NA <- rm_data$row_rm
  groups <- rep(NA, scd$getncells)
  groups[to_train_on] <- as.vector(
    scd$select(b_cells = to_train_on)$getfeature(feature)
  )
  groups[to_predict_on][!rm_data$row_rm] <- as.vector(
    classification$model$groups
  )
  pgroups <- rep(NA, scd$getncells)
  pgroups[to_train_on][!rm_data_train$row_rm] <- classification$model$proba
  pgroups[to_predict_on][!rm_data$row_rm] <- classification$model$proba.test
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
  data <- apply(data, 2, scRNAtools::ascb)
  row_rm <- scRNAtools::clean_data(data = data, group_by = group_by)
  if (!any(b_cells[!row_rm])) {
    if (v) {
      print("warning: predicted set empty. Using predicting set.")
    }
    b_cells <- rep(TRUE, nrow(data))
  }
  print(dim(data[b_cells & !row_rm, ]))
  return(
    list(
      data = data[b_cells & !row_rm, ],
      group_by = group_by[b_cells & !row_rm],
      row_rm = row_rm[b_cells],
      b_cells = !b_cells[!row_rm]
    )
  )
}

get_features <- function(scd, features, v = TRUE) {
  data <- c()
  if (!missing(features)) {
    if (length(features) == 1) {
      data <- scd$getfeature(features)
    } else {
      data <- scd$select(features = features)$getfeatures
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
          which(scd$getgenes %in% genes)
        ]
    }
  }
  return(data)
}

clean_data <- function(data, group_by) {
  row_rm <- apply(
    X = data,
    MARGIN = 1,
    FUN = function(x) {
      any(is.na(x))
    }
  )
  return(row_rm)
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

#' @importFrom plsgenomics multinomial.spls.stab
multinomial_spls_stab_training <- function(by, data, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_msplsstab.Rdata"))) {
    print(paste0(file, "_msplsstab.Rdata found. skipping training step..."))
    load(paste0(file, "_msplsstab.Rdata"))
  } else {
    fit <- plsgenomics::multinom.spls.stab(
      X = data,
      Y = by$by,
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(ncol(data), 3)),
      svd.decompose = FALSE,
      center.X = TRUE,
      scale.X = TRUE,
      weighted.center = TRUE,
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

#' @importFrom plsgenomics multinomial.spls.stab
multinomial_spls_cv_training <- function(by, data, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_msplscv.Rdata"))) {
    print(paste0(file, "_msplscv.Rdata found. skipping training step..."))
    load(paste0(file, "_msplscv.Rdata"))
  } else {
    fit <- plsgenomics::multinom.spls.cv(
      X = data,
      Y = by$by,
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(ncol(data), 3)),
      svd.decompose = FALSE,
      center.X = TRUE,
      scale.X = TRUE,
      weighted.center = TRUE,
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

#' @importFrom plsgenomics multinomial.spls
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
      X.center = TRUE,
      X.scale = TRUE,
      weighted.center = TRUE
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

#' @importFrom plsgenomics logistic.spls.stab
logistic_spls_stab_training <- function(by, data, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_lsplsstab.Rdata"))) {
    print(paste0(file, "_lsplsstab.Rdata found. skipping training step..."))
    load(paste0(file, "_lsplsstab.Rdata"))
  } else {
    fit <- plsgenomics::logit.spls.stab(
      X = data,
      Y = by$by,
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(ncol(data), 3)),
      svd.decompose = FALSE,
      center.X = TRUE,
      scale.X = TRUE,
      weighted.center = TRUE,
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

#' @importFrom plsgenomics logistic.spls.stab
logistic_spls_cv_training <- function(by, data, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_lsplscv.Rdata"))) {
    print(paste0(file, "_lsplscv.Rdata found. skipping training step..."))
    load(paste0(file, "_lsplscv.Rdata"))
  } else {
    fit <- plsgenomics::logit.spls.cv(
      X = data,
      Y = by$by,
      lambda.ridge.range = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      lambda.l1.range = seq(0.05, 0.95, by = 0.05),
      ncomp.range = c(1:min(ncol(data), 3)),
      svd.decompose = FALSE,
      center.X = TRUE,
      scale.X = TRUE,
      weighted.center = TRUE,
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
  fit_pls <- logistic_pls_cv_training(
    by = fit$fit$group_by,
    data = data_train[, fit$fit$selected],
    ncores = ncores,
    file = file,
    force = force
  )
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

#' @importFrom plsgenomics logistic.spls
logistic_spls_classification <- function(
    fit, data, data_train, ncores, file, force) {
  if (!missing(file) & file.exists(paste0(file, "_lspls.Rdata"))) {
    print(paste0(file, "_lspls.Rdata found. skipping classification step..."))
    load(paste0(file, "_lspls.Rdata"))
  } else {
    model <- plsgenomics::logit.spls(
      Xtrain = data_train,
      Ytrain = fit$by,
      lambda.ridge = fit$lambda.ridge.opt,
      lambda.l1 = fit$lambda.l1.opt,
      ncomp = fit$ncomp.opt,
      Xtest = data,
      svd.decompose = FALSE,
      X.center = TRUE,
      X.scale = TRUE,
      weighted.center = TRUE
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

#' @importFrom plsgenomics mrpls.cv
logistic_pls_cv_training <- function(by, data, ncores, file, force){
  if (!missing(file) & file.exists(paste0(file, "_lplscv.Rdata"))) {
    print(paste0(file, "_lplscv.Rdata found. skipping training step..."))
    load(paste0(file, "_lplscv.Rdata"))
  } else {
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
