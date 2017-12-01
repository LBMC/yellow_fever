#' classification for scRNASeq data
#'
#' @param scd is a scRNASeq data object
#' @param feature is the name of the feature to classifiy on
#' @param training is a vector of cell names to train on
#' @param algo is the algorithm to use
#' @param output_file if non empty file to save the classification results
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
  output_file = "") {

  to_train_on <- !is.na(scd$getfeature(feature))
  print("training on:")
  group_by <- group_by_transform(
    scd$select(b_cells = to_train_on)$getfeature(feature)
  )
  print("building training set...")
  data_train <- scRNAtools::get_data(scd, to_train_on, features, genes)

  rm_data_train <- clean_data(data_train, group_by)
  data_train <- rm_data_train$data
  group_by <- rm_data_train$group_by

  model_type <- scRNAtools::pls_type(by)
  print("training PLS...")
  algo_training <- get(paste0(model_type, "_", algo , "_training"))
  training <- algo_training(
    by = group_by,
    data = data_train,
    ncores = ncores,
    file = paste0(output_file, "_training")
  )

  print("building set to predict...")
  data <- scRNAtools::get_data(scd, !to_train_on, features, genes)
  rm_data <- clean_data(data, group_by)
  data <- rm_data$data
  print("predicting from PLS model...")
  algo_classification <- get(paste0(model_type, "_", algo , "_classification"))
  classification <- algo_classification(
    fit = training,
    data = data,
    data_train = data_train,
    ncores = ncores,
    file = paste0(output_file, "_classification")
  )
  classification$rm_NA_training <- rm_data_train$row_rm
  classification$rm_NA <- rm_data$row_rm
  groups <- rep(NA, scd$getncells)
  groups[to_train_on] <- as.vector(
    scd$select(b_cells = to_train_on)$getfeature(feature)
  )
  groups[!to_train_on][!rm_data$row_rm] <- as.vector(
    classification$model$groups
  )
  return(list(
    classification = classification,
    scd = scd,
    feature = feature,
    groups = groups
  ))
}

################################################################################
# helper functions

logspace <- function( d1, d2, n) {
  return(exp(log(10) * seq(d1, d2, length.out = n)))
}

get_data <- function(scd, b_cells, features, genes) {
  if (missing(features) & missing(genes)) {
    stop("error: missing features and / or genes to classify on")
  }
  data <- c()
  if (!missing(features)) {
    if (length(features) == 1) {
      data <- scd$select(
        b_cells = b_cells)$getfeature(features)
    } else {
      data <- scd$select(b_cells = b_cells, features = features)$getfeatures
    }
  }
  if (!missing(genes)) {
    if (length(genes) == 1) {
      data <- cbind(
        data,
        scd$select(b_cells = b_cells)$getgene(genes)
      )
    } else {
      data <- cbind(
        data,
        scd$select(b_cells = b_cells)$getcounts
      )
    }
  }
  data <- apply(as.matrix(data), c(1, 2), as.numeric)
  data <- apply(data, 2, scRNAtools::ascb)
  print(dim(data))
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
  if (missing(group_by)) {
    return(list(
      data = data[!row_rm, ], row_rm = row_rm
    ))
  } else {
    return(list(
      data = data[!row_rm, ], group_by = group_by[!row_rm], row_rm = row_rm
    ))
  }
}

group_by_transform <- function(by) {
  by <- as.factor(as.vector(by))
  print(summary(by))
  return(list(by = as.numeric(by) - 1, names = levels(by)))
}

pls_type <- function(by){
  if (length(levels(by)) == 2) {
    return("logistic")
  }
  return("multinomial")
}

################################################################################
# multinomial functions

#' @importFrom plsgenomics multinomial.spls.stab
multinomial_spls_stab_training <- function(by, data, ncores, file) {
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
multinomial_spls_cv_training <- function(by, data, ncores, file) {
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
    fit, data, data_train, ncores, file) {
  fit$fit$selected <- plsgenomics::stability.selection(fit$fit)$selected.predictor
  print("training PLS with selected genes...")
  fit_pls <- multinomial_pls_cv_training(
    by = fit$fit$group_by,
    data = data_train[, fit$fit$selected],
    ncores = ncores,
    file = file
  )
  model <- multinomial_pls_cv_classification(
    fit = fit_pls$fit,
    data = data[, fit$fit$selected],
    data_train = data_train[, fit$fit$selected],
    file = file
  )
  return(list(
    fit_spls = fit,
    fit_pls = fit_pls,
    model = model$model
  ))
}

#' @importFrom plsgenomics multinomial.spls
multinomial_spls_classification <- function(
    fit, data, data_train, ncores, file) {
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
      Xtest = data
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
    fit, data, data_train, ncores, file) {
  model_spls <- multinomial_spls_classification(
    fit = fit$fit,
    data = data,
    data_train = data_train,
    ncores = ncores,
    file = file)
  fit_pls <- multinomial_pls_cv_training(
    by = model_spls$fit$group_by,
    data = data_train[, model_spls$model$A.full],
    ncores = ncores,
    file = file
  )
  model <- multinomial_pls_cv_classification(
    fit = fit_pls$fit,
    data = data[, model_spls$model$A.full],
    data_train = data[, model_spls$model$A.full],
    file = file)
  return(list(
    fit_spls = model_spls$fit,
    model_spls = model_spls$model,
    fit_pls = fit_pls,
    model = model$model
  ))
}

multinomial_pls_stab_training <- function(by, data, ncores, file) {
  print("no stability selection for multinomial pls switching to CV")
  return(multinomial_pls_cv_training(data, by, ncores, file))
}

#' @importFrom plsgenomics mrpls.cv
multinomial_pls_cv_training <- function(by, data, ncores, file){
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
      ncompMax = 3,
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
    fit, data, data_train, ncores, file) {
  multinomial_pls_cv_classification(fit, data, data_train, ncores, file)
}

#' @importFrom plsgenomics mrpls
multinomial_pls_cv_classification <- function(
    fit, data, data_train, ncores, file) {
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
logistic_spls_stab_training <- function(by, data, ncores, file) {
  if (!missing(file) & file.exists(paste0(file, "_msplsstab.Rdata"))) {
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
logistic_spls_cv_training <- function(by, data, ncores, file) {
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
    fit, data, data_train, ncores, file) {
  fit$fit$selected <- plsgenomics::stability.selection(fit$fit)$selected.predictor
  print("training PLS with selected genes...")
  fit_pls <- logistic_pls_cv_training(
    by = fit$fit$group_by,
    data = data_train[, fit$fit$selected],
    ncores = ncores,
    file = file
  )
  model <- logistic_pls_cv_classification(
    fit = fit_pls$fit,
    data = data[, fit$fit$selected],
    data_train = data_train[, fit$fit$selected],
    file = file
  )
  return(list(
    fit_spls = fit,
    fit_pls = fit_pls,
    model = model$model
  ))
}

#' @importFrom plsgenomics logistic.spls
logistic_spls_classification <- function(
    fit, data, data_train, ncores, file) {
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
      Xtest = data
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
    fit, data, data_train, ncores, file) {
  model_spls <- logistic_spls_classification(
    fit = fit$fit,
    data = data,
    data_train = data_train,
    ncores = ncores,
    file = file)
  fit_pls <- logistic_pls_cv_training(
    by = model_spls$fit$group_by,
    data = data_train[, model_spls$model$A.full],
    ncores = ncores,
    file = file
  )
  model <- logistic_pls_cv_classification(
    fit = fit_pls$fit,
    data = data[, model_spls$model$A.full],
    data_train = data[, model_spls$model$A.full],
    file = file)
  return(list(
    fit_spls = model_spls$fit,
    model_spls = model_spls$model,
    fit_pls = fit_pls,
    model = model$model
  ))
}

logistic_pls_stab_training <- function(by, data, ncores, file) {
  print("no stability selection for logistic pls switching to CV")
  return(logistic_pls_cv_training(data, by, ncores, file))
}

#' @importFrom plsgenomics mrpls.cv
logistic_pls_cv_training <- function(by, data, ncores, file){
  if (!missing(file) & file.exists(paste0(file, "lplscv.Rdata"))) {
    print(paste0(file, "_lplscv.Rdata found. skipping training step..."))
    load(paste0(file, "_lplscv.Rdata"))
  } else {
    fit <- plsgenomics::rpls.cv(
      Xtrain = data,
      Ytrain = by$by,
      LambdaRange = signif(
        logspace(d1 = -2, d2 = 3, n = 21), digits = 3
      ),
      ncompMax = 3,
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
    fit, data, data_train, ncores, file) {
  logistic_pls_cv_classification(fit, data, data_train, ncores, file)
}

#' @importFrom plsgenomics mrpls
logistic_pls_cv_classification <- function(
    fit, data, data_train, ncores, file) {
  if (!missing(file) & file.exists(paste0(file, "_lpls.Rdata"))) {
    print(paste0(file, "_lpls.Rdata found. skipping classification step..."))
    load(paste0(file, "_lpls.Rdata"))
  } else {
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


#
# #' @import plsgenomics
# spls_training <- function(
#   by, data, ncores = 1, file, force, skip, randomize) {
#   data <- as.matrix(data)
#   is_group <- factorize(by)
#   groups_names <- levels(as.factor(by))
#   is_group <- as.numeric(is_group) - 1
#   print("sparse pls tuning")
#   print(table(is_group))
#   model_tunning <- logit.spls.cv(
#     X = data,
#     Y = is_group,
#     lambda.ridge.range = signif(logspace(d1 = -2, d2 = 3, n = 21), digits = 3),
#     lambda.l1.range = seq(0.05,0.95,by = 0.1),
#     ncomp.range = c(1:min(4,ncol(data) - 1)),
#     adapt = FALSE,
#     maxIter = 100,
#     svd.decompose = FALSE,
#     return.grid = TRUE,
#     ncores = ncores,
#     nfolds = 5,
#     center.X = TRUE,
#     scale.X = TRUE,
#     weighted.center = TRUE,
#     nrun = 50,
#     verbose = TRUE)
#   parameters <- list(
#     is_group = is_group,
#     lambda.ridge = model_tunning$lambda.ridge.opt,
#     lambda.l1 = model_tunning$lambda.l1.opt,
#     ncomp = model_tunning$ncomp.opt,
#     groups_names = groups_names
#   )
#   model <- list(
#     model = NA,
#     parameters = parameters,
#     cv_error = model_tunning$cv.grid,
#     fit = model_tunning
#   )
#   if (!missing(file)) {
#     save(model, file = paste0(file, ".Rdata"))
#   }
#   return(model)
# }
#
# #' @import plsgenomics
# stab_spls_training <- function(
#   by, data, ncores = 1, file, force, skip, randomize) {
#   data <- as.matrix(data)
#   is_group <- factorize(by)
#   groups_names <- levels(as.factor(by))
#   is_group <- as.numeric(is_group) - 1
#   print("sparse pls tuning")
#   print(table(is_group))
#   model_tunning <- logit.spls.stab(
#     X = data,
#     Y = is_group,
#     lambda.ridge.range = signif(
#       logspace(d1 = -2, d2 = 3, n = 21), digits = 3
#     ),
#     lambda.l1.range = seq(0.05,0.95,by = 0.05),
#     ncomp = 1,
#     adapt = FALSE,
#     maxIter = 100,
#     svd.decompose = FALSE,
#     ncores = ncores,
#     nfolds = 5,
#     center.X = TRUE,
#     scale.X = TRUE,
#     weighted.center = TRUE,
#     nresamp = 100,
#     verbose = TRUE
#   )
#   if (!missing(file)) {
#     save(model_tunning, file = paste0(file, ".Rdata"))
#   }
#   return(model_tunning)
# }
#
# #' @import plsgenomics
# spls_classification <- function(
#   fit, data, data_train, file) {
#   data <- as.matrix(data)
#   data_train <- as.matrix(data_train)
#   print("sparse pls classification")
#   print("min error")
#   print(min(fit$cv_error$error))
#   model <- logit.spls(
#     Xtrain = data_train,
#     Ytrain = fit$parameters$is_group,
#     lambda.ridge = fit$parameters$lambda.ridge,
#     lambda.l1 = fit$parameters$lambda.l1,
#     ncomp = fit$parameters$ncomp,
#     Xtest = data,
#     adapt = TRUE,
#     maxIter = 100,
#     svd.decompose = FALSE,
#     center.X = TRUE,
#     scale.X = TRUE,
#     weighted.center = TRUE
#   )
#   model$groups_names <- fit$parameters$groups_names
#   model$groups <- model$groups_names[model$hatYtest + 1]
#   if (!missing(file)) {
#     save(model, file = paste0(file, ".Rdata"))
#   }
#   return(model)
# }
#
# ################################################################################
# # sparse logistic pls
#
# #' @import plsgenomics
# mspls_training <- function(
#   by, data, ncores = 1, file, force, skip, randomize) {
#   data <- as.matrix(data)
#   is_group <- factorize(by)
#   groups_names <- levels(as.factor(by))
#   is_group <- as.numeric(is_group) - 1
#   print("sparse multinomial pls tuning")
#   print(table(is_group))
#   model_tunning <- multinom.spls.cv(
#     X = data,
#     Y = is_group,
#     lambda.ridge.range = signif(
#       logspace(d1 = -2, d2 = 3, n = 21), digits = 3
#     ),
#     lambda.l1.range = seq(0.05,0.95,by = 0.1),
#     ncomp.range = c(1:min(4,ncol(data) - 1)),
#     adapt = FALSE,
#     maxIter = 100,
#     svd.decompose = FALSE,
#     return.grid = TRUE,
#     ncores = ncores,
#     nfolds = 5,
#     center.X = TRUE,
#     scale.X = TRUE,
#     weighted.center = TRUE,
#     nrun = 50,
#     verbose = TRUE)
#   parameters <- list(
#     is_group = is_group,
#     lambda.ridge = model_tunning$lambda.ridge.opt,
#     lambda.l1 = model_tunning$lambda.l1.opt,
#     ncomp = model_tunning$ncomp.opt,
#     groups_names = groups_names
#   )
#   model <- list(
#     model = NA,
#     parameters = parameters,
#     cv_error = model_tunning$cv.grid,
#     fit = model_tunning
#   )
#   if (!missing(file)) {
#     save(model, file = paste0(file, ".Rdata"))
#   }
#   return(model)
# }
#
# #' @import plsgenomics
# stab_mspls_training <- function(by, data, ncores = 1, file, force, skip,
#     randomize) {
#   data <- as.matrix(data)
#   is_group <- factorize(by)
#   is_group <- as.numeric(is_group) - 1
#   print("sparse multinomial pls tuning")
#   print(table(is_group))
#   model_tunning <- plsgenomics::multinom.spls.stab(
#     X = data,
#     Y = is_group,
#     lambda.ridge.range = signif(
#       logspace(d1 = -2, d2 = 3, n = 21), digits = 3
#     ),
#     lambda.l1.range = seq(0.05, 0.95, by = 0.05),
#     ncomp.range = c(1:3),
#     svd.decompose = FALSE,
#     ncores = ncores,
#     center.X = TRUE,
#     scale.X = TRUE,
#     weighted.center = TRUE,
#     verbose = TRUE
#   )
#   if (!missing(file)) {
#     save(model_tunning, file = paste0(file, ".Rdata"))
#   }
#   return(model_tunning)
# }
#
# #' @import plsgenomics
# mspls_classification <- function(fit, data, data_train, file) {
#   data <- as.matrix(data)
#   data_train <- as.matrix(data_train)
#   print("sparse multinomial pls classification")
#   print("min error")
#   print(min(fit$cv_error$error))
#   model <- multinom.spls(
#     Xtrain = data_train,
#     Ytrain = fit$parameters$is_group,
#     lambda.ridge = fit$parameters$lambda.ridge,
#     lambda.l1 = fit$parameters$lambda.l1,
#     ncomp = fit$parameters$ncomp,
#     Xtest = data,
#     adapt = TRUE,
#     maxIter = 100,
#     svd.decompose = FALSE,
#     center.X = TRUE,
#     scale.X = TRUE,
#     weighted.center = TRUE
#   )
#   model$groups_names <- fit$parameters$groups_names
#   model$groups <- model$groups_names[model$hatYtest + 1]
#   if (!missing(file)) {
#     save(model, file = paste0(file, ".Rdata"))
#   }
#   return(model)
# }
#
# ################################################################################
# # sparse binomial pls + binomial pls
#
# #' @import plsgenomics
# spls_pls_training <- function(
#   by, data, ncores = 1, file, force, skip = 0, randomize = FALSE) {
#   if (randomize) {
#     classification_sparse_r <- vector("list", 100)
#     classification_sparse_r <- lapply(
#       classification_sparse_r,
#       FUN = function(x, by, data, file, ncores) {
#           r_select <- sample(
#             rownames(data),
#             round(nrow(data) * 0.20)
#           )
#           training_sparse <- stab_spls_training(
#             factorize(by[r_select]),
#             data[r_select,],
#             ncores,
#             paste0(file, "spls_pls_training_a_r"))
#           classification_sparse <- spls_classification(
#             training_sparse,
#             data,
#             data[r_select,],
#             paste0(file, "spls_pls_training_b_r")
#           )
#           retun(
#             list(
#               training_sparse = training_sparse,
#               classification_sparse = classification_sparse,
#               r_select = r_select
#             )
#           )
#         },
#         by,
#         data,
#         file,
#         ncores
#       )
#     return(classification_sparse_r)
#   }
#   if (skip <= 0) {
#     training_sparse <- stab_spls_training(
#       by = by,
#       data = data,
#       ncores = ncores,
#       file = paste0(file, "spls_pls_training_a")
#     )
#     skip <- 1
#   }
#   if (skip <= 1) {
#     if (!exists("training_sparse")) {
#       load(paste0(file, "spls_pls_training_a.Rdata"))
#       training_sparse <- model
#     }
#     classification_sparse <- spls_classification(
#       fit = training_sparse,
#       data = data,
#       data_train = data,
#       file = paste0(file, "spls_pls_training_b")
#     )
#     skip <- 2
#   }
#   if (!exists("classification_sparse")) {
#     load(paste0(file, "spls_pls_training_b.Rdata"))
#     classification_sparse <- model
#   }
#   c_select <- classification_sparse$A
#   print("selected:")
#   print(colnames(data)[c_select])
#   if (!missing(force)) {
#     print("forced:")
#     print(force)
#     if (length(force) >= 1) {
#       if (!is.numeric(force[1])) {
#         force <- which(colnames(data) %in% force)
#       }
#     }
#     c_select <- unique(c(c_select, force))
#     print("selected + forced:")
#     print(colnames(data)[c_select])
#   }
#   if (skip <= 2) {
#     training <- pls_training(
#       by = by,
#       data = data[ ,c_select],
#       ncores = ncores,
#       file = paste0(file, "spls_pls_training_c")
#     )
#   }else {
#     if (!exists("training")) {
#       load(paste0(file, "mspls_mpls_training_c.Rdata"))
#       training <- model
#     }
#   }
#   training$A.full <- c_select
#   training$training_sparse <- training_sparse
#   training$classification_sparse <- classification_sparse
#   save(training, training_sparse, classification_sparse,
#     file = paste0(file, "mspls_mpls_training_all"))
#   return(training)
# }
#
# spls_pls_classification <- function(fit, data, data_train, file) {
#   return(pls_classification(
#     fit = fit,
#     data = data[ ,fit$A.full],
#     data_train = data_train[ ,fit$A.full],
#     file = file)
#   )
# }
#
# ################################################################################
# # sparse multinomial pls + multinomial pls
#
# mspls_mpls_training <- function(by, data, ncores = 1, file, force = c(), skip = 0, randomize = FALSE) {
#   if (skip <= 0) {
#     training_sparse <- stab_mspls_training(
#       by = by,
#       data = data,
#       ncores = ncores,
#       file = paste0(file, "mspls_mpls_training_a")
#     )
#     skip <- 1
#   }
#   if (skip <= 1) {
#     if (!exists("training_sparse")) {
#       load(paste0(file, "mspls_mpls_training_a.Rdata"))
#       training_sparse <- model
#     }
#     classification_sparse <- mspls_classification(
#       fit = training_sparse,
#       data = data,
#       data_train = data,
#       file = paste0(file, "mspls_mpls_training_b")
#     )
#     skip <- 2
#   }
#   if (!exists("classification_sparse")) {
#     load(paste0(file, "mspls_mpls_training_b.Rdata"))
#     classification_sparse <- model
#   }
#   c_select <- classification_sparse$A.full
#   print("selected:")
#   print(colnames(data)[c_select])
#   if (!missing(force)) {
#     print("forced:")
#     print(force)
#     if (length(force) >= 1) {
#       if (!is.numeric(force[1])) {
#         force <- which(colnames(data) %in% force)
#       }
#     }
#     c_select <- unique(c(c_select, force))
#     print("selected + forced:")
#     print(colnames(data)[c_select])
#   }
#   if (skip <= 2) {
#     training <- mpls_training(
#       by, data[,c_select], ncores, paste0(file, "mspls_mpls_training_c")
#     )
#     skip <- 3
#   }else {
#     if (!exists("training")) {
#       load(paste0(file, "mspls_mpls_training_c.Rdata"))
#       training <- model
#     }
#   }
#   training$A.full <- c_select
#   training$training_sparse <- training_sparse
#   training$classification_sparse <- classification_sparse
#   save(training, training_sparse, classification_sparse,
#     file = paste0(file, "mspls_mpls_training_all"))
#   return(training)
# }
#
# mspls_mpls_classification <- function(fit, data, data_train, file) {
#   return(mpls_classification(
#     fit = fit,
#     data = data[ ,fit$A.full],
#     data_train = data_train[ ,fit$A.full],
#     file = file)
#   )
# }
