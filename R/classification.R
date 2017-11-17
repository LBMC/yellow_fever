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
  training,
  ncores = 4,
  algo = "mspls_mpls",
  output_file = "") {

  algo_training <- paste0(algo, "_training")
  to_train_on <- !is.na(scd$getfeature(feature))
  training <- algo_training(
    by = scd$select(b_cells = to_train_on)$getfeature(feature),
    data = scd$select(b_cells = to_train_on)$getcounts,
    ncores = ncores,
    file = paste0(output_file, "_training")
  )

  algo_classification <- paste0(algo, "_classification")
  classification <- algo_classification(
    fit = training,
    data = scd$select(b_cells = !to_train_on)$getcounts,
    data_train = scd$select(b_cells = to_train_on)$getcounts,
    file = paste0(output_file, "_classification")
  )
  return(list(
    training = training,
    classification = classification,
    scd = scd,
    feature = feature
  ))
}

################################################################################
# sparse pls
logspace <- function( d1, d2, n) exp(log(10) * seq(d1, d2, length.out = n))

#' @import plsgenomics
spls_training <- function(
  by, data, ncores = 1, file, force, skip, randomize) {
  data <- as.matrix(data)
  is_group <- factorize(by)
  groups_names <- levels(as.factor(by))
  is_group <- as.numeric(is_group) - 1
  print("sparse pls tuning")
  print(table(is_group))
  model_tunning <- rirls.spls.tune(
    X = data,
    Y = is_group,
    lambda.ridge.range = signif(logspace(d1 <- -2, d2 <- 3, n = 10), digits = 3),
    lambda.l1.range = seq(0.05,0.95,by = 0.1),
    ncomp.range = c(1:min(4,ncol(data) - 1)),
    adapt = FALSE,
    maxIter = 100,
    svd.decompose = FALSE,
    return.grid = TRUE,
    ncores = ncores,
    nfolds = 5,
    center.X = TRUE,
    scale.X = TRUE,
    weighted.center = TRUE,
    nrun = 50,
    verbose = TRUE)
  parameters <- list(
    is_group = is_group,
    lambda.ridge = model_tunning$lambda.ridge.opt,
    lambda.l1 = model_tunning$lambda.l1.opt,
    ncomp = model_tunning$ncomp.opt,
    groups_names = groups_names
  )
  model <- list(
    model = NA,
    parameters = parameters,
    cv_error = model_tunning$cv.grid,
    fit = model_tunning
  )
  if (!missing(file)) {
    save(model, file = paste0(file, ".RData"))
  }
  return(model)
}

#' @import plsgenomics
stab_spls_training <- function(
  by, data, ncores = 1, file, force, skip, randomize) {
  data <- as.matrix(data)
  is_group <- factorize(by)
  groups_names <- levels(as.factor(by))
  is_group <- as.numeric(is_group) - 1
  print("sparse pls tuning")
  print(table(is_group))
  model_tunning <- rirls.spls.stab(
    X = data,
    Y = is_group,
    lambda.ridge.range = signif(
      logspace(d1 <- -2, d2 <- 3, n = 10), digits = 3
    ),
    lambda.l1.range = seq(0.05,0.95,by = 0.05),
    ncomp = 1,
    adapt = FALSE,
    maxIter = 100,
    svd.decompose = FALSE,
    ncores = ncores,
    nfolds = 5,
    center.X = TRUE,
    scale.X = TRUE,
    weighted.center = TRUE,
    nresamp = 100,
    verbose = TRUE
  )
  if (!missing(file)) {
    save(model_tunning, file = paste0(file, ".RData"))
  }
  return(model_tunning)
}

#' @import plsgenomics
spls_classification <- function(
  fit, data, data_train, file) {
  data <- as.matrix(data)
  data_train <- as.matrix(data_train)
  print("sparse pls classification")
  print("min error")
  print(min(fit$cv_error$error))
  model <- rirls.spls(
    Xtrain = data_train,
    Ytrain = fit$parameters$is_group,
    lambda.ridge = fit$parameters$lambda.ridge,
    lambda.l1 = fit$parameters$lambda.l1,
    ncomp = fit$parameters$ncomp,
    Xtest = data,
    adapt = TRUE,
    maxIter = 100,
    svd.decompose = FALSE,
    center.X = TRUE,
    scale.X = TRUE,
    weighted.center = TRUE
  )
  model$groups_names <- fit$parameters$groups_names
  model$groups <- model$groups_names[model$hatYtest + 1]
  if (!missing(file)) {
    save(model, file = paste0(file, ".RData"))
  }
  return(model)
}

################################################################################
# sparse multinomial pls

#' @import plsgenomics
mspls_training <- function(
  by, data, ncores = 1, file, force, skip, randomize) {
  data <- as.matrix(data)
  is_group <- factorize(by)
  groups_names <- levels(as.factor(by))
  is_group <- as.numeric(is_group) - 1
  print("sparse multinomial pls tuning")
  print(table(is_group))
  model_tunning <- m.rirls.spls.tune(
    X = data,
    Y = is_group,
    lambda.ridge.range = signif(
      logspace(d1 <- -2, d2 <- 3, n = 10), digits = 3
    ),
    lambda.l1.range = seq(0.05,0.95,by = 0.1),
    ncomp.range = c(1:min(4,ncol(data) - 1)),
    adapt = FALSE,
    maxIter = 100,
    svd.decompose = FALSE,
    return.grid = TRUE,
    ncores = ncores,
    nfolds = 5,
    center.X = TRUE,
    scale.X = TRUE,
    weighted.center = TRUE,
    nrun = 50,
    verbose = TRUE)
  parameters <- list(
    is_group = is_group,
    lambda.ridge = model_tunning$lambda.ridge.opt,
    lambda.l1 = model_tunning$lambda.l1.opt,
    ncomp = model_tunning$ncomp.opt,
    groups_names = groups_names
  )
  model <- list(
    model = NA,
    parameters = parameters,
    cv_error = model_tunning$cv.grid,
    fit = model_tunning
  )
  if (!missing(file)) {
    save(model, file = paste0(file, ".RData"))
  }
  return(model)
}

#' @import plsgenomics
stab_mspls_training <- function(by, data, ncores = 1, file, force, skip, randomize) {
  data <- as.matrix(data)
  is_group <- factorize(by)
  groups_names <- levels(as.factor(by))
  is_group <- as.numeric(is_group) - 1
  print("sparse multinomial pls tuning")
  print(table(is_group))
  model_tunning <- m.rirls.spls.stab(
    X = data,
    Y = is_group,
    lambda.ridge.range = signif(
      logspace(d1 <- -2, d2 <- 3, n = 10), digits = 3
    ),
    lambda.l1.range = seq(0.05,0.95,by = 0.05),
    ncomp = 1,
    adapt = FALSE,
    maxIter = 100,
    svd.decompose = FALSE,
    ncores = ncores,
    nfolds = 5,
    center.X = TRUE,
    scale.X = TRUE,
    weighted.center = TRUE,
    nresamp = 100,
    verbose = TRUE
  )
  if (!missing(file)) {
    save(model_tunning, file = paste0(file, ".RData"))
  }
  return(model_tunning)
}

#' @import plsgenomics
mspls_classification <- function(fit, data, data_train, file) {
  data <- as.matrix(data)
  data_train <- as.matrix(data_train)
  print("sparse multinomial pls classification")
  print("min error")
  print(min(fit$cv_error$error))
  model <- m.rirls.spls(
    Xtrain = data_train,
    Ytrain = fit$parameters$is_group,
    lambda.ridge = fit$parameters$lambda.ridge,
    lambda.l1 = fit$parameters$lambda.l1,
    ncomp = fit$parameters$ncomp,
    Xtest = data,
    adapt = TRUE,
    maxIter = 100,
    svd.decompose = FALSE,
    center.X = TRUE,
    scale.X = TRUE,
    weighted.center = TRUE
  )
  model$groups_names <- fit$parameters$groups_names
  model$groups <- model$groups_names[model$hatYtest + 1]
  if (!missing(file)) {
    save(model, file = paste0(file, ".RData"))
  }
  return(model)
}

################################################################################
# sparse binomial pls + binomial pls

#' @import plsgenomics
spls_pls_training <- function(
  by, data, ncores = 1, file, force, skip = 0, randomize = FALSE) {
  if (randomize) {
    classification_sparse_r <- vector("list", 100)
    classification_sparse_r <- lapply(
      classification_sparse_r,
      FUN = function(x, by, data, file, ncores) {
          r_select <- sample(
            rownames(data),
            round(nrow(data) * 0.20)
          )
          training_sparse <- spls_training(
            factorize(by[r_select]),
            data[r_select,],
            ncores,
            paste0(file, "spls_pls_training_a_r"))
          classification_sparse <- spls_classification(
            training_sparse,
            data,
            data[r_select,],
            paste0(file, "spls_pls_training_b_r")
          )
          retun(
            list(
              training_sparse = training_sparse,
              classification_sparse = classification_sparse,
              r_select = r_select
            )
          )
        },
        by,
        data,
        file,
        ncores
      )
    return(classification_sparse_r)
  }
  if (skip <= 0) {
    training_sparse <- spls_training(
      by = by,
      data = data,
      ncores = ncores,
      file = paste0(file, "spls_pls_training_a")
    )
    skip <- 1
  }
  if (skip <= 1) {
    if (!exists("training_sparse")) {
      load(paste0(file, "spls_pls_training_a.RData"))
      training_sparse <- model
    }
    classification_sparse <- spls_classification(
      fit = training_sparse,
      data = data,
      data_train = data,
      file = paste0(file, "spls_pls_training_b")
    )
    skip <- 2
  }
  if (!exists("classification_sparse")) {
    load(paste0(file, "spls_pls_training_b.RData"))
    classification_sparse <- model
  }
  c_select <- classification_sparse$A
  print("selected:")
  print(colnames(data)[c_select])
  if (!missing(force)) {
    print("forced:")
    print(force)
    if (length(force) >= 1) {
      if (!is.numeric(force[1])) {
        force <- which(colnames(data) %in% force)
      }
    }
    c_select <- unique(c(c_select, force))
    print("selected + forced:")
    print(colnames(data)[c_select])
  }
  if (skip <= 2) {
    training <- pls_training(
      by = by,
      data = data[,c_select],
      ncores = ncores,
      file = paste0(file, "spls_pls_training_c")
    )
  }else {
    if (!exists("training")) {
      load(paste0(file, "mspls_mpls_training_c.RData"))
      training <- model
    }
  }
  training$A.full <- c_select
  training$training_sparse <- training_sparse
  training$classification_sparse <- classification_sparse
  save(training, training_sparse, classification_sparse,
    file = paste0(file, "mspls_mpls_training_all"))
  return(training)
}

spls_pls_classification <- function(fit, data, data_train, file) {
  return(pls_classification(
    fit = fit,
    data = data[,fit$A.full],
    data_train = data_train[,fit$A.full],
    file = file)
  )
}

################################################################################
# sparse multinomial pls + multinomial pls

mspls_mpls_training <- function(by, data, ncores = 1, file, force = c(), skip = 0, randomize = FALSE) {
  if (skip <= 0) {
    training_sparse <- mspls_training(
      by = by,
      data = data,
      ncores = ncores,
      file = paste0(file, "mspls_mpls_training_a")
    )
    skip <- 1
  }
  if (skip <= 1) {
    if (!exists("training_sparse")) {
      load(paste0(file, "mspls_mpls_training_a.RData"))
      training_sparse <- model
    }
    classification_sparse <- mspls_classification(
      fit = training_sparse,
      data = data,
      data_train = data,
      file = paste0(file, "mspls_mpls_training_b")
    )
    skip <- 2
  }
  if (!exists("classification_sparse")) {
    load(paste0(file, "mspls_mpls_training_b.RData"))
    classification_sparse <- model
  }
  c_select <- classification_sparse$A.full
  print("selected:")
  print(colnames(data)[c_select])
  if (!missing(force)) {
    print("forced:")
    print(force)
    if (length(force) >= 1) {
      if (!is.numeric(force[1])) {
        force <- which(colnames(data) %in% force)
      }
    }
    c_select <- unique(c(c_select, force))
    print("selected + forced:")
    print(colnames(data)[c_select])
  }
  if (skip <= 2) {
    training <- mpls_training(by, data[,c_select], ncores, paste0(file, "mspls_mpls_training_c"))
    skip <- 3
  }else {
    if (!exists("training")) {
      load(paste0(file, "mspls_mpls_training_c.RData"))
      training <- model
    }
  }
  training$A.full <- c_select
  training$training_sparse <- training_sparse
  training$classification_sparse <- classification_sparse
  save(training, training_sparse, classification_sparse,
    file = paste0(file, "mspls_mpls_training_all"))
  return(training)
}

mspls_mpls_classification <- function(fit, data, data_train, file) {
  return(mpls_classification(
    fit = fit,
    data = data[,fit$A.full],
    data_train = data_train[,fit$A.full],
    file = file)
  )
}
