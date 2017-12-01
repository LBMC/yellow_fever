#' DEA for scRNASeq data
#'
#' @param scd is a scRNASeq data object
#' @param formula_null null model for the LR test
#' @param formula_full full model for the LR test
#' @param output_file if non empty file to save the classification results
#' @return a list() with the result of the DEA
#' @examples
#' \dontrun {
#' classif <- classification(scd, 500000)
#' }
#' @export DEA
DEA <- function(scd, formula_null, formula_full, b_cells,
  cpus = 4, v = F, folder_name) {
  expressed <- scd$getgenes[colSums(scd$getcounts) > 0]
  scd_DEA <- scd$select(genes = expressed)
  genes_list <- as.list(scd_DEA$select(b_cells = b_cells)$getgenes)
  names(genes_list) <- scd_DEA$select(b_cells = b_cells)$getgenes
  features <- formula_to_features(
    scd = scd_DEA$select(b_cells = b_cells),
    formula_full = formula_full
  )
  results <- lapply_parallel(
    genes_list = genes_list,
    counts = scd$select(b_cells = b_cells)$getcounts,
    features = features,
    formula_null = formula_null,
    formula_full = formula_full,
    cpus = cpus,
    v = v,
    folder_name = folder_name
  )
  results_unlised <- as.data.frame(do.call(rbind, results))
  results_unlised$gene <- names(results)
  passed <- !(is.na(results_unlised$pval))
  results_unlised$padj <- NA
  results_unlised$padj[passed] <- stats::p.adjust(results_unlised$pval[passed])
  return(results_unlised)
}

#' @importFrom parallel mclapply
lapply_parallel <- function(genes_list, counts, features,
    formula_null, formula_full,
    cpus = 1, v, folder_name) {
  results <- list()
  if (cpus > 1) {
    results <- parallel::mclapply(
      X = genes_list,
      FUN  = function(x, counts, features, formula_null, formula_full,
          v, folder_name){
        data <- data.frame(y = round(counts[ ,colnames(counts) %in% x]))
        data <- cbind(features, data)
        DEA_gene(
          data = data,
          formula_null,
          formula_full,
          gene_name = x,
          v = v,
          folder_name = folder_name
        )
      },
      mc.cores = cpus,
      counts = counts,
      features = features,
      formula_null = formula_null,
      formula_full = formula_full,
      v = v,
      folder_name = folder_name
    )
  } else {
    results <- lapply(
      X = genes_list,
      FUN  = function(x, counts, features, formula_null, formula_full,
          v, folder_name){
        data <- data.frame(y = round(counts[ ,colnames(counts) %in% x]))
        data <- cbind(features, data)
        DEA_gene(
          data = data,
          formula_null,
          formula_full,
          gene_name = x,
          v = v,
          folder_name = folder_name
        )
      },
      counts = counts,
      features = features,
      formula_null = formula_null,
      formula_full = formula_full,
      v = v,
      folder_name = folder_name
    )
  }
}

DEA_gene <- function(data, formula_null, formula_full, gene_name,
    family = "nbinom1", link = "log",
    v, folder_name) {
  print(summary(data))
  models_result <- DEA_fit(
    data = data,
    formula_null = formula_null,
    formula_full = formula_full,
    gene_name = gene_name,
    family = family,
    link = link,
    v = v,
    folder_name = folder_name
  )
  LRT_result <- DEA_LRT(
    models_result,
    gene_name = gene_name,
    v = v,
    folder_name = folder_name
  )
  DEA_result <- DEA_format(
    LRT_result = LRT_result,
    models_result = models_result)
  return(DEA_result)
}

DEA_fit <- function(data, formula_null, formula_full, gene_name,
    family = "nbinom1", link = "log",
    v, folder_name) {
  tmp_file <- ""
  if (!missing(folder_name)) {
    system(paste0("mkdir -p ", folder_name, "/DEA_fit/"))
    tmp_file <- paste0(folder_name, "/DEA_fit/fit_", gene_name, ".Rdata")
  }
  if (file.exists(tmp_file)) {
    load(tmp_file)
  } else {
    formulas <- list(
      formula_null = formula_null,
      formula_full = formula_full
    )
    models_result <- list()
    models_result[["file"]] <- tmp_file
    models_result[["is_zi"]] <- zi_test(
      data = data,
      formula_full = formula_full,
      gene_name = gene_name,
      family = family,
      link = link,
      threshold = 0.05,
      v = v)
    for (formula in names(formulas)) {
      models_result[[formula]] <- ziNB_fit(
        data = data,
        formula = formulas[[formula]],
        gene_name = gene_name,
        family = family,
        link = link,
        zi = models_result[["is_zi"]],
        v = v
      )
    }
    if (!missing(folder_name)) {
      save(
        models_result, gene_name, formulas,
        file = tmp_file
      )
    }
  }
  return(models_result)
}

#' importFrom lmtest lrtest
DEA_LRT <- function(models_result, gene_name, v, folder_name) {
  tmp_file <- ""
  LRT_result <- NA
  if (!missing(folder_name)) {
    system(paste0("mkdir -p ", folder_name, "/DEA_LRT/"))
    tmp_file <- paste0(folder_name, "/DEA_LRT/LRT_", gene_name, ".Rdata")
  }
  if (file.exists(tmp_file)) {
    load(tmp_file)
  } else {
    LRT_result <- tryCatch({
      lmtest::lrtest(
        models_result[["formula_full"]],
        models_result[["formula_null"]]
      )
    }, error = function(e){
      if (v) {
        print(paste0("error: DEA_LRT for gene ", gene_name))
        print(e)
      }
      return(NA)
    })
    if (!missing(folder_name)) {
      save(
        LRT_result, gene_name,
        file = tmp_file
      )
    }
  }
  return(list(LRT = LRT_result, file = tmp_file))
}

DEA_format <- function(LRT_result, models_result) {
  results <- tryCatch({
    data.frame(
      zi = models_result[["is_zi"]],
      pvalue = LRT_result[["LRT"]][["Pr(>Chisq)"]][2],
      null_loglik = models_result[["formula_null"]]$loglik,
      null_df = models_result[["formula_null"]]$npar,
      full_loglik = models_result[["formula_full"]]$loglik,
      full_df = models_result[["formula_full"]]$npar,
      LRT_file = LRT_result[["file"]],
      models_file = models_result[["file"]]
    )
  }, error = function(e){
    data.frame(
      zi = models_result[["is_zi"]],
      pvalue = NA,
      null_loglik = NA,
      null_df = NA,
      full_loglik = NA,
      full_df = NA,
      LRT_file = NA,
      models_file = NA
    )
  })
  if (results$zi & !is.na(results$pvalue)) {
    results$null_zi <- models_result[["formula_null"]]$pz
    results$full_zi <- models_result[["formula_full"]]$pz
  } else {
    results$null_zi <- NA
    results$full_zi <- NA
  }
  return(as.vector(results))
}

#' importFrom glmmADMB glmmadmb
ziNB_fit <- function(data, formula, gene_name,
    family = "nbinom1", link = "log", zi = TRUE,
    v) {
  if (v) {
    if (zi) {
      print(paste0(gene_name, " : ziNB : ", formula))
    } else {
      print(paste0(gene_name, " : NB : ", formula))
    }
  }
  model <- tryCatch({
    glmmADMB::glmmadmb(
      as.formula(formula),
      data = data,
      zeroInflation = zi,
      family = family,
      link = link,
      mcmc = FALSE
    )
  }, error = function(e){
    if (v) {
      if (zi) {
        print(paste0(
            "error: ziNB_fit for ", gene_name,
            " with zero-inflation"
        ))
      } else {
        print(paste0("error: ziNB_fit for ", gene_name))
      }
      print(e)
    }
    return(list(residuals = NA))
  })
  return(model)
}

formula_to_features <- function(scd, formula_full){
  features <- c()
  for(feature in colnames(scd$getfeatures)) {
    features <- c(
      features,
      grepl(
        pattern = feature,
        x = formula_full
      )
    )
  }
  features[1] <- TRUE
  features <- apply(
    X = scd$getfeatures[, features],
    MARGIN = 2,
    FUN = as.factor
  )
  return(features)
}

#' @importFrom MASS glm.nb
#' @importFrom pscl zeroinfl
#' @importFROM pscl vuong
zi_test <- function(data, formula_full, gene_name,
    family = "nbinom1", link = "log", threshold = 0.05, v = F){
  if(max(data$y) == 0){
    stop(paste0("error : ", gene_name, " contains only zeros"))
  }
  m1 <- ziNB_fit(
    data = data,
    formula = formula_full,
    gene_name = gene_name,
    family = family,
    link = link,
    zi = FALSE,
    v = v
  )
  m2 <- ziNB_fit(
    data = data,
    formula = formula_full,
    gene_name = gene_name,
    family = family,
    link = link,
    zi = TRUE,
    v = v
  )
  if (is.na(m1$residuals[1]) & is.na(m2$residuals[1])) {
    if(sum(data$y == 0) < 3){
      if (v) {
        print(paste0(
          "no zeroinfl detected with error for",
          gene_name, " (#y == 0 < 3)"
        ))
      }
      return(FALSE)
    } else {
      if (v) {
        print(paste0(
          "zeroinfl detected with error for",
          gene_name, " (#y == 0 >= 3)"
        ))
      }
      return(TRUE)
    }
  }
  if (is.na(m1$residuals[1])) {
    if (v) {
      print(paste0("zeroinfl detected with error for ",
        gene_name, " (no NB fit)"))
    }
    return(TRUE)
  }
  if (is.na(m2$residuals[1])) {
    if (v) {
      print(paste0("no zeroinfl detected with error for ",
        gene_name, " (no ziNB fit)"))
    }
    return(FALSE)
  }
  test <- vuong_test(m1 = m1, m2 = m2, v = T)
  if(test < threshold){
    if (v) {
      print(paste0("zeroinfl detected for ", gene_name, " (vuong test)"))
    }
    return(TRUE)
  }
  if (v) {
    print(paste0("no zeroinfl detected for ", gene_name, " (vuong test)"))
  }
  return(FALSE)
}

vuong_test <- function (m1, m2, digits = getOption("digits"), v){
  m1y <- m1$frame$y
  m2y <- m2$frame$y
  m1n <- length(m1y)
  m2n <- length(m2y)
  if (m1n == 0 | m2n == 0)
    stop("Could not extract dependent variables from models.")
  if (m1n != m2n)
    stop(paste("Models appear to have different numbers of observations.\n",
        "Model 1 has ", m1n, " observations.\n", "Model 2 has ",
        m2n, " observations.\n", sep = ""))
  if (any(m1y != m2y))
    stop(paste("Models appear to have different values on dependent variables.\n"))
  p1 <- predprob(m1, zi = FALSE)
  p2 <- predprob(m2, zi = TRUE)
  if (!all(colnames(p1) == colnames(p2)))
    stop("Models appear to have different values on dependent variables.\n")
  whichCol <- match(m1y, colnames(p1))
  whichCol2 <- match(m2y, colnames(p2))
  if (!all(whichCol == whichCol2))
    stop("Models appear to have different values on dependent variables.\n")
  m1p <- rep(NA, m1n)
  m2p <- rep(NA, m2n)
  for (i in 1:m1n) {
      m1p[i] <- p1[i, whichCol[i]]
      m2p[i] <- p2[i, whichCol[i]]
  }
  k1 <- length(coef(m1))
  k2 <- length(coef(m2))
  lm1p <- log(m1p)
  lm2p <- log(m2p)
  m <- lm1p - lm2p
  bad1 <- is.na(lm1p) | is.nan(lm1p) | is.infinite(lm1p)
  bad2 <- is.na(lm2p) | is.nan(lm2p) | is.infinite(lm2p)
  bad3 <- is.na(m) | is.nan(m) | is.infinite(m)
  bad <- bad1 | bad2 | bad3
  neff <- sum(!bad)
  if (any(bad)) {
    cat("NA or numerical zeros or ones encountered in fitted probabilities\n")
    cat(paste("dropping these", sum(bad), "cases, but proceed with caution\n"))
  }
  aic.factor <- (k1 - k2)/neff
  bic.factor <- (k1 - k2)/(2 * neff) * log(neff)
  v <- rep(NA, 3)
  arg1 <- matrix(m[!bad], nrow = neff, ncol = 3, byrow = FALSE)
  arg2 <- matrix(c(0, aic.factor, bic.factor), nrow = neff,
      ncol = 3, byrow = TRUE)
  num <- arg1 - arg2
  s <- apply(num, 2, sd)
  numsum <- apply(num, 2, sum)
  v <- numsum/(s * sqrt(neff))
  names(v) <- c("Raw", "AIC-corrected", "BIC-corrected")
  pval <- rep(NA, 3)
  msg <- rep("", 3)
  for (j in 1:3) {
    if (v[j] > 0) {
      pval[j] <- 1
      msg[j] <- "model1 < model2"
    }
    else {
      pval[j] <- pnorm(v[j])
      msg[j] <- "model1 < model2"
    }
  }
  out <- data.frame(v, msg, format.pval(pval))
  names(out) <- c("Vuong z-statistic", "H_A", "p-value")
  return(as.numeric(as.vector(pval[3])))
}

#' importFrom extraDistr dzinb
predprob <- function(obj, zi = F){
  yhat <- predict(obj, type="response")
  y <- obj$frame$y
  yUnique <- c(0:max(y))
  nUnique <- length(yUnique)
  p <- matrix(NA,length(yhat),nUnique)
  dimnames(p) <- list(NULL,yUnique)
  for(i in 1:nUnique){
    p[, i] <- dnbinom(
      mu = yhat,
      size = obj$alpha,
      x = yUnique[i]
    )
    if (zi) {
      p[yhat == 0, i] <- obj$pz + (1 - obj$pz) * p[yhat == 0, i]
    }
  }
  return(p)
}
