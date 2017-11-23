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
#' @importFrom parallel mclapply
#' @export DEA
DEA <- function(scd, formula_null, formula_full, b_cells,
  cpus = 4, v = F, tmp_folder) {
  genes_list <- as.list(scd$select(b_cells = b_cells)$getgenes)
  names(genes_list) <- scd$select(b_cells = b_cells)$getgenes
  features <- formula_to_features(
    scd = scd$select(b_cells = b_cells),
    formula_full = formula_full
  )
  # parallel::mclapply(
  lapply(
    X = genes_list,
    FUN  = function(x, counts, features, formula_null, formula_full, v, tmp_folder){
      data <- data.frame(y = round(counts[ ,colnames(counts) %in% x]))
      data <- cbind(features, data)
      DEA_gene(
        data = data,
        formula_null,
        formula_full,
        gene_name = x,
        v = v,
        tmp_folder = tmp_folder
      )
    },
    # mc.cores = cpus,
    counts = scd$select(b_cells = b_cells)$getcounts,
    features = features,
    formula_null = formula_null,
    formula_full = formula_full,
    v = v,
    tmp_folder = tmp_folder
  )
}

DEA_gene <- function(data, formula_null, formula_full, gene_name,
    family = "nbinom1", link = "log",
    v, tmp_folder) {
  hist(data$y, breaks = length(data$y), main = gene_name)
  models_result <- DEA_fit(
    data = data,
    formula_null = formula_null,
    formula_full = formula_full,
    gene_name = gene_name,
    family = family,
    link = link,
    v = v,
    tmp_folder = tmp_folder
  )
  LRT_result <- DEA_LRT(
    models_result,
    gene_name = gene_name,
    v = v,
    tmp_folder = tmp_folder
  )
  DEA_results <- DEA_format(
    LRT_result = LRT_result,
    models_result = models_result)
  return(DEA_results)
}

DEA_fit <- function(data, formula_null, formula_full, gene_name,
    family = "nbinom1", link = "log",
    v, tmp_folder) {
  tmp_file <- ""
  if (!missing(tmp_folder)) {
    tmp_file <- paste0(tmp_folder, "/DEA_fit_", gene_name, ".Rdata")
  }
  if (file.exists(tmp_file)) {
    load(tmp_file)
  } else {
    formulas <- list(
      formula_null = formula_null,
      formula_full = formula_full
    )
    models_result <- list()
    if (zi_test(data, gene_name = gene_name, v = v)) {
      for (formula in names(formulas)) {
        models_result[[formula]] <- ziNB_fit(
          data = data,
          formula = formulas[[formula]],
          gene_name = gene_name,
          family = family,
          link = link,
          v = v
        )
      }
    } else {
      for (formula in names(formulas)) {
        models_result[[formula]] <- NB_fit(
          data = data,
          formula = formulas[[formula]],
          gene_name = gene_name,
          v = v
        )
      }
    }
    if (!missing(tmp_folder)) {
      save(
        models_result, gene_name, formulas,
        file = tmp_file
      )
    }
  }
  return(models_result)
}

#' importFrom lmtest lrtest
DEA_LRT <- function(models_result, gene_name, v, tmp_folder) {
  tmp_file <- ""
  LRT_result <- NULL
  if (!missing(tmp_folder)) {
    tmp_file <- paste0(tmp_folder, "/DEA_LRT_", gene_name, ".Rdata")
  }
  if (file.exists(tmp_file)) {
    load(tmp_file)
  } else {
    LRT_results <- tryCatch({
      if (class(models_result[["formula_null"]]) == "glmerMod") {
        stats::anova(
          models_result[["formula_null"]],
          models_result[["formula_full"]]
        )
      } else {
        lmtest::lrtest(
          models_result[["formula_full"]],
          models_result[["formula_null"]]
        )
      }
    }, error = function(e){
      if (v) {
        print(paste0("error: DEA_LRT for gene ", gene_name))
        print(e)
      }
      return(rep(NA, 8))
    })
    if (!missing(tmp_folder)) {
      save(
        LRT_results, gene_name,
        file = tmp_file
      )
    }
  }
  return(LRT_result)
}

DEA_format <- function(LRT_result, models_result) {

}

#' importFrom glmmADMB glmmadmb
ziNB_fit <- function(data, formula, gene_name,
    family = "nbinom1", link = "log",
    v) {
  model <- tryCatch({
    glmmADMB::glmmadmb(
      as.formula(formula),
      data = data,
      zeroInflation = TRUE,
      family = family,
      link = link,
      mcmc = FALSE
    )
  }, error = function(e){
    if (v) {
      print(paste0("error: ziNB_fit for gene ", gene_name))
      print(e)
    }
    tryCatch({
      glmmADMB::glmmadmb(
        as.formula(fornula),
        data = data,
        zeroInflation = TRUE,
        family = family,
        link = link,
        mcmc = TRUE)
    }, error = function(e){
      if (v) {
        print(paste0("error: ziNB_fit with mcmc for ", gene_name))
        print(e)
      }
      return(list(residuals = NA))
    })
  })
  return(model)
}

#' importFrom lme4 glmer.nb
NB_fit <- function(data, formula, gene_name,
    v) {
  print(formula)
  model <- tryCatch({
    glmer.nb(
      as.formula(formula),
      data = data,
      control = glmerControl(optCtrl = list(maxfun = 10000)))
  }, error = function(e){
    if (v) {
      print(paste0("error: NB_fit for gene ", gene_name))
      print(e)
    }
      tryCatch({
        glmer.nb(
          as.formula(formula),
          data = data,
          control = glmerControl(
            optimizer="bobyqa",
            optCtrl = list(maxfun = 10000)
          )
        )
      }, error = function(e){
        if (v) {
          print(paste0("error: NB_fit with bobyqa for gene ", gene_name))
          print(e)
        }
        return(list(residuals = NA))
      })
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
zi_test <- function(data, threshold=0.05, v = F, gene_name){
  if(max(data$y) == 0){
    exit(paste0("error : ", gene_name, " contains only zeros"))
  }
  if(sum(data$y == 0) < 3){
    if (v) {
      print("no zeroinfl detected (#y == 0 < 3)")
    }
    return(FALSE)
  }
  tryCatch({
    m1 <- tryCatch({
        MASS::glm.nb(y ~ 1,
               data = data)
      }, error = function(e){
        return(TRUE)
      })
    m2 <- tryCatch({
        pscl::zeroinfl(y ~ 1,
                 data = data,
                 dist = "negbin")
      }, error = function(e){
        return(TRUE)
      })
    if(is.atomic(m1) | is.atomic(m2)){
      if (v) {
        print("zeroinfl detected with error")
      }
      return(TRUE)
    }
    test <- vuong_test(m1 = m1, m2 = m2, v = v)
    if(test < threshold){
      if (v) {
        print("zeroinfl detected (vuong test)")
      }
      return(TRUE)
    }
    if (v) {
      print("no zeroinfl detected (vuong test)")
    }
    return(FALSE)
  }, error = function(e){
    if (v) {
      print("zeroinfl detected with error (vuong test)")
    }
    return(TRUE)
  })
}

vuong_test <- function (m1, m2, digits = getOption("digits"), v){
  m1y <- m1$y
  m2y <- m2$y
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
  p1 <- predprob(m1)
  p2 <- predprob(m2)
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
