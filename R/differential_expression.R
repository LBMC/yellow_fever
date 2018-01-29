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
DEA <- function(scd, formula_null, formula_full, b_cells, zi_threshold = 0.9,
    cpus = 4, v = F, folder_name) {
  scd_DEA <- scd$select(b_cells = b_cells)
  scd_DEA <- scd_DEA$select(
    genes = expressed(scd = scd_DEA, zi_threshold = zi_threshold)
  )
  genes_list <- as.list(ERCC(scd_DEA, minus = TRUE))
  names(genes_list) <- ERCC(scd_DEA, minus = TRUE)
  features <- formula_to_features(
    scd = scd_DEA,
    formula_full = formula_full
  )
  results <- lapply_parallel(
    genes_list = genes_list,
    counts = scd_DEA$getcounts,
    features = features,
    formula_null = formula_null,
    formula_full = formula_full,
    cpus = cpus,
    v = v,
    folder_name = folder_name
  )
  results_unlisted <- unlist_results(results)
  return(results_unlisted)
}

#' helper function for DEA for scRNASeq data
#'
#' @param paraload_file paraload table file
#' @param scd object to analyses
#' @param job_DEA_number number of DEA per paraload jobs
#' @return write a paraload ready table with the parameters (as line) on which
#' to run every DEA_paraload() commands
#' @examples
#' \dontrun{
#' DEA_paraload_parameters('results/cell_type/DEA/paraload_file.txt')
#' }
#' @export DEA_paraload_parameters
DEA_paraload_parameters <- function(
  paraload_file,
  scd,
  job_DEA_number = 5,
  formula_null, formula_full, b_cells, zi_threshold = 0.9, cpus, folder_name) {
  job_number <- scd$getngenes / job_DEA_number
  print(job_number)
  scd_DEA <- scd$select(b_cells = b_cells)
  scd_DEA <- scd_DEA$select(genes = expressed(scd_DEA, zi_threshold))
  save(scd_DEA, file = paste0(paraload_file, ".Rdata"))
  parameters <- data.frame(
    analysis = rep("DEA", job_number),
    job_number = rep(job_DEA_number, job_number),
    range_from = seq(
      from = 0,
      to = scd$getngenes - job_DEA_number,
      by = job_DEA_number) + 1,
    range_to = seq(
      from = job_DEA_number,
      to = scd$getngenes,
      by = job_DEA_number),
    formula_null = gsub(" ", "", formula_null),
    formula_full = gsub(" ", "", formula_full),
    folder_name = folder_name,
    scd_path = paste0(paraload_file, ".Rdata"),
    cpus = cpus
  )
  utils::write.table(
    parameters,
    file = paraload_file,
    append = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

#' helper function to run DEA in a script
#'
#' @param scd_file RData file with an scTools data object 'data' saved in it
#' @param Args agument passsed when used within a script
#' @return write a paraload ready table with the parameters (as line) on which
#' to run every DEA() commands
#' @examples
#' \dontrun{
#' DEA_pbs()
#' }
#' @export DEA_pbs
DEA_pbs <- function(Args = commandArgs()) {
  print(getwd())
  if (length(Args) > 4) {
    args <- utils::read.table(Args[6])
    job_number <- as.numeric(args[, 2])
    range_from <- as.vector(args[, 3])
    range_to <- as.vector(args[, 4])
    formula_null <- as.vector(args[, 5])
    formula_full <- as.vector(args[, 6])
    folder_name <- as.vector(args[, 7])
    scd_path <- as.vector(args[, 8])
    cpus <- as.numeric(args[, 9])
    print(date())
    load(scd_path)
    genes <- scd_DEA$getgenes[range_from:range_to]
    DEA <- DEA(
      scd = scd_DEA$select(genes = genes),
      formula_null = formula_null,
      formula_full = formula_full,
      b_cells = rep(TRUE, scd_DEA$getncells),
      cpus = cpus,
      v = T,
      folder_name = folder_name
    )
  }
  utils::write.table(
    genes,
    file = Args[7],
    append = F, row.names = F, col.names = F, quote = F
  )
  print("DEA done for genes")
  print(genes)
}

unlist_results <- function(results){
  results_unlisted <- as.data.frame(do.call(rbind, results))
  results_unlisted$gene <- names(results)
  passed <- !is.na(results_unlisted$pval) & results_unlisted$pval != ""
  results_unlisted[
    !passed, -c(1:3, ncol(results_unlisted))
    ] <- NA
  results_unlisted$pvalue <- as.vector(results_unlisted$pvalue)
  results_unlisted$padj <- NA
  results_unlisted$padj[passed] <- stats::p.adjust(
    results_unlisted$pvalue[passed],
    method = "BH"
  )
  return(results_unlisted)
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
    models_result = models_result,
    v = v)
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
    if (formula_null == formula_full) {
      models_result[["formula_null"]] <- ziNB_fit(
        data = data,
        formula = formulas[["formula_null"]],
        gene_name = gene_name,
        family = family,
        link = link,
        zi = models_result[["is_zi"]],
        v = v
      )
      models_result[["formula_full"]] <- models_result[["formula_null"]]
    } else {
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

#' importFrom glmmADMB glmmadmb
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
      anova(
        models_result[["formula_null"]],
        models_result[["formula_full"]]
      )
    }, error = function(e){
      if (v) {
        print(paste0("error: DEA_LRT for gene ", gene_name))
        print(e)
      }
      return(
        data.frame(
          NoPar = c(NA, NA),
          LogLik = c(NA, NA),
          Df = c(NA, NA),
          Deviance = c(NA, NA),
          "Pr(>Chi)" = c(NA, NA),
          stringsAsFactors = FALSE
        )
      )
    })
    if (!missing(folder_name)) {
      save(
        LRT_result, gene_name,
        file = tmp_file
      )
    }
  }
  if (v) {
    print(paste0(
      "LRT : ", LRT_result[["Pr(>Chi)"]][2], " (", gene_name, ")"
    ))
  }
  return(list(LRT = LRT_result, file = tmp_file))
}

DEA_format <- function(LRT_result, models_result, v) {
  results <- tryCatch({
    LRT <- data.frame(
      null_loglik = LRT_result[["LRT"]][["LogLik"]][1],
      full_loglik = LRT_result[["LRT"]][["LogLik"]][2],
      null_Deviance = LRT_result[["LRT"]][["Deviance"]][1],
      full_Deviance = LRT_result[["LRT"]][["Deviance"]][2],
      null_df = LRT_result[["LRT"]][["Df"]][1],
      full_df = LRT_result[["LRT"]][["Df"]][2],
      null_npar = LRT_result[["LRT"]][["NoPar"]][1],
      full_npar = LRT_result[["LRT"]][["NoPar"]][2],
      pvalue = LRT_result[["LRT"]][["Pr(>Chi)"]][2],
      stringsAsFactors = FALSE
    )
    model_null <- DEA_format_ziNB(
      models_result[["formula_null"]], models_result[["is_zi"]]
    )
    names(model_null) <- paste0(
      "null_",
      names(model_null)
    )
    model_full <- DEA_format_ziNB(
      models_result[["formula_full"]], models_result[["is_zi"]]
    )
    names(model_full) <- paste0(
      "full_",
      names(model_full)
    )
    c(model_null, model_full, unlist(LRT))
  }, error = function(e){
    if (v) {
      print("error: DEA_format")
      print(e)
    }
    return(NA)
  })
  results <- c(
    unlist(data.frame(
      zi = as.vector(models_result[["is_zi"]]),
      LRT_file = as.vector(LRT_result[["file"]]),
      models_file = as.vector(models_result[["file"]]),
      stringsAsFactors = FALSE
    )),
    results
  )
  return(results)
}

DEA_format_ziNB <- function(model, zi) {
  result <- tryCatch({
    unlist(data.frame(
      zi = ifelse(zi, model$pz, NA),
      alpha = model$alpha
    ))
  }, error = function(e){
    return(NA)
  })
  results_b <- tryCatch({
    if ("b" %in% names(model)) {
      b_estimate <- as.data.frame(model[["b"]])
      b_estimate <- unlist(t(model[["b"]]))
      names(b_estimate) <- paste0(
        "fixed_", names(model$b))
      return(b_estimate)
    }
    return(unlist(data.frame(fixed = NA)))
  }, error = function(e){
    return(NA)
  })
  results_S <- tryCatch({
    if ("S" %in% names(model)) {
      S_estimate <- as.data.frame(model[["S"]])
      S_estimate <- unlist(S_estimate)
      names(S_estimate) <- paste0(
        "mixed_", names(model$S))
      return(S_estimate)
    }
    return(unlist(data.frame(mixed = NA)))
  }, error = function(e){
    return(NA)
  })
  result <- c(result, results_b, results_S)
  return(result)
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
    glmmADMB::admbControl(
      impSamp = 1000,
      maxfn = 10000,
      imaxf = 10000,
    )
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
zi_test <- function(data, formula_full, gene_name,
    family = "nbinom1", link = "log", threshold = 0.05, v = F){
  if (v) {
    print(paste0(
      "zi test: zi : ", sum(data$y == 0) / length(data$y), " for ", gene_name
    ))
  }
  if (max(data$y) == 0){
    stop(paste0("error : ", gene_name, " contains only zeros"))
  }
  m1 <- tryCatch({
      glm.nb(y ~ 1,
        data = data,
        control = glm.control(maxit = 10000)
      )
    }, error = function(e){
      if (v) {
        print("zi test: error in NB model")
      }
      return(list(residuals = NA))
    })
  m2 <- tryCatch({
      zeroinfl(y ~ 1,
        data = data,
        dist = "negbin",
        control = zeroinfl.control(maxit = 10000)
      )
    }, error = function(e){
      if (v) {
        print("zi test: error in ziNB model")
      }
      return(list(residuals = NA))
    })
  if (is.na(m1$residuals[1]) & is.na(m2$residuals[1])) {
    if(sum(data$y == 0) < 3){
      if (v) {
        print(paste0(
          "zi test: no zeroinfl detected with error for",
          gene_name, " (#y == 0 < 3)"
        ))
      }
      return(FALSE)
    } else {
      if (v) {
        print(paste0(
          "zi test: zeroinfl detected with error for",
          gene_name, " (#y == 0 >= 3)"
        ))
      }
      return(TRUE)
    }
  }
  if (is.na(m1$residuals[1])) {
    if (v) {
      print(paste0("zi test: zeroinfl detected with error for ",
        gene_name, " (no NB fit)"))
    }
    return(TRUE)
  }
  if (is.na(m2$residuals[1])) {
    if (v) {
      print(paste0("zi test: no zeroinfl detected with error for ",
        gene_name, " (no ziNB fit)"))
    }
    return(FALSE)
  }
  test <- vuong_test(m1 = m1, m2 = m2, v = T)
  if (test < threshold){
    if (v) {
      print(paste0("zi test: zeroinfl detected for ", gene_name, " (vuong test)"))
    }
    return(TRUE)
  }
  if (v) {
    print(paste0("zi test: no zeroinfl detected for ", gene_name, " (vuong test)"))
  }
  return(FALSE)
}

#' importFrom pscl predprob
vuong_test <- function (m1, m2, digits = getOption("digits"), v = F){
    m1y <- m1$y
    m2y <- m2$y
    m1n <- length(m1y)
    m2n <- length(m2y)
    if (m1n == 0 | m2n == 0){
      if (v) {
        stop("Could not extract dependent variables from models.")
      }
    }
    if (m1n != m2n) {
      if (v) {
        stop(paste("Models appear to have different numbers of observations.\n",
            "Model 1 has ", m1n, " observations.\n", "Model 2 has ",
            m2n, " observations.\n", sep = ""))
        }
    }
    if (any(m1y != m2y)) {
      if (v){
        stop(paste("Models appear to have different values on dependent variables.\n"))
      }
    }
    p1 <- predprob(m1)
    p2 <- predprob(m2)
    if (!all(colnames(p1) == colnames(p2))) {
      if (v){
        stop("Models appear to have different values on dependent variables.\n")
      }
    }
    whichCol <- match(m1y, colnames(p1))
    whichCol2 <- match(m2y, colnames(p2))
    if (!all(whichCol == whichCol2)) {
      if (v){
        stop("Models appear to have different values on dependent variables.\n")
      }
    }
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
      if (v){
        cat("NA or numerical zeros or ones encountered in fitted probabilities\n")
        cat(paste("dropping these", sum(bad), "cases, but proceed with caution\n"))
      }
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
            pval[j] <- 1 - pnorm(v[j])
            msg[j] <- "model1 > model2"
        } else {
            pval[j] <- pnorm(v[j])
            msg[j] <- "model2 > model1"
        }
    }
    out <- data.frame(v, msg, format.pval(pval))
    names(out) <- c("Vuong z-statistic", "H_A", "p-value")
    return(as.numeric(as.vector(pval[3])))
}
