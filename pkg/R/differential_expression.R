#' DEA for scRNASeq data
#'
#' @param scd is a scRNASeq data object
#' @param formula_null null model for the LR test
#' @param formula_full full model for the LR test
#' @param folder_name folder to write down intermendiate DEA results
#' @param continuous vector of feature names that should be keep continuous (not
#' converted into factors
#' @param b_cells boolean vector of cells to keep in scd
#' @param zi_threshold (default: 0.9) work with genes with less that 0.9 zeros
#' @param zero_threshold (default: 5) counts lesser than this threshold will be
#' considered as zeros for the binomial regression
#' @param model_family (default: "nbinom1") model to use for the DEA in
#' c("nbinom1", "binomial")
#' @param cpus (default: 4) number of cpus to use
#' @param v (default: F) print stuff
#' @return a list() with the result of the DEA
#' @examples
#' \dontrun {
#'
#'   b_cells <- scd$getfeature("QC_good") %in% T &
#'     !is.na(scd$getfeature("DEA_cell_type")) &
#'     scd$getfeature("day") %in% day
#'   mbatch_DEA_cell_type_DEA <- DEA(
#'     scd = scd,
#'     formula_null = "y ~ (1|batch)",
#'     formula_full = "y ~ (1|batch) + DEA_cell_type",
#'     b_cells = b_cells,
#'     cpus = 16,
#'     v = T,
#'     folder_name = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA")
#'   )
#'
#' }
#' @export DEA
DEA <- function(scd, formula_null, formula_full, b_cells, zi_threshold = 0.9,
    model_family = "nbinom1", cpus = 4, v = F, folder_name,
    zero_threshold = 5, continuous = c()) {
  scd_DEA <- scd$select(b_cells = b_cells)
  scd_DEA <- scd_DEA$select(
    genes = expressed(scd = scd_DEA, zi_threshold = zi_threshold)
  )
  genes_list <- as.list(ERCC(scd_DEA, minus = TRUE))
  names(genes_list) <- ERCC(scd_DEA, minus = TRUE)
  features <- formula_to_features(
    scd = scd_DEA,
    formula_full = formula_full,
    continuous = continuous
  )
  results <- lapply_parallel(
    genes_list = genes_list,
    counts = scd_DEA$getcounts,
    features = features,
    formula_null = formula_null,
    formula_full = formula_full,
    model_family = model_family,
    zero_threshold = zero_threshold,
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
  expected_size <- max( unlist(lapply(results, FUN = length)) )
  results_expected <- lapply(results, FUN = function(x, expected_size) {
    if (length(x) == expected_size) {
      return(x)
    }
  },
  expected_size)[[1]]
  results_expected[names(results_expected)] <- NA
  results <- lapply(results, FUN = function(x, expected_size, expected_res) {
      x_size <- length(x)
      if (expected_size > x_size) {
        expected_res[names(x)] <- x
        x <- expected_res
      }
      return(x)
    },
    expected_size,
    expected_res = results_expected
  )
  results_unlisted <- as.data.frame(do.call(rbind, results))
  results_unlisted$gene <- names(results)
  results_unlisted$pvalue <- as.numeric(as.vector(results_unlisted$pvalue))
  results_unlisted$padj <- stats::p.adjust(
    results_unlisted$pvalue,
    method = "BH"
  )
  return(results_unlisted)
}

#' @importFrom parallel mclapply
lapply_parallel <- function(genes_list, counts, features,
    formula_null, formula_full, model_family = "nbinom1", zero_threshold = 5,
    cpus = 1, v, folder_name) {
  results <- list()
  if (cpus > 1) {
    results <- parallel::mclapply(
      X = genes_list,
      FUN  = function(x, counts, features, formula_null, formula_full,
          model_family, zero_threshold, v, folder_name){
        data <- data.frame(y = round(counts[ ,colnames(counts) %in% x]))
        data <- cbind(features, data)
        if (model_family %in% "binomial") {
          data$y <- as.numeric(data$y > zero_threshold)
        }
        DEA_gene(
          data = data,
          formula_null = formula_null,
          formula_full = formula_full,
          model_family = model_family,
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
      model_family = model_family,
      zero_threshold = zero_threshold,
      v = v,
      folder_name = folder_name
    )
  } else {
    results <- lapply(
      X = genes_list,
      FUN  = function(x, counts, features, formula_null, formula_full,
          v, model_family, zero_threshold, folder_name){
        data <- data.frame(y = round(counts[ ,colnames(counts) %in% x]))
        data <- cbind(features, data)
        if (model_family %in% "binomial") {
          data$y <- as.numeric(data$y > zero_threshold)
        }
        DEA_gene(
          data = data,
          formula_null,
          formula_full,
          gene_name = x,
          model_family = model_family,
          v = v,
          folder_name = folder_name
        )
      },
      counts = counts,
      features = features,
      formula_null = formula_null,
      formula_full = formula_full,
      model_family = model_family,
      zero_threshold = zero_threshold,
      v = v,
      folder_name = folder_name
    )
  }
}

DEA_gene <- function(data, formula_null, formula_full, gene_name,
    model_family = "nbinom1", link = "log",
    v, folder_name) {
  models_result <- DEA_fit(
    data = data,
    formula_null = formula_null,
    formula_full = formula_full,
    gene_name = gene_name,
    model_family = model_family,
    link = link,
    v = v,
    folder_name = folder_name
  )
  LRT_result <- DEA_LRT(
    models_result,
    gene_name = gene_name,
    v = v,
    folder_name = folder_name,
    model_family = model_family
  )
  DEA_result <- DEA_format(
    LRT_result = LRT_result,
    models_result = models_result,
    model_family = model_family,
    v = v)
  return(DEA_result)
}

DEA_fit <- function(data, formula_null, formula_full, gene_name,
    model_family = "nbinom1", link = "log",
    v, folder_name) {
  tmp_file <- ""
  if (!missing(folder_name)) {
    system(paste0("mkdir -p ", folder_name, "/DEA_fit/"))
    tmp_file <- paste0(folder_name, "/DEA_fit/fit_", gene_name, ".Rdata")
    system(paste0("mkdir -p ", folder_name, "/DEA_LRT/"))
    tmp_file_LRT <- paste0(folder_name, "/DEA_LRT/LRT_", gene_name, ".Rdata")
  }
  formulas <- list(
    formula_null = formula_null,
    formula_full = formula_full
  )
  models_result <- list()
  for (formula in names( formulas )) {
    models_result[[formula]] <- list(residuals = NA)
  }
  if (file.exists(tmp_file)) {
    if (v) {
      print(paste0("fit cache found for ", gene_name))
    }
    load(tmp_file)
  }
  if (!("file" %in% names(models_result))) {
    models_result[["file"]] <- tmp_file
  }
  if (!("is_zi" %in% names(models_result)) & model_family %in% "nbinom1") {
    models_result[["is_zi"]] <- zi_test(
      data = data,
      gene_name = gene_name,
      model_family = model_family,
      link = link,
      threshold = 0.05,
      v = v)
  }
  try_left <- 1
  FUN <- ziNB_fit
  simplified <- FALSE
  if (model_family %in% "binomial") {
    FUN <- binomial_fit
    simplified <- TRUE
  }
  while(try_left > 0) {
    if (formula_null == formula_full) {
      if ( is.na( residuals(models_result[["formula_null"]])[1]) ) {
        models_result[["formula_null"]] <- FUN(
          data = data,
          formula = formulas[["formula_null"]],
          gene_name = gene_name,
          model_family = model_family,
          link = link,
          zi = models_result[["is_zi"]],
          v = v
        )
      }
      models_result[["formula_full"]] <- models_result[["formula_null"]]
    } else {
      for (formula in names(formulas)) {
        if ( is.na(residuals(models_result[[formula]])[1]) ) {
          models_result[[formula]] <- FUN(
            data = data,
            formula = formulas[[formula]],
            gene_name = gene_name,
            model_family = model_family,
            link = link,
            zi = models_result[["is_zi"]],
            v = v
          )
          if (file.exists(tmp_file_LRT)) {
            system(paste0("rm ", tmp_file_LRT))
          }
        }
      }
    }
    try_left <- 0
    if ((is.na( residuals(models_result[["formula_null"]])[1] )|
        is.na( residuals(models_result[["formula_full"]])[1] )) &
        !simplified)
    {
      if (v) {
        print(paste0("error: in gene ", gene_name, " simplifying formula..."))
      }
      formulas <- simplify_formula(formulas)
      models_result[["formula_null"]]$residuals[1] <- NA
      models_result[["formula_full"]]$residuals[1] <- NA
      try_left <- 1
      if (!base::grepl("\\(1\\|(.*)\\)",
                       formulas[[names(formulas)[2]]], perl = T) &
          !base::grepl("batch",
                       formulas[[names(formulas)[2]]], perl = T)
          ) {
        simplified <- TRUE
      }
    }
  }
  models_result[["formulas"]] <- formulas
  if (!missing(folder_name)) {
    save(
      models_result, gene_name, formulas,
      file = tmp_file
    )
  }
  return(models_result)
}

#' importFrom glmmADMB glmmadmb
DEA_LRT <- function(models_result, gene_name, v, folder_name,
                    model_family = "nbinom1") {
  tmp_file <- ""
  LRT_result <- NA
  if (!missing(folder_name)) {
    system(paste0("mkdir -p ", folder_name, "/DEA_LRT/"))
    tmp_file <- paste0(folder_name, "/DEA_LRT/LRT_", gene_name, ".Rdata")
  }
  if (file.exists(tmp_file)) {
    if (v) {
      print(paste0("LRT cache found for ", gene_name))
    }
    load(tmp_file)
  } else {
    LRT_result <- tryCatch({
      anova(
        models_result[[ "formula_null" ]],
        models_result[[ "formula_full" ]]
      )
    }, error = function(e){
      if (v) {
        print(paste0("error: DEA_LRT for gene ", gene_name))
        print(e)
      }
      if (model_family %in% "nbinom1") {
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
      }
      if (model_family %in% "binomial") {
        return(
          data.frame(
            Df = c(NA, NA),
            AIC = c(NA, NA),
            BIC = c(NA, NA),
            logLik = c(NA, NA),
            deviance = c(NA, NA),
            Chisq = c(NA, NA),
            "Chisq Df" = c(NA, NA),
            "Pr(>Chisq)" = c(NA, NA),
            stringsAsFactors = FALSE
          )
        )
      }
    })
    if (model_family %in% "nbinom1") {
        if (!models_result[["is_zi"]] &
          !models_result[[ "formula_null" ]][["admb"]] &
          !base::grepl("\\(1\\|(.*)\\)",
                      models_result[["formulas"]][["formula_null"]],
                      perl = T)) {
        LRT_result <- data.frame(
          NoPar = c(NA, NA),
          LogLik = LRT_result[, 4],
          Df = LRT_result[, 3],
          Deviance = c(NA, NA),
          "Pr(>Chi)" = LRT_result[, min(8, ncol(LRT_result))],
          stringsAsFactors = FALSE
        )
      }
    }
    if (!missing(folder_name)) {
      save(
        LRT_result, gene_name,
        file = tmp_file
      )
    }
  }
  if (v) {
    if (model_family %in% "nbinom1") {
      print(paste0(
        "LRT : ", LRT_result[["Pr(>Chi)"]][2], " (", gene_name, ")"
      ))
    }
    if (model_family %in% "binomial") {
      print(paste0(
        "LRT : ", LRT_result[["Pr(>Chisq)"]][2], " (", gene_name, ")"
      ))
    }
  }
  return(list(LRT = LRT_result, file = tmp_file))
}

DEA_format <- function(LRT_result, models_result, v, model_family = "nbinom1") {
  if (model_family %in% "nbinom1") {
    return(DEA_format_ziNB(LRT_result, models_result, v))
  }
  if (model_family %in% "binomial") {
    return(DEA_format_binomial(LRT_result, models_result, v))
  }
}

DEA_format_binomial <- function(LRT_result, models_result, v) {
  results <- tryCatch({
    data.frame(
      null_loglik = LRT_result[["LRT"]][1, 4],
      full_loglik = LRT_result[["LRT"]][2, 4],
      null_Deviance = LRT_result[["LRT"]][1, 5],
      full_Deviance = LRT_result[["LRT"]][2, 5],
      null_df = LRT_result[["LRT"]][1, 1],
      full_df = LRT_result[["LRT"]][2, 1],
      pvalue = LRT_result[["LRT"]][2, 8],
      stringsAsFactors = FALSE
    )
  }, error = function(e){
    if (v) {
      print("error: DEA_format")
      print(e)
    }
    return(NA)
  })
  results <- c(
    unlist(data.frame(
      LRT_file = as.vector(LRT_result[["file"]]),
      models_file = as.vector(models_result[["file"]]),
      stringsAsFactors = FALSE
    )),
    results
  )
  return(results)
}


DEA_format_ziNB <- function(LRT_result, models_result, v) {
  results <- tryCatch({
    LRT <- data.frame(
      null_loglik = LRT_result[["LRT"]][1, 2],
      full_loglik = LRT_result[["LRT"]][2, 2],
      null_Deviance = LRT_result[["LRT"]][1, 4],
      full_Deviance = LRT_result[["LRT"]][2, 4],
      null_df = LRT_result[["LRT"]][1, 3],
      full_df = LRT_result[["LRT"]][2, 3],
      null_npar = LRT_result[["LRT"]][1, 1],
      full_npar = LRT_result[["LRT"]][2, 1],
      pvalue = LRT_result[["LRT"]][2, 5],
      stringsAsFactors = FALSE
    )
    model_null <- DEA_format_ziNB_model(
      models_result[["formula_null"]], models_result[["is_zi"]]
    )
    names(model_null) <- paste0(
      "null_",
      names(model_null)
    )
    model_full <- DEA_format_ziNB_model(
      models_result[["formula_full"]], models_result[["is_zi"]]
    )
    names(model_full) <- paste0(
      "full_",
      names(model_full)
    )
    return( c(model_null, model_full, unlist(LRT)) )
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

DEA_format_ziNB_model <- function(model, zi) {
  result <- tryCatch({
    unlist(data.frame(
      zi = ifelse(zi, model$pz, NA),
      alpha = model$alpha
    ))
  }, error = function(e){
    return(c(NA, NA))
  })
  return(result)
}

#' importFrom glmmADMB glmmadmb
#' importFrom MASS glm.nb
ziNB_fit <- function(data, formula, gene_name,
    model_family = "nbinom1", link = "log", zi = TRUE,
    v) {
  if (v) {
    if (length(zi) != 0 & zi) {
      print(paste0(gene_name, " : ziNB : ", formula))
    } else {
      print(paste0(gene_name, " : NB : ", formula))
    }
  }
  model <- tryCatch({
      model <- glmmADMB::glmmadmb(
        as.formula(formula),
        data = data,
        zeroInflation = zi,
        family = model_family,
        link = link,
        mcmc = FALSE,
        admb.opts = glmmADMB::admbControl(
          maxfn = 10000,
          imaxf = 10000,
          maxph = 100
        )
      )
      model[["admb"]] <- TRUE
      return(model)
  }, error = function(e){
    if (v) {
      if (zi) {
        print(paste0(
            "error: ziNB_fit for ", gene_name
        ))
      } else {
        print(paste0("error: NB_fit for ", gene_name))
      }
      print(e)
    }
    if (!zi &
        !base::grepl("\\(1\\|(.*)\\)", formula, perl = T) &
        !base::grepl("batch", formula, perl = T) ) {
      model <- tryCatch({
        if (v) {
          print(paste0("error: ADMB:NB_fit for ", gene_name,
                       " trying with MASS : ", formula))
        }
        MASS::glm.nb(
          as.formula(formula),
          data = data
        )
      }, error = function(e){
        if (v) {
          print(paste0("error: NB_fit for ", gene_name))
          print(e)
        }
        return(list(residuals = NA))
      })
      model[["admb"]] <- FALSE
      return(model)
    } else {
      return(list(residuals = NA, admb = TRUE))
    }
  })
  return(model)
}

#' importFrom lme4 glmer
binomial_fit <- function(data, formula, gene_name,
    model_family = "binomial", link = "log", zi = TRUE,
    v) {
  if (v) {
    print(paste0(gene_name, " : binomial : ", formula))
  }
  model <- tryCatch({
    lme4::glmer(as.formula(formula),
                data = data,
                family=binomial(link=log), nAGQ=0,
                control=glmerControl(optimizer = "nloptwrap")
                )
  }, error = function(e){
    if (v) {
      print(paste0("error: binomial_fit for ", gene_name))
      print(e)
    }
    return(list(residuals = NA))
  })
  return(model)
}

formula_to_features <- function(scd, formula_full, continuous = c()){
  b_features <- c()
  for(feature in colnames(scd$getfeatures)) {
    b_features <- c(
      b_features,
      grepl(
        pattern = feature,
        x = formula_full
      )
    )
  }
  b_features
  b_features[1] <- TRUE
  features <- apply(
    X = scd$getfeatures[, b_features],
    MARGIN = 2,
    FUN = as.factor
  )
  features <- as.data.frame(features)
  if (length(continuous) > 0) {
    if (any(!(continuous %in% colnames(features)))) {
      stop(paste0("error: continuous variable in ",
                  continuous, " not in formula"))
    }
    b_features <- which(colnames(features) %in% continuous)
    for (i in b_features){
      features[, i] <- as.numeric(as.vector(features[, i]))
    }
  }
  return(features)
}

#' @importFrom MASS glm.nb
#' @importFrom pscl zeroinfl
zi_test <- function(data, formula = y ~ 1, gene_name,
    model_family = "nbinom1", link = "log", threshold = 0.05, v = F){
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

simplify_formula <- function (formulas, v = F){
  if (!base::grepl("\\(1\\|(.*)\\)", formulas[[names( formulas )[2]]], perl = T) ) {
    for (formula in names(formulas)) {
      print( formulas[[formula]] )
      formulas[[formula]] <- gsub(
        "batch", "1", formulas[[formula]], perl = T
      )
      print( formulas[[formula]] )
    }
  }
  if (base::grepl("\\(1\\|(.*)\\)", formulas[[names( formulas )[2]]], perl = T) ) {
    for (formula in names(formulas)) {
      print( formulas[[formula]] )
      formulas[[formula]] <- gsub(
        "\\(1\\|([^+]*)\\)", "\\1", formulas[[formula]], perl = T
      )
      print( formulas[[formula]] )
    }
  }
  return(formulas)
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


#' up_down param information from DEA results
#'
#' @param folder_name name of the folder where to find the DEA results
#' @return a data.frame() with the farameters of the DEA
#' @examples
#' \dontrun {
#'
#'   b_cells <- scd$getfeature("QC_good") %in% T &
#'     !is.na(scd$getfeature("DEA_cell_type")) &
#'     scd$getfeature("day") %in% day
#'   dea_parameters <- up_down(
#'     folder_name = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA")
#'   )
#'
#' }
#' @export up_down
up_down <- function(folder_name) {
  f_DEA_fit <- paste0(folder_name, "/DEA_fit/")
  f_DEA_LRT <- paste0(folder_name, "/DEA_LRT/")
  DEA_fits <- list.files(f_DEA_fit)
  DEA_LRTs <- list.files(f_DEA_LRT)
  fits_coefs <- NA
  for ( DEA_fit in DEA_fits ) {
    fits_coefs <- add_coef(data = fits_coefs,
                           file = paste0(f_DEA_fit, DEA_fit))
  }
  fits_coefs$pvalue <- NA
  for ( DEA_LRT in DEA_LRTs ) {
    fits_coefs <- add_pvalue(data = fits_coefs,
                           file = paste0(f_DEA_LRT, DEA_LRT))
  }
  fits_coefs$padj <- stats::p.adjust(
    fits_coefs$pvalue,
    method = "BH"
  )
  return(fits_coefs)
}

add_coef <- function(data, file) {
  load(file)
  if ( length(models_result[["formula_full"]][["residuals"]]) > 1 ) {
    coef <- NA
    if ( "b" %in% names(models_result[["formula_full"]])) {
      coefs <- as.data.frame(
        t(models_result[["formula_full"]][["b"]])
      )
    } else {
      coefs <- as.data.frame(
        t(models_result[["formula_full"]][["coefficients"]])
      )
    }
    rownames(coefs) <- gene_name
    if (is.data.frame(data)) {
      data <- rbind(coefs, data)
    } else {
      data <- coefs
    }
  }
  return(data)
}

add_pvalue <- function(data, file) {
  load(file)
  if (gene_name %in% rownames(data)) {
    if ("Pr(>Chi)" %in% names(LRT_result)) {
      data[gene_name, "pvalue"] <- LRT_result[["Pr(>Chi)"]][2]
    } else {
      if ("Pr..Chi." %in% names(LRT_result)) {
        data[gene_name, "pvalue"] <- LRT_result[["Pr..Chi."]][2]
      } else {
        data[gene_name, "pvalue"] <- NA
      }
    }
  }
  return(data)
}
