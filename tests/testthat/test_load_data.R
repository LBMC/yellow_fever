library(scRNAtools)
context("load scRNASeq data")

data_test <- scdata$new(
  infos = data.frame(
    size = c(1:10),
    sex = c(rep("F", 5), rep("M", 5))
  ),
  counts = matrix(
    1:1000,
    nrow = 10,
    dimname = list(
      paste0("cell", 1:10),
      paste0("gene", 1:100)
      )
  ))

test_that("loading data", {
  expect_identical(
    data_test$getgenes,
    paste0("gene", 1:100))
  expect_identical(
    data_test$getcells,
    paste0("cell", 1:10))
  expect_equal(
    as.numeric(data_test$getgene("gene2")),
    11:20)
  expect_equal(
    data_test$getfeature("sex"),
    as.factor(c(rep("F", 5), rep("M", 5))))
  expect_identical(
    data_test$getcounts,
    matrix(
      1:1000,
      nrow = 10,
      dimname = list(
        paste0("cell", 1:10),
        paste0("gene", 1:100)
        )
      )
  )
  expect_identical(
    data_test$getfeatures,
    data.frame(
      size = c(1:10),
      sex = c(rep("F", 5), rep("M", 5))
    )
  )
})

test_that("selecting data", {
  expect_identical(
    data_test$getcountsw(
      cells = paste0("cell", 2:4),
      genes = paste0("gene", 66:80)
    ),
    matrix(
      1:1000,
      nrow = 10,
      dimname = list(
        paste0("cell", 1:10),
        paste0("gene", 1:100)
        )
      )[2:4, 66:80]
  )
  expect_identical(
    data_test$getcountso(
      cells = paste0("cell", 2:4),
      genes = paste0("gene", 66:80)
    ),
    matrix(
      1:1000,
      nrow = 10,
      dimname = list(
        paste0("cell", 1:10),
        paste0("gene", 1:100)
        )
      )[-c(2:4), -c(66:80)]
  )
  expect_identical(
    data_test$getfeaturesw(
      cells = paste0("cell", 2:4),
      features = "sex"
    ),
    data.frame(
      size = c(1:10),
      sex = c(rep("F", 5), rep("M", 5))
    )[c(2:4), 2]
  )
  expect_identical(
    data_test$getfeatureso(
      cells = paste0("cell", 2:4),
      features = "sex"
    ),
    data.frame(
      size = c(1:10),
      sex = c(rep("F", 5), rep("M", 5))
    )[-c(2:4), -2]
  )
})
