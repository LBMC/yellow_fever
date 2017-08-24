library(scRNAtools)
context("load scRNASeq data")

data_test <- scdata$new(
  infos = data.frame(
    size = c(5:1),
    sex = c(rep("F", 5)),
    id = paste0(5:1, "cell")
  ),
  counts = matrix(
    1:500,
    nrow = 5,
    dimname = list(
      paste0(5:1,"cell"),
      paste0("gene", 1:100)
      )
  ))

test_that("loading data", {
  expect_identical(
    data_test$getgenes,
    paste0("gene", 1:100))
  expect_identical(
    data_test$getcells,
    paste0(1:5, "cell"))
  expect_equal(
    as.numeric(data_test$getgene("gene2")),
    10:6)
  expect_equal(
    data_test$getfeature("sex"),
    as.factor(c(rep("F", 5))))
  expect_identical(
    data_test$getcounts,
    matrix(
      1:500,
      nrow = 5,
      dimname = list(
        paste0(5:1, "cell"),
        paste0("gene", 1:100)
        )
      )[order(5:1), ]
  )
  expect_identical(
    data_test$getfeatures,
    data.frame(
      size = c(5:1),
      sex = c(rep("F", 5)),
      id = paste0(5:1, "cell")
    )[order(5:1), ]
)
})

test_that("selecting data", {
  expect_identical(
    data_test$getcountsw(
      cells = paste0(2:4, "cell"),
      genes = paste0("gene", 66:80)
    ),
    matrix(
      1:500,
      nrow = 5,
      dimname = list(
        paste0(5:1, "cell"),
        paste0("gene", 1:100)
        )
      )[order(5:1), ][2:4, 66:80]
  )
  expect_identical(
    data_test$getcountso(
      cells = paste0(2:4, "cell"),
      genes = paste0("gene", 66:80)
    ),
    matrix(
      1:500,
      nrow = 5,
      dimname = list(
        paste0(5:1, "cell"),
        paste0("gene", 1:100)
        )
      )[order(5:1), ][-c(2:4), -c(66:80)]
  )
  expect_identical(
    data_test$getfeaturesw(
      cells = paste0(2:4, "cell"),
      features = "sex"
    ),
    data.frame(
      size = c(5:1),
      sex = c(rep("F", 5)),
      id = paste0(5:1, "cell")
    )[order(5:1), ][c(2:4), 2]
  )
  expect_identical(
    data_test$getfeatureso(
      cells = paste0(2:4, "cell"),
      features = "sex"
    ),
    data.frame(
      size = c(5:1),
      sex = c(rep("F", 5)),
      id = paste0(5:1, "cell")
    )[order(5:1), ][-c(2:4), -2]
  )
})

data_test$add(
  infos = data.frame(
    size = c(6:9),
    sex = c(rep("M", 4)),
    id = paste0(6:9, "cell")
  ),
  counts = matrix(
    501:900,
    nrow = 4,
    dimname = list(
      paste0(6:9, "cell"),
      paste0("gene", 1:100)
      )
  )
)

test_that("adding data", {
  expect_identical(
    data_test$getgenes,
    paste0("gene", 1:100))
  expect_identical(
    data_test$getcells,
    paste0(1:9, "cell"))
  expect_equal(
    as.numeric(data_test$getgene("gene2")),
    c(10:6,505:508))
  expect_equal(
    data_test$getfeature("sex"),
    as.factor(c(rep("F", 5), rep("M", 4)))
    )
  expect_identical(
    data_test$getcounts,
    rbind(matrix(
      1:500,
      nrow = 5,
      dimname = list(
        paste0(5:1, "cell"),
        paste0("gene", 1:100)
        )
      )[order(5:1), ],
      matrix(
        501:900,
        nrow = 4,
        dimname = list(
          paste0(6:9, "cell"),
          paste0("gene", 1:100)
          )
      )
    )
  )
  expect_identical(
    data_test$getfeatures,
    rbind(
      data.frame(
          size = c(5:1),
          sex = c(rep("F", 5)),
          id = paste0(5:1, "cell")
        )[order(5:1), ],
      data.frame(
        size = c(6:9),
        sex = c(rep("M", 4)),
        id = paste0(6:9, "cell")
      )
    )
  )
})
