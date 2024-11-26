test_that("matrices dimensions", {
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = 5, nlfB = 4,
                                    fAeffect = 2, fBeffect = 3, plot = FALSE)[[1]]
  expect_equal(dim(mean_mat), c(5, 4))
})

test_that("factor A stepwise effect ratio", {
  faeff <- 2
  fA <- 5
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = 4,
                                    fAeffect = faeff, fBeffect = 3, plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  faend.1 <- mean_mat[fA,1]
  expect_equal(faend.1/f1.1, faeff)
})

test_that("factor B stepwise effect", {
  faeff <- 2
  fA <- 5
  fbeff <- 3
  fB <- 4
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff, plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  fbend.1 <- mean_mat[1,fB]
  expect_equal(fbend.1/f1.1, fbeff)
})

test_that("factor A end effect ratio", {
  faeff <- 2
  fA <- 5
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = 4,
                                    fAeffect = faeff, fBeffect = 3, endincrement = TRUE,
                                    plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  fanext.1 <- mean_mat[2,1]
  expect_equal(fanext.1/f1.1, faeff)
})

test_that("factor B end effect ration", {
  faeff <- 2
  fA <- 5
  fbeff <- 3
  fB <- 4
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff, endincrement = TRUE,
                                    plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  fbnext.1 <- mean_mat[1,2]
  expect_equal(fbnext.1/f1.1, fbeff)
})
