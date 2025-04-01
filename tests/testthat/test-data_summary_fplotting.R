faeff <- 2
fA <- 3
fbeff <- 0.5
fB <- 2

test_that("data missingness in checked", {
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  meansd_mats$mean.mat[2,2] <- NA
  expect_warning(graph_twoway_assumptions(matrices_obj = meansd_mats))
  })
