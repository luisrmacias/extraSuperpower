faeff <- 2
fA <- 3
fbeff <- 0.5
fB <- 2

test_that("repeated data input", {
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       rho=0.8, withinf="fB", plot = FALSE)
  p <- graph_twoway_assumptions(group_size = 10, matrices_obj = meansd_mats)
  expect_no_error(p)
})
