test_that("graph input check", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  sampsizes <- seq(6,8,2)
  iterations <- 100
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  expect_error(plot_powercurves(simindep))
})
