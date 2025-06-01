test_that("workflow", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  set.seed(15440804)
  sampsizes <- seq(6,15,2)
  iterations <- 100
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  res <- test_power_overkn(simindep)
  expect_true(all(res$power_table[seq(1,13,3),2] > res$power_table[seq(2,14,3),2]))
})
