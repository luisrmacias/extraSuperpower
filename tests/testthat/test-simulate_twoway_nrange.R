test_that("sample sizes match design", {
  faeff <- 10
  fA <- 2
  fbeff <- 0.5
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  set.seed(15440804)
  sampsizes <- seq(4,10,2)
  iterations <- 10
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  obspercond <- sapply(simindep, function(x) table(x$n[x$cond=="V1"]))
  expect_equal(as.vector(obspercond), sampsizes*iterations)
})
