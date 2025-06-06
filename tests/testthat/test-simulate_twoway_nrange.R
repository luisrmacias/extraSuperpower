test_that("balanced sample sizes match design", {
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

test_that("unbalanced sample sizes match design", {
  faeff <- 10
  fA <- 2
  fbeff <- 0.5
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)
  gsize <- matrix(c(6,4,6,4), 2, 2, byrow = TRUE)
  sampsizes <- 0:2
  iterations <- 10
  simindep <- simulate_twoway_nrange(nset = sampsizes, matrices_obj = mean_mat,
                                     group_size = gsize, balanced = FALSE, nsims = iterations)
  obspercond <- sapply(simindep, function(x) table(x$n, x$cond))
  nsim <- unlist(sapply(simindep, function(x) table(x$cond)/iterations))
  expect_true(all(nsim==sapply(seq(sampsizes), function(x) t(gsize)+sampsizes[x])))
})

test_that("type of design input check", {
  faeff <- 10
  fA <- 2
  fbeff <- 0.5
  fB <- 2
  sampsizes <- seq(4,10,2)
  iterations <- 100
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)
  expect_error(simindep <- simulate_twoway_nrange(nset = sampsizes, matrices_obj = mean_mat, nsims = iterations,
                                     repeated_measurements = TRUE))
  expect_no_error(simindep <- simulate_twoway_nrange(nset = sampsizes, matrices_obj = mean_mat, nsims = iterations))
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff, plot = FALSE)
  expect_error(simindep <- simulate_twoway_nrange(nset = sampsizes, matrices_obj = mean_mat, nsims = iterations,
                                                  repeated_measurements = TRUE))
  expect_no_error(simindep <- simulate_twoway_nrange(nset = sampsizes, matrices_obj = mean_mat, nsims = iterations))
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = 0.5, withinf = "both")
  expect_error(simindep <- simulate_twoway_nrange(nset = sampsizes, matrices_obj = mean_mat,
                                                  group_size = gsize, nsims = iterations))
  expect_no_error(simrep <- simulate_twoway_nrange(nset = sampsizes, matrices_obj = mean_mat, nsims = iterations,
                                                  repeated_measurements = TRUE))
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = 0.5, withinf = "both", plot = FALSE)
  expect_error(simindep <- simulate_twoway_nrange(nset = sampsizes, matrices_obj = mean_mat,
                                                  group_size = gsize, nsims = iterations))
  expect_no_error(simrep <- simulate_twoway_nrange(nset = sampsizes, matrices_obj = mean_mat, nsims = iterations,
                                                   repeated_measurements = TRUE))
})
