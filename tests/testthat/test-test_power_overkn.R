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

test_that("default output class", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  sampsizes <- seq(6,8,2)
  iterations <- 100
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  res <- test_power_overkn(data = simindep)
  expect_s3_class(object = res[[1]], class = "data.frame")
  expect_s3_class(object = res[[2]], class = "gg")
})


test_that("no graph option output class", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  sampsizes <- seq(6,8,2)
  iterations <- 100
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  res <- test_power_overkn(data = simindep, plot = FALSE)
  expect_s3_class(object = res, class = "data.frame")
})


test_that("graph options", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  sampsizes <- seq(6,8,2)
  iterations <- 100
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  expect_error(test_power_overkn(data = simindep, target_power = 80))
  expect_no_error(test_power_overkn(data = simindep, target_power = 0.9))
  expect_error(test_power_overkn(data = simindep, title = 2))
  expect_no_error(test_power_overkn(data = simindep, title = "Test Title"))
})


test_that("permutation is performed", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  sampsizes <- seq(6,8,2)
  iterations <- 100
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  tt <- system.time(test_power_overkn(data = simindep, plot = FALSE, test = "permutation"))[3]
  expect_gt(object = tt, expected = 14)
})

test_that("permutation is performed", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  iterations <- 100
  sampsizes <- seq(6,8,2)
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  tt <- system.time(test_power_overkn(data = simindep, plot = FALSE, test = "permutation"))[3]
  expect_gt(object = tt, expected = 14)

  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = 0.7, withinf = "both")

  simrep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations, repeated_measurements = TRUE)
  tt <- system.time(test_power_overkn(data = simrep, plot = FALSE, test = "permutation"))[3]
  expect_gt(object = tt, expected = 14)
})


test_that("rank is performed", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  iterations <- 100
  sampsizes <- seq(6,8,2)
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  tt <- system.time(test_power_overkn(data = simindep, plot = FALSE, test = "rank"))[3]
  expect_gt(object = tt, expected = 2)

  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = 0.7, withinf = "both")

  simrep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations, repeated_measurements = TRUE)
  tt <- system.time(test_power_overkn(data = simrep, plot = FALSE, test = "rank"))[3]
  expect_gt(object = tt, expected = 1)
})

test_that("unbalanced designs are tested", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 3
  iterations <- 100
  gsize <- matrix(rep(6:8, 2), 2, 3, byrow = TRUE)
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat,
                                     nset = 0:3, balanced = FALSE, group_size = gsize, nsims = iterations)
  res <- test_power_overkn(data = simindep)
  expect_s3_class(object = res[[1]], class = "data.frame")
  expect_s3_class(object = res[[2]], class = "gg")
})

