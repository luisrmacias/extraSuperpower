test_that("workflow", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff, plot = FALSE)

  set.seed(15440804)
  sampsizes <- seq(9,15,2)
  iterations <- 30
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  res <- test_power_overkn(simindep)
  expect_true(all(res$power_table$power[res$power_table$effect=="fA"] >
                  res$power_table$power[res$power_table$effect=="fB"]))
})

test_that("default output class", {
  faeff <- 1.3
  fA <- 2
  fbeff <- 1.2
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff, plot = FALSE)

  sampsizes <- seq(6,8,2)
  iterations <- 5
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
                                    fAeffect = faeff, fBeffect = fbeff, plot = FALSE)

  sampsizes <- seq(6,8,2)
  iterations <- 5
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)
  res <- test_power_overkn(data = simindep, plot = FALSE)
  expect_s3_class(object = res, class = "data.frame")
})

faeff <- 1.3
fA <- 2
fbeff <- 1.2
fB <- 2
sampsizes <- c(6, 8)
iterations <- 20
mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                  fAeffect = faeff, fBeffect = fbeff)

set.seed(15440804)
simindep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes, nsims = iterations)

mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                  fAeffect = faeff, fBeffect = fbeff,
                                  rho = 0.7, withinf = "both", plot = FALSE)

simrep <- simulate_twoway_nrange(matrices_obj = mean_mat, nset = sampsizes,
                                 nsims = iterations, repeated_measurements = TRUE)

test_that("graph options", {
  expect_error(test_power_overkn(data = simindep, target_power = 80))
  expect_no_error(test_power_overkn(data = simindep, target_power = 0.9))
  expect_error(test_power_overkn(data = simindep, title = 2))
  expect_no_error(test_power_overkn(data = simindep, title = "Test Title"))
})


test_that("permutation is performed", {
  expect_message(test_power_overkn(data = simindep, plot = FALSE, test = "permutation"), regexp = "permutation")

  expect_message(test_power_overkn(data = simrep, plot = FALSE, test = "permutation"), regexp = "permutation")
})


test_that("rank is performed", {
  expect_message(test_power_overkn(data = simindep, plot = FALSE, test = "rank"), regexp = "rank")

  expect_message(test_power_overkn(data = simrep, plot = FALSE, test = "rank"), regexp = "rank")
})

test_that("unbalanced designs are tested", {
  gsize <- matrix(rep(6:7, 2), 2, 2, byrow = TRUE)
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)
  simindep <- simulate_twoway_nrange(matrices_obj = mean_mat,
                                     nset = 0:3, balanced = FALSE, group_size = gsize, nsims = iterations)
  res <- test_power_overkn(data = simindep)
  expect_s3_class(object = res[[1]], class = "data.frame")
  expect_s3_class(object = res[[2]], class = "gg")
})

