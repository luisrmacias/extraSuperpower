test_that("format check works", {
  faeff <- 1
  fA <- 2
  fbeff <- 3
  fB <- 4
  group_size <- 5
  rho <- 0.75
  fwithin <- "fB"
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)

  expect_error(twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 3))

  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = rho, withinf = fwithin)

  expect_no_error(twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 3))
})

test_that("design check works", {
  faeff <- 1
  fA <- 2
  fbeff <- 3
  fB <- 4
  rho <- 0.75
  fwithin <- "fB"
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = rho, withinf = fwithin)
  group_size <- 5.4
  expect_error(twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 3))

  group_size <- 10
  expect_error(twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 3,
                                             balanced = FALSE))

  group_size <- c(5, 6)
  expect_error(twoway_simulation_independent(group_size = group_size, mean_mat, nsims = 3,
                                             balanced = TRUE))

})

test_that("factor A as within factor", {
  faeff <- 1
  fA <- 2
  fbeff <- 3
  fB <- 2
  rho <- 0.9
  fwithin <- "fA"
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = rho, withinf = fwithin)

  group_size <- 10
  set.seed(15440804)
  sim <- twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 1)
  expect_gt(cor.test(sim$simulated_data$y[sim$simulated_data$fA=="fA_A"], sim$simulated_data$y[sim$simulated_data$fA=="fA_B"])$estimate,
            rho)
  expect_lt(cor.test(sim$simulated_data$y[sim$simulated_data$fB=="fB_a"], sim$simulated_data$y[sim$simulated_data$fB=="fB_b"])$estimate,
            0.2)
})


test_that("factor B as within factor", {
  faeff <- 3
  fA <- 4
  fbeff <- 3
  fB <- 4
  rho <- 0.9
  fwithin <- "fB"
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = rho, withinf = fwithin)

  group_size <- 10
  set.seed(15440804)
  sim <- twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 1)$simulated_data
  expect_gt(cor.test(sim$simulated_data$y[sim$simulated_data$fB=="fB_a"], sim$simulated_data$y[sim$simulated_data$fB=="fB_b"])$estimate,
            rho)
  expect_lt(cor.test(sim$simulated_data$y[sim$simulated_data$fA=="fA_A"], sim$simulated_data$y[sim$simulated_data$fA=="fA_B"])$estimate,
            0.2)
})
