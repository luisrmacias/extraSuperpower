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
  expect_error(twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 3,
                                            balanced = FALSE))

  group_size <- 10
  expect_error(twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 3,
                                             balanced = FALSE))

  group_size <- c(5, 6)
  expect_error(twoway_simulation_independent(group_size = group_size, mean_mat, nsims = 3,
                                             balanced = TRUE))

})

test_that("factor A as within factor", {
  faeff <- 0.5
  fA <- 2
  fbeff <- 10
  fB <- 2
  rho <- 0.9
  fwithin <- "fA"
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = rho, withinf = fwithin)

  group_size <- 100
  set.seed(15440804)
  sim <- twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 1)$simulated_data
  fBcor1 <- cor.test(sim$y[sim$cond=="A_a"], sim$y[sim$cond=="A_b"])$estimate
  expect_lt(abs(fBcor1),
            0.2)
  fBcor2 <- cor.test(sim$y[sim$cond=="B_a"], sim$y[sim$cond=="B_b"])$estimate
  expect_lt(abs(fBcor2),
            0.2)
  fAcor1 <- cor.test(sim$y[sim$cond=="A_a"], sim$y[sim$cond=="B_a"])$estimate
  expect_gt(abs(fAcor1),
            rho-0.1)
  fAcor2 <- cor.test(sim$y[sim$cond=="A_b"], sim$y[sim$cond=="B_b"])$estimate
  expect_gt(abs(fAcor2),
            rho-0.1)
})


test_that("factor B as within factor", {
  faeff <- 10
  fA <- 2
  fbeff <- 0.5
  fB <- 2
  rho <- 0.9
  fwithin <- "fB"
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = rho, withinf = fwithin)

  group_size <- 100
  set.seed(15440804)
  sim <- twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 1)$simulated_data
  fBcor1 <- cor.test(sim$y[sim$cond=="A_a"], sim$y[sim$cond=="A_b"])$estimate
  expect_gt(abs(fBcor1),
            rho-0.1)
  fBcor2 <- cor.test(sim$y[sim$cond=="B_a"], sim$y[sim$cond=="B_b"])$estimate
  expect_gt(abs(fBcor2),
            rho-0.1)
  fAcor1 <- cor.test(sim$y[sim$cond=="A_a"], sim$y[sim$cond=="B_a"])$estimate
  expect_lt(abs(fAcor1),
            0.2)
  fAcor2 <- cor.test(sim$y[sim$cond=="A_b"], sim$y[sim$cond=="B_b"])$estimate
  expect_lt(abs(fAcor2),
           0.2)
})
