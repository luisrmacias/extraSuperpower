faeff <- 2
fA <- 3
fbeff <- 0.5
fB <- 2

test_that("data missingness indirectly", {
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  meansd_mats$mean.mat[2,2] <- NA
  expect_warning(graph_twoway_assumptions(matrices_obj = meansd_mats))
  })

test_that("data missingness directly", {
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  set.seed(15020404)
  sim <- twoway_simulation_independent(group_size = 5, matrices_obj = meansd_mats, nsims = 1)
  sim$y[c(1,26)] <- NA
  res <- summarySE(data = sim, measurevar = "y", groupvars = c("fA", "fB"), na.rm = TRUE)
  expect_true(all(complete.cases(res)))
})

