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
                                            balanced = TRUE))

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

test_that("simulated values are normally distributed", {
  nlevfA <- 2
  nlevfB <- 2
  label_list <- list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB])
  group_size <- 100
  iterations <- 1
  rho <- -0.9
  fwithin <- "fB"
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 5, fBeffect = 0.2, plot = FALSE,
                                   sdratio = 0.1, sdproportional = FALSE,
                                   rho = rho, withinf = fwithin,
                                   label_list = label_list)
  set.seed(160724)
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs, nsims = iterations)$simulated_data
  distest <- tapply(simdat$y, simdat$cond, shapiro.test)
  expect_true(all(sapply(distest, "[", "p.value")>0.05))

  nlevfA <- 3
  nlevfB <- 6
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 5, fBeffect = 0.2, plot = FALSE, sdratio = 0.1, sdproportional = FALSE,
                                   groupswinteraction = matrix(c(2,2,3,3,1,4,2,5,3,6), 5, 2, byrow = TRUE), interact = 50,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]),
                                   rho = rho, withinf = fwithin)
  refs$sigmat <- Matrix::nearPD(refs$sigmat)$mat
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs, nsims = iterations)$simulated_data
  distest <- tapply(simdat$y, simdat$cond, shapiro.test)
  pvals <- unlist(sapply(distest, "[", "p.value"))
  expect_true(all(p.adjust(pvals)>0.05))
})

test_that("skewed and truncated input checks work", {
  nlevfA <- 2
  nlevfB <- 4
  group_size <- 100
  iterations <- 1
  suplim <- 12
  rho <- 0.8
  fwithin <- "both"
  matlist <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 2, fBeffect = 2, plot = FALSE, sdratio = 0.3,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]),
                                   rho = rho, withinf = fwithin)

  expect_warning(twoway_simulation_correlated(group_size = group_size, matrices_obj = matlist, nsims = iterations,
                                               distribution = "skewed", superior_limit = suplim))

  expect_warning(twoway_simulation_correlated(group_size = group_size, matrices_obj = matlist, nsims = iterations,
                                               distribution = "truncated.normal", shape.parameter = 2))
})


# test_that("simulated values are skewed", {
#   nlevfA <- 2
#   nlevfB <- 4
#   label_list <- list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB])
#   group_size <- 100
#   iterations <- 1
#   rho <- 0.6
#   fwithin <- "fB"
#   refs <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
#                                    fAeffect = 2, fBeffect = 0.5, plot = FALSE,
#                                    label_list = label_list, rho = rho, withinf = fwithin)
#   set.seed(160724)
#   simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs,
#                                          nsims = iterations,
#                                          distribution = "skewed", shape.parameter = 0.1)$simulated_data
#   nsubabmean <- min(tapply(simdat$y, simdat$cond, function(x)sum(x>mean(x))))
#   expect_gte(nsubabmean, group_size/2)
#
#   nlevfA <- 3
#   nlevfB <- 6
#   group_size <- 200
#   matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
#                                    fAeffect = 5, fBeffect = 0.2, plot = FALSE, sdratio = 0.1, sdproportional = FALSE)
#   simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist,
#                                           distribution = "skewed", skewness = 0.1, nsims = iterations)
#   distpvals <- ks.test(simdat$y[simdat$cond=="V1"], "pnorm", matlist$mean.mat[1,1], matlist$sd.mat)$p.value
#   distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V6"], "pnorm", matlist$mean.mat[1,6], matlist$sd.mat)$p.value)
#   distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V7"], "pnorm", matlist$mean.mat[2,1], matlist$sd.mat)$p.value)
#   distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V12"], "pnorm", matlist$mean.mat[2,6], matlist$sd.mat)$p.value)
#   distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V13"], "pnorm", matlist$mean.mat[3,1], matlist$sd.mat)$p.value)
#   distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V18"], "pnorm", matlist$mean.mat[3,6], matlist$sd.mat)$p.value)
#   expect_true(all(distpvals<0.07))
# })
