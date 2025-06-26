
set.seed(1552104)

test_that("format check works", {
  faeff <- 1
  fA <- 2
  fbeff <- 3
  fB <- 2
  group_size <- 5
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                     fAeffect = faeff, fBeffect = fbeff)
  dimnames(mean_mat$matrices_obj$mean.mat) <- NULL
  expect_error(twoway_simulation_independent(group_size = group_size, mean_mat, nsims = 3))
})

test_that("design check works", {
  faeff <- 1
  fA <- 2
  fbeff <- 3
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff)
  group_size <- 5.4
  expect_error(twoway_simulation_independent(group_size = group_size, mean_mat, nsims = 3))

  group_size <- c(5, 6, 4)
  expect_error(twoway_simulation_independent(group_size = group_size, mean_mat, nsims = 3,
                                             balanced = FALSE))

  group_size <- 10
  expect_error(twoway_simulation_independent(group_size = group_size, mean_mat, nsims = 3,
                                             balanced = FALSE))

  group_size <- c(5, 6)
  expect_error(twoway_simulation_independent(group_size = group_size, mean_mat, nsims = 3,
                                             balanced = TRUE))
})

test_that("input format check works", {
  faeff <- 1
  fA <- 2
  fbeff <- 3
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = 0.6, withinf = "both")
  expect_error(twoway_simulation_independent(10, mean_mat))
})

test_that("number of levels generated correspond to model", {
  nlevfA <- 2
  nlevfB <- 2
  group_size <- 5
  mean.mat <- matrix(c(1, 0, 0, 0), nlevfA, nlevfB, dimnames = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  sd.mat <- 2
  matlist <- list(mean.mat=mean.mat, sd.mat=sd.mat)
  simdat <- twoway_simulation_independent(group_size = group_size, matlist, nsims = 3)
  expect_equal(c(nlevfA, nlevfB, 3), dim(table(simdat$treatment, simdat$groups, simdat$iteration)))

  nlevfA <- 4
  nlevfB <- 2
  group_size <- 5
  mean.mat <- matrix(1:8, nlevfA, nlevfB, dimnames = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  sd.mat <- 2
  matlist <- list(mean.mat=mean.mat, sd.mat=sd.mat)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = 3)
  expect_equal(c(nlevfA, nlevfB, 3), dim(table(simdat$groups, simdat$treatment, simdat$iteration)))
})

test_that("number of subject generated correspond to sample size", {
  nlevfA <- 2
  nlevfB <- 2
  group_size <- 5
  iterations <- 3
  mean.mat <- matrix(c(1, 0, 0, 0), nlevfA, nlevfB, dimnames = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  sd.mat <- 2
  matlist <- list(mean.mat=mean.mat, sd.mat=sd.mat)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  expect_equal(nrow(simdat), nlevfA*nlevfB*group_size*iterations)

  nlevfA <- 3
  nlevfB <- 6
  mean.mat <- matrix(1:9, nlevfA, nlevfB, dimnames = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  matlist <- list(mean.mat=mean.mat, sd.mat=sd.mat)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  expect_equal(nrow(simdat), nlevfA*nlevfB*group_size*iterations)
})

test_that("mean simulated values match mean matrix", {
  nlevfA <- 2
  nlevfB <- 2
  group_size <- 400
  iterations <- 1
  mean.mat <- matrix(c(1, 0.001, 0.001, 0.001), nlevfA, nlevfB,
                     dimnames = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  sd.mat <- 0.001
  matlist <- list(mean.mat=mean.mat, sd.mat=sd.mat)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gmeans_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = mean)
  expect_true(all(abs(gmeans_simdat$y - mean.mat)<1e-4))

  nlevfA <- 3
  nlevfB <- 6
  mean.mat <- matrix(1:9, nlevfA, nlevfB, dimnames = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  matlist <- list(mean.mat=mean.mat, sd.mat=sd.mat)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gmeans_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = mean)
  expect_true(all(abs(gmeans_simdat$y - mean.mat)<1e-4))
})

test_that("mean simulated values match mean matrix with interaction", {
  nlevfA <- 2
  nlevfB <- 2
  group_size <- 400
  iterations <- 1
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 0.01, fBeffect = 0.01, plot = FALSE, sdratio = 0.1, sdproportional = FALSE,
                                   groupswinteraction = c(2, 2), interact = -50,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gmeans_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = mean)
  expect_true(all(abs(gmeans_simdat$y - matlist$mean.mat)<1e-1))

  nlevfA <- 3
  nlevfB <- 6
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 0.01, fBeffect = 0.01, plot = FALSE, sdratio = 0.1, sdproportional = FALSE,
                                   groupswinteraction = matrix(c(2,2,3,3,1,4,2,5,3,6), 5, 2, byrow = TRUE), interact = 50,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gmeans_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = mean)
  expect_true(all(abs(gmeans_simdat$y - matlist$mean.mat)<1e-1))
})

test_that("sd of simulated values match sd matrix", {
  nlevfA <- 2
  nlevfB <- 2
  group_size <- 400
  iterations <- 1
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 100, fBeffect = 100, plot = FALSE, sdratio = 0.1,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gsds_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = sd)
  expect_true(all((abs(gsds_simdat$y) - matlist$sd.mat)/matlist$sd.mat<0.1))

  nlevfA <- 3
  nlevfB <- 6
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 100, fBeffect = 100, plot = FALSE, sdratio = 0.1,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gsds_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = sd)
  expect_true(all((abs(gsds_simdat$y) - matlist$sd.mat)/matlist$sd.mat<0.1))
})

test_that("simulated values are normally distributed", {
  nlevfA <- 2
  nlevfB <- 2
  label_list <- list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB])
  group_size <- 100
  iterations <- 1
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 0.01, fBeffect = 0.01, plot = FALSE, sdratio = 0.1, sdproportional = FALSE,
                                   label_list = label_list)
  set.seed(160724)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  distest <- tapply(simdat$y, simdat$cond, shapiro.test)
  expect_true(all(sapply(distest, "[", "p.value")>0.05))

  nlevfA <- 3
  nlevfB <- 6
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 0.01, fBeffect = 0.01, plot = FALSE, sdratio = 0.1, sdproportional = FALSE,
                                   groupswinteraction = matrix(c(2,2,3,3,1,4,2,5,3,6), 5, 2, byrow = TRUE), interact = 50,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
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
  matlist <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 2, fBeffect = 2, plot = FALSE, sdratio = 0.3,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))

  expect_warning(twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations,
                                          distribution = "skewed", superior_limit = suplim))

  expect_warning(twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations,
                                               distribution = "truncated.normal", skewness = 5))
})

test_that("simulated values are skewed", {
  nlevfA <- 2
  nlevfB <- 4
  label_list <- list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB])
  group_size <- 100
  iterations <- 1
  matlist <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 2, fBeffect = 0.5, plot = FALSE,
                                   label_list = label_list)
  set.seed(160724)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist,
                                          nsims = iterations, distribution = "skewed", skewness = 0.01)
  nsubabmean <- min(tapply(simdat$y, simdat$cond, function(x)sum(x>mean(x))))
  expect_gte(nsubabmean, group_size/2)

  nlevfA <- 3
  nlevfB <- 6
  group_size <- 200
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 5, fBeffect = 0.2, plot = FALSE, sdratio = 0.1, sdproportional = FALSE)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist,
                                          distribution = "skewed", skewness = 0.1, nsims = iterations)
  distpvals <- ks.test(simdat$y[simdat$cond=="V1"], "pnorm", matlist$mean.mat[1,1], matlist$sd.mat)$p.value
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V6"], "pnorm", matlist$mean.mat[1,6], matlist$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V7"], "pnorm", matlist$mean.mat[2,1], matlist$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V12"], "pnorm", matlist$mean.mat[2,6], matlist$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V13"], "pnorm", matlist$mean.mat[3,1], matlist$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="V18"], "pnorm", matlist$mean.mat[3,6], matlist$sd.mat)$p.value)
  expect_true(all(distpvals<0.07))
})

test_that("simulated values respect truncation limit", {
  nlevfA <- 2
  nlevfB <- 4
  group_size <- 100
  iterations <- 1
  suplim <- 12
  matlist <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 2, fBeffect = 2, plot = FALSE, sdratio = 0.3,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  set.seed(160724)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations,
                                          distribution = "truncated.normal", superior_limit = suplim)

  expect_lt(max(simdat$y[simdat$cond=="V5"]), suplim)
})


test_that("unbalanced designs are simulated", {
  nlevfA <- 2
  nlevfB <- 4
  iterations <- 1
  matlist <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 2, fBeffect = 2, plot = FALSE, sdratio = 0.3,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  group_size <- matrix(seq(12, 6, -2), nlevfA, nlevfB, byrow = TRUE)
  set.seed(160724)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations,
                                          balanced = FALSE)
  expect_equal(table(simdat$cond), as.vector(t(group_size)), ignore_attr=TRUE)
})
