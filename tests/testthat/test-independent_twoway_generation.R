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
  group_size <- 1000
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
  group_size <- 1000
  iterations <- 1
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 0.01, fBeffect = 0.01, plot = FALSE, sdratio = 0.1, sdproportional = FALSE,
                                   groupswinteraction = c(2, 2), interact = -50,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gmeans_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = mean)
  expect_true(all(abs(gmeans_simdat$y - matlist$mean.mat)<1e-4))

  nlevfA <- 3
  nlevfB <- 6
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 0.01, fBeffect = 0.01, plot = FALSE, sdratio = 0.1, sdproportional = FALSE,
                                   groupswinteraction = matrix(c(2,2,3,3,1,4,2,5,3,6), 5, 2, byrow = TRUE), interact = 50,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gmeans_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = mean)
  expect_true(all(abs(gmeans_simdat$y - matlist$mean.mat)<1e-4))
})

test_that("sd of simulated values match sd matrix", {
  nlevfA <- 2
  nlevfB <- 2
  group_size <- 1000
  iterations <- 1
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 100, fBeffect = 100, plot = FALSE, sdratio = 0.1,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gsds_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = sd)
  expect_true(all((abs(gsds_simdat$y) - matlist$sd.mat)/matlist$sd.mat<0.05))

  nlevfA <- 3
  nlevfB <- 6
  matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 100, fBeffect = 100, plot = FALSE, sdratio = 0.1,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
  gsds_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = sd)
  expect_true(all((abs(gsds_simdat$y) - matlist$sd.mat)/matlist$sd.mat<0.05))
})

# test_that("simulated values are normally distributed", {
#   nlevfA <- 2
#   nlevfB <- 2
#   group_size <- 200
#   iterations <- 5
#   matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
#                                    fAeffect = 0.01, fBeffect = 0.01, plot = FALSE, sdratio = 0.1, sdproportional = FALSE,
#                                    label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
#   simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
#   gmeans_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = mean)
#   expect_true(all(abs(gmeans_simdat$y - matlist$mean.mat)<1e-4))
#
#   nlevfA <- 3
#   nlevfB <- 6
#   matlist <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
#                                    fAeffect = 0.01, fBeffect = 0.01, plot = FALSE, sdratio = 0.1, sdproportional = FALSE,
#                                    groupswinteraction = matrix(c(2,2,3,3,1,4,2,5,3,6), 5, 2, byrow = TRUE), interact = 50,
#                                    label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
#   simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = iterations)
#   gmeans_simdat <- aggregate(y ~ groups + treatment, data = simdat, FUN = mean)
#   expect_true(all(abs(gmeans_simdat$y - matlist$mean.mat)<1e-4))
# })
