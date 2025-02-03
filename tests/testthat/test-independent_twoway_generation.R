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

test_that("simulated values match expectations", {
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
