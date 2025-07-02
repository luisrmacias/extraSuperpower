set.seed(25089)
test_that("tests are done according to input", {
  nlevfA <- 4
  nlevfB <- 4
  group_size <- 5
  gwint <- matrix(c(1,3,1,4,2,4), 3, 2, byrow = TRUE)
  nlist <- list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB])
  matlist <- calculate_mean_matrix(refmean = 3, nlfA = nlevfA, nlfB = nlevfB, fAeffect = 2, fBeffect = 1.2,
                                   groupswinteraction = gwint, interact = 1.4, sdproportional = FALSE,
                                   label_list = nlist)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = 3)
  expect_message(twoway_simulation_testing(data = simdat), "independent observations")

  matlist <- calculate_mean_matrix(refmean = 3, nlfA = nlevfA, nlfB = nlevfB, fAeffect = 2, fBeffect = 1.2,
                                   groupswinteraction = gwint, interact = 1.4, sdproportional = FALSE,
                                   label_list = nlist, plot = FALSE)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = 3)
  expect_message(twoway_simulation_testing(data = simdat), "independent observations")

  matlist <- calculate_mean_matrix(refmean = 3, nlfA = nlevfA, nlfB = nlevfB, fAeffect = 2, fBeffect = 1.2,
                                   groupswinteraction = gwint, interact = 1.4, sdproportional = FALSE,
                                   label_list = nlist, rho = 0.8, withinf = "fB")
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = matlist, nsims = 3)
  expect_message(twoway_simulation_testing(data = simdat), "repeated observations")

  matlist <- calculate_mean_matrix(refmean = 3, nlfA = nlevfA, nlfB = nlevfB, fAeffect = 2, fBeffect = 1.2,
                                   groupswinteraction = gwint, interact = 1.4, sdproportional = FALSE,
                                   label_list = nlist, rho = 0.8, withinf = "fB", plot = FALSE)
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = matlist, nsims = 3)
  expect_message(twoway_simulation_testing(data = simdat), "repeated observations")
})

test_that("correct independent sample test is performed", {
  nlevfA <- 4
  nlevfB <- 4
  group_size <- 5
  gwint <- matrix(c(1,3,1,4,2,4), 3, 2, byrow = TRUE)
  nlist <- list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB])
  matlist <- calculate_mean_matrix(refmean = 3, nlfA = nlevfA, nlfB = nlevfB, fAeffect = 2, fBeffect = 1.2,
                                   groupswinteraction = gwint, interact = 1.4, sdproportional = FALSE,
                                   label_list = nlist)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = 3)
  aov_time <- system.time(twoway_simulation_testing(simdat, test = "ANOVA"))[3]
  per_time <- system.time(twoway_simulation_testing(simdat, test = "permutation"))[3]
  rank_time <- system.time(twoway_simulation_testing(simdat, test = "rank"))[3]
  expect_true(all(c(per_time>rank_time), c(rank_time>aov_time)))
})

test_that("correct repeated sample test is performed", {
  nlevfA <- 3
  nlevfB <- 4
  label_list <- list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB])
  group_size <- 5
  iterations <- 3
  rho <- -0.9
  fwithin <- "fB"
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 1.1, fBeffect = 0.8, plot = FALSE, sdratio = 0.2, sdproportional = FALSE,
                                groupswinteraction = matrix(c(2,2,3,3,1,4), 3, 2, byrow = TRUE), interact = 0.8,
                                label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]),
                                rho = rho, withinf = fwithin)
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs, nsims = iterations)
  aov_time <- system.time(twoway_simulation_testing(simdat, test = "ANOVA"))[3]
  per_time <- system.time(twoway_simulation_testing(simdat, test = "permutation"))[3]
  rank_time <- system.time(twoway_simulation_testing(simdat, test = "rank"))[3]
  expect_true(all(c(per_time>rank_time), c(rank_time>aov_time)))

  fwithin <- "fA"
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 1.1, fBeffect = 0.8, plot = FALSE, sdratio = 0.2, sdproportional = FALSE,
                                groupswinteraction = matrix(c(2,2,3,3,1,4), 3, 2, byrow = TRUE), interact = 0.8,
                                label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]),
                                rho = rho, withinf = fwithin)
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs, nsims = iterations)
  aov_time <- system.time(twoway_simulation_testing(simdat, test = "ANOVA"))[3]
  per_time <- system.time(twoway_simulation_testing(simdat, test = "permutation"))[3]
  rank_time <- system.time(twoway_simulation_testing(simdat, test = "rank"))[3]
  expect_true(all(c(per_time>rank_time), c(rank_time>aov_time)))

  fwithin <- "both"
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 1.1, fBeffect = 0.8, sdratio = 0.2, sdproportional = FALSE,
                                groupswinteraction = matrix(c(2,2,3,3,1,4), 3, 2, byrow = TRUE), interact = 0.8,
                                label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]),
                                rho = rho, withinf = fwithin)
  simdat <- twoway_simulation_correlated(group_size = 10, matrices_obj = refs, nsims = iterations)
  aov_time <- system.time(twoway_simulation_testing(simdat, test = "ANOVA"))[3]
  per_time <- system.time(twoway_simulation_testing(simdat, test = "permutation"))[3]
  rank_time <- system.time(twoway_simulation_testing(simdat, test = "rank"))[3]
  expect_true(all(c(per_time>rank_time), c(rank_time>aov_time)))
})

test_that("correct within factor is used", {
  nlevfA <- 3
  nlevfB <- 3
  nlist <- list(groups=LETTERS[1:nlevfA], time=letters[1:nlevfB])
  gwint <- matrix(c(2,3, 3, 3), 2, 2, byrow = TRUE)
  iterations <- 50
  rho <- -0.9

  fwithin <- "fB"
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 2, fBeffect = 1, sdratio = 0.2,
                                ##groupswinteraction = gwint, interact = 1.25,
                                label_list = nlist,
                                rho = rho, withinf = fwithin)
  group_size <- 5
  simdatB <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs, nsims = iterations)

  fwithin <- "fA"
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 2, fBeffect = 1, sdratio = 0.2,
                                ##groupswinteraction = gwint, interact = 1.25,
                                label_list = nlist,
                                rho = rho, withinf = fwithin)
  simdatA <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs, nsims = iterations)

  fwithin <- "both"
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 2, fBeffect = 1, sdratio = 0.2,
                                ##groupswinteraction = gwint, interact = 1.25,
                                label_list = nlist,
                                rho = rho, withinf = fwithin)
  group_size <- 15
  simdatboth <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs, nsims = iterations)

  res_fB <- twoway_simulation_testing(simdatB, test = "rank")
  res_fA <- twoway_simulation_testing(simdatA, test = "rank")
  res_both <- twoway_simulation_testing(simdatboth, test = "rank")
  expect_gt(res_fA[1,2], 0.9)
  expect_lt(res_fA[2,2], 0.15)
  expect_gt(res_fB[1,2], 0.9)
  expect_lt(res_fB[2,2], 0.15)
  expect_gt(res_both[1,2], 0.9)
  expect_lt(res_both[2,2], 0.15)

  res_fB <- twoway_simulation_testing(simdatB, test = "ANOVA")
  res_fA <- twoway_simulation_testing(simdatA, test = "ANOVA")
  res_both <- twoway_simulation_testing(simdatboth, test = "ANOVA")
  expect_gt(res_fA[1,2], 0.9)
  expect_lt(res_fA[2,2], 0.15)
  expect_gt(res_fB[1,2], 0.9)
  expect_lt(res_fB[2,2], 0.15)
  expect_gt(res_both[1,2], 0.9)
  expect_lt(res_both[2,2], 0.15)

  res_fB <- twoway_simulation_testing(simdatB, test = "permutation")
  res_fA <- twoway_simulation_testing(simdatA, test = "permutation")
  res_both <- twoway_simulation_testing(simdatboth, test = "permutation")
  expect_gt(res_fA[1,2], 0.9)
  expect_lt(res_fA[2,2], 0.15)
  expect_gt(res_fB[1,2], 0.9)
  expect_lt(res_fB[2,2], 0.15)
  expect_gt(res_both[1,2], 0.9)
  expect_lt(res_both[2,2], 0.15)
})


test_that("unbalanced designs", {
  nlevfA <- 3
  nlevfB <- 3
  nlist <- list(groups=LETTERS[1:nlevfA], time=letters[1:nlevfB])
  gwint <- matrix(c(2,3, 3, 3), 2, 2, byrow = TRUE)
  iterations <- 20
  rho <- -0.9
  group_size <- matrix(rep(8:6, 3), 3, 3, byrow = TRUE)
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 1.3, fBeffect = 1, sdratio = 0.2,
                                groupswinteraction = gwint, interact = 1.1,
                                label_list = nlist)
  sindatind <- twoway_simulation_independent(group_size = group_size, matrices_obj = refs, balanced = FALSE)
  resind <- twoway_simulation_testing(sindatind)
  expect_equal(dim(resind), c(3,6))
})
