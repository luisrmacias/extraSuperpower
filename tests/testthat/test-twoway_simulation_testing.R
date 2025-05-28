test_that("correct independent sample test is performed", {
  nlevfA <- 4
  nlevfB <- 4
  group_size <- 5
  gwint <- matrix(c(1,3,1,4,2,4), 3, 2, byrow = TRUE)
  nlist <- list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB])
  matlist <- calculate_mean_matrix(refmean = 3, nlfA = nlevfA, nlfB = nlevfB, fAeffect = 2, fBeffect = 1.2,
                                   groupswinteraction = gwint, interact = 1.4, sdproportional = FALSE,
                                   label_list = nlist)
  #matrix(1:8, nlevfA, nlevfB, dimnames = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]))
  #sd.mat <- 2
  #matlist <- list(mean.mat=mean.mat, sd.mat=sd.mat)
  #group_size <- matrix(c(6,3), nlevfA, nlevfB, byrow = TRUE)
  simdat <- twoway_simulation_independent(group_size = group_size, matrices_obj = matlist, nsims = 10)
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
  iterations <- 100
  rho <- -0.9
  fwithin <- "fB"
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 1.1, fBeffect = 0.8, plot = FALSE, sdratio = 0.2, sdproportional = FALSE,
                                groupswinteraction = matrix(c(2,2,3,3,1,4), 3, 2, byrow = TRUE), interact = 0.8,
                                label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]),
                                rho = rho, withinf = fwithin)
  refs$sigmat <- Matrix::nearPD(refs$sigmat)$mat
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
  refs$sigmat <- Matrix::nearPD(refs$sigmat)$mat
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs, nsims = iterations)
  aov_time <- system.time(twoway_simulation_testing(simdat, test = "ANOVA"))[3]
  per_time <- system.time(twoway_simulation_testing(simdat, test = "permutation"))[3]
  rank_time <- system.time(twoway_simulation_testing(simdat, test = "rank"))[3]
  expect_true(all(c(per_time>rank_time), c(rank_time>aov_time)))

  fwithin <- "both"
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 1.1, fBeffect = 0.8, plot = FALSE, sdratio = 0.2, sdproportional = FALSE,
                                groupswinteraction = matrix(c(2,2,3,3,1,4), 3, 2, byrow = TRUE), interact = 0.8,
                                label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]),
                                rho = rho, withinf = fwithin)
  refs$sigmat <- Matrix::nearPD(refs$sigmat)$mat
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs, nsims = iterations)
  aov_time <- system.time(twoway_simulation_testing(simdat, test = "ANOVA"))[3]
  per_time <- system.time(twoway_simulation_testing(simdat, test = "permutation"))[3]
  rank_time <- system.time(twoway_simulation_testing(simdat, test = "rank"))[3]
  expect_true(all(c(per_time>rank_time), c(rank_time>aov_time)))
})

