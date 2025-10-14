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
  expect_message(twoway_simulation_testing(simdat),regexp = "ANOVA")
  expect_message(twoway_simulation_testing(simdat, test = "permutation"), regexp = "permutation")
  expect_message(twoway_simulation_testing(simdat, test = "rank"), regexp = "rank")
})

test_that("independent unbalanced testing", {
  nlevfA <- 3
  nlevfB <- 3
  nlist <- list(groups=LETTERS[1:nlevfA], time=letters[1:nlevfB])
  gwint <- matrix(c(2,3, 3, 3), 2, 2, byrow = TRUE)
  iterations <- 20
  rho <- -0.9
  group_size <- matrix(rep(8:6, 3), 3, 3, byrow = TRUE)
  set.seed(14534)
  refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 1.3, fBeffect = 1, sdratio = 0.2,
                                groupswinteraction = gwint, interact = 1.1,
                                label_list = nlist)
  sindatind <- twoway_simulation_independent(group_size = group_size, matrices_obj = refs, balanced = FALSE, nsims = iterations)
  resind_ANOVA <- twoway_simulation_testing(sindatind)
  resind_perm <- twoway_simulation_testing(sindatind, test = "permutation")
  resind_rank <- twoway_simulation_testing(sindatind, test = "rank")
  expect_gt(resind_ANOVA[1,4], resind_ANOVA[2,4])
  expect_gt(resind_ANOVA[2,4], resind_ANOVA[3,4])
  expect_gt(resind_perm[1,4], resind_perm[2,4])
  expect_gt(resind_perm[2,4], resind_perm[3,4])
  expect_gt(resind_rank[1,4], resind_rank[2,4])
  expect_gt(resind_rank[2,4], resind_rank[3,4])
})

nlevfA <- 3
nlevfB <- 3
nlist <- list(groups=LETTERS[1:nlevfA], time=letters[1:nlevfB])
gwint <- matrix(c(2,3, 3, 3), 2, 2, byrow = TRUE)
iterations <- 35
rho <- -0.9

fwithin <- "fB"
refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                              fAeffect = 2, fBeffect = 1, sdratio = 0.2,
                              groupswinteraction = gwint, interact = 1.25,
                              label_list = nlist,
                              rho = rho, withinf = fwithin)
group_size <- 5
simdatB <- twoway_simulation_correlated(group_size = group_size,
                                        matrices_obj = refs, nsims = iterations)

unbal_group_size <- matrix(c(6,6,5,6,5,4,6,4,3), 3, 3, byrow = TRUE)
unbal_simdatB <- twoway_simulation_correlated(group_size = unbal_group_size,
                                              matrices_obj = refs, balanced = FALSE,
                                              loss = "sequential", nsims = iterations)


fwithin <- "fA"
refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                              fAeffect = 2, fBeffect = 1, sdratio = 0.2,
                              groupswinteraction = gwint, interact = 1.25,
                              label_list = nlist,
                              rho = rho, withinf = fwithin)
simdatA <- twoway_simulation_correlated(group_size = group_size,
                                        matrices_obj = refs, nsims = iterations)
unbal_simdatA <- twoway_simulation_correlated(group_size = unbal_group_size,
                                              matrices_obj = refs, balanced = FALSE,
                                              loss = "sequential", nsims = iterations)


fwithin <- "both"
refs <- calculate_mean_matrix(refmean = 1, nlfA = nlevfA, nlfB = nlevfB,
                              fAeffect = 2, fBeffect = 1, sdratio = 0.2,
                              groupswinteraction = gwint, interact = 1.25,
                              label_list = nlist,
                              rho = rho, withinf = fwithin)
group_size <- 15
simdatboth <- twoway_simulation_correlated(group_size = group_size,
                                           matrices_obj = refs, nsims = iterations)

unbal_simdatboth <- twoway_simulation_correlated(group_size = unbal_group_size,
                                              matrices_obj = refs, balanced = FALSE,
                                              loss = "sequential", nsims = iterations)

test_that("correct repeated measures test is performed", {

  expect_message(fB_AOV <- twoway_simulation_testing(simdatB),regexp = "ANOVA")
  expect_message(fB_per <- twoway_simulation_testing(simdatB, test = "permutation"), regexp = "permutation")
  expect_message(fB_rank <- twoway_simulation_testing(simdatB, test = "rank"), regexp = "rank")

  expect_message(fA_AOV <- twoway_simulation_testing(simdatA),regexp = "ANOVA")
  expect_message(fA_per <- twoway_simulation_testing(simdatA, test = "permutation"), regexp = "permutation")
  expect_message(fA_rank <- twoway_simulation_testing(simdatA, test = "rank"), regexp = "rank")

  expect_message(fAB_AOV <- twoway_simulation_testing(simdatboth),regexp = "ANOVA")
  expect_message(fAB_per <- twoway_simulation_testing(simdatboth, test = "permutation"), regexp = "permutation")
  expect_message(fAB_rank <- twoway_simulation_testing(simdatboth, test = "rank"), regexp = "rank")
})



test_that("correct within factor is used", {
  res_fB <- twoway_simulation_testing(simdatB, test = "rank")
  res_fA <- twoway_simulation_testing(simdatA, test = "rank")
  res_both <- twoway_simulation_testing(simdatboth, test = "rank")

  expect_gt(res_fA[1,2], 0.9)
  expect_gt(res_fA[2,2], 0.9)
  expect_gt(res_fB[1,2], 0.9)
  expect_lt(res_fB[2,2], 0.7)
  expect_gt(res_both[1,2], 0.9)
  expect_gt(res_both[2,2], 0.7)

  res_fB <- twoway_simulation_testing(simdatB, test = "ANOVA")
  res_fA <- twoway_simulation_testing(simdatA, test = "ANOVA")
  res_both <- twoway_simulation_testing(simdatboth, test = "ANOVA")
  expect_gt(res_fA[1,2], 0.9)
  expect_gt(res_fA[2,2], 0.9)
  expect_gt(res_fB[1,2], 0.9)
  expect_lt(res_fB[2,2], 0.5)
  expect_gt(res_both[1,2], 0.9)
  expect_gt(res_both[2,2], 0.7)

  res_fB <- twoway_simulation_testing(simdatB, test = "permutation")
  res_fA <- twoway_simulation_testing(simdatA, test = "permutation")
  res_both <- twoway_simulation_testing(simdatboth, test = "permutation")
  expect_gt(res_fA[1,2], 0.9)
  expect_gt(res_fA[2,2], 0.9)
  expect_gt(res_fB[1,2], 0.9)
  expect_lt(res_fB[2,2], 0.5)
  expect_gt(res_both[1,2], 0.9)
  expect_gt(res_both[2,2], 0.7)
})


test_that("repeated unbalanced testing", {
  res_fB <- twoway_simulation_testing(unbal_simdatB, test = "rank")
  res_fA <- twoway_simulation_testing(unbal_simdatA, test = "rank")
  res_both <- twoway_simulation_testing(unbal_simdatboth, test = "rank")
  expect_gt(res_fA[1,4], 0.9)
  expect_gt(res_fA[2,4], 0.9)
  expect_gt(res_fB[1,4], 0.9)
  expect_lt(res_fB[2,4], 0.5)
  expect_gt(res_both[1,4], 0.9)
  expect_lt(res_both[2,4], 0.5)

  res_fB <- twoway_simulation_testing(unbal_simdatB, test = "ANOVA")
  res_fA <- twoway_simulation_testing(unbal_simdatA, test = "ANOVA")
  res_both <- twoway_simulation_testing(unbal_simdatboth, test = "ANOVA")
  expect_gt(res_fA[1,4], 0.9)
  expect_gt(res_fA[2,4], 0.8)
  expect_gt(res_fB[1,4], 0.9)
  expect_lt(res_fB[2,4], 0.5)
  expect_gt(res_both[1,4], 0.5)
  expect_lt(res_both[2,4], 0.3)

  res_fB <- twoway_simulation_testing(unbal_simdatB, test = "permutation")
  res_fA <- twoway_simulation_testing(unbal_simdatA, test = "permutation")
  res_both <- twoway_simulation_testing(unbal_simdatboth, test = "permutation")
  expect_gt(res_fA[1,4], 0.9)
  expect_gt(res_fA[2,4], 0.8)
  expect_gt(res_fB[1,4], 0.9)
  expect_lt(res_fB[2,4], 0.5)
  expect_gt(res_both[1,4], 0.9)
  expect_lt(res_both[2,4], 0.5)

})
