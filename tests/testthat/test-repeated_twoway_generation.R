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
                                    fAeffect = faeff, fBeffect = fbeff, plot = FALSE)

  expect_error(twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 3))

  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = rho, withinf = fwithin)

  expect_no_error(twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 3))

  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    rho = rho, withinf = fwithin, plot = FALSE)

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
  expect_error(twoway_simulation_correlated(group_size = group_size, mean_mat, nsims = 3,
                                             balanced = TRUE))
})

test_that("factor A as within factor", {
  faeff <- 2
  fA <- 2
  fbeff <- 2
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


test_that("both factors within", {
  faeff <- 10
  fA <- 2
  fbeff <- 0.5
  fB <- 2
  rho <- 0.9
  fwithin <- "both"
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
  expect_gt(abs(fAcor1),
            rho-0.1)
  fAcor2 <- cor.test(sim$y[sim$cond=="A_b"], sim$y[sim$cond=="B_b"])$estimate
  expect_gt(abs(fAcor2),
            rho-0.1)
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
  expect_true(sum(p.adjust(pvals)>0.05)>15)
})

test_that("skewed and truncated input checks work", {
  nlevfA <- 2
  nlevfB <- 4
  group_size <- 30
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
                                               distribution = "truncated.normal", shape = 2))

  expect_warning(twoway_simulation_correlated(group_size = group_size, matrices_obj = matlist, nsims = iterations,
                                              shape = 2))

  expect_error(twoway_simulation_correlated(group_size = group_size, matrices_obj = matlist, nsims = iterations,
                                              distribution = "skewed", shape = "fB"))

  expect_error(twoway_simulation_correlated(group_size = group_size, matrices_obj = matlist, nsims = iterations,
                                            distribution = "skewed", shape = seq(2, 6, 2)))


  fwithin="fB"
  group_size <- rep(seq(25, 10, -5), 2)
  group_size <- matrix(group_size, nlevfA, nlevfB, byrow = TRUE)
  matlist <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                   fAeffect = 2, fBeffect = 2, plot = FALSE, sdratio = 0.3,
                                   label_list = list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB]),
                                   rho = rho, withinf = fwithin)

  expect_no_warning(twoway_simulation_correlated(group_size = group_size, matrices_obj = matlist,
                                                 nsims = iterations, balanced = FALSE, loss="sequential",
                                                 distribution = "truncated.normal", inferior_limit = 3))

  expect_warning(twoway_simulation_correlated(group_size = group_size, matrices_obj = matlist,
                                                 nsims = iterations, balanced = FALSE, loss="sequential",
                                                 distribution = "truncated.normal", superior_limit = suplim))
  })


test_that("simulated values are skewed", {
  nlevfA <- 2
  nlevfB <- 4
  label_list <- list(groups=LETTERS[1:nlevfA], treatment=letters[1:nlevfB])
  group_size <- 150
  iterations <- 1
  rho <- 0.3
  fwithin <- "fB"
  refs <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 2, fBeffect = 0.5, plot = FALSE,
                                sdproportional = FALSE, sdratio = 0.1,
                                label_list = label_list, rho = rho, withinf = fwithin)
  set.seed(160724)
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs,
                                         nsims = iterations,
                                         distribution = "skewed", shape = rep(0.8, nlevfB))$simulated_data

  distpvals <- ks.test(simdat$y[simdat$cond=="A_a"], "pnorm", refs$mean.mat[1,1], refs$sd.mat)$p.value
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="A_b"], "pnorm", refs$mean.mat[1,2], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="A_c"], "pnorm", refs$mean.mat[1,3], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="B_a"], "pnorm", refs$mean.mat[2,1], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="B_b"], "pnorm", refs$mean.mat[2,2], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="B_d"], "pnorm", refs$mean.mat[2,4], refs$sd.mat)$p.value)
  expect_true(any(distpvals<0.05))

  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs,
                                         nsims = iterations,
                                         distribution = "skewed", shape = rep(4, nlevfA))$simulated_data
  distpvals <- ks.test(simdat$y[simdat$cond=="A_a"], "pnorm", refs$mean.mat[1,1], refs$sd.mat)$p.value
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="A_b"], "pnorm", refs$mean.mat[1,2], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="A_c"], "pnorm", refs$mean.mat[1,3], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="B_a"], "pnorm", refs$mean.mat[2,1], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="B_b"], "pnorm", refs$mean.mat[2,2], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="B_d"], "pnorm", refs$mean.mat[2,4], refs$sd.mat)$p.value)
  expect_true(any(distpvals>0.05))

  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs,
                                         nsims = iterations,
                                         distribution = "skewed", shape = rep(3, prod(nlevfA, nlevfB)))$simulated_data
  distpvals <- ks.test(simdat$y[simdat$cond=="A_a"], "pnorm", refs$mean.mat[1,1], refs$sd.mat)$p.value
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="A_b"], "pnorm", refs$mean.mat[1,2], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="A_c"], "pnorm", refs$mean.mat[1,3], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="B_a"], "pnorm", refs$mean.mat[2,1], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="B_b"], "pnorm", refs$mean.mat[2,2], refs$sd.mat)$p.value)
  distpvals <- c(distpvals, ks.test(simdat$y[simdat$cond=="B_d"], "pnorm", refs$mean.mat[2,4], refs$sd.mat)$p.value)
  expect_true(any(distpvals>0.05))
})

test_that("loss is random", {
  nlevfA <- 3
  nlevfB <- 4
  label_list <- list(groups=LETTERS[1:nlevfA], time=letters[1:nlevfB])
  n <- 10
  group_size <- rep(10,prod(nlevfA, nlevfB))
  group_size <- matrix(group_size, nlevfA, nlevfB, byrow = TRUE)
  group_size[1,] <- group_size[1,]-c(0,2,4,5)
  group_size[2,] <- group_size[2,]-c(0,2,3,4)
  group_size[3,] <- group_size[3,]-c(0,0,4,5)

  iterations <- 1
  rho <- 0.3
  fwithin <- "fB"
  refs <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 2, fBeffect = 0.5, plot = FALSE,
                                sdproportional = TRUE, sdratio = 0.1,
                                label_list = label_list, rho = rho, withinf = fwithin)

  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs,
                                         balanced = FALSE, loss = "random", nsims = iterations)$simulated_data
  seq_check <- all(apply(table(simdat$time, simdat$subject), 2,
                         function(x) all(x[1]>=x[2] & x[x]>=x[3] & x[3]>=x[4])))
  res <- table(simdat$groups, simdat$time)
  expect_true(all(!seq_check & all(res==group_size)))
})


test_that("loss is sequential", {
  nlevfA <- 3
  nlevfB <- 4
  label_list <- list(groups=LETTERS[1:nlevfA], time=letters[1:nlevfB])
  n <- 10
  group_size <- rep(10,prod(nlevfA, nlevfB))
  group_size <- matrix(group_size, nlevfA, nlevfB, byrow = TRUE)
  group_size[1,] <- group_size[1,]-c(0,2,4,5)
  group_size[2,] <- group_size[2,]-c(0,2,3,4)
  group_size[3,] <- group_size[3,]-c(0,0,4,5)

  iterations <- 1
  rho <- 0.3
  fwithin <- "fB"
  refs <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 2, fBeffect = 0.5, plot = FALSE,
                                sdproportional = TRUE, sdratio = 0.1,
                                label_list = label_list, rho = rho, withinf = fwithin)

  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs,
                                         balanced = FALSE, loss = "sequential", nsims = iterations)$simulated_data
  seq_check <- all(apply(table(simdat$time, simdat$subject), 2,
                         function(x) all(x[1]>=x[2] & x[x]>=x[3] & x[3]>=x[4])))
  res <- table(simdat$groups, simdat$time)
  expect_true(all(seq_check & all(res==group_size)))

  fwithin <- "fA"
  group_size <- rep(10,prod(nlevfA, nlevfB))
  group_size <- matrix(group_size, nlevfA, nlevfB, byrow = TRUE)
  group_size[,2] <- group_size[,2]-c(1,2,3)
  group_size[,3] <- group_size[,3]-c(2,2,4)
  group_size[,4] <- group_size[,4]-c(2,2,4)
  refs <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 2, fBeffect = 0.5, plot = FALSE,
                                sdproportional = TRUE, sdratio = 0.1,
                                label_list = label_list, rho = rho, withinf = fwithin)

  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs,
                                         balanced = FALSE, loss = "sequential", nsims = iterations)$simulated_data
  seq_check <- all(apply(table(simdat$groups, simdat$subject), 2,
                         function(x) all(x[1]>=x[2] & x[2]>=x[3])))
  res <- table(simdat$groups, simdat$time)
  expect_true(all(seq_check & all(res==group_size)))


  fwithin <- "both"
  group_size <- rep(10,prod(nlevfA, nlevfB))
  group_size <- matrix(group_size, nlevfA, nlevfB, byrow = TRUE)
  group_size[1,] <- group_size[1,]-c(0,0,1,1)
  group_size[2,] <- group_size[2,]-c(1,1,1,1)
  group_size[3,] <- group_size[3,]-c(1,1,2,2)

  refs <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 2, fBeffect = 0.5, plot = FALSE,
                                sdproportional = TRUE, sdratio = 0.1,
                                label_list = label_list, rho = rho, withinf = fwithin)
  simdat <- twoway_simulation_correlated(group_size = group_size, matrices_obj = refs,
                               nsims = iterations, balanced = FALSE, loss="sequential",
                               distribution = "truncated.normal", inferior_limit = 3)$simulated_data
  seq_check <- all(apply(table(simdat$time, simdat$subject), 2,
                         function(x) all(x[1]>=x[2] & x[x]>=x[3] & x[3]>=x[4])))
  res <- table(simdat$groups, simdat$time)
  expect_true(all(seq_check & all(res==group_size)))
})

test_that("loss input check works", {
  nlevfA <- 3
  nlevfB <- 4
  label_list <- list(groups=LETTERS[1:nlevfA], time=letters[1:nlevfB])
  n <- 10
  group_size <- rep(10,prod(nlevfA, nlevfB))
  group_size <- matrix(group_size, nlevfA, nlevfB, byrow = TRUE)
  group_size[1,] <- group_size[1,]-c(0,2,4,5)
  group_size[2,] <- group_size[2,]-c(0,2,3,4)
  group_size[3,] <- group_size[3,]-c(0,0,4,5)

  iterations <- 1
  rho <- 0.3
  fwithin <- "fB"
  refs <- calculate_mean_matrix(refmean = 10, nlfA = nlevfA, nlfB = nlevfB,
                                fAeffect = 2, fBeffect = 0.5, plot = FALSE,
                                sdproportional = TRUE, sdratio = 0.1,
                                label_list = label_list, rho = rho, withinf = fwithin)

  expect_error(twoway_simulation_correlated(group_size = group_size, matrices_obj = refs,
                                         balanced = FALSE, nsims = iterations))
})
