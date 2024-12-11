faeff <- 2
fA <- 3
fbeff <- 0.5
fB <- 2
rho <- 0.8

test_that("factor A covariance", {
  fwithin <- "fA"
  ##constant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       sdproportional = FALSE,
                                       plot = FALSE)
  sd <- meansd_mats[[2]]
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])

  cov_mat <- gencovariancemat(cor_mat, sd, withinf =fwithin,
                       nlfA = facdim[1], nlfB = facdim[2])
  levacov <- cov_mat[grep("_a", names(cov_mat[,1])),grep("_a", names(cov_mat[1,]))]
  levacov <- c(levacov[upper.tri(levacov)], levacov[lower.tri(levacov)])
  levbcov <- cov_mat[grep("_b", names(cov_mat[,1])),grep("_b", names(cov_mat[1,]))]
  levbcov <- c(levbcov[upper.tri(levbcov)], levbcov[lower.tri(levbcov)])
  expect_true(all(all(c(levacov, levbcov)==rho*sd^2) & all(diag(cov_mat)==sd^2)))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(identical(linpk::cor2cov(cor_mat, sd), strippedres))

  ##cell variant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                        fAeffect = faeff, fBeffect = fbeff,
                                        plot = FALSE)
  sd_mat <- meansd_mats[[2]]
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(correlation_matrix = cor_mat, sd_matrix = sd_mat,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  expect_true(identical(cov2cor(cov_mat), cor_mat))
})


test_that("factor B covariance", {
  fwithin <- "fB"
  #constant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       sdproportional = FALSE,
                                       plot = FALSE)
  sd <- meansd_mats[[2]]
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                           rho = rho, withinf =fwithin,
                           nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(cor_mat, sd_matrix = sd,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  levacov <- cov_mat[grep("A_", names(cov_mat[,1])),grep("A_", names(cov_mat[1,]))]
  levacov <- c(levacov[upper.tri(levacov)], levacov[lower.tri(levacov)])
  levbcov <- cov_mat[grep("B_", names(cov_mat[,1])),grep("B_", names(cov_mat[1,]))]
  levbcov <- c(levbcov[upper.tri(levbcov)], levbcov[lower.tri(levbcov)])
  expect_true(all(all(c(levacov, levbcov)==rho*sd^2) & all(diag(cov_mat)==sd^2)))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(identical(linpk::cor2cov(cor_mat, sd), strippedres))

  ##cell variant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  sd_mat <- meansd_mats[[2]]
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(correlation_matrix = cor_mat, sd_matrix = sd_mat,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  expect_true(identical(cov2cor(cov_mat), cor_mat))
})

test_that("both factor covariance with constant correlation", {
  fwithin <- "both"
  #constant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       sdproportional = FALSE,
                                       plot = FALSE)
  sd <- meansd_mats[[2]]
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                           rho = rho, withinf =fwithin,
                           nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(cor_mat, sd_matrix = sd,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  levacov <- cov_mat[grep("A_", names(cov_mat[,1])),grep("A_", names(cov_mat[1,]))]
  levacov <- c(levacov[upper.tri(levacov)], levacov[lower.tri(levacov)])
  levbcov <- cov_mat[grep("B_", names(cov_mat[,1])),grep("B_", names(cov_mat[1,]))]
  levbcov <- c(levbcov[upper.tri(levbcov)], levbcov[lower.tri(levbcov)])
  expect_true(all(all(c(levacov, levbcov)==rho*sd^2) & all(diag(cov_mat)==sd^2)))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(identical(linpk::cor2cov(cor_mat, sd), strippedres))
  ##cell variant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  sd_mat <- meansd_mats[[2]]
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(correlation_matrix = cor_mat, sd_matrix = sd_mat,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  expect_true(all.equal(cov2cor(cov_mat), cor_mat))
})
