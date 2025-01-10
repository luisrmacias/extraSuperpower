faeff <- 2
fA <- 3
fbeff <- 0.5
fB <- 2
rho <- 0.8

separate_covariance_level <- function(cmat, level)
{
  levcov <- cmat[grep(level, names(cmat[,1])),grep(level, names(cmat[1,]))]
  levcov <- c(levcov[upper.tri(levcov)], levcov[lower.tri(levcov)])
}

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
  levelcolumns <- paste0("_", letters[1:fB])
  levcov <- sapply(levelcolumns, function(x) separate_covariance_level(cov_mat, x))
  expect_true(all(all(as.vector(levcov)==rho*sd^2) & all(diag(cov_mat)==sd^2)))
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
  levelcolumns <- paste0(LETTERS[1:fA], "_")
  levcov <- sapply(levelcolumns, function(x) separate_covariance_level(cov_mat, x))
  expect_true(all(all(as.vector(levcov)==rho*sd^2) & all(diag(cov_mat)==sd^2)))
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
  levcov <- cov_mat[upper.tri(cov_mat)|lower.tri(cov_mat)]
  expect_true(all(all(levcov==rho*sd^2) & all(diag(cov_mat)==sd^2)))
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

##covariance with correlation gradient
##factor A is within factor

rho <- c(0.4,0.8)
test_that("factor A covariance with correlation gradient", {
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
  levelcolumns <- paste0("_", letters[1:fB])
  levcov <- sapply(levelcolumns, function(x) separate_covariance_level(cov_mat, x))
  covariances_expected <- rep(c(seq(rho[1], rho[2], length.out=fA-1), seq(rho[2], rho[1], length.out=fA-1)[-1])*sd^2, 2)
  covariances_expected <- rep(covariances_expected, fB)
  expect_true(all(all.equal(as.vector(levcov), covariances_expected) & all.equal(as.vector(diag(cov_mat)), rep(sd^2, prod(facdim)))))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(all.equal(linpk::cor2cov(cor_mat, sd), strippedres))

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

test_that("both factor covariance with correlation gradient", {
  fwithin <- "both"
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
  levelcolumns <- paste0(LETTERS[1:fA], "_")
  levcov <- sapply(levelcolumns, function(x) separate_covariance_level(cov_mat, x))
  covariances_expected <- rep(c(seq(rho[1], rho[2], length.out=fB-1), seq(rho[2], rho[1], length.out=fB-1)[-1])*sd^2, 2)
  covariances_expected <- rep(covariances_expected, fA)
  expect_true(all(all.equal(as.vector(levcov), covariances_expected) & all.equal(as.vector(diag(cov_mat)), rep(sd^2, prod(fA, fB)))))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(all.equal(linpk::cor2cov(cor_mat, sd), strippedres))

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

test_that("factor B covariance with correlation gradient", {
  fwithin <- "fB"
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
  levcov <- cov_mat[upper.tri(cov_mat)|lower.tri(cov_mat)]
  covariances_expected <- rep(c(seq(rho[1], rho[2], length.out=fB-1), seq(rho[2], rho[1], length.out=fB-1)[-1])*sd^2, 2)
  covariances_expected <- rep(covariances_expected, fA)
  expect_true(all(all.equal(as.vector(levcov), covariances_expected) & all.equal(as.vector(diag(cov_mat)), rep(sd^2, prod(fA, fB)))))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(all.equal(linpk::cor2cov(cor_mat, sd), strippedres))

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
