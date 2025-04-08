faeff <- 2
fA <- 3
fbeff <- 0.5
fB <- 4
rho <- 0.8

separate_covariance_level <- function(cmat, level)
{
  levcov <- cmat[grep(level, names(cmat[,1])),grep(level, names(cmat[1,]))]
  c(levcov[upper.tri(levcov)], levcov[lower.tri(levcov)])
}


test_that("sd matrix dimension check works", {
  fwithin <- "fA"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  sd <- meansd_mats[[2]]
  sd <- sd[,-4]
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  expect_error(gencovariancemat(cor_mat, sd, withinf =fwithin,
                                nlfA = facdim[1], nlfB = facdim[2]))
  sd <- meansd_mats[[2]]
  sd <- sd[-3,]
  expect_error(gencovariancemat(cor_mat, sd, withinf =fwithin,
                                nlfA = facdim[1], nlfB = facdim[2]))
})

test_that("sd matrix name assignment works", {
  fwithin <- "fA"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])

  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])

  sd <- meansd_mats[[2]]
  dimnames(sd) <- list(treatment=c("Ctrl", "MedA", "MedB"),
                       time=paste("t", 1:4, sep="_"))
  expect_message(cov_mat <- gencovariancemat(cor_mat, sd, withinf =fwithin,
                              nlfA = facdim[1], nlfB = facdim[2]))
  sdnames <- dimnames(sd)
  sdnames <- expand.grid(sdnames[[2]], sdnames[[1]])
  sdnames <- paste(sdnames$Var2, sdnames$Var1, sep = "_")

  expect_equal(dimnames(cov_mat)[[1]], sdnames)
})

test_that("sd matrix name check works", {
  fwithin <- "fA"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])

  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])

  sd <- meansd_mats[[2]]
  dimnames(sd) <- list(treatment=c("Ctrl", "MedA", "MedB"),
                       time=paste("t", 1:4, sep="_"))
  expect_error(gencovariancemat(cor_mat, sd, withinf =fwithin,
                                nlfA = facdim[1], nlfB = facdim[2],
                                label_list = list(intervention=c("Placebo", "Tiritin", "Toroton"),
                                                  time = paste("day", 1:4, sep = "_"))))
})

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
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(correlation_matrix = cor_mat, sd_matrix = sd_mat,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  expect_true(all.equal(cov2cor(cov_mat), cor_mat))
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
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(correlation_matrix = cor_mat, sd_matrix = sd_mat,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  expect_true(all.equal(cov2cor(cov_mat), cor_mat))
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
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(correlation_matrix = cor_mat, sd_matrix = sd_mat,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  expect_true(all.equal(cov2cor(cov_mat), cor_mat))
})

####################################################################################
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
  expect_true(all(all.equal(cor_mat*sd^2, cov_mat) & all.equal(as.vector(diag(cov_mat)), rep(sd^2, prod(fA, fB)))))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(all.equal(linpk::cor2cov(cor_mat, sd), strippedres))

  ##cell variant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  sd_mat <- meansd_mats[[2]]
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(correlation_matrix = cor_mat, sd_matrix = sd_mat,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  expect_true(all.equal(cov2cor(cov_mat), cor_mat))
})

##both factors have a different, constant correlation
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
  fastcovmat <- cor_mat*tcrossprod(as.vector(matrix(sd, nrow = fA, ncol = fB)))
  levelcolumns <- paste0("_", letters[1:fB])
  levcov <- sapply(levelcolumns, function(x) separate_covariance_level(fastcovmat, x))
  expect_true(all(all(as.vector(levcov)==rho[1]*sd^2) & all(diag(fastcovmat)==sd^2)))

  levelcolumns <- paste0(LETTERS[1:fA], "_")
  levcov <- sapply(levelcolumns, function(x) separate_covariance_level(fastcovmat, x))
  expect_true(all(as.vector(levcov)==rho[2]*sd^2))

  strippedres <- matrix(fastcovmat, prod(facdim), prod(facdim))
  expect_true(all.equal(linpk::cor2cov(cor_mat, sd), strippedres))

  ##cell variant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  sd_mat <- meansd_mats[[2]]
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(correlation_matrix = cor_mat, sd_matrix = sd_mat,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])

  fastcovmat <- cor_mat*tcrossprod(as.vector(sd_mat))
  expect_true(all.equal(cov2cor(fastcovmat), cor_mat))
})

rho <- matrix(c(0.8, 0.7, 0.4, 0.3), 2, 2)

##both factors have a different, varying correlations
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
  fastcovmat <- cor_mat*tcrossprod(as.vector(matrix(sd, nrow = fA, ncol = fB)))
  levelcolumns <- paste0("_", letters[1:fB])
  levcov <- sapply(levelcolumns, function(x) separate_covariance_level(fastcovmat, x))
  if(all(sapply(2:fB, function(x) identical(levcov[,1], levcov[,x]))))
  {
    levcov <- levcov[,1]
  }
  covvector <- rho[1,]*sd^2
  covariances_expected <- rep(c(covvector, rev(covvector[1:(length(covvector)-1)])),2)
  expect_true(all(all.equal(levcov,covariances_expected) & all(diag(fastcovmat)==sd^2)))

  levelcolumns <- paste0(LETTERS[1:fA], "_")
  levcov <- sapply(levelcolumns, function(x) separate_covariance_level(fastcovmat, x))
  if(all(sapply(2:fA, function(x) identical(levcov[,1], levcov[,x]))))
  {
    levcov <- levcov[,1]
  }
  # covvector <- rho[2,]*sd^2
  # covariances_expected <- rep(c(covvector, rev(covvector[1:(length(covvector)-1)])),2)
  # expect_true(all(all.equal(levcov,covariances_expected) & all(diag(fastcovmat)==sd^2)))
  strippedres <- matrix(fastcovmat, prod(facdim), prod(facdim))
  expect_true(all.equal(linpk::cor2cov(cor_mat, sd), strippedres))

  ##cell variant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  sd_mat <- meansd_mats[[2]]
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(correlation_matrix = cor_mat, sd_matrix = sd_mat,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])

  fastcovmat <- cor_mat*tcrossprod(as.vector(sd_mat))
  expect_true(all.equal(cov2cor(fastcovmat), cor_mat))
})



test_that("within factor consistency check works", {
  rho <- 0.8
  fwithin <- "fA"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  sd <- meansd_mats[[2]]
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])

  expect_error(gencovariancemat(cor_mat, sd, withinf ="fB",
                                nlfA = facdim[1], nlfB = facdim[2]))
})

##test correlation and covariance concordance when standard deviation is 0

test_that("covariance matrix with 0 standard deviation", {
  fwithin <- "fB"
  rho <- 0.8
  #constant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  sd <- meansd_mats[[2]]
  sd[1,3:4] <- 0
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- gencovariancemat(cor_mat, sd_matrix = sd,
                              withinf = fwithin, nlfA = facdim[1], nlfB = facdim[2])
  expect_equal(cor_mat*tcrossprod(as.vector(t(sd))), cov_mat)
})

