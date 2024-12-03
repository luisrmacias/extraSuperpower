test_that("factor A correlation", {
  faeff <- 2
  fA <- 3
  fbeff <- 3
  fB <- 2
  rho <- 0.9
  fwhithin <-"fA"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                   fAeffect = faeff, fBeffect = fbeff,
                                   plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencovmat(meansd_mats[[1]], meansd_mats[[2]],
                           rho = rho, withinf =fwhithin,
                           nlfA = facdim[1], nlfB = facdim[2])[[1]]
  levacor <- cor_mat[grep("_a", names(cor_mat[,1])),grep("_a", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("_b", names(cor_mat[,1])),grep("_b", names(cor_mat[1,]))]
  expect_true(all(c(levacor, levbcor)==rho | c(levacor, levbcor)==1))
})

test_that("factor B correlation", {
  faeff <- 2
  fA <- 2
  fbeff <- 3
  fB <- 3
  rho <- 0.9
  fwhithin <- "fB"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencovmat(meansd_mats[[1]], meansd_mats[[2]],
                       rho = rho, withinf =fwhithin,
                       nlfA = facdim[1], nlfB = facdim[2])[[1]]
  levacor <- cor_mat[grep("A_", names(cor_mat[,1])),grep("A_", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("B_", names(cor_mat[,1])),grep("B_", names(cor_mat[1,]))]
  expect_true(all(c(levacor, levbcor)==rho | c(levacor, levbcor)==1))
})

test_that("both factor correlation", {
  faeff <- 2
  fA <- 2
  fbeff <- 3
  fB <- 3
  rho <- 0.9
  fwhithin <- "both"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencovmat(meansd_mats[[1]], meansd_mats[[2]],
                       rho = rho, withinf =fwhithin,
                       nlfA = facdim[1], nlfB = facdim[2])[[1]]
  triupper <- cor_mat[upper.tri(cor_mat)]
  trilower <- cor_mat[lower.tri(cor_mat)]
  expect_true(all(c(triupper, trilower)==rho & all(diag(cor_mat)==1)))
})

test_that("factor A covariance", {
  faeff <- 2
  fA <- 2
  fbeff <- 3
  fB <- 2
  rho <- 0.9
  fwhithin <- "fA"
  ##constant standard deviation
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff,
                                    sdproportional = FALSE,
                                    plot = FALSE)
  sd <- meansd_mats[[2]]
  facdim <- dim(meansd_mats[[1]])
  covcor_mats <- gencovmat(meansd_mats[[1]], sd,
                       rho = rho, withinf =fwhithin,
                       nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- covcor_mats[[2]]
  levacov <- cov_mat[grep("_a", names(cov_mat[,1])),grep("_a", names(cov_mat[1,]))]
  levacov <- c(levacov[upper.tri(levacov)], levacov[lower.tri(levacov)])
  levbcov <- cov_mat[grep("_b", names(cov_mat[,1])),grep("_b", names(cov_mat[1,]))]
  levbcov <- c(levbcov[upper.tri(levbcov)], levbcov[lower.tri(levbcov)])
  expect_true(all(c(levacov, levbcov)==rho*sd^2 & diag(cov_mat)==sd^2))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(identical(linpk::cor2cov(covcor_mats[[1]], sd), strippedres))

  ##cell variant standard deviation
  ## meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
  ##                                      fAeffect = faeff, fBeffect = fbeff,
  ##                                      plot = FALSE)
  ## sd_mat <- meansd_mats[[2]]
  ## facdim <- dim(meansd_mats[[1]])
  ## covcor_mats <- gencovmat(meansd_mats[[1]], sd_mat,
  ##                          rho = rho, withinf =fwhithin,
  ##                          nlfA = facdim[1], nlfB = facdim[2])
  ## cov_mat <- covcor_mats[[2]]
  ## levacov <- cov_mat[grep("_a", names(cov_mat[,1])),grep("_a", names(cov_mat[1,]))]
  ## levacov <- c(levacov[upper.tri(levacov)], levacov[lower.tri(levacov)])
  ## levbcov <- cov_mat[grep("_b", names(cov_mat[,1])),grep("_b", names(cov_mat[1,]))]
  ## levbcov <- c(levbcov[upper.tri(levbcov)], levbcov[lower.tri(levbcov)])
})


test_that("factor B covariance", {
  faeff <- 2
  fA <- 2
  fbeff <- 3
  fB <- 3
  rho <- 0.9
  fwhithin <- "fB"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       sdproportional = FALSE,
                                       plot = FALSE)
  sd <- meansd_mats[[2]]
  facdim <- dim(meansd_mats[[1]])
  covcor_mats <- gencovmat(meansd_mats[[1]], sd,
                           rho = rho, withinf =fwhithin,
                           nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- covcor_mats[[2]]
  levacov <- cov_mat[grep("A_", names(cov_mat[,1])),grep("A_", names(cov_mat[1,]))]
  levacov <- c(levacov[upper.tri(levacov)], levacov[lower.tri(levacov)])
  levbcov <- cov_mat[grep("B_", names(cov_mat[,1])),grep("B_", names(cov_mat[1,]))]
  levbcov <- c(levbcov[upper.tri(levbcov)], levbcov[lower.tri(levbcov)])
  expect_true(all(c(levacov, levbcov)==rho*sd^2 & diag(cov_mat)==sd^2))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(identical(linpk::cor2cov(covcor_mats[[1]], sd), strippedres))
})

test_that("both factor covariance", {
  faeff <- 2
  fA <- 4
  fbeff <- 3
  fB <- 5
  rho <- 0.9
  fwhithin <- "both"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       sdproportional = FALSE,
                                       plot = FALSE)
  sd <- meansd_mats[[2]]
  facdim <- dim(meansd_mats[[1]])
  covcor_mats <- gencovmat(meansd_mats[[1]], sd,
                           rho = rho, withinf =fwhithin,
                           nlfA = facdim[1], nlfB = facdim[2])
  cov_mat <- covcor_mats[[2]]
  levacov <- cov_mat[grep("A_", names(cov_mat[,1])),grep("A_", names(cov_mat[1,]))]
  levacov <- c(levacov[upper.tri(levacov)], levacov[lower.tri(levacov)])
  levbcov <- cov_mat[grep("B_", names(cov_mat[,1])),grep("B_", names(cov_mat[1,]))]
  levbcov <- c(levbcov[upper.tri(levbcov)], levbcov[lower.tri(levbcov)])
  expect_true(all(c(levacov, levbcov)==rho*sd^2 & diag(cov_mat)==sd^2))
  strippedres <- matrix(cov_mat, prod(facdim), prod(facdim))
  expect_true(identical(linpk::cor2cov(covcor_mats[[1]], sd), strippedres))
})
