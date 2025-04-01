faeff <- 2
fA <- 3
fbeff <- 0.5
fB <- 2

test_that("different dimensions of mean and correlation throws error", {
  fwithin <-"fA"
  rho <- 0.8
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  expect_error(gencorrelationmat(meansd_mats[[1]],
                    rho = rho, withinf =fwithin,
                    nlfA = facdim[2], nlfB = facdim[2]))

  expect_error(gencorrelationmat(meansd_mats[[1]],
                                 rho = rho, withinf =fwithin,
                                 nlfA = facdim[1], nlfB = facdim[1]))
  })

test_that("transferability of level labels", {
  fwithin <-"fA"
  rho <- 0.8
  ##test with 'small' number of levels
  factors_levels <- list(treatment=letters[1:fA], time=c("before", "after"))
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       label_list = factors_levels,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  vec1 <- rep(dimnames(meansd_mats[[1]])[[1]], each=fB)
  vec2 <- rep(dimnames(meansd_mats[[1]])[[2]], fA)
  expect_identical(paste(vec1, vec2, sep = "_"), colnames(cor_mat))

  ##test with 'large' number of levels
  fA <- 5
  fB <- 3
  factors_levels <- list(treatment=letters[1:fA], time=c("before", "during", "after"))
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       label_list = factors_levels,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  vec1 <- rep(dimnames(meansd_mats[[1]])[[1]], each=fB)
  vec2 <- rep(dimnames(meansd_mats[[1]])[[2]], fA)
  expect_identical(paste(vec1, vec2, sep = "_"), colnames(cor_mat))

  ##test with same number of levels for both factors
  fA <- 4
  fB <- 4
  factors_levels <- list(treatment=letters[1:fA], time=c("before", "early", "late", "after"))
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       label_list = factors_levels,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  vec1 <- rep(dimnames(meansd_mats[[1]])[[1]], each=fB)
  vec2 <- rep(dimnames(meansd_mats[[1]])[[2]], fA)
  expect_identical(paste(vec1, vec2, sep = "_"), colnames(cor_mat))
})

test_that("label incompatibilities are detected", {
  ##factors are inverted
  fA <- 3
  fB <- 3
  fwithin <-"fA"
  rho <- 0.8
  factors_levels <- list(treatment=letters[1:fA], time=c("before", "during", "after"))
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       label_list = factors_levels,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  expect_error(gencorrelationmat(meansd_mats[[1]],
                                 rho = rho, withinf =fwithin,
                                 nlfA = facdim[1], nlfB = facdim[2],
                                 label_list = rev(factors_levels)))

  ##different classes
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2],
                               label_list = list(treatment=as.character(factors_levels$treatment), time=factors_levels$time))
  expect_true(is(cor_mat, "matrix"))
})

test_that("error is shown when rho is out of format", {
  fwithin <-"fA"
  rho <- c(0.4, 0.2, 0.3)
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  expect_error(gencorrelationmat(meansd_mats[[1]], rho=rho, withinf = fwithin, nlfA = 3, nlfB = 2))
})

test_that("correlation input format and design are consistent", {
  fwithin <-"fA"
  rho <- matrix(c(0.4, 0.2, 0.3, 0.1), 2, 2)
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  expect_error(gencorrelationmat(meansd_mats[[1]], rho=rho, withinf = fwithin, nlfA = 3, nlfB = 2))
})

test_that("warning appears when late labelling is used", {
  fwithin <-"fA"
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  levels_list <- list(treatment=c("control", "MedA", "MedB"), time=c("before", "after"))
  expect_warning(gencorrelationmat(meansd_mats[[1]], rho=0.8, withinf = fwithin, nlfA = 3, nlfB = 2,
                                    label_list = levels_list))
})

##fixed correlation
test_that("correlation when factor A is the repeated factor", {
  fwithin <-"fA"
  #positive correlation
  rho <- 0.8
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  levacor <- cor_mat[grep("_a", names(cor_mat[,1])),grep("_a", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("_b", names(cor_mat[,1])),grep("_b", names(cor_mat[1,]))]
  expect_identical(unique(c(levbcor, levacor)), c(1, rho))
  #negative correlation
  rho <- -0.8
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  levacor <- cor_mat[grep("_a", names(cor_mat[,1])),grep("_a", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("_b", names(cor_mat[,1])),grep("_b", names(cor_mat[1,]))]
  expect_identical(unique(c(levbcor, levacor)), c(1, rho))
})

test_that("correlation when factor B is the repeated factor", {
  faeff <- 2
  fA <- 2
  fbeff <- 3
  fB <- 3

  fwithin <- "fB"
  ##positive correlation
  rho <- 0.9
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  levacor <- cor_mat[grep("A_", names(cor_mat[,1])),grep("A_", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("B_", names(cor_mat[,1])),grep("B_", names(cor_mat[1,]))]
  expect_identical(unique(c(levbcor, levacor)), c(1, rho))

  ##negative correlation
  rho <- -0.9
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  levacor <- cor_mat[grep("A_", names(cor_mat[,1])),grep("A_", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("B_", names(cor_mat[,1])),grep("B_", names(cor_mat[1,]))]
  expect_identical(unique(c(levbcor, levacor)), c(1, rho))
})


test_that("correlation when both factors are repeated", {
  faeff <- 2
  fA <- 2
  fbeff <- 3
  fB <- 3
  fwithin <- "both"
  #positive correlation
  rho <- 0.9
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  triupper <- cor_mat[upper.tri(cor_mat)]
  trilower <- cor_mat[lower.tri(cor_mat)]
  expect_true(all(c(triupper, trilower)==rho & all(diag(cor_mat)==1)))

  #negative correlation
  rho <- -0.9
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  triupper <- cor_mat[upper.tri(cor_mat)]
  trilower <- cor_mat[lower.tri(cor_mat)]
  expect_true(all(c(triupper, trilower)==rho & all(diag(cor_mat)==1)))
})


##correlation gradient

test_that("correlation gradient when factor A is the repeated factor", {
  fwithin <-"fA"
  fA <- 4
  fB <- 3
  #positive correlation
  rho <- c(-0.8, 0)
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  levacor <- cor_mat[grep("_a", names(cor_mat[,1])),grep("_a", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("_b", names(cor_mat[,1])),grep("_b", names(cor_mat[1,]))]
  topvector <- c(1, seq(rho[1], rho[2], length.out=fA-1))
  expect_identical(as.vector(levacor[1,]), topvector)
  expect_identical(as.vector(levacor[fA,]), rev(topvector))
  expect_identical(as.vector(levbcor[,1]), topvector)
  expect_identical(as.vector(levbcor[,fA]), rev(topvector))
  #negative correlation
  rho <- c(-0.8, 0)
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  levacor <- cor_mat[grep("_a", names(cor_mat[,1])),grep("_a", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("_b", names(cor_mat[,1])),grep("_b", names(cor_mat[1,]))]
  topvector <- c(1, seq(rho[1], rho[2], length.out=fA-1))
  expect_identical(as.vector(levacor[1,]), topvector)
  expect_identical(as.vector(levacor[fA,]), rev(topvector))
  expect_identical(as.vector(levbcor[,1]), topvector)
  expect_identical(as.vector(levbcor[,fA]), rev(topvector))
})

test_that("correlation gradient when factor B is the repeated factor", {
  fwithin <-"fB"
  fA <- 4
  fB <- 6
  #positive correlation gradient
  rho <- c(0.8, 0.2)
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  levacor <- cor_mat[grep("A_", names(cor_mat[,1])),grep("A_", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("B_", names(cor_mat[,1])),grep("B_", names(cor_mat[1,]))]
  topvector <- c(1, seq(rho[1], rho[2], length.out=fB-1))
  expect_identical(as.vector(levacor[1,]), topvector)
  expect_identical(as.vector(levacor[fB,]), rev(topvector))
  expect_identical(as.vector(levbcor[,1]), topvector)
  expect_identical(as.vector(levbcor[,fB]), rev(topvector))
  #negative correlation
  rho <- c(-0.8, 0)
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  levacor <- cor_mat[grep("A_", names(cor_mat[,1])),grep("A_", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("B_", names(cor_mat[,1])),grep("B_", names(cor_mat[1,]))]
  topvector <- c(1, seq(rho[1], rho[2], length.out=fB-1))
  expect_identical(as.vector(levacor[1,]), topvector)
  expect_identical(as.vector(levacor[fB,]), rev(topvector))
  expect_identical(as.vector(levbcor[,1]), topvector)
  expect_identical(as.vector(levbcor[,fB]), rev(topvector))
})

test_that("correlation gradient when both factors are repeated", {
  fwithin <-"both"
  fA <- 4
  fB <- 6
  #positive correlation gradient
  rho <- matrix(c(0.8, -0.8, 0.2, -0.1), 2, 2)
  meansd_mats <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                       fAeffect = faeff, fBeffect = fbeff,
                                       plot = FALSE)
  facdim <- dim(meansd_mats[[1]])
  cor_mat <- gencorrelationmat(meansd_mats[[1]],
                               rho = rho, withinf =fwithin,
                               nlfA = facdim[1], nlfB = facdim[2])
  levacor <- cor_mat[grep("A_", names(cor_mat[,1])),grep("A_", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("B_", names(cor_mat[,1])),grep("B_", names(cor_mat[1,]))]
  expect_identical(as.vector(levacor), as.vector(levbcor))
  corleva <- cor_mat[grep("_a", names(cor_mat[,1])),grep("_a", names(cor_mat[1,]))]
  corlevb <- cor_mat[grep("_b", names(cor_mat[,1])),grep("_b", names(cor_mat[1,]))]
  expect_identical(as.vector(corleva), as.vector(corlevb))
  topvector <- c(1, seq(rho[2,1], rho[2,2], length.out=fB-1))
  expect_identical(as.vector(levacor[1,]), topvector)
  expect_identical(as.vector(levacor[fB,]), rev(topvector))
  topvector <- c(1, seq(rho[1,1], rho[1,2], length.out=fA-1))
  expect_identical(as.vector(corleva[1,]), topvector)
  expect_identical(as.vector(corleva[fA,]), rev(topvector))
})
