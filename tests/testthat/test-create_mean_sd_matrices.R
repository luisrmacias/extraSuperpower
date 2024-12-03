test_that("matrices dimensions", {
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = 5, nlfB = 4,
                                    fAeffect = 2, fBeffect = 3, plot = FALSE)[[1]]
  expect_equal(dim(mean_mat), c(5, 4))
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = 4, nlfB = 5,
                                    fAeffect = 2, fBeffect = 3, plot = FALSE)[[1]]
  expect_equal(dim(mean_mat), c(4, 5))
})

test_that("standard deviation proportionality", {
  mean_mats <- calculate_mean_matrix(refmean = 10, nlfA = 5, nlfB = 4,
                                     fAeffect = 2, fBeffect = 3,
                                     plot = FALSE)
  expect_true(all(mean_mats[[2]]/mean_mats[[1]]==0.2))
  sdcoef <- 0.1
  mean_mats <- calculate_mean_matrix(refmean = 10, nlfA = 5, nlfB = 4,
                                    fAeffect = 2, fBeffect = 3, sdratio = sdcoef,
                                    plot = FALSE)
  expect_true(all(mean_mats[[2]]/mean_mats[[1]]==sdcoef))
})

test_that("constant standard deviation", {
  mean_mats <- calculate_mean_matrix(refmean = 10, nlfA = 5, nlfB = 4,
                                     fAeffect = 2, fBeffect = 3,
                                     sdproportional = FALSE,
                                     plot = FALSE)
  expect_equal(mean(mean_mats[[1]])*0.2, mean_mats[[2]])
  sdratio <- 0.1
  mean_mats <- calculate_mean_matrix(refmean = 10, nlfA = 5, nlfB = 4,
                                     fAeffect = 2, fBeffect = 3,
                                     sdratio= sdratio, sdproportional = FALSE,
                                     plot = FALSE)
  expect_equal(mean(mean_mats[[1]])*sdratio, mean_mats[[2]])
})

test_that("factor A stepwise effect ratio", {
  faeff <- 2
  fA <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = 2,
                                    fAeffect = faeff, fBeffect = 1, plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  faend.1 <- mean_mat[fA,1]
  expect_equal(faend.1/f1.1, faeff)
  fA <- 5
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = 4,
                                    fAeffect = faeff, fBeffect = 3, plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  faend.1 <- mean_mat[fA,1]
  expect_equal(faend.1/f1.1, faeff)
})


test_that("factor B stepwise effect", {
  faeff <- 1
  fA <- 2
  fbeff <- 3
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff, plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  fbend.1 <- mean_mat[1,fB]
  expect_equal(fbend.1/f1.1, fbeff)

  faeff <- 2
  fA <- 5
  fbeff <- 3
  fB <- 4
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff, plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  fbend.1 <- mean_mat[1,fB]
  expect_equal(fbend.1/f1.1, fbeff)
})

test_that("factor A end effect ratio", {
  faeff <- 2
  fA <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = 4,
                                    fAeffect = faeff, fBeffect = 1, endincrement = TRUE,
                                    plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  fanext.1 <- mean_mat[2,1]
  expect_equal(fanext.1/f1.1, faeff)

  fA <- 5
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = 4,
                                    fAeffect = faeff, fBeffect = 3, endincrement = TRUE,
                                    plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  fanext.1 <- mean_mat[2,1]
  expect_equal(fanext.1/f1.1, faeff)
})

test_that("factor B end effect ratio", {
  faeff <- 1
  fA <- 4
  fbeff <- 3
  fB <- 2
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff, endincrement = TRUE,
                                    plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  fbnext.1 <- mean_mat[1,2]
  expect_equal(fbnext.1/f1.1, fbeff)

  faeff <- 2
  fA <- 5
  fbeff <- 3
  fB <- 4
  mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                    fAeffect = faeff, fBeffect = fbeff, endincrement = TRUE,
                                    plot = FALSE)[[1]]
  f1.1 <- mean_mat[1,1]
  fbnext.1 <- mean_mat[1,2]
  expect_equal(fbnext.1/f1.1, fbeff)
})

test_that("interaction modelled as defined with stepwise effect", {
  faeff <- 2
  fA <- 5
  fbeff <- 3
  fB <- 4
  ginteract <- expand.grid(1:2,3:4)
  intereff <- 1.5
  int_mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                        fAeffect = faeff, fBeffect = fbeff,
                                        groupswinteraction = ginteract, interact = intereff,
                                        plot = FALSE)[[1]]

  noint_mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                          fAeffect = faeff, fBeffect = fbeff,
                                          plot = FALSE)[[1]]
  rowindices <- unique(ginteract[,1])
  colindices <- unique(ginteract[,2])
  ratio <- int_mean_mat[rowindices,colindices]/noint_mean_mat[rowindices, colindices]
  expect_true(all(ratio==intereff))
})

test_that("interaction modelled as defined with end effect", {
  faeff <- 2
  fA <- 5
  fbeff <- 3
  fB <- 4
  ginteract <- expand.grid(1:2,3:4)
  intereff <- 1.5
  int_mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                        fAeffect = faeff, fBeffect = fbeff,
                                        endincrement = TRUE,
                                        groupswinteraction = ginteract, interact = intereff,
                                        plot = FALSE)[[1]]

  noint_mean_mat <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                          fAeffect = faeff, fBeffect = fbeff,
                                          endincrement = TRUE,
                                          plot = FALSE)[[1]]
  rowindices <- unique(ginteract[,1])
  colindices <- unique(ginteract[,2])
  ratio <- int_mean_mat[rowindices,colindices]/noint_mean_mat[rowindices, colindices]
  expect_true(all(ratio==intereff))
})
