test_that("matrices dimensions", {
  iA <- 2
  iB <- 3

  a <- 5
  b <- 4
  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a)
  fBvec <- genvecs(change = iB, reps = b)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames)
  expect_equal(dim(mean_mat), c(a, b))

  fAvec <- genvecs(change = iA, reps = a, bystart = FALSE, scaler = 10)
  fBvec <- genvecs(change = iB, reps = b, bystart = FALSE, scaler = 10)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames, bystart = FALSE)
  expect_equal(dim(mean_mat), c(a, b))

  a <- 4
  b <- 5
  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a)
  fBvec <- genvecs(change = iB, reps = b)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames)
  expect_equal(dim(mean_mat), c(a, b))

  fAvec <- genvecs(change = iA, reps = a, bystart = FALSE, scaler = 10)
  fBvec <- genvecs(change = iB, reps = b, bystart = FALSE, scaler = 10)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames, bystart = FALSE)
  expect_equal(dim(mean_mat), c(a, b))
})

test_that("factor A stepwise effect ratio", {
  iA <- 2
  iB <- 1

  a <- 2
  b <- 2

  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a)
  fBvec <- genvecs(change = iB, reps = b)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames)
  f1.1 <- mean_mat[1,1]
  faend.1 <- mean_mat[a,1]
  expect_equal(faend.1/f1.1, iA)

  iA <- 5
  a  <- 5
  fAvec <- genvecs(change = iA, reps = a)
  fnames <- list(fA=letters[1:a], fB=1:b)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames)
  f1.1 <- mean_mat[1,1]
  faend.1 <- mean_mat[a,1]
  expect_equal(faend.1/f1.1, iA)
})


test_that("factor B stepwise effect", {
  iA <- 1
  iB <- 2

  a <- 2
  b <- 2

  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a)
  fBvec <- genvecs(change = iB, reps = b)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames)
  f1.1 <- mean_mat[1,1]
  fbend.1 <- mean_mat[1,b]
  expect_equal(fbend.1/f1.1, iB)

  iA <- 2
  iB <- 3

  a <- 5
  b <- 4

  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a)
  fBvec <- genvecs(change = iB, reps = b)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames)
  f1.1 <- mean_mat[1,1]
  fbend.1 <- mean_mat[1,b]
  expect_equal(fbend.1/f1.1, iB)
})

test_that("combined A and B stepwise effect", {
  iA <- 2
  iB <- 3

  a <- 5
  b <- 4

  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a)
  fBvec <- genvecs(change = iB, reps = b)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames)
  betaA <- unique(diff(fAvec))
  betaB <- diff(fBvec)[1]
  cellab <-  1+(betaA*(a-1))+(betaB*(b-1))
  aend.bend <- mean_mat[a,b]
  expect_equal(cellab, aend.bend)

  iA <- 2
  iB <- 3

  a <- 3
  b <- 3

  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a)
  fBvec <- genvecs(change = iB, reps = b)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames)
  betaA <- unique(diff(fAvec))
  betaB <- diff(fBvec)[1]
  cellab <-  1+(betaA*(a-1))+(betaB*(b-1))
  aend.bend <- mean_mat[a,b]
  expect_equal(cellab, aend.bend)
})


test_that("factor A end effect ratio", {
  iA <- 2
  iB <- 1
  a  <- 2
  b  <- 4
  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a, bystart = FALSE, scaler = 10)
  fBvec <- genvecs(change = iB, reps = b, bystart = FALSE, scaler = 10)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames, bystart = FALSE)
  f1.1 <- mean_mat[1,1]
  fanext.1 <- mean_mat[2,1]
  expect_equal(fanext.1/f1.1, iA)

  a <- 5
  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a, bystart = FALSE, scaler = 10)
  fBvec <- genvecs(change = iB, reps = b, bystart = FALSE, scaler = 10)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames, bystart = FALSE)
  f1.1 <- mean_mat[1,1]
  fanext.1 <- mean_mat[2,1]
  expect_equal(fanext.1/f1.1, iA)
})

test_that("factor B end effect ratio", {

  iA <- 1
  iB <- 3
  a <- 4
  b <- 2
  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a, bystart = FALSE, scaler = 10)
  fBvec <- genvecs(change = iB, reps = b, bystart = FALSE, scaler = 10)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames, bystart = FALSE)
  f1.1 <- mean_mat[1,1]
  fbnext.1 <- mean_mat[1,2]
  expect_equal(fbnext.1/f1.1, iB)

  iA <- 2
  iB <- 3
  a  <- 5
  b  <- 4
  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a, bystart = FALSE, scaler = 10)
  fBvec <- genvecs(change = iB, reps = b, bystart = FALSE, scaler = 10)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames, bystart = FALSE)
  f1.1 <- mean_mat[1,1]
  fbnext.1 <- mean_mat[1,2]
  expect_equal(fbnext.1/f1.1, iB)
})

test_that("combined A and B end effect", {
  iA <- 2
  iB <- 3

  a <- 5
  b <- 4

  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a, bystart = FALSE, scaler = 10)
  fBvec <- genvecs(change = iB, reps = b, bystart = FALSE, scaler = 10)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames, bystart = FALSE)
  betaA <- unique(diff(fAvec))
  betaB <- diff(fBvec)[1]
  cellab <-  1+(betaA*(a-1))+(betaB*(b-1))
  aend.bend <- mean_mat[a,b]
  expect_equal(cellab, aend.bend)

  iA <- 2
  iB <- 3

  a <- 3
  b <- 3

  fnames <- list(fA=letters[1:a], fB=1:b)
  fAvec <- genvecs(change = iA, reps = a, bystart = FALSE, scaler = 10)
  fBvec <- genvecs(change = iB, reps = b, bystart = FALSE, scaler = 10)
  mean_mat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = iA, a, b, label_list = fnames, bystart = FALSE)
  betaA <- unique(diff(fAvec))
  betaB <- diff(fBvec)[1]
  cellab <-  1+(betaA*(a-1))+(betaB*(b-1))
  aend.bend <- mean_mat[a,b]
  expect_equal(cellab, aend.bend)
})
