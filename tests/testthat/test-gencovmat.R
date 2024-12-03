test_that("factor A correlation", {
  faeff <- 2
  fA <- 3
  fbeff <- 3
  fB <- 2
  rho <- 0.9
  modobj <- calculate_mean_matrix(refmean = 10, nlfA = fA, nlfB = fB,
                                   fAeffect = faeff, fBeffect = fbeff,
                                   rho = rho, withinf = "fA",
                                   plot = FALSE)
  levacor <- cor_mat[grep("_a", names(cor_mat[,1])),grep("_a", names(cor_mat[1,]))]
  levbcor <- cor_mat[grep("_b", names(cor_mat[,1])),grep("_b", names(cor_mat[1,]))]
  expect_true(all(c(levacor, levbcor)==rho | c(levacor, levbcor)==1))
})
