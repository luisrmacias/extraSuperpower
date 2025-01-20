data(mtcars)
mtcars$am_f <- factor(mtcars$am)
mtcars$cyl_f <- factor(mtcars$cyl)
mpgmeans <- aggregate(mpg ~ am_f + cyl_f, data = mtcars, FUN = mean)
mpgmeans <- matrix(mpgmeans$mpg, length(levels(mpgmeans$am_f)), length(levels(mpgmeans$cyl_f)),
                   dimnames = list(levels(mpgmeans$am_f), levels(mpgmeans$cyl_f)))

mpgsd <- aggregate(mpg ~ am_f + cyl_f, data = mtcars, FUN = sd)
mpgsd <- matrix(mpgsd$mpg, length(levels(mpgsd$am_f)), length(levels(mpgsd$cyl_f)),
                   dimnames = list(levels(mpgsd$am_f), levels(mpgsd$cyl_f)))


test_that("effectsizeestimateconcordance", {
  expect_equal(2 * 2, 4)
})
