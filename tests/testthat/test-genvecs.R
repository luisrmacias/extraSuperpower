test_that("vectors change accordingly", {
  expect_equal(extraSuperpower:::genvecs(change = 2, reps = 4), c(3/3, 4/3, 5/3, 6/3))
  expect_equal(extraSuperpower:::genvecs(change = 2, reps = 4, bystart = TRUE, scaler = 10), 1:4)
})

