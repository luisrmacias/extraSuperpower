test_that("vectors change accordingly", {
  expect_equal(genvecs(change = 2, reps = 4), c(3/3, 4/3, 5/3, 6/3))
})
