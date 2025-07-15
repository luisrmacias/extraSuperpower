test_that("vectors change accordingly", {
  expect_equal(genvecs(change = 2, reps = 4), c(3/3, 4/3, 5/3, 6/3))
  expect_equal(genvecs(change = 2, reps = 4, bystart = FALSE, scaler = 10), 1:4)
})

test_that("data entry type", {
  expect_error(genvecs(change = 1.2, reps = 4.2))
  expect_error(genvecs(change = "A", reps = 4))
  expect_error(genvecs(change = -2, reps = 4))
})
