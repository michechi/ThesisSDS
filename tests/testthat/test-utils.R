test_that("find_mode returns the correct mode for a unimodal vector", {
  expect_equal(find_mode(c(1, 1, 1, 2, 3)), 1)
  expect_equal(find_mode(c(3, 3, 2, 1)), 3)
})

test_that("find_mode handles ties", {
  set.seed(42)
  result <- find_mode(c(1, 1, 2, 2))
  expect_true(result %in% c(1, 2))
})

test_that("find_mode works with a single element", {
  expect_equal(find_mode(7), 7)
})

test_that(".symmetrize_matrix produces a symmetric matrix", {
  m <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), 3, 3)
  m[1, 2] <- 5
  m[1, 3] <- 7
  m[2, 3] <- 9
  result <- mfmsbm:::.symmetrize_matrix(m)
  expect_true(isSymmetric(result))
  expect_equal(result[2, 1], 5)
  expect_equal(result[3, 1], 7)
  expect_equal(result[3, 2], 9)
})
