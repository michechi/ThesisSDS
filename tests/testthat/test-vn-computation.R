test_that("logsumexp handles normal values", {
  expect_equal(logsumexp(log(2), log(3)), log(5), tolerance = 1e-10)
})

test_that("logsumexp handles -Inf", {
  expect_equal(logsumexp(-Inf, 0), 0)
  expect_equal(logsumexp(-Inf, -Inf), -Inf)
})

test_that("falling_factorial computes correctly", {
  expect_equal(falling_factorial(5, 3), 60)  # 5 * 4 * 3
  expect_equal(falling_factorial(4, 1), 4)
  expect_equal(falling_factorial(4, 0), 1)
})

test_that("rising_factorial computes correctly", {
  expect_equal(rising_factorial(2, 3), 24)  # 2 * 3 * 4
  expect_equal(rising_factorial(1, 4), 24)  # 1 * 2 * 3 * 4
  expect_equal(rising_factorial(5, 0), 1)
  expect_equal(rising_factorial(5, 1), 5)
})

test_that("log_vn_miller returns a vector of correct length", {
  result <- log_vn_miller(1, 10, 15, lambda = 1)
  expect_length(result, 15)
  expect_true(all(is.finite(result[1:10])))
  expect_true(all(result[11:15] == -Inf))
})

test_that("log_vn_gnedin returns finite values for valid inputs", {
  result <- log_vn_gnedin(50, 3, 0.5)
  expect_true(is.finite(result))
})

test_that("log_vn_approx returns finite values", {
  result <- log_vn_approx(1000, 3, 1)
  expect_true(is.finite(result))
})
