test_that("labels_to_matrix produces correct binary matrix", {
  m <- labels_to_matrix(c(1, 1, 2, 3))
  expect_equal(dim(m), c(4, 3))
  expect_equal(m[1, ], c(1, 0, 0))
  expect_equal(m[3, ], c(0, 1, 0))
  expect_equal(m[4, ], c(0, 0, 1))
})

test_that("coclustering_matrix produces correct dimensions and range", {
  z_post <- matrix(c(1, 1, 2, 1, 2, 2, 1, 1, 2), nrow = 3)
  cc <- coclustering_matrix(z_post)
  expect_equal(dim(cc), c(3, 3))
  expect_true(all(cc >= 0 & cc <= 1))
  expect_true(all(diag(cc) == 1))
})

test_that("edge_probability produces correct dimensions and range", {
  Y <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
  ep <- edge_probability(c(1, 1, 2), Y, 1, 1)
  expect_equal(dim(ep), c(3, 3))
  expect_true(all(diag(ep) == 0))
  expect_true(all(ep >= 0 & ep <= 1))
})

test_that("posterior_k returns correct dimensions", {
  pk <- posterior_k(1, 20, 5, 5)
  expect_equal(dim(pk), c(5, 5))
  expect_true(all(pk >= 0))
})
