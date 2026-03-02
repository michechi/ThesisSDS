test_that("cogibbs_poisson runs without error on a small case", {
  set.seed(123)
  sbm <- generate_sbm(10, 2, 0.5, 0.1, "balanced")
  lv <- log_vn_miller(1, 10, 20)
  z0 <- c(1, 2, 3, sample(1:3, 7, replace = TRUE))
  fit <- cogibbs_poisson(5, 3, sbm$A0, 10, 1, 1, 1, lv, z0)

  expect_type(fit, "list")
  expect_equal(dim(fit$z_post), c(10, 5))
  expect_length(fit$num_k, 5)
  expect_true(all(fit$num_k >= 1))
})

test_that("cogibbs_gnedin runs without error on a small case", {
  set.seed(456)
  sbm <- generate_sbm(10, 2, 0.5, 0.1, "balanced")
  z0 <- c(1, 2, 3, sample(1:3, 7, replace = TRUE))
  fit <- cogibbs_gnedin(5, 3, sbm$A0, 10, 1, 1, 1, 0.5, z0)

  expect_type(fit, "list")
  expect_equal(dim(fit$z_post), c(10, 5))
  expect_length(fit$num_k, 5)
  expect_true(all(fit$num_k >= 1))
})
