test_that("generate_sbm produces correct dimensions", {
  set.seed(1)
  sbm <- generate_sbm(30, 3, 0.5, 0.1, "balanced")
  expect_equal(dim(sbm$A0), c(30, 30))
  expect_length(sbm$z0, 30)
  expect_equal(dim(sbm$Q), c(3, 3))
})

test_that("generate_sbm adjacency matrix is symmetric and binary", {
  set.seed(2)
  sbm <- generate_sbm(20, 2, 0.5, 0.1, "balanced")
  expect_true(isSymmetric(sbm$A0))
  expect_true(all(sbm$A0 %in% c(0, 1)))
  expect_true(all(diag(sbm$A0) == 0))
})

test_that("generate_sbm balanced produces roughly equal clusters", {
  set.seed(3)
  sbm <- generate_sbm(30, 3, 0.5, 0.1, "balanced")
  tab <- table(sbm$z0)
  expect_equal(length(tab), 3)
  expect_true(all(tab == 10))
})

test_that("generate_sbm unbalanced produces unequal clusters", {
  set.seed(4)
  sbm <- generate_sbm(40, 2, 0.5, 0.1, "unbalanced")
  tab <- table(sbm$z0)
  expect_equal(length(tab), 2)
  expect_true(tab[1] != tab[2])
})

test_that("generate_sbm validates type_network", {
  expect_error(generate_sbm(10, 2, 0.5, 0.1, "invalid"))
})
