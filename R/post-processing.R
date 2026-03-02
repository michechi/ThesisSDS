#' Convert Cluster Labels to a Binary Membership Matrix
#'
#' Transforms a vector of cluster labels into a binary `V x H` matrix where
#' `M[v, h] = 1` if node `v` belongs to cluster `h`.
#'
#' @param clust_lab Integer vector of length `V` with cluster labels.
#' @return A binary matrix of dimension `V x H` where `H = max(clust_lab)`.
#' @export
#' @examples
#' labels_to_matrix(c(1, 1, 2, 2, 3))
labels_to_matrix <- function(clust_lab) {
  V <- length(clust_lab)
  H <- max(clust_lab)
  M <- matrix(0L, nrow = V, ncol = H)
  for (v in seq_len(V)) {
    M[v, clust_lab[v]] <- 1L
  }
  M
}

#' Compute the Posterior Co-clustering Matrix
#'
#' Computes the `V x V` posterior co-clustering matrix from MCMC samples of
#' cluster assignments. Element `[v, u]` gives the fraction of iterations in
#' which nodes `v` and `u` were assigned to the same cluster.
#'
#' @param z_post A `V x N_iter` integer matrix of posterior cluster assignments.
#' @return A `V x V` numeric matrix with co-clustering probabilities.
#' @export
#' @examples
#' z_post <- matrix(c(1,1,2, 1,2,2, 1,1,2), nrow = 3)
#' coclustering_matrix(z_post)
coclustering_matrix <- function(z_post) {
  V <- nrow(z_post)
  N_iter <- ncol(z_post)
  coclust <- matrix(0, nrow = V, ncol = V)
  for (t in seq_len(N_iter)) {
    Z <- labels_to_matrix(z_post[, t])
    coclust <- coclust + Z %*% t(Z)
  }
  coclust / N_iter
}

#' Compute the Estimated Edge Probability Matrix
#'
#' Estimates the matrix of edge probabilities given a point estimate of
#' cluster memberships and the adjacency matrix. Uses a Beta-Binomial
#' conjugate update.
#'
#' @param memb Integer vector of cluster memberships (length `V`).
#' @param Y A `V x V` binary adjacency matrix.
#' @param a Shape1 hyperparameter of the Beta prior on edge probabilities.
#' @param b Shape2 hyperparameter of the Beta prior on edge probabilities.
#' @return A `V x V` numeric matrix of estimated edge probabilities (zero diagonal).
#' @export
#' @examples
#' Y <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3, 3)
#' edge_probability(c(1, 1, 2), Y, 1, 1)
edge_probability <- function(memb, Y, a, b) {
  z <- stats::model.matrix(~ factor(memb) - 1)
  V <- nrow(Y)
  Abs_Freq <- t(z) %*% Y %*% z
  diag(Abs_Freq) <- diag(Abs_Freq) / 2
  Tot <- t(z) %*% matrix(1, V, V) %*% z
  diag(Tot) <- (diag(Tot) - table(memb)) / 2
  Rel_Freq <- (a + Abs_Freq) / (a + b + Tot)
  edge_matr <- z %*% Rel_Freq %*% t(z)
  diag(edge_matr) <- 0
  edge_matr
}

#' Compute the Posterior Distribution on K
#'
#' Computes `p(k | t)` for `k = 1, ..., upto_k` and `t = 1, ..., upto_t`
#' under the MFM model with a zero-truncated Poisson prior.
#'
#' @param gamma_ Dirichlet concentration parameter.
#' @param n Number of observations.
#' @param upto_k Maximum number of components.
#' @param upto_t Maximum number of clusters.
#' @return A `upto_k x upto_t` matrix with `p(k | t)` values.
#' @export
#' @importFrom extraDistr dtpois
#' @examples
#' posterior_k(1, 50, 10, 10)
posterior_k <- function(gamma_, n, upto_k, upto_t) {
  if (upto_t > n) stop("p(k|t) is undefined for t > n.")
  log_v <- log_vn_miller(gamma_, n, upto_t, lambda = 1)
  p <- matrix(0, nrow = upto_k, ncol = upto_t)
  for (k in seq_len(upto_k)) {
    for (t in seq_len(min(k, upto_t))) {
      b <- lgamma(k + 1) - lgamma(k - t + 1) - lgamma(k * gamma_ + n) +
        lgamma(k * gamma_) +
        extraDistr::dtpois(k, lambda = 1, a = 0, b = Inf, log = TRUE)
      p[k, t] <- exp(b - log_v[t])
    }
  }
  p
}
