#' Convergence Test for MFM-SBM
#'
#' Evaluates the convergence behaviour of the collapsed Gibbs sampler by
#' computing the average Rand Index across `N` replicated datasets at
#' increasing sample sizes.
#'
#' @param K True number of communities.
#' @param a Shape1 hyperparameter of the Beta prior.
#' @param b Shape2 hyperparameter of the Beta prior.
#' @param gamma_ Dirichlet concentration parameter.
#' @param lambda Rate parameter for the zero-truncated Poisson prior.
#' @param p Within-community edge probability.
#' @param q Between-community edge probability.
#' @param n_vec Integer vector of sample sizes to test.
#' @param N Number of replicated datasets at each sample size.
#' @param K_start Initial number of clusters.
#' @param M Number of MCMC iterations.
#' @param burn_in Number of burn-in iterations.
#' @param type_network Type of network: `"balanced"` or `"unbalanced"`.
#' @return A numeric vector of average Rand Index values, one per element of `n_vec`.
#' @export
#' @examples
#' \dontrun{
#' ri <- convergence_test(
#'   K = 2, a = 1, b = 1, gamma_ = 1, lambda = 1,
#'   p = 0.5, q = 0.1, n_vec = c(20, 40, 60),
#'   N = 5, K_start = 5, M = 100, burn_in = 50
#' )
#' }
convergence_test <- function(K, a, b, gamma_, lambda, p, q, n_vec, N, K_start,
                              M, burn_in, type_network = "unbalanced") {
  ri_avg <- rep(NA_real_, length(n_vec))

  for (num in seq_along(n_vec)) {
    n <- n_vec[num]
    message("n = ", n)
    ri_i <- rep(NA_real_, N)

    for (i in seq_len(N)) {
      data_gen <- generate_sbm(n, K, p, q, type_network = type_network)
      A <- data_gen$A0
      Z0 <- data_gen$z0
      z_start <- c(sample(seq_len(K_start), K_start, replace = FALSE),
                   sample(seq_len(K_start), n - K_start, replace = TRUE))
      logVn <- log_vn_miller(gamma_, n, n + 10, lambda)
      fit_po <- cogibbs_poisson(M, K_start, A, n, a, b, gamma_, logVn, z_start)
      coclust <- coclustering_matrix(fit_po$z_post[, (burn_in + 1):M])

      if (requireNamespace("mcclust.ext", quietly = TRUE) &&
          requireNamespace("fossil", quietly = TRUE)) {
        memb_vi <- mcclust.ext::minVI(coclust, method = "avg", max.k = 20)$cl
        ri_i[i] <- fossil::rand.index(memb_vi, Z0)
      }
    }

    ri_avg[num] <- mean(ri_i, na.rm = TRUE)
  }

  ri_avg
}
