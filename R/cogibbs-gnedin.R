#' Collapsed Gibbs Sampler with Gnedin Prior
#'
#' Runs a collapsed Gibbs sampler for the MFM-SBM model using a
#' Gnedin prior on the number of components.
#'
#' @param M Number of MCMC iterations (positive integer).
#' @param K_start Initial number of clusters (positive integer).
#' @param A A symmetric `n x n` binary adjacency matrix.
#' @param n Number of nodes.
#' @param a Shape1 hyperparameter of the Beta prior on edge probabilities.
#' @param b Shape2 hyperparameter of the Beta prior on edge probabilities.
#' @param gamma_ Dirichlet concentration parameter (positive scalar).
#' @param gamma_gn Gnedin's gamma parameter (in `(0, 1)`).
#' @param z_start Integer vector of length `n` with initial cluster assignments.
#' @return A named list with components:
#'   \describe{
#'     \item{z_post}{An `n x M` matrix of cluster assignments at each iteration.}
#'     \item{num_k}{Integer vector of length `M` with the number of clusters
#'       at each iteration.}
#'     \item{z_start}{The initial cluster assignment vector.}
#'   }
#' @export
#' @importFrom stats rbeta
#' @examples
#' \dontrun{
#' sbm <- generate_sbm(20, 2, 0.5, 0.1, "balanced")
#' z0 <- c(sample(1:3, 3, replace = FALSE), sample(1:3, 17, replace = TRUE))
#' fit <- cogibbs_gnedin(50, 3, sbm$A0, 20, 1, 1, 1, 0.5, z0)
#' }
cogibbs_gnedin <- function(M, K_start, A, n, a, b, gamma_, gamma_gn, z_start) {
  num_k <- rep(NA_integer_, M)
  z_post <- matrix(NA_integer_, nrow = n, ncol = M)
  z_update <- z_start
  Q <- matrix(NA_real_, nrow = K_start, ncol = K_start)

  for (iter in seq_len(M)) {
    k <- length(unique(z_update))
    num_k[iter] <- k

    for (r in seq_len(k)) {
      for (s in r:k) {
        a_sum <- sum(A[which(z_update == r), which(z_update == s)])
        a_rs <- if (r == s) a_sum / 2 else a_sum
        n_r <- sum(z_update == r)
        n_s <- sum(z_update == s)
        n_rs <- if (r == s) n_r * (n_r - 1) / 2 else n_r * n_s
        Q[r, s] <- stats::rbeta(1, a_rs + a, n_rs - a_rs + b)
      }
    }
    Q <- .symmetrize_matrix(Q)

    for (i in seq_len(n)) {
      pr_tbl <- rep(NA_real_, k)
      for (tbl in seq_len(k)) {
        indx <- (seq_len(n))[-i]
        l_prod <- 0
        for (ind in indx) {
          par_sum <- log(Q[tbl, z_update[ind]]) * A[i, ind] +
            log(1 - Q[tbl, z_update[ind]]) * (1 - A[i, ind])
          l_prod <- l_prod + par_sum
        }
        pr_tbl[tbl] <- (sum(z_update == tbl) + gamma_) * exp(l_prod)
      }

      c_i <- length(unique(z_update[-i]))

      mAi <- 1
      for (tt in seq_len(c_i)) {
        j_ct <- which(z_update[-i] == tt)
        somma <- sum(A[j_ct, i])
        mAi <- mAi * beta(a, b)^(-1) * beta(somma + a, length(j_ct) - somma + b)
      }

      pr_tbl <- c(pr_tbl,
                   exp(log_vn_gnedin(n, c_i + 1, gamma_gn) -
                         log_vn_gnedin(n, c_i, gamma_gn)) * gamma_ * mAi)

      k_old <- k
      zi_old <- z_update[i]
      z_update[i] <- sample(seq_len(k + 1), 1, replace = FALSE, prob = pr_tbl)
      k <- length(unique(z_update))

      if (k_old > k) {
        Q <- as.matrix(Q[-zi_old, -zi_old])
        z_update[z_update > zi_old] <- z_update[z_update > zi_old] - 1L
      }

      if (k_old < k) {
        new_prb <- rep(NA_real_, k)
        for (group in seq_len(k - 1)) {
          a_sum <- sum(A[which(z_update == group), which(z_update == k)])
          n_group <- sum(z_update == group)
          n_k <- sum(z_update == k)
          n_gk <- n_group * n_k
          new_prb[group] <- stats::rbeta(1, a_sum + a, n_gk - a_sum + b)
        }
        new_prb[k] <- stats::rbeta(1, a, b)
        Q <- as.matrix(rbind(
          cbind(Q, new_prb[-k]),
          new_prb
        ))
      }

      if (max(z_update) > k) {
        new_prb <- rep(NA_real_, max(z_update))
        for (group in (seq_len(max(z_update)))[-zi_old]) {
          a_sum <- sum(A[which(z_update == group), which(z_update == max(z_update))])
          n_r <- sum(z_update == group)
          n_s <- sum(z_update == max(z_update))
          n_rs <- n_r * n_s
          new_prb[group] <- stats::rbeta(1, a_sum + a, n_rs - a_sum + b)
        }
        new_prb[max(z_update)] <- stats::rbeta(1, a, b)
        Q_par <- cbind(Q, new_prb[-max(z_update)])
        Q <- as.matrix(rbind(Q_par, new_prb))
        Q <- as.matrix(Q[-zi_old, -zi_old])
        z_update[z_update > zi_old] <- z_update[z_update > zi_old] - 1L
      }
    }

    z_post[, iter] <- z_update
  }

  list(z_post = z_post, num_k = num_k, z_start = z_start)
}
