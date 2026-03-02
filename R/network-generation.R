#' Generate a Stochastic Block Model Network
#'
#' Generates a random undirected network from a Stochastic Block Model (SBM)
#' with `K` communities and `n` nodes.
#'
#' @param n Number of nodes (positive integer).
#' @param K Number of communities (positive integer, `K <= n`).
#' @param p Within-community edge probability.
#' @param q Between-community edge probability.
#' @param type_network Character string: `"balanced"` for equal-sized communities,
#'   `"unbalanced"` for progressively unequal sizes.
#' @return A named list with components:
#'   \describe{
#'     \item{z0}{Integer vector of length `n` with the true community labels.}
#'     \item{Q}{The `K x K` block probability matrix.}
#'     \item{A0}{The `n x n` symmetric binary adjacency matrix (zero diagonal).}
#'   }
#' @export
#' @examples
#' set.seed(42)
#' sbm <- generate_sbm(n = 50, K = 2, p = 0.5, q = 0.1, type_network = "balanced")
#' dim(sbm$A0)
#' table(sbm$z0)
generate_sbm <- function(n, K, p, q, type_network = c("balanced", "unbalanced")) {
  type_network <- match.arg(tolower(type_network),
                            choices = c("balanced", "unbalanced"))

  if (type_network == "balanced") {
    community_size <- rep(trunc(n / K), K)
    if (trunc(n / K) * K != n) {
      community_size[K] <- community_size[K] + (n - sum(community_size))
    }
    z0 <- rep(NA_integer_, n)
    for (group in seq_len(K)) {
      for (i in seq_len(community_size[group])) {
        z0[trunc(n / K) * (group - 1L) + i] <- group
      }
    }
  } else {
    community_size <- rep(trunc(n / K), K)
    for (i in seq_len(length(community_size) - 1L)) {
      half_i <- floor(community_size[i] / 2)
      community_size[i] <- half_i
      community_size[i + 1L] <- community_size[i + 1L] + half_i
    }
    community_size[K] <- community_size[K] + (n - sum(community_size))
    z0 <- rep(NA_integer_, n)
    for (group in seq_len(K)) {
      start_i <- if (group == 1L) 1L else sum(community_size[1:(group - 1L)]) + 1L
      for (i in start_i:(community_size[group] + start_i - 1L)) {
        z0[i] <- group
      }
    }
  }

  Q <- matrix(q, nrow = K, ncol = K)
  diag(Q) <- p

  A0 <- matrix(0L, nrow = n, ncol = n)
  for (row in seq_len(n)) {
    for (col in row:n) {
      if (row != col) {
        A0[row, col] <- stats::rbinom(1, 1, Q[z0[row], z0[col]])
      }
    }
  }
  A0 <- .symmetrize_matrix(A0)

  list(z0 = z0, Q = Q, A0 = A0)
}
