#' Log-Sum-Exp Trick
#'
#' Computes `log(exp(a) + exp(b))` in a numerically stable way.
#'
#' @param a A numeric scalar.
#' @param b A numeric scalar.
#' @return The log-sum-exp of `a` and `b`.
#' @export
#' @examples
#' logsumexp(100, 101)
#' logsumexp(-Inf, 0)
logsumexp <- function(a, b) {
  m <- max(a, b)
  if (m == -Inf) {
    -Inf
  } else {
    log(exp(a - m) + exp(b - m)) + m
  }
}

#' Compute log V_n(t) Using Miller's Method
#'
#' Computes `log(V_n(t))` for `t = 1, ..., upto` under a zero-truncated
#' Poisson prior on the number of components, using the method described
#' in Miller and Harrison (2018).
#'
#' @param gamma_ Dirichlet concentration parameter (positive scalar).
#' @param n Number of observations (positive integer).
#' @param upto Maximum value of `t` for which to compute log V_n(t).
#' @param lambda Rate parameter for the zero-truncated Poisson prior
#'   (positive scalar, default 1).
#' @return A numeric vector of length `upto` with `log(V_n(t))` values.
#' @export
#' @importFrom extraDistr dtpois
#' @examples
#' log_vn_miller(1, 50, 60, lambda = 1)
log_vn_miller <- function(gamma_, n, upto, lambda = 1) {
  tolerance <- 1e-100
  log_v <- rep(0, upto)
  for (t in seq_len(upto)) {
    if (t > n) {
      log_v[t] <- -Inf
    } else {
      a <- 0
      cc <- -Inf
      k <- 1
      p <- 0
      dif <- a - cc
      while ((abs(dif) > tolerance) || (p < (1 - tolerance))) {
        if (k >= t) {
          a <- cc
          b <- lgamma(k + 1) - lgamma(k - t + 1) - lgamma(k * gamma_ + n) +
            lgamma(k * gamma_) + extraDistr::dtpois(k, lambda, 0, log = TRUE)
          cc <- logsumexp(a, b)
        }
        p <- p + exp(log(extraDistr::dtpois(k, 1, 0)))
        k <- k + 1
        dif <- if (a == cc) 0 else a - cc
      }
      log_v[t] <- cc
    }
  }
  log_v
}

#' Compute log V_n(t) Using the Gnedin Prior
#'
#' Computes `log(V_n(t))` under the Gnedin prior on the number of components.
#'
#' @param n Number of observations (positive integer).
#' @param t Number of clusters (positive integer, `t <= n`).
#' @param gamma_gn Gnedin's gamma parameter (in `(0, 1)`).
#' @return A numeric scalar: `log(V_n(t))` under the Gnedin prior.
#' @export
#' @examples
#' log_vn_gnedin(50, 3, 0.5)
log_vn_gnedin <- function(n, t, gamma_gn) {
  lfactorial(t - 1) + lgamma(2 - gamma_gn) - lgamma(3 - gamma_gn - t) +
    lgamma(gamma_gn + 1) - lgamma(gamma_gn - n + t + 1) -
    (lfactorial(n - 1) + lgamma(gamma_gn + 2) - lgamma(gamma_gn - n + 3))
}

#' Compute log V_n(t) Using the Large-n Approximation
#'
#' Approximation of `log(V_n(t))` valid for large `n` (Theorem 5.1).
#'
#' @param n Number of observations (positive integer, ideally `n > 150`).
#' @param t Number of clusters (positive integer).
#' @param gamma_ Dirichlet concentration parameter.
#' @return A numeric scalar: the approximated `log(V_n(t))`.
#' @export
#' @importFrom extraDistr dtpois
#' @examples
#' log_vn_approx(1000, 3, 1)
log_vn_approx <- function(n, t, gamma_) {
  gam_t <- gamma_ * t
  lfactorial(t) - lfactorial(n) + lgamma(gam_t) - log(n^(gam_t - 1)) +
    log(extraDistr::dtpois(t, 1, 0))
}

#' Compute V_n(t) via Direct Summation
#'
#' Computes `V_n(t)` by direct summation (not on the log scale). Uses
#' [falling_factorial()] and [rising_factorial()] internally. Suitable only
#' for small `n` due to numerical overflow.
#'
#' @param n Number of observations (positive integer).
#' @param t Number of clusters (positive integer).
#' @param gamma_ Dirichlet concentration parameter.
#' @return A numeric scalar: `V_n(t)`.
#' @export
#' @importFrom extraDistr dtpois
#' @examples
#' vn_summation(10, 2, 1)
vn_summation <- function(n, t, gamma_) {
  k_start <- if (t != 0) t else 1
  tol <- 1e-100
  ans <- 0
  k <- k_start
  repeat {
    num <- falling_factorial(k, t)
    den <- rising_factorial(gamma_ * k, n)
    poi <- extraDistr::dtpois(k, 1, 0)
    next_term <- (num / den) * poi
    ans <- ans + next_term
    if (abs(next_term) < tol) break
    k <- k + 1
  }
  ans
}

#' Falling Factorial
#'
#' Computes the falling factorial `k^{(t)} = k * (k-1) * ... * (k-t+1)`.
#'
#' @param k A non-negative integer.
#' @param t A non-negative integer (`t <= k`).
#' @return The falling factorial as a numeric scalar.
#' @export
#' @examples
#' falling_factorial(5, 3)  # 5 * 4 * 3 = 60
falling_factorial <- function(k, t) {
  if (t == 0 || t == 1) {
    return(max(t * k, 1))
  }
  result <- k
  for (ind in seq_len(t - 1)) {
    result <- as.numeric(result * (k - ind))
  }
  result
}

#' Rising Factorial (Pochhammer Symbol)
#'
#' Computes the rising factorial `(x)_n = x * (x+1) * ... * (x+n-1)`.
#'
#' @param x A positive numeric scalar.
#' @param n A non-negative integer.
#' @return The rising factorial as a numeric scalar.
#' @export
#' @examples
#' rising_factorial(2, 3)  # 2 * 3 * 4 = 24
rising_factorial <- function(x, n) {
  if (n == 0) return(1)
  result <- x
  if (n == 1) return(result)
  for (j in seq_len(n - 1)) {
    result <- result * (x + j)
  }
  result
}
