#' Sensitivity Analysis for MFM-SBM
#'
#' Runs a sensitivity analysis over a grid of hyperparameters for both the
#' zero-truncated Poisson and Gnedin priors. For each simulated dataset,
#' the Gibbs sampler is run with each hyperparameter value and performance
#' metrics are recorded.
#'
#' @param n Number of nodes.
#' @param N Number of simulated datasets.
#' @param M Number of MCMC iterations per run.
#' @param K True number of communities.
#' @param K_start Initial number of clusters for the Gibbs sampler.
#' @param p Within-community edge probability.
#' @param q Between-community edge probability.
#' @param a Shape1 hyperparameter of the Beta prior.
#' @param b Shape2 hyperparameter of the Beta prior.
#' @param gamma_ Dirichlet concentration parameter.
#' @param BI Number of burn-in iterations to discard.
#' @param par_po Numeric vector of Poisson lambda values to test
#'   (default `seq(1, 10, 2)`).
#' @param par_gn Numeric vector of Gnedin gamma values to test
#'   (default `seq(0.1, 0.9, 0.2)`).
#' @param type_network Type of network: `"balanced"` or `"unbalanced"`.
#' @return A named list with matrices of results for each metric and prior type:
#'   `TIME.gn`, `TIME.po`, `prob.gn`, `prob.po`, `k_es.gn`, `k_es.po`,
#'   `RaIn.gn`, `RaIn.po`, `peVI.gn`, `peVI.po`.
#' @export
#' @examples
#' \dontrun{
#' sa <- sensitivity_analysis(50, 5, 100, 2, 5, 0.5, 0.1, 1, 1, 1, 50)
#' }
sensitivity_analysis <- function(n, N, M, K, K_start, p, q, a, b, gamma_, BI,
                                  par_po = seq(1, 10, 2),
                                  par_gn = seq(0.1, 0.9, 0.2),
                                  type_network = "unbalanced") {
  n_par <- length(par_po)

  TIME_po <- prob_po <- k_es_po <- RaIn_po <- peVI_po <-
    matrix(NA_real_, nrow = N, ncol = n_par)
  TIME_gn <- prob_gn <- k_es_gn <- RaIn_gn <- peVI_gn <-
    matrix(NA_real_, nrow = N, ncol = n_par)

  for (i in seq_len(N)) {
    message("Dataset ", i, " / ", N)
    data_gen <- generate_sbm(n, K, p, q, type_network = type_network)
    A <- data_gen$A0
    Z0 <- data_gen$z0
    z_start <- c(sample(seq_len(K_start), K_start, replace = FALSE),
                 sample(seq_len(K_start), n - K_start, replace = TRUE))

    for (ii in seq_along(par_po)) {
      time_start <- Sys.time()
      logVn <- log_vn_miller(gamma_, n, n + 10, par_po[ii])
      fit_po <- cogibbs_poisson(M, K_start, A, n, a, b, gamma_, logVn, z_start)
      time_end <- Sys.time()

      TIME_po[i, ii] <- as.numeric(difftime(time_end, time_start, units = "secs"))
      prob_po[i, ii] <- mean(fit_po$num_k[-seq_len(BI)] == K)
      k_es_po[i, ii] <- find_mode(fit_po$num_k[-seq_len(BI)])

      coclust <- coclustering_matrix(fit_po$z_post[, (BI + 1):M])
      if (requireNamespace("mcclust.ext", quietly = TRUE)) {
        memb_vi <- mcclust.ext::minVI(coclust, method = "avg", max.k = 20)$cl
        if (requireNamespace("fossil", quietly = TRUE)) {
          peVI_po[i, ii] <- fossil::compare(memb_vi, Z0, method = "vi")
          RaIn_po[i, ii] <- fossil::compare(memb_vi, Z0, method = "rand")
        }
      }
    }

    for (jj in seq_along(par_gn)) {
      time_start <- Sys.time()
      fit_gn <- cogibbs_gnedin(M, K_start, A, n, a, b, gamma_, par_gn[jj], z_start)
      time_end <- Sys.time()

      TIME_gn[i, jj] <- as.numeric(difftime(time_end, time_start, units = "secs"))
      prob_gn[i, jj] <- mean(fit_gn$num_k[-seq_len(BI)] == K)
      k_es_gn[i, jj] <- find_mode(fit_gn$num_k[-seq_len(BI)])

      coclust <- coclustering_matrix(fit_gn$z_post[, (BI + 1):M])
      if (requireNamespace("mcclust.ext", quietly = TRUE)) {
        memb_vi <- mcclust.ext::minVI(coclust, method = "avg", max.k = 20)$cl
        if (requireNamespace("fossil", quietly = TRUE)) {
          peVI_gn[i, jj] <- fossil::compare(memb_vi, Z0, method = "vi")
          RaIn_gn[i, jj] <- fossil::compare(memb_vi, Z0, method = "rand")
        }
      }
    }
  }

  list(TIME.gn = TIME_gn, TIME.po = TIME_po,
       prob.gn = prob_gn, prob.po = prob_po,
       k_es.gn = k_es_gn, k_es.po = k_es_po,
       RaIn.gn = RaIn_gn, RaIn.po = RaIn_po,
       peVI.gn = peVI_gn, peVI.po = peVI_po)
}
