#' Backward-Compatible Aliases
#'
#' These aliases map the original function names used in the thesis scripts
#' to the new package function names. They are provided so that existing
#' analysis code continues to work without modification.
#'
#' @inheritParams cogibbs_poisson
#' @name compat-aliases
#' @keywords internal
NULL

#' @rdname compat-aliases
#' @export
fcogibbs_po <- cogibbs_poisson

#' @rdname compat-aliases
#' @export
fcogibbs_gn <- cogibbs_gnedin

#' @rdname compat-aliases
#' @export
fgenera_data <- generate_sbm

#' @rdname compat-aliases
#' @export
log_Vn.M <- log_vn_miller

#' @rdname compat-aliases
#' @export
l.Vn.Gn <- log_vn_gnedin

#' @rdname compat-aliases
#' @export
flv_n.app <- log_vn_approx

#' @rdname compat-aliases
#' @export
fv_nt <- vn_summation

#' @rdname compat-aliases
#' @export
fk_t <- falling_factorial

#' @rdname compat-aliases
#' @export
fgamk_i <- rising_factorial

#' @rdname compat-aliases
#' @export
vec2mat <- labels_to_matrix

#' @rdname compat-aliases
#' @export
pr_cc <- coclustering_matrix

#' @rdname compat-aliases
#' @export
edge_est <- edge_probability

#' @rdname compat-aliases
#' @export
post_k <- posterior_k

#' @rdname compat-aliases
#' @export
sens_anal <- sensitivity_analysis
