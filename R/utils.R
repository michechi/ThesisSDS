#' Find the Mode of a Vector
#'
#' Computes the mode (most frequent value) of a numeric or integer vector.
#' When there are multiple modes, one is selected at random.
#'
#' @param vect A numeric or integer vector.
#' @return The modal value. If multiple modes exist, one is chosen at random.
#' @export
#' @examples
#' find_mode(c(1, 2, 2, 3, 3))
#' find_mode(c(1, 1, 1, 2, 3))
find_mode <- function(vect) {
  tab <- table(vect)
  who_max <- as.numeric(names(which(tab == max(tab))))
  n_max <- length(who_max)
  if (n_max != 1L) {
    who_max[sample.int(n_max, 1L)]
  } else {
    who_max
  }
}

#' Symmetrize a Square Matrix
#'
#' Copies the upper triangle of a square matrix to the lower triangle,
#' making the matrix symmetric. Replaces the dependency on
#' `gdata::lowerTriangle()` / `gdata::upperTriangle()`.
#'
#' @param mat A square numeric matrix.
#' @return The symmetrized matrix (lower triangle mirrors the upper triangle).
#' @keywords internal
.symmetrize_matrix <- function(mat) {
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  mat
}
