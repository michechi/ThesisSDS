#' Dolphin Social Network
#'
#' An undirected social network of frequent associations between 62 bottlenose
#' dolphins in a community living off Doubtful Sound, New Zealand.
#' Compiled by Lusseau et al. (2003).
#'
#' @format A list with two elements:
#' \describe{
#'   \item{A}{A 62 x 62 binary symmetric adjacency matrix.}
#'   \item{Z0}{An integer vector of length 62 with the known community
#'     memberships (two communities).}
#' }
#' @source Lusseau D., Schneider K., Boisseau O.J., Haase P., Slooten E.,
#'   Dawson S.M. (2003). The bottlenose dolphin community of Doubtful Sound
#'   features a large proportion of long-lasting associations.
#'   \emph{Behavioral Ecology and Sociobiology}, 54, 396--405.
#' @examples
#' data(dolphin)
#' dim(dolphin$A)
#' table(dolphin$Z0)
"dolphin"
