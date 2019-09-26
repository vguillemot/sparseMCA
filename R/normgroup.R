#' Compute the l1,2 group norm of a vector.
#'
#' @param u A vector of numerics
#' @param g A group factor
#' @return The \eqn{l_{1,2}} group norm of \eqn{u}: \eqn{\sum_g ||u_g||_2}.
#' @examples
#' normgroup(c(0, 0,-0.5, 0.5), gl(2, 2))
#' @export
normgroup <- function(u, g, w = 1) {
  # sum(tapply(X = u, INDEX = g, FUN = norm2))
  if (length(w) == 1) {
    return(sum(tapply(X = u, INDEX = g, FUN = norm2, w = w)))
  } else if (length(w) == length(u)) {
    res <- 0
    for (gk in unique(g)) {
      res <- res + norm2(u[g == gk], w[g == gk])
    }
    return(res)
  } else {
    stop("Mis-specified inputs.")
  }
}
