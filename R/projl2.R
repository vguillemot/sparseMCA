#' Normalize a (weighted) vector.
#'
#' @param u A vector of numerics
#' @param w A (positive) scalar or a vector of (positive) weights that has the same length as u. By default equal to 1.
#' @return The normalized version of \eqn{u}.
#' @examples
#' normalize(c(0, 0,-0.5, 0.5))
#' @export
projl2 <- function(u, w = 1) {
  normu <- norm2(u, w)
  if (normu == 0) {
    return(u)
  }
  return(u / normu)
}
