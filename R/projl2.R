#' Normalize a vector.
#'
#' @param u A vector of numerics
#' @return The normalized version of \eqn{u}.
#' @examples
#' projl2(c(0, 0,-0.5, 0.5))
#' @export
projl2 <- function(u) {
  normu <- norm2(u)
  if (normu == 0) {
    return(u)
  }
  return(u / normu)
}
