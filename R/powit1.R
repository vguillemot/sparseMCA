#' Power iteration engine, without orthogonality constraints
#'
#' @param X
#' @param U0
#' @param V0
#' @param au
#' @param av
#' @param Gu
#' @param Gv
#' @param eps.pi
#' @param itermax.pi
#'
#' @return
#' @export
#'
#' @examples
powit1 <- function(X, U0, V0, au, av, Gu, Gv, eps.pi, itermax.pi) {
  uold <- unew <- U0
  vold <- vnew <- V0
  for (iter in 1:itermax.pi) {
    vnew <- projgroup(t(X) %*% uold, a = av, Gv)$x
    unew <- projgroup(X %*% vnew, a = au, Gu)$x
    if (norm2(vnew - vold) < eps.pi &&
        norm2(unew - uold) < eps.pi)
      break
    vold <- vnew
    uold <- unew
  }

  return(list(U = unew, V = vnew, iter = iter))
}
