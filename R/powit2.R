#' Power iteration engine, with orthogonality constraints
#'
#' @param X
#' @param U0
#' @param V0
#' @param Uorth
#' @param Vorth
#' @param au
#' @param av
#' @param Gu
#' @param Gv
#' @param eps.pi
#' @param eps.pocs
#' @param itermax.pi
#' @param itermax.pocs
#'
#' @return
#' @export
#'
#' @examples
powit2 <- function(X, U0, V0, Uorth, Vorth, au, av, Gu, Gv, eps.pi, eps.pocs, itermax.pi, itermax.pocs) {
  uold <- unew <- U0
  vold <- vnew <- V0
  for (iter in 1:itermax.pi) {
    vnew <- projgrouporth(x=t(X) %*% uold, r=av, g=Gv, M=Vorth,
                       itermax = itermax.pocs, eps = eps.pocs)$x
    unew <- projgrouporth(x=X %*% vnew, r=au, g = Gu, M=Uorth,
                       itermax = itermax.pocs, eps = eps.pocs)$x
    if ( norm2(vnew - vold) < eps.pi && norm2(unew - uold) < eps.pi ) break
    vold <- vnew
    uold <- unew
  }

  return(list(U=unew, V=vnew, iter=iter))
}
