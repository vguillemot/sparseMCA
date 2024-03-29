#' Compute the projection of a vector onto the intersection of the l1 ball,
#' the l2 ball and the space orthogonal to M.
#'
#' @param x A vector of numerics
#' @param a The radius (>0) of the l1 ball
#' @param M A matrix of vectors
#' @param itermax The maximum number of iterations
#' @param eps Precision
#' @return The projection of \eqn{x} onto \eqn{B_1(a) \cap B_2 \cap M^\perp}.
#' @examples
#' projl1l2orth(1:10, a=10, M=projl1l2(rnorm(10))$x)
#' @export
projl1l2orth <- function(x, a=1, M, itermax=5000, eps=1e-16) {
  xold <- xnew <- x
  for (k in 1:itermax) {
    if (is.null(M)) {
      xnew <- projl1l2(xold, a=a)$x
    } else {
      xnew <- projl1l2(projorth(xold, M), a=a)$x
    }
    if ( norm2(xnew - xold) < eps ) break
    xold <- xnew
  }
  list(x=xnew, k=k)
}
