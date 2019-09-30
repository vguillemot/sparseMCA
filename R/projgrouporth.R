#' Compute the projection of a vector onto the intersection of the
#' group l1-l2 ball, the l2 ball and the space orthogonal to M.
#'
#' @param x A vector of numerics
#' @param r The radius (>0) of the l1 ball
#' @param g Vector of groups.
#' @param M A matrix of vectors
#' @param itermax The maximum number of iterations
#' @param eps Precision
#' @return The projection of \eqn{x} onto \eqn{B_1(a) \cap B_2 \cap M^\perp}.
#' @examples
#' projgrouporth(x=1:10, r=10, g = 1:10, M=projl2(rnorm(10)))
#' @export
projgrouporth <- function(x=x, r=1, g=g, M=M, itermax=5000, eps=1e-16) {
  xold <- xnew <- x
  for (k in 1:itermax) {
    if (is.null(M)) {
      xnew <- projgroup(x=xold, r=r, g=g)
    } else {
      xnew <- projgroup(x=projorth(x=xold, M), r=r, g=g)
    }
    if ( norm2(xnew - xold) < eps ) break
    xold <- xnew
  }
  list(x=xnew, k=k)
}
