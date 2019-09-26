#' Compute the weighted l2 norm of a vector.
#'
#' @param u A vector of numerics.
#' @param w A (positive) scalar or a vector of (positive) weights that has the same length as u. By default equal to 1.
#' @return The weighted l2 norm of \eqn{u}: \eqn{\sqrt{\sum_i w_i u_i^2}}.
#' @examples
#' norm2(c(0, 0,-0.5, 0.5))
#' @export
norm2 <- function(u, w = 1) sqrt(sum(w * u**2))
