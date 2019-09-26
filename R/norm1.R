#' Compute the weighted l1 norm of a vector.
#'
#' @param u A vector of numerics
#' @param w A (positive) scalar or a vector of (positive) weights that has the same length as u. By default equal to 1.
#' @return The weighted \eqn{l_1} norm of \eqn{u}: \eqn{\sum_i w_i |u_i|}.
#' @examples
#' norm1(c(0, 0,-0.5, 0.5))
#' @export
norm1 <- function(u, w = 1) sum(w * abs(u))
