#' Proximal operator of the group l1l2 norm
#'
#' @param x vector of numerics
#' @param g vector describing the non-overlapping groups
#' @param lambda regularization factor
#'
#' @return The transformation of the input vector by the proximal operator of the group l1l2 norm
#' @export proxgroup
#'
#' @examples
#' x <- rnorm(50)
#' g <- gl(10, 5)
#' y <- proxgroup(x, g, 2)
#' layout(t(1:2))
#' plot(x, col = g)
#' plot(y, col = g)
#' 
proxgroup <- function(x, g, lambda) {
  ave(x, g, FUN = function(xsub) proxl2(xsub, lambda))
}
