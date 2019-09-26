#' Compute the sum of squares.
#'
#' @param x A vector of numeric values
#' @return The sum of the squared coefficients of vector x.
#' @examples 
#' ssq(1:10)
#' @export
ssq <- function(x) sum(x**2)
