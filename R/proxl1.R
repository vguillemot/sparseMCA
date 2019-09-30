#' Proximal operator of the l1 norm
#'
#' @param x a vetor of numeric values
#' @param lambda 
#'
#' @return
#' @export
#' x <- c(-1, 0.2, 0.2, 1)
#' proxl1(x, 0) # Does not change the input vector
#' proxl1(x, 1) 
#' proxl2(x, 2) # Outputs a null vector
#' @examples
proxl1 <- function(x, lambda) {
  return(sign(x) * pmax(0, abs(x) - lambda))
}
