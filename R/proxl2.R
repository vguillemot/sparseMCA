#' Proximal operator of the l2 norm
#'
#' @param x a vetor of numeric values
#' @param lambda 
#'
#' @return
#' @export
#' x <- c(-1, 0.2, 0.2, 1)
#' proxl2(x, 0) # Does not change the input vector
#' proxl2(x, 1) # Outputs a vector linearly linked to the input vector but with a lower l2 norm
#' proxl2(x, 2) # Outputs a null vector
#' @examples
proxl2 <- function(x, lambda) {
  norm2x <- norm2(x)
  if (norm2x == 0) return(0*x)
  return(x * max(0, (norm2x - lambda) / norm2x))
}
