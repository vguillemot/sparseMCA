#' Compute the projection of a vector onto the l1 ball.
#'
#' @param x A vector of numerics
#' @param r The radius (>0) of the l1 ball
#' @param g A vector of groups
#' @return The projection of \eqn{v} onto the l1 ball of radius \eqn{a}.
#' @examples
#' projgroup(x=1:10, r=3, g=1:10)
#' @export
projgroup <- function(x=x, r=r, g=g) {
  y <- tapply(x, g, norm2)
  if (r == 1) {
    gmax <- unique(g)[which.max(y)]
    res <- x
    res[g != gmax] <- 0
    return(projl2(res))
  }
  res.proj <- projl1l2(y, r)
  return(projl2(proxgroup(x, g, res.proj$lambda)))
}

