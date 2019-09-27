#' Compute the projection of a vector onto the l1 ball.
#'
#' @param v A vector of numerics
#' @param a The radius (>0) of the l1 ball
#' @param g A vector of groups
#' @return The projection of \eqn{v} onto the l1 ball of radius \eqn{a}.
#' @examples
#' projgroup(1:10, rep(1:2, e=5), 3)
#' @export
projgroup <- function(v, a, g) {
  y <- tapply(x, g, norm2)
  if (r == 1) {
    gmax <- unique(g)[which.max(y)]
    res <- x
    res[g != gmax] <- 0
    return(normalize(res))
  }
  res.proj <- projl1l2(y, r)
  return(projl2(proxg(x, g, res.proj$lambda)))
}

