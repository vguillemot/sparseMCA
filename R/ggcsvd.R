#' Groups sparse Generalized SVD of a matrix
#'
#' @param X a (data) matrix;
#' @param R the desired rank of the singular decomposition;
#' @param au The radiuses (>0) of the
#' $L_1$ ball for each left vector.
#' Default to the maximum possible radius, such that the result
#' is the same as the result of a regular SVD.
#' @param av The radiuses
#' (>0) of the $L_1$ balls for each right vector.
#' Default to the maximum possible radius, such that the result
#' is the same as the result of a regular SVD.
#' @param Gu a vector describing the groups for the lines.
#' @param Gv a vector describing the groups for the columns.
#' @param itermax.pi The maximum number of iterations for the power iteration.
#' @param itermax.pocs The maximum number of iterations for POCS.
#' @param eps.pi Precision for the power iteration.
#' @param eps.pocs Precision for POCS.
#' @param init How to initialize the algorithm.
#' Either "svd" (default) or "rand" to intialize with,
#' respectively, the results of a regular SVD or random
#' vectors.
#' @param order_sv Boolean. Should the singular values
#' be artificially ordered. Default to TRUE.
#' @return Pseudo-singular vectors (U and V) and values (D), and the number of iterations.
#' @examples
#' X <- matrix(rnorm(20), 5, 4)
#' ggcsvd(X) # By default, no sparsity constraint is imposed on the decomposition
#' svd(X) # Should give roughly the same results
#' @author Vincent Guillemot
#' @export
ggcsvd <- function(X,
                   R = 2,
                   au = rep(sqrt(nrow(X)), R),
                   av = rep(sqrt(ncol(X)), R),
                   Gu = seq(NROW(X)),
                   Gv = seq(NCOL(X)),
                   itermax.pi = 1000,
                   itermax.pocs = 1000,
                   eps.pi = 1e-16,
                   eps.pocs = 1e-16,
                   init = "svd",
                   order_sv = TRUE) {
  I <- nrow(X)
  J <- ncol(X)

  if (I == 1 & J == 1)
    stop("Are you sure you want to perform the SVD of a scalar?")
  if (R <= 1)
    stop("R should be > 1")

  nas <- is.na(X)
  if (sum(nas) > 0) {
    X[nas] <- mean(X[!nas])
    warning(sprintf("Number of NAs:%i", sum(nas)))
  }

  # Build initialization matrices either with SVD
  # (prefered method) or randomly
  if (init=="svd") {
    svdx <- svd(X, nu=R, nv=R)
    U0 <- svdx$u
    V0 <- svdx$v
  } else if( init=="rand") {
    U0 <- 1/(I-1) * mvrnorm(n = I, mu = rep(0,R),
                            Sigma = diag(R), empirical = TRUE)
    V0 <- 1/(J-1) * mvrnorm(n = J, mu = rep(0,R),
                            Sigma = diag(R), empirical = TRUE)
  } else {
    stop("init should be either svd or rand.")
  }
  U <- matrix(0, I, R)
  V <- matrix(0, J, R)

  iter <- rep(NA, R)

  ## Power iteration without orth projection
  res.powit1 <- powit1(X,
                       U0[,1,drop=FALSE], V0[,1,drop=FALSE],
                       au[1], av[1],
                       Gu, Gv,
                       eps.pi,
                       itermax.pi)
  U[,1] <- res.powit1$U
  V[,1] <- res.powit1$V
  iter[1] <- res.powit1$iter

  if (R > 1) {
    for (r in 2:R) {
      ## Power Iteration with orth projection
      res.powit2 <- powit2(X,                        # original matrix
                           U0[,r], V0[,r],           # initialization vectors
                           U[,1:(r-1),drop=FALSE], V[,1:(r-1),drop=FALSE], # orth constraint
                           au[r], av[r],             # sparsity constrain
                           Gu, Gv,                   # groups
                           eps.pi, eps.pocs,         # precision
                           itermax.pi, itermax.pocs) # max iteration

      U[,r] <- res.powit2$U
      V[,r] <- res.powit2$V
      iter[r] <- res.powit2$iter

    }
  }

  D <- diag(t(U) %*% X %*% V)
  if (order_sv) {
    oD <- order(D, decreasing = TRUE)
  } else {
    oD <- 1:R
  }
  res <- list(U=U[,oD], V=V[,oD], D=D[oD], iter=iter)
  return(res)
}
