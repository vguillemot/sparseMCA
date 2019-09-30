#' Compute the projection of a vector onto the intersection of the l1 ball
#' and the l2 ball.
#'
#' @param x A vector of numerics
#' @param a The radius (>0) of the l1 ball
#' @return The "projection" of \eqn{x} onto \eqn{B_1(a) \cap B_2}.
#' @examples
#' projl1l2(1:10, 21)
#' @export
projl1l2 <- function(x, a = 1) {
  # set.seed(1)
  # Check if constraints are already satisfied
  norm2_x <- norm2(x)
  if (a == 1) {
    res <- x
    res[x == max(abs(x))] <- 0
    return(list(x=x/norm2_x, lambda = 1, k=NaN))
  }
  if ( norm2_x < 1e-32 ) return(list(x=x, lambda = NaN))
  if ( sum(abs(x/norm2_x)) <= a ) return(list(x=x/norm2_x, lambda = 0, k=NaN))
  uneq <- x != 0
  L <- sum(!uneq)
  p <- abs(x[uneq])
  # Check if multiple maximum
  MAX <- max(p)
  bMAX <- p == MAX
  nMAX <- sum(bMAX)
  if (a < sqrt(nMAX)){
    print(a)
    print(sqrt(nMAX))
    stop("Impossible to project, minimum ratio is : ", sqrt(nMAX))
  } else if (a == sqrt(nMAX)){
    warning("a == sqrt(nMAX)")
    x_soft        <- rep(0, length(x))
    x_soft[bMAX]  <- 1/sqrt(nMAX)
    return(list(x=x_soft, k=NaN))
  }
  #Initialize parameters
  s_1 <- s_2 <- nb <- 0
  iter <- 0
  while (T) {
    iter <- iter + 1
    N <- length(p)
    if (N==0) {
      warning("length(p) = 0")
      break
    }
    #Choose next a_k
    a_k     <- p[round(runif(1, 1, N), 0)]
    while(a_k == MAX){
      a_k     <- p[round(runif(1, 1, N), 0)]
    }
    # print(a_k)
    #Make a partition of list p
    p_inf_ak <- p < a_k
    p_sup_ak <- p > a_k
    p_high  <- p[p_inf_ak]
    p_low   <- p[p_sup_ak]
    #Evaluation decreasing rank of a_k
    nb_a_k  <- sum(p == a_k)
    k       <- nb + sum(p_sup_ak) + nb_a_k
    #Compute value of the constraint
    aksq <- a_k^2
    s_low_1 <- sum(p_low) + nb_a_k*a_k
    s_low_2 <- ssq(p_low) + nb_a_k*aksq
    psi_a_k <- (s_1 + s_low_1 - k*a_k)/sqrt(s_2 + s_low_2 - 2*a_k*(s_1 + s_low_1) + k*aksq)
    #Choose partition depending on the constraint
    if (psi_a_k > a){
      if ( length(p_low) == 0 ) break
      p         <- p_low
    }else{
      if (length(p_high) == 0){
        break
      }else{
        a_k_1     <- max(p_high)
        psi_a_k_1 <- (s_1 + s_low_1 - k*a_k_1)/sqrt(s_2 + s_low_2 - 2*a_k_1*(s_1 + s_low_1) + k*a_k_1^2)
        if (psi_a_k_1 > a){
          break
        }
        p   <- p_high
        nb  <- k
        s_1 <- s_1 + s_low_1
        s_2 <- s_2 + s_low_2
      }
    }
  }
  #Compute lambda and thus the soft tresholded vector to return
  lambda <- a_k - (a*sqrt((k - psi_a_k^2)/(k - a^2))-psi_a_k)*(s_1 + s_low_1 - k*a_k)/(psi_a_k*(k))
  x_soft <- sign(x)*pmax(0, abs(x) - lambda)
  return( list(x=x_soft / norm2(x_soft) , lambda = lambda, k=iter) )
}
