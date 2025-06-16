kummer <- function(a, b, z, eps = 1e-06) {
  #' Confluent \eqn{D}-Hypergeometric Function
  #'
  #' Computes the Kummer's function, or confluent hypergeometric function.
  #'
  #' @aliases kummer
  #'
  #' @usage kummer(a, b, z, eps = 1e-06)
  #' @param a numeric.
  #' @param b numeric
  #' @param z numeric vector.
  #' @param eps numeric. Precision for the sum (default 1e-06).
  #' @return A numeric value: the value of the Kummer's function,
  #' with two attributes \code{attr(, "epsilon")} (precision of the result) and \code{attr(, "k")} (number of iterations).
  #'
  #' @details The Kummer's confluent hypergeometric function is given by:
  #' \deqn{\displaystyle{_1 F_1\left(a, b; z\right) = \sum_{n = 0}^{+\infty}{ \frac{ (a)_n }{ (b)_n } \frac{z^n}{n!} }}}
  #'
  #' where \eqn{(z)_p} is the Pochhammer symbol (see \code{\link{pochhammer}}).
  #'
  #' The \code{eps} argument gives the required precision for its computation.
  #' It is the \code{attr(, "epsilon")} attribute of the returned value.
  #'
  #' @author Pierre Santagostini, Angélina El Ghaziri, Nizar Bouhlel
  #' 
  #' @references El Ghaziri, A., Bouhlel, N., Sapoukhina, N., Rousseau, D.,
  #' On the importance of non-Gaussianity in chlorophyll fluorescence imaging.
  #' Remote Sensing 15(2), 528 (2023).
  #' \doi{10.3390/rs15020528}
  #' 
  #' @export
  
  d <- 1
  
  n <- 0
  res <- 0
  
  while ((max(Re(d)) > eps) & (all(is.finite(d)))) {
    res <- res + d
    # cat(d, "\t", res, "\n")
    n <- n + 1
    # d <- ( pochhammer(a, n) / pochhammer(b, n) ) * ( z^n / factorial(n) )
    d <- exp(
      lnpochhammer(a, n) - lnpochhammer(b, n) + n*log(z) - lfactorial(n) #sapply(z, function(x) { n*log(x) })
    )
  }
  # cat("\n")
  if (n == 0) {
    res <- rep(1, length(z))
  }
  
  attr(res, "k") <- n
  attr(res, "epsilon") <- Re(d)
  return(Re(res))
}