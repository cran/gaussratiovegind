pnormratio <- function (z, bet, rho, delta) {
  #' Cumulative Distribution of a Normal Ratio Distribution
  #'
  #' Cumulative distribution of the ratio of two independent Gaussian distributions<!-- with strictly positive means and variances -->.
  #'
  #' @aliases pnormratio
  #'
  #' @usage pnormratio(z, bet, rho, delta)
  #' @param z length \eqn{p} vector of quantiles.
  #' @param bet,rho,delta numeric values. The parameters \eqn{(\beta, \rho, \delta_y)} of the distribution, see Details.
  #'
  #' @details Let two independant random variables
  #' \eqn{X \sim N(\mu_x, \sigma_x)} and \eqn{Y \sim N(\mu_y, \sigma_y)}.
  #' <!-- with \eqn{\mu_x > 0} and \eqn{\mu_y > 0}. -->
  #' 
  #' If we denote \eqn{\displaystyle{ f_Z(z; \beta, \rho, \delta_y)}}
  #' the probability distribution function of the ratio
  #' \eqn{\displaystyle{Z = \frac{X}{Y}}},
  #' with \eqn{\beta = \frac{\mu_x}{\mu_y}},
  #' \eqn{\displaystyle{\rho = \frac{\sigma_y}{\sigma_x}}}
  #' and \eqn{\displaystyle{\delta_y = \frac{\sigma_y}{\mu_y}}}
  #' (see [dnormratio()], Details section).
  #' 
  #' The probability distribution for \eqn{Z} is given by:
  #' \deqn{\displaystyle{F(z; \beta, \rho, \delta_y) = \int_{-\infty}^z{f_Z(z; \beta, \rho, \delta_y)}}}
  #' 
  #' This integral is computed using numerical integration.
  #' 
  #' @return Numeric: the value of density.
  #' 
  #' @seealso [dnormratio()]: density function.
  #' 
  #' [rnormratio()]: sample simulation.
  #' 
  #' [estparnormratio()]: parameter estimation.
  #' 
  #' @author Pierre Santagostini, Angélina El Ghaziri, Nizar Bouhlel
  #'
  #' @references El Ghaziri, A., Bouhlel, N., Sapoukhina, N., Rousseau, D.,
  #' On the importance of non-Gaussianity in chlorophyll fluorescence imaging.
  #' Remote Sensing 15(2), 528 (2023).
  #' \doi{10.3390/rs15020528}
  #' 
  #' Marsaglia, G. 2006. Ratios of Normal Variables.
  #' Journal of Statistical Software 16.
  #' \doi{10.18637/jss.v016.i04}
  #'
  #' Díaz-Francés, E., Rubio, F.J.,
  #' On the existence of a normal approximation to the distribution of the ratio of two independent normal random variables.
  #' Stat Papers 54, 309–323 (2013).
  #' \doi{10.1007/s00362-012-0429-2}
  #'
  #' @examples
  #' # First example
  #' beta1 <- 0.15
  #' rho1 <- 5.75
  #' delta1 <- 0.22
  #' pnormratio(0, bet = beta1, rho = rho1, delta = delta1)
  #' pnormratio(0.5, bet = beta1, rho = rho1, delta = delta1)
  #' curve(pnormratio(x, bet = beta1, rho = rho1, delta = delta1), from = -0.1, to = 0.7)
  #'
  #' # Second example
  #' beta2 <- 2
  #' rho2 <- 2
  #' delta2 <- 2
  #' pnormratio(0, bet = beta2, rho = rho2, delta = delta2)
  #' pnormratio(0.5, bet = beta2, rho = rho2, delta = delta2)
  #' curve(pnormratio(x, bet = beta2, rho = rho2, delta = delta2), from = -0.1, to = 0.7)
  #'
  #' @importFrom stats integrate
  #' @export

  f <- function(z) {
    dnormratio(z, bet = bet, rho = rho, delta = delta)
  }
  res <- sapply(z, function(z) integrate(f, lower = -Inf, upper = z)$value)
  return(res)
}
