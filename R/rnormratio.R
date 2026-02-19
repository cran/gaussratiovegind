rnormratio <- function (n, bet, rho, delta, r = 0) {
  #' Simulate from a Normal Ratio Distribution
  #'
  #' Simulate data from a ratio of two <!-- independent --> Gaussian distributions<!-- with strictly positive means and variances -->.
  #'
  #' @aliases rnormratio
  #'
  #' @usage rnormratio(n, bet, rho, delta, r = 0)
  #' @param n integer. Number of observations. If \code{length(n) > 1}, the length is taken to be the nmber required.
  #' @param bet,rho,delta numeric values. The parameters \eqn{(\beta, \rho, \delta_y)} of the distribution, see Details.
  #' @param r numeric. The correlation coefficient. Default \code{r=0} (the two distributions are considered independent).
  #'
  #' @details Let two random variables
  #' \eqn{X \sim N(\mu_x, \sigma_x)} and \eqn{Y \sim N(\mu_y, \sigma_y)}
  #' <!-- (\eqn{\mu_x > 0}, \eqn{\mu_y > 0}) -->
  #' with probability densities \eqn{f_X} and \eqn{f_Y}.
  #' 
  #' The parameters of the distribution of the ratio \eqn{Z = \frac{X}{Y}} are:
  #' \eqn{\displaystyle{\beta = \frac{\mu_x}{\mu_y}}},
  #' \eqn{\displaystyle{\rho = \frac{\sigma_y}{\sigma_x}}},
  #' \eqn{\displaystyle{\delta_y = \frac{\sigma_y}{\mu_y}}}
  #' and \eqn{\displaystyle{r = Cor(X, Y) = \frac{Cov(X, Y)}{\sigma_x \sigma_y}}}.
  #' 
  #' \eqn{\mu_x}, \eqn{\sigma_x}, \eqn{\mu_y} and \eqn{\sigma_y} are computed from
  #' \eqn{\beta}, \eqn{\rho} and \eqn{\delta_y} (by fixing arbitrarily \eqn{\mu_x = 1}).
  #' 
  #' If \eqn{X} and \eqn{Y} are independent, i.e. \eqn{r=0}, two random samples
  #' \eqn{\left( x_1, \dots, x_n \right)} and
  #' \eqn{\left( y_1, \dots, y_n \right)} are simulated.
  #' 
  #' If \eqn{X} and \eqn{Y} are not independent,
  #' a sample \eqn{\left( (x_1, y_1), \dots, (x_n, y_n) \right)}
  #' is simpulated using [MASS::mvrnorm()].
  #' 
  #' Then \eqn{\displaystyle{\left( \frac{x_1}{y_1}, \dots, \frac{x_n}{y_n} \right)}} is returned.
  #' 
  #' @return A numeric vector: the produced sample.
  #' 
  #' @seealso [dnormratio()]: probability density of a normal ratio.
  #' 
  #' [pnormratio()]: probability distribution function.
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
  #' Pham-Gia, T., Turkkan, N., Marchand, E. (2006)
  #' Density of the Ratio of Two Normal Random Variables and Applications,
  #' Communications in Statistics - Theory and Methods, 35:9, 1569-1591.
  #' \doi{10.1080/03610920600683689}
  #'
  #' @examples
  #' # First example: ratio of independent variables
  #' beta1 <- 0.15
  #' rho1 <- 5.75
  #' delta1 <- 0.22
  #' rnormratio(20, bet = beta1, rho = rho1, delta = delta1)
  #'
  #' # Second example: ratio of correlated variables
  #' beta2 <- 0.24
  #' rho2 <- 4.21
  #' delta2 <- 0.25
  #' r2 <- 0.8
  #' rnormratio(20, bet = beta2, rho = rho2, delta = delta2, r = r2)
  #'
  #' @importFrom stats rnorm
  #' @importFrom MASS mvrnorm
  #' @export

  muy <- 1
  mux <- bet
  sigmay <- delta
  sigmax <- delta/rho
  
  if (r == 0) {
    x <- rnorm(n, mean = mux, sd = sigmax)
    y <- rnorm(n, mean = muy, sd = sigmay)
  } else {
    Mu <- c(mux, muy)
    Sigma <- matrix(c(sigmax^2, r*sigmax*sigmay, r*sigmax*sigmay, sigmay^2),
                    nrow = 2)
    X <- MASS::mvrnorm(n = n, mu = Mu, Sigma = Sigma)
      x <- X[, 1]; y <- X[, 2]
  }
  
  return(x/y)
}
