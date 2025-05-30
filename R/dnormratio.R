dnormratio <- function (z, bet, rho, delta) {
  #' Density Function of a Normal Ratio Distribution
  #'
  #' Density of the ratio of two independent Gaussian distributions<!-- with strictly positive means and variances -->.
  #'
  #' @aliases dnormratio
  #'
  #' @usage dnormratio(z, bet, rho, delta)
  #' @param z length \eqn{p} numeric vector.
  #' @param bet,rho,delta numeric values. The parameters \eqn{(\beta, \rho, \delta_y)} of the distribution, see Details.
  #'
  #' @details Let two independant random variables
  #' \eqn{X \sim N(\mu_x, \sigma_x)} and \eqn{Y \sim N(\mu_y, \sigma_y)}.
  #' <!-- with \eqn{\mu_x > 0} and \eqn{\mu_y > 0}. -->
  #' 
  #' If we denote \eqn{\beta = \frac{\mu_x}{\mu_y}}, \eqn{\displaystyle{\rho = \frac{\sigma_y}{\sigma_x}}}
  #' and \eqn{\displaystyle{\delta_y = \frac{\sigma_y}{\mu_y}}},
  #' the probability distribution function of the ratio \eqn{\displaystyle{Z = \frac{X}{Y}}}
  #' is given by:
  #' \deqn{\displaystyle{ f_Z(z; \beta, \rho, \delta_y) = \frac{\rho}{\pi (1 + \rho^2 z^2)} \left[ \exp{\left(-\frac{\rho^2 \beta^2 + 1}{2\delta_y^2}\right)} + \sqrt{\frac{\pi}{2}} \ q \ \text{erf}\left(\frac{q}{\sqrt{2}}\right) \exp\left(-\frac{\rho^2 (z-\beta)^2}{2 \delta_y^2 (1 + \rho^2 z^2)}\right) \right] }}
  #' with \eqn{\displaystyle{ q = \frac{1 + \beta \rho^2 z}{\delta_y \sqrt{1 + \rho^2 z^2}} }}
  #' and \eqn{\displaystyle{ \text{erf}\left(\frac{q}{\sqrt{2}}\right) = \frac{2}{\sqrt{\pi}} \int_0^{\frac{q}{\sqrt{2}}}{\exp{(-t^2)}\ dt} }}
  #' 
  #' Another expression of this density, used by the [estparnormratio()] function, is:
  #' \deqn{\displaystyle{
  #' f_Z(z; \beta, \rho, \delta_y) = \frac{\rho}{\pi (1 + \rho^2 z^2)} \ \exp{\left(-\frac{\rho^2 \beta^2 + 1}{2\delta_y^2}\right)} \ {}_1 F_1\left( 1, \frac{1}{2}; \frac{1}{2 \delta_y^2} \frac{(1 + \beta \rho^2 z)^2}{1 + \rho^2 z^2} \right)
  #' }}
  #' 
  #' where \eqn{_1 F_1\left(a, b; x\right)} is the confluent hypergeometric function
  #' (Kummer's function):
  #' \deqn{\displaystyle{
  #' _1 F_1\left(a, b; x\right) = \sum_{n = 0}^{+\infty}{ \frac{ (a)_n }{ (b)_n } \frac{x^n}{n!} }
  #' }}
  #' 
  #' @return Numeric: the value of density.
  #' 
  #' @seealso [pnormratio()]: probability distribution function.
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
  #' dnormratio(0, bet = beta1, rho = rho1, delta = delta1)
  #' dnormratio(0.5, bet = beta1, rho = rho1, delta = delta1)
  #' curve(dnormratio(x, bet = beta1, rho = rho1, delta = delta1), from = -0.1, to = 0.7)
  #'
  #' # Second example
  #' beta2 <- 2
  #' rho2 <- 2
  #' delta2 <- 2
  #' dnormratio(0, bet = beta2, rho = rho2, delta = delta2)
  #' dnormratio(0.5, bet = beta2, rho = rho2, delta = delta2)
  #' curve(dnormratio(x, bet = beta2, rho = rho2, delta = delta2), from = -0.1, to = 0.7)
  #'
  #' @export

  q <- (1+bet*rho^2*z)/(delta*sqrt(1+rho^2*z^2))
  denom = rho/(pi*(1+rho^2*z^2))*(
    exp(-(rho^2*bet^2+1)/(2*delta^2)) +
      sqrt(pi/2)*q*.erf(q/sqrt(2))*exp(-0.5*(rho^2*(z-bet)^2)/(delta^2*(1+rho^2*z^2)))
  )
  return(denom)
}

.erf <- Vectorize(function(x) 2*pnorm(sqrt(2)*x) - 1)
