estparnormratio <- function(z, eps = 1e-06) {
  #' Estimation of the Parameters of a Normal Ratio Distribution
  #'
  #' Estimation of the parameters of <!-- the distribution of a ratio of two distributions -->
  #' a ratio \eqn{\displaystyle{Z = \frac{X}{Y}}},
  #' \eqn{X} and \eqn{Y} being two independant random variables distributed
  #' according to Gaussian distributions,
  #' using the EM (estimation-maximization) algorithm.
  #'
  #' @aliases estparnormratio
  #'
  #' @usage estparnormratio(z, eps = 1e-6)
  #' @param z numeric matrix or data frame.
  #' @param eps numeric. Precision for the estimation of the parameters.
  #' @return A list of 3 elements \code{beta}, \code{rho}, \code{delta}:
  #' the estimated parameters of the \eqn{Z} distribution
  #' \eqn{\hat{\beta}}, \eqn{\hat{\rho}}, \eqn{\hat{\delta}_y},
  #' with two attributes \code{attr(, "epsilon")} (precision of the result) and \code{attr(, "k")} (number of iterations).
  #'
  #' @details Let a random variable: \eqn{\displaystyle{Z = \frac{X}{Y}}},
  #' 
  #' \eqn{X} and \eqn{Y} being normally distributed:
  #' \eqn{X \sim N(\mu_x, \sigma_x)} and \eqn{Y \sim N(\mu_y, \sigma_y)}.
  #' 
  #' The density probability of \eqn{Z} is:
  #' \deqn{\displaystyle{
  #' f_Z(z; \beta, \rho, \delta_y) = \frac{\rho}{\pi (1 + \rho^2 z^2)} \ \exp{\left(-\frac{\rho^2 \beta^2 + 1}{2\delta_y^2}\right)} \ {}_1 F_1\left( 1, \frac{1}{2}; \frac{1}{2 \delta_y^2} \frac{(1 + \beta \rho^2 z)^2}{1 + \rho^2 z^2} \right)
  #' }}
  #' 
  #' with: \eqn{\displaystyle{\beta = \frac{\mu_x}{\mu}_y}},
  #' \eqn{\displaystyle{\rho = \frac{\sigma_y}{\sigma_x}}},
  #' \eqn{\displaystyle{\delta_y = \frac{\sigma_y}{\mu_y}}}.
  #' 
  #' and \eqn{_1 F_1\left(a, b; x\right)} is the confluent hypergeometric function
  #' (Kummer's function):
  #' \deqn{\displaystyle{
  #' _1 F_1\left(a, b; x\right) = \sum_{n = 0}^{+\infty}{ \frac{ (a)_n }{ (b)_n } \frac{x^n}{n!} }
  #' }}
  #' 
  #' The parameters \eqn{\beta}, \eqn{\rho}, \eqn{\delta_y} of the \eqn{Z} distribution
  #' are estimated with the EM algorithm, as presented in El Ghaziri et al.
  #' The computation uses the \code{\link{kummerM}} function.
  #' 
  #' This uses an iterative algorithm.
  #'
  #' The precision for the estimation of the parameters is given by the \code{eps} parameter.
  #'
  #' @author Pierre Santagostini, Angélina El Ghaziri, Nizar Bouhlel
  #'
  #' @references El Ghaziri, A., Bouhlel, N., Sapoukhina, N., Rousseau, D.,
  #' On the importance of non-Gaussianity in chlorophyll fluorescence imaging.
  #' Remote Sensing 15(2), 528 (2023).
  #' \doi{10.3390/rs15020528}
  #'
  #' @seealso [dnormratio()]: probability density of a normal ratio.
  #' 
  #' [rnormratio()]: sample simulation.
  #' 
  #' @examples
  #' \donttest{
  #' # First example
  #' beta1 <- 0.15
  #' rho1 <- 5.75
  #' delta1 <- 0.22
  #' 
  #' set.seed(1234)
  #' z1 <- rnormratio(800, bet = beta1, rho = rho1, delta = delta1)
  #' 
  #' estparnormratio(z1)
  #' 
  #' # Second example
  #' beta2 <- 0.24
  #' rho2 <- 4.21
  #' delta2 <- 0.25
  #' 
  #' set.seed(1234)
  #' z2 <- rnormratio(800, bet = beta2, rho = rho2, delta = delta2)
  #' 
  #' estparnormratio(z2)
  #' }
  #' 
  #' @importFrom graphics plot
  #' @export

  kummA <- function(x) {
    Re(kummerM(2, 1.5, x))/Re(kummerM(1, 0.5, x, eps = eps))
  }

  kummB <- function(x) {
    Re(kummerM(2, 0.5, x))/Re(kummerM(1, 0.5, x, eps = eps))
  }

  # Number of observations
  n <- length(z)

  #initialisation
  mux <- 1; muy <- 1
  variancex <- 1; variancey <- 1
  iter <- 0
  diff <- 100

  while (diff > eps) {
    
    # Current values: mux, variancex; muy, variancey
    # Updates: mux1, variancex1; muy1, variancey1

    mu <- 0.5*( (z^2/variancex) + 1/variancey )
    gam <- muy/variancey + (mux/variancex)*z

    A <- sapply(gam^2/(4*mu), kummA)

    B <- sapply(gam^2/(4*mu), kummB)

    mux1 <- (1/n)*sum(z*(gam/mu)*A)
    variancex1 <- (1/n)*sum((z^2)*(1/mu)*B) - mux1^2

    muy1 <- (1/n)*sum((gam/mu)*A)
    variancey1 <- (1/n)*sum((1/mu)*B )-muy1^2 # 1

    diff <- max(abs(c(mux1 - mux,
                      variancex1 - variancex,
                      muy1 - muy,
                      variancey1 - variancey)
                    ))
    mux <- mux1; variancex <- variancex1
    muy <- muy1; variancey <- variancey1
    
    iter <- iter + 1
    
  }

  sigmax <- sqrt(variancex); sigmay <- sqrt(variancey)

  res <- list()
  # res$x <- list(mu = mux, sigma = sigmax)
  # res$y <- list(mu = muy, sigma = sigmay)

  res$beta <- mux/muy
  res$rho <- sigmay/sigmax
  res$delta <- sigmay/muy

  attr(res, "k") <- iter
  attr(res, "epsilon") <- eps

  return(res)
}
