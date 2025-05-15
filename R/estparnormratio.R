#' Estimation of the Parameters of a Normal Ratio Distribution
#'
#' Estimation of the parameters of a ratio \eqn{\displaystyle{Z = \frac{X}{Y}}},
#' \eqn{X} and \eqn{Y} being two independent random variables distributed
#' according to Gaussian distributions,
#' using the EM (estimation-maximization) algorithm or variational inference.
#' Depending on the estimation method, the \code{estparnormatio} function calls
#' \code{estparEM} (EM algorithm) or \code{estparVB} (variational Bayes).
#'
#' @name estparnormratio
#' @aliases estparnormratio,estparEM,estparVB
#'
#' @usage estparnormratio(z, method = c("EM", "VB"), eps = 1e-06,
#'                        display = FALSE, mux0 = 1, sigmax0 = 1,
#'                        alphax0 = NULL, betax0 = NULL, muy0 = 1, sigmay0 = 1,
#'                        alphay0 = NULL, betay0 = NULL)
#' @usage estparEM(z, eps = 1e-06,  display = FALSE, #plot = display,
#'                        mux0 = 1, sigmax0 = 1, muy0 = 1, sigmay0 = 1)
#' @usage estparVB(z, eps = 1e-06, display = FALSE, mux0 = 1, sigmax0 = 1,
#'                        alphax0 = 1, betax0 = 1, muy0 = 1, sigmay0 = 1,
#'                        alphay0 = 1, betay0 = 1)
#' @param z numeric.
#' @param method the method used to estimate the parameters of the distribution.
#' It can be \code{"EM"} (expectation-maximization) or \code{"VB"} (Variational Bayes).
#' @param eps numeric. Precision for the estimation of the parameters (see Details).
#' @param display logical. When \code{TRUE} the successive values of the stop criterion
#' (distance between successive values) is printed.
#' @param mux0,sigmax0,muy0,sigmay0 initial values of the means and
#' standard deviations of the \eqn{X} and \eqn{Y} variables. Default:
#' \code{mux0 = 1, sigmax0 = 1, muy0 = 1, sigmay0 = 1}.
#' @param alphax0,betax0,alphay0,betay0 initial values for the variational
#' Bayes method. Omitted if \code{method="EM"}.
#' If \code{method="VB"}, if omitted, they are set to 1.
#' @return A list of 3 elements \code{beta}, \code{rho}, \code{delta}:
#' the estimated parameters of the \eqn{Z} distribution
#' \eqn{\hat{\beta}}, \eqn{\hat{\rho}}, \eqn{\hat{\delta}_y},
#' with three attributes \code{attr(, "epsilon")} (precision of the result),
#' \code{attr(, "k")} (number of iterations) and \code{attr(, "method")} (estimation method).
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
#' with: \eqn{\displaystyle{\beta = \frac{\mu_x}{\mu_y}}},
#' \eqn{\displaystyle{\rho = \frac{\sigma_y}{\sigma_x}}},
#' \eqn{\displaystyle{\delta_y = \frac{\sigma_y}{\mu_y}}}.
#' 
#' and \eqn{_1 F_1\left(a, b; x\right)} is the confluent hypergeometric function
#' (Kummer's function):
#' \deqn{\displaystyle{
#' _1 F_1\left(a, b; x\right) = \sum_{n = 0}^{+\infty}{ \frac{ (a)_n }{ (b)_n } \frac{x^n}{n!} }
#' }}
#' 
#' If \code{method = "EM"}, the means and standard deviations
#' \eqn{\mu_x}, \eqn{\sigma_x}, \eqn{\mu_y} and \eqn{\sigma_y}
#' are estimated with the EM algorithm, as presented in El Ghaziri et al.
#' If \code{method = "VB"}, they are estimated with the variational Bayes method
#' as presented in Bouhlel et al.
#' 
#' Then the parameters \eqn{\beta}, \eqn{\rho}, \eqn{\delta_y} of the \eqn{Z} distribution
#' are computed from these means and standard deviations.
#' 
#' The estimation of \eqn{\mu_x}, \eqn{\sigma_x}, \eqn{\mu_y} and \eqn{\sigma_y}
#' uses an iterative algorithm.
#' The precision for their estimation is given by the \code{eps} parameter.
#' 
#' The computation uses the \code{\link{kummer}} function.
#' 
#' If there are ties in the `z` vector, it generates a warning,
#' as there should be no ties in data distributed among a continuous distribution.
#'
#' @author Pierre Santagostini, Angélina El Ghaziri, Nizar Bouhlel
#'
#' @references El Ghaziri, A., Bouhlel, N., Sapoukhina, N., Rousseau, D.,
#' On the importance of non-Gaussianity in chlorophyll fluorescence imaging.
#' Remote Sensing 15(2), 528 (2023).
#' \doi{10.3390/rs15020528}
#' 
#' Bouhlel, N., Mercier, F., El Ghaziri, A., Rousseau, D., 
#' Parameter Estimation of the Normal Ratio Distribution with Variational Inference.
#' 2023 31st European Signal Processing Conference (EUSIPCO), Helsinki, Finland, 2023, pp. 1823-1827.
#' \doi{10.23919/EUSIPCO58844.2023.10290111}
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
#' # With the EM method:
#' estparnormratio(z1, method = "EM")
#' 
#' # With the variational method:
#' estparnormratio(z1, method = "VB")
#' 
#' # Second example
#' beta2 <- 0.24
#' rho2 <- 4.21
#' delta2 <- 0.25
#' 
#' set.seed(1234)
#' z2 <- rnormratio(800, bet = beta2, rho = rho2, delta = delta2)
#' 
#' # With the EM method:
#' estparnormratio(z2, method = "EM")
#' 
#' # With the variational method:
#' estparnormratio(z2, method = "VB")
#' }
#' 
#' @export

estparnormratio <- function(z, method = c("EM", "VB"), eps = 1e-06,
                            display = FALSE, mux0 = 1, sigmax0 = 1,
                            alphax0 = NULL, betax0 = NULL, muy0 = 1,
                            sigmay0 = 1, alphay0 = NULL, betay0 = NULL) {
  
  method <- match.arg(method)
  
  # Continuous distribution => there should be no ties.
  if(length(unique(z)) < length(z))
    warning("Presence of ties in the data")
  
  res <- switch(method,
         EM = {
           estparEM(z, eps = eps,  display = display, #plot = plot,
                    mux0 = mux0, sigmax0 = sigmax0, muy0 = muy0,
                    sigmay0 = sigmay0)
         },
         VB = {
           if (is.null(alphax0)) alphax0 <- 1
           if (is.null(betax0)) betax0 <- 1
           if (is.null(alphay0)) alphay0 <- 1
           if (is.null(betay0)) betay0 <- 1
           estparVB(z, eps = eps, display = display, #plot = plot,
                    mux0 = mux0, sigmax0 = sigmax0,
                    alphax0 = alphax0, betax0 = betax0,
                    muy0 = muy0, sigmay0 = sigmay0,
                    alphay0 = alphay0, betay0 = betay0)
         })
  
  attr(res, "method") <- method

  return(res)
}

estparEM <- function(z, eps = 1e-06,  display = FALSE, #plot = display,
                     mux0 = 1, sigmax0 = 1, muy0 = 1, sigmay0 = 1) {
  #' @rdname estparnormratio
  #' @export
  
  # kummA <- function(x) {
  #   Re(kummer(2, 1.5, x))/Re(kummer(1, 0.5, x, eps = eps))
  # }
  # 
  # kummB <- function(x) {
  #   Re(kummer(2, 0.5, x))/Re(kummer(1, 0.5, x, eps = eps))
  # }
  
  variancex0 <- sigmax0^2; variancey0 <- sigmay0^2
  
  # Number of observations
  n <- length(z)
  
  #initialisation
  mux <- mux0; muy <- muy0
  variancex <- variancex0; variancey <- variancey0
  
  iter <- 0
  diff <- 100
  
  while ((!is.nan(diff)) & (diff > eps)) {
    
    # print(iter)
    
    # Current values: mux, variancex; muy, variancey
    # Updates: mux1, variancex1; muy1, variancey1
    
    mu <- 0.5*( (z^2/variancex) + 1/variancey )
    invmu <- 1/mu
    gam <- muy/variancey + (mux/variancex)*z
    gmu <- gam^2*invmu/4
    
    # cat("\n"); print(range(gmu))
    
    # if (any(gmu > 50)) {
    #   A <- 0.5
    #   B <- gmu
    # } else {
    num1 <- kummer(2, 1.5, gmu, eps = eps)
    num2 <- kummer(2, 0.5, gmu, eps = eps)
    denom <- kummer(1, 0.5, gmu, eps = eps)
    A <- num1/denom
    B <- num2/denom
    # }
    
    # A <- sapply(gam^2/(4*mu), kummA)
    # 
    # B <- sapply(gam^2/(4*mu), kummB)
    
    mux1 <- (1/n)*sum(z*(gam*invmu)*A)
    # mux1 <- 1
    variancex1 <- (1/n)*sum((z^2)*invmu*B) - mux1^2
    # variancex1 <- 1
    
    muy1 <- (1/n)*sum((gam/mu)*A)
    # muy1 <- 1
    variancey1 <- (1/n)*sum(invmu*B )-muy1^2
    # variancey1 <- 1
    
    diff <- max(abs(c(mux1 - mux,
                      variancex1 - variancex,
                      muy1 - muy,
                      variancey1 - variancey)
    ))
    
    if (display) cat(diff, "  ")
    
    mux <- mux1; variancex <- variancex1
    muy <- muy1; variancey <- variancey1
    
    iter <- iter + 1
    
    # print(c(
    #   gmu.min = min(gmu), gmu.max = max(gmu),
    #   mux = mux, variancex = variancex,
    #   muy = muy, variancey = variancey,
    #   beta = mux/muy, rho = sqrt(variancey)/sqrt(variancex), delta = sqrt(variancey)/muy,
    #   diff = diff
    # ))
    # cat("\n")
  }
  if (display) cat("\n")
  
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

estparVB <- function(z, eps = 1e-06, display = FALSE, mux0 = 1, sigmax0 = 1,
                     alphax0 = 1, betax0 = 1, muy0 = 1, sigmay0 = 1,
                     alphay0 = 1, betay0 = 1) {
  #' @rdname estparnormratio
  #' @export
  
  # print(c(mux0 = mux0, variancex0 = variancex0, alphax0 = alphax0,
  #         betax0 = betax0, muy0 = muy0, variancey0 = variancey0,
  #         alphay0 = alphay0, betay0 = betay0))
  
  variancex0 <- sigmax0^2; variancey0 <- sigmay0^2
  
  # Sample size
  n <- length(z)
  
  # # Initialization
  # mux0 <- 1; muy0 <- 2
  # variancex0 <- 2; variancey0 <- 1
  # alphax0 <- betax0 <- alphay0 <- betay0 <- 1
  
  mux <- mux0; muy <- muy0
  variancex <- variancex0; variancey <- variancey0
  alphax <- alphax0 + n/2; alphay <- alphay0 + n/2
  betax <- betax0; betay <- betay0
  
  iter <- 0
  delta2 <- 100
  # Delta <- numeric(0)
  
  # while (delta2 >= eps^2 & iter <= 1e4) {
  # while (iter <= 1e4) {
  while (delta2 >= eps^2) {
      
    iter <- iter + 1
    
    # print(c(iter = iter, delta2 = delta2,
    #         mux = mux, variancex = variancex, alphax = alphax,
    #         betax = betax, muy = muy, variancey = variancey,
    #         alphay = alphay, betay = betay))
    
    mux1 <- mux; variancex1 <- variancex
    muy1 <- muy; variancey1 <- variancey
    # alphax0 <- alphax; alphay0 <- alphay
    betax1 <- betax; betay1 <- betay
    
    # print(summary(z))

    abx <- alphax/betax
    aby <- alphay/betay
    
    A <- 0.5*(z^2 * abx + aby)
    B <- mux * z * abx + muy * aby
    BA <- B^2/(4*A)
    
    # print(summary(BA))
    
    if (any(BA > 50) & iter < 10) {
      # cat(range(BA), "> 50\n")
      rapport1 <- 0.5
      rapport2 <- BA
    } else {
      # cat(range(BA), "<= 50\n")
      num1 <- kummer(2, 1.5, BA, eps = eps^2)
      num2 <- kummer(2, 0.5, BA, eps = eps^2)
      denom <- kummer(1, 0.5, BA, eps = eps^2)
      
      rapport1 <- num1/denom
      rapport2 <- num2/denom
      
      # cat("\n", iter, "  "); print(c(range(BA)))
    }
    
    invA <- 1/A
    Ey <- B*invA * rapport1
    sumEy <- sum(Ey)
    sumzEy <- sum(z*Ey)
    
    Ey2 <- invA * rapport2
    sumEy2 <- sum(Ey2)
    sumz2Ey2 <- sum(z^2*Ey2)
    
    betax <- betax0 + 0.5 * sumz2Ey2 + n/2 * (variancex + mux^2) - mux * sumzEy
    
    variancex <- 1/(n * abx + 1/variancex0)
    # variancex <- 1
    
    mux <- variancex * (mux0/variancex0 + abx * sumzEy)
    # mux <- 1
    
    betay <- betay0 + 0.5 * sumEy2 + n/2 * (variancey + muy^2) - muy * sumEy
    
    variancey <- 1/(n*aby + 1/variancey0)
    # variancey <- 1
    
    muy <- variancey * (muy0/variancey0 + aby * sumEy)
    # muy <- 1
    
    # print(c(iter = iter,
    #   mux = mux, variancex = variancex, alphax = alphax, betax = betax, # muy = muy,
    #   variancey = variancey, alphay = alphay, betay = betay,
    #   beta = mux/muy, rho = sqrt((betay * (alphax-1)) / (betax * (alphay-1))), delta = 1/muy * sqrt(betay / (alphay-1))
    # ))
    
    delta2 <- sum( c(mux - mux1, variancex - variancex1,
                     muy - muy1, variancey - variancey1,
                     betax - betax1, betay - betay1)^2 )
    
    if (display) cat(delta2, "  ")
    
    # cat(delta2 > 0.71); if (all(BA <= 50)) {print(c(range(num1), range(num2), range(denom)))} else {print(c(range(BA)))}
    
    # print(c(iter, delta2))
    # Delta <- c(Delta, delta2)
    
    # print(c(
    #   iter = iter,
    #   BA.min = min(BA), BA.max = max(BA),
    #   mux = mux, variancex = variancex, betax = betax,
    #   muy = muy, variancey = variancey, betay = betay,
    #   beta = mux/muy,
    #   rho = sqrt((betay * (alphax-1)) / (betax * (alphay-1))),
    #   delta = 1/muy * sqrt(betay / (alphay-1)),
    #   delta2 = delta2
    # ))
    # cat("\n")
  }
  # print(length(Delta)); print(iter)
  # plot(1:iter, Delta, type = "l", ylim = c(0, 0.8))
  if (display) cat("\n")
  
  res <- list()
  
  res$beta <- mux/muy
  res$rho <- sqrt((betay * (alphax-1)) / (betax * (alphay-1)))
  res$delta <- 1/muy * sqrt(betay / (alphay-1))
  
  attr(res, "k") <- iter
  attr(res, "epsilon") <- eps
  
  return(res)
}
