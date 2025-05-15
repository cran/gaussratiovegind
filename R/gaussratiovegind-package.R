#' @description
#' This package provides tools for the distribution of the ratio of Gaussian distributions:
#' \describe{
#' \item{\link{dnormratio}}{Probability density function}
#' \item{\link{pnormratio}}{Cumulative distribution function}
#' \item{\link{rnormratio}}{Sample simulation}
#' \item{\link{estparnormratio}}{Parameter estimation using the EM (expectation-maximization) algorithm or the variational Bayes method}
#' }
#' @details
#' It is well known that the distribution of Gaussian ratios does not follow a Gaussian distribution.
#' 
#' This issue arises in many vegetation indices based on the ratio of spectral reflectance values
#' from different bands (e.g., near-infrared, red, or chlorophyll content bands) commonly used
#' in remote sensing and agricultural studies. However, the lack of awareness among users of
#' vegetation indices about this non-Gaussian nature often leads to incorrect statistical modeling
#' and interpretation.
#' 
#' This package provides tools to accurately handle and analyze the distribution of such ratios.
#' @references El Ghaziri, A., Bouhlel, N., Sapoukhina, N., Rousseau, D.,
#' On the importance of non-Gaussianity in chlorophyll fluorescence imaging.
#' Remote Sensing 15(2), 528 (2023).
#' \doi{10.3390/rs15020528}
#' 
#' Bouhlel, N., Mercier, F., El Ghaziri, A., Rousseau, D., 
#' Parameter Estimation of the Normal Ratio Distribution with Variational Inference.
#' 2023 31st European Signal Processing Conference (EUSIPCO), Helsinki, Finland, 2023, pp. 1823-1827.
#' \doi{10.23919/EUSIPCO58844.2023.10290111}
#' @keywords internal
"_PACKAGE"

NULL
