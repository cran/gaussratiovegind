#' Statistics on Chlorophyll Fluorescence Parameters
#' 
#' Mean and standard deviation values on healthy and diseased tissues
#' of chlorophyll fluorescence parameters \eqn{F_0} (minimum fluorescence)
#' and \eqn{F_m} (maximum fluorescence) for a dataset of *Arabidopsis thaliana*
#' plants infected with fungal pathogen data;
#' parameters of the distribution of the ratio
#' \eqn{\displaystyle{\frac{F_v}{F_m} = \frac{F_m - F_0}{F_m}}}.
#' 
#' @aliases arabidopsis
#' @format A data frame with 10 rows and 6 columns:
#' \describe{
#'    \item{time}{times of the acquisition of chlorophyll fluorescence images}
#'    \item{condition}{indicates if the plant was inoculated: \code{healthy} (inoculated with water) or \code{diseased} (inoculated with the pathogen)}
#'    \item{mF0, sF0}{Mean and standard deviation values of the chlorophyll parameter \eqn{F_0}}
#'    \item{mFm, sFm}{Mean and standard deviation values of the chlorophyll parameter \eqn{F_m}}
#'    \item{beta, rho, delta}{the \eqn{\beta}, \eqn{\rho} and \eqn{\delta_y} parameters of the distribution of \eqn{\displaystyle{\frac{F_v}{F_m} = \frac{F_m - F_0}{F_m}}} (distributed according to a normal ratio distribution, see Details)}
#' }
#' @details
#' On each leaf picture, the \eqn{F_0} and \eqn{F_m} values are normally distributed.
#' Hence, \eqn{\displaystyle{\frac{F_0}{F_m}}} is a ratio of two normal distributions.
#' 
#' Let \eqn{\mu_{F_0}} and \eqn{\sigma_{F_0}} the mean and standard deviation of \eqn{F_0}
#' and \eqn{\mu_{F_m}} and \eqn{\sigma_{F_m}} the mean and standard deviation of \eqn{F_m}.
#' The parameters \eqn{\beta}, \eqn{\rho} and \eqn{\delta_y} are given by:
#' \deqn{\beta = \frac{\mu_{F_0}}{\mu_{F_m}}}
#' \deqn{\rho = \frac{\sigma_{F_m}}{\sigma_{F_0}}}
#' \deqn{\delta_y = \frac{\sigma_{F_m}}{\mu_{F_m}}}
#' 
#' @references El Ghaziri, A., Bouhlel, N., Sapoukhina, N., Rousseau, D.,
#' On the importance of non-Gaussianity in chlorophyll fluorescence imaging.
#' Remote Sensing 15(2), 528 (2023).
#' \doi{10.3390/rs15020528}
#' 
#' Pavicic, M., Overmyer, K., Rehman, A.u., Jones, P., Jacobson, D., Himanen, K.
#' Image-Based Methods to Score Fungal Pathogen Symptom Progression and Severity in Excised Arabidopsis Leaves.
#' Plants, 10, 158 (2021).
#' \doi{10.3390/plants10010158}
"arabidopsis"