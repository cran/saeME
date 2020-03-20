#' saeME: Small Area Estimation with Measurement Error
#'
#' The sae with measurement error provides function for small area estimation when auxiliary variable is measured with error,
#' and function for mean squared error estimation using jackknife method.
#' This package implement model of Fay Herriot with Measurement Error developed by Ybarra and Lohr (2008).
#'
#' @section Authors:
#' Muhammad Rifqi Mubarak, Azka Ubaidillah
#'
#' @section Email:
#' Muhammad Rifqi Mubarak \email{16.9304@@stis.ac.id}
#'
#' @section Functions:
#' \describe{
#'     \item{\code{\link{FHme}}}{Gives the EBLUP of SAE with Measurement Error}
#'     \item{\code{\link{mse_FHme}}}{Gives the MSE of SAE with Measurement Error}
#'     }
#'
#' @references Ybarra, L.M. and Lohr, S. L. (2008). Small area estimation when auxiliary information is measured with error. Biometrika 95, 919-931.
#' @docType package
#' @name saeME-package
#' @import expm
#' @import MASS
#' @import stats
NULL
