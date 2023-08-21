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
#'     \item{\code{\link{FHme}}}{Gives the EBLUP for each domain based on Fay-Herriot with measurement error model.}
#'     \item{\code{\link{mse_FHme}}}{Gives the MSE for each domain using the jackknife method.}
#'     }
#'
#' @references Ybarra, L.M. and Lohr, S. L. (2008). Small area estimation when auxiliary information is measured with error. Biometrika 95, 919-931.
#' @importFrom MASS ginv
#' @import stats
"_PACKAGE"


