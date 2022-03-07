#' Design matrix for w of the integrated emulator.
#'
#' A dataset containing the values of design matrix w.
#'
#' @format A tibble with 120 rows and 3 variables:
#' \describe{
#'   \item{x1}{input of w matrix, quarterly heating demand}
#'   \item{x2}{input of w matrix, a quarterly gas price.}
#'   \item{x3}{input of w matrix, a quarterly electricity price.}
#' }
#' @source latin hypercube
"w_design"

#' Design matrix for z of the integrated emulator.
#'
#' A dataset containing the values of design matrix z.
#'
#' @format A tibble with 120 rows and 2 variables:
#' \describe{
#'   \item{z1}{input of z matrix, efficiency of gas boiler.}
#'   \item{z2}{input of z matrix, efficiency of heat pump.}
#' }
#' @source latin hypercube
"z_design"

#' Output vector of the energy system model.
#'
#' A vector containing the simulations of energy system model.
#'
#' @format A vector of size 120:
#' \describe{
#'   \item{y}{a vector containing the simulations.}
#' }
#' @source output of simulations
"y_design"

#' Validation (test) matrix for w of the integrated emulator.
#'
#' A dataset containing the values of validation matrix w.
#'
#' @format A tibble with 36 rows and 3 variables:
#' \describe{
#'   \item{x1}{input of w matrix, quarterly heating demand}
#'   \item{x2}{input of w matrix, a quarterly gas price.}
#'   \item{x3}{input of w matrix, a quarterly electricity price.}
#' }
#' @source latin hypercube
"w_valid"

#' Validation (test) matrix for z of the integrated emulator.
#'
#' A dataset containing the values of validation matrix z.
#'
#' @format A tibble with 36 rows and 2 variables:
#' \describe{
#'   \item{z1}{input of z matrix, efficiency of gas boiler.}
#'   \item{z2}{input of z matrix, efficiency of heat pump.}
#' }
#' @source latin hypercube
"z_valid"


#' Output vector of the energy system model used for validation.
#'
#' A vector containing the simulations of energy system model.
#'
#' @format A vector of size 36:
#' \describe{
#'   \item{y}{a vector containing the simulations.}
#' }
#' @source output of simulations
"y_valid"

#' Predictive means for w of the integrated emulator.
#'
#' A dataset containing the values of predictive means for w.
#'
#' @format A tibble with 36 rows and 3 variables:
#' \describe{
#'   \item{mu1}{predictive mean for quarterly heating demand}
#'   \item{mu2}{predictive mean for a quarterly gas price.}
#'   \item{mu3}{predictive mean for a quarterly electricity price.}
#' }
#' @source outputs of GP emulation and two DLMs.
"mean_pred"

#' Predictive variance-covariance for w of the integrated emulator.
#'
#' A dataset containing the values of variance and covariances for w.
#'
#' @format A tibble with 36 rows and 6 variables:
#' \describe{
#'   \item{var1}{predictive variance for quarterly heating demand.}
#'   \item{var2}{predictive variance for a quarterly gas price.}
#'   \item{var3}{predictive variance for a quarterly electricity price.}
#'   \item{covar12}{predictive covariance between x1 and x2 inputs.}
#'   \item{covar13}{predictive covariance between x1 and x3 inputs.}
#'   \item{covar23}{predictive covariance between x2 and x3 inputs.}
#' }
#' @source outputs of GP emulation and two DLMs.
"covar_pred"


