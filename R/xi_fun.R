#' Function to compute \eqn{\xi_{ik}} (independent case) and \eqn{\xi_{i}}
#' (dependent case) of integrated emulator with a squared exponential kernel.
#'
#' @param w_i a vector of size k.
#' @param mu a vector of predictive mean values of length k.
#' @param Sigma a variance-covariance matrix.
#' @param params a list with parameter values.
#' @param ind logical indicating if we adopt independence assumption.
#' @param det_comp a determinant component in the calculations.
#' @param LS.inv an inverse of the sum of covariance matrix and Lambda.
#'
#' @return a vector (independent case) or a single value (dependent case).
#' @export
#'
#' @examples
xi_fun <- function(w_i, mu, Sigma, params, ind = "TRUE",
                   det_comp = NULL, LS.inv = NULL) {

  d <- length(w_i)
  if(ind == "TRUE") {
    #sigma_sq <- diag(Sigma)
    xi_val <- sapply(1:d, function(k)
      1/(sqrt(1+2*Sigma[k, k]/params$delta_par[k]^2))*
        exp(-(mu[k]-w_i[k])^2/(2*Sigma[k, k]+params$delta_par[k]^2)))
  } else {
    if(is.null(det_comp)) {
      Lambda.mat <- diag(params$delta_par^2/2, nrow = d, ncol = d)
      LS.inv <- solve(Lambda.mat + Sigma)
      det_comp <- 1/sqrt(det((Lambda.mat+Sigma)%*%solve(Lambda.mat)))
    }
    xi_val <- det_comp*exp(-1/2*t(w_i-mu)%*%LS.inv%*%(w_i-mu))
  }
  return(xi_val)
}
