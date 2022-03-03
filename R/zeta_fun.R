#' Function to compute \eqn{\zeta_{ijk}} (independent case) and \eqn{\zeta_{ij}}
#' (dependent case) of integrated emulator with a squared exponential kernel.
#'
#' @param w_i a vector of size k.
#' @param w_j a vector of size k.
#' @param mu a vector of predictive mean values of lenght k.
#' @param Sigma a variance-covariance matrix.
#' @param params a list with parameter values.
#' @param ind logical indicating if we adopt independece assumption.
#' @param det_comp a determinant component in the calculations.
#' @param Gamma.inv an inverse of Gamma matrix.
#' @param GS.inv an inverse of the sum of covariance matrix and Gamma.
#'
#' @return a vector (independent case) or a single value.
#' @export
#'
#' @examples
zeta_fun <- function(w_i, w_j, mu, Sigma, params, ind = "TRUE",
                     det_comp = NULL, Gamma.inv = NULL, GS.inv = NULL) {

  d <- length(w_i)

  if(ind == "TRUE") {
    sigma_sq <- diag(Sigma)
    zeta_val <- sapply(1:d, function(k) 1/sqrt(1+4*sigma_sq[k]/params$delta_par[k]^2)*
                         exp(-(((w_i[k]+w_j[k])/2 - mu[k])^2)/(params$delta_par[k]^2/2+
                                                                 2*sigma_sq[k])-
                               (w_i[k]-w_j[k])^2/(2*params$delta_par[k]^2)))
  } else {
    if(is.null(det_comp)) {
      Gamma.mat <- diag(params$delta_par^2/4, nrow = d, ncol = d)
      Gamma.inv <- solve(Gamma.mat)
      GS.inv <- solve(Gamma.mat + Sigma)
      det_comp <- 1/sqrt(det((Gamma.mat+Sigma)%*%Gamma.inv))
    }

    diff_w <- w_i - w_j
    sum_w <- w_i + w_j

    zeta_val <- det_comp*exp(-1/8*t(diff_w)%*%Gamma.inv%*%diff_w)*
      exp(-1/2*t(sum_w/2-mu)%*%GS.inv%*%(sum_w/2-mu))
  }
  return(zeta_val)
}
