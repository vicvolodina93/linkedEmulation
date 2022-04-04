#' Fuction to compute \eqn{\psi_{jk}} for d inputs of integrated emulator
#' for (independent case) and (dependent case) with a squared exponential kernel.
#'
#' @param w_j a vector of size k.
#' @param mu a vector of predictive mean values of length k.
#' @param Sigma a variance-covariance matrix.
#' @param params a list with parameter values.
#' @param ind logical indicating if we adopt independence assumption.
#' @param Lambda.mat Lambda matrix.
#' @param LS.inv logical indicating if we adopt independence assumption.
#'
#' @return a vector.
#' @export
#'
#' @examples
psi_fun <- function(w_j, mu, Sigma, params, ind = "TRUE",
                    Lambda.mat = NULL, LS.inv = NULL) {

  d <- length(w_j)
  if(ind == "TRUE") {
    #sigma_sq <- diag(Sigma)
    psi_val <- sapply(1:d, function(k)
      1/sqrt(1+Sigma[k, k]/params$delta_par[k]^2)*
        exp(-(mu[k]-w_j[k])^2/(2*Sigma[k, k]+params$delta_par[k]^2))*
        (2*Sigma[k, k]*w_j[k] + params$delta_par[k]^2*mu[k])/
        (2*Sigma[k, k]+params$delta_par[k]^2))

  } else {
    if(is.null(Lambda.mat)) {
      Lambda.mat <- diag(params$delta_par^2/2, nrow = d, ncol = d)
      LS.inv <- solve(Lambda.mat + Sigma)
    }
    det_comp <- 1/sqrt(det((Lambda.mat+Sigma)%*%solve(Lambda.mat)))
    xi_val <- xi_fun(w_j, mu, Sigma, params, ind = "FALSE",
                     det_comp, LS.inv)
    psi_val <- c()
    for(k in 1:d) {
      e_k <- matrix(rep(0, d), ncol = d)
      e_k[, k] <- 1
      psi_val[k] <- e_k%*%(Lambda.mat%*%LS.inv%*%mu+Sigma%*%LS.inv%*%w_j)*xi_val
    }
  }
  return(psi_val)
}
