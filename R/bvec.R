#' Function to compute m by 1 column vector \eqn{I} (independent and dependent
#' cases) of integrated emulator with squared exponential kernel.
#'
#' @param w an \eqn{m} (number of design points) by \eqn{d} (dimensions)
#' matrix of design points.
#' @param z a vector of length d1 of exogenous inputs.
#' @param mu a vector of predictive mean values of length k.
#' @param Sigma a variance-covariance matrix.
#' @param z_design a design matrix of exogenous inputs.
#' @param params a list with parameter values.
#' @param ind logical argument indicating if we adopt independence assumption.
#' @param det_comp a determinant component in the calculations.
#' @param LS.inv an inverse of the sum of covariance matrix and Lambda.
#'
#' @return a vector of covariances between input of interest and design point.
#' @export
#'
#' @examples
bvec <- function(w, z =  NULL, mu, Sigma, z_design = NULL, params,
                 ind = "TRUE", det_comp = NULL, LS.inv = NULL) {

  m <- dim(w)[1] # number of training points
  d <- dim(w)[2] # number of dimensions
  d1 <- dim(z_design)[2]
  I <- matrix(nrow = m, ncol = 1)

  if(ind == "TRUE") {
    for(i in 1:m) {
      b <- xi_fun(w[i, ], mu, Sigma, params)
      if(is.null(z_design)) {
        I[i, 1] <- prod(b)
      } else {
        b1 <- sapply(1:d1, function(x) exp(-(z_design[i, x] - z[x])^2/
                                             (params$delta_par[d+x])^2))
        I[i, 1] <- b*prod(b1)
      }

    }
  } else {
    for(i in 1:m) {
      b <- xi_fun(w[i, ], mu, Sigma, params, ind, det_comp, LS.inv)
      if(is.null(z_design)) {
        I[i, 1] <- b
      } else {
        b1 <- sapply(1:d1, function(x) exp(-(z_design[i, x] - z[x])^2/
                                             (params$delta_par[d+x])^2))
        I[i, 1] <- b*prod(b1)
      }
    }
  }
  return(I)
}
