#' Function to compute \eqn{B} matrix of Linked Emulator with squared
#' exponential correlation function.
#'
#' @param w an \eqn{m} (number of design points) by \eqn{d} (dimensions)
#' matrix of design points
#' @param z a vector of length d1 of exogenous inputs.
#' @param mu a vector of length \eqn{d} of predictive mean.
#' @param Sigma a \eqn{d} by \eqn{d} variance-covariance matrix.
#' @param z_design a design matrix of exogenous inputs.
#' @param params a list with parameter values.
#' @param ind logical argument indicating if we adopt independence
#' assumption for \eqn{W}.
#'
#' @return an \eqn{d} by \eqn{m} matrix \eqn{B}.
#' @export
#'
#' @examples
Bmatrix <- function(w, z = NULL, mu, Sigma, z_design = NULL, params,
                    ind = "TRUE") {

  m <- dim(w)[1] # number of training points
  d <- dim(w)[2] # number of dimensions
  d1 <- dim(z_design)[2] # number of dimensions for z

  B <- matrix(nrow = d, ncol = m)
  for(j in 1:m) {
    if(!is.null(z_design)) {
      z_vec <- sapply(1:d1, function(x)
        exp(-(z_design[j, x] - z[x])^2/(params$delta_par[d+x])^2))
    }
    if(ind == "TRUE") {
      xi_vec <- xi_fun(w[j, ], mu, Sigma, params, ind)
      psi_vec <- psi_fun(w[j, ], mu, Sigma, params, ind)
      for(l in 1:d) {
        if(!is.null(z_design)) {
          B[l, j] = psi_vec[l]*prod(xi_vec[-l])*prod(z_vec)
        } else {
          B[l, j] = psi_vec[l]*prod(xi_vec[-l])
        }
      }
    } else {
      psi_vec <- psi_fun(w[j, ], mu, Sigma, params, ind)
      for(l in 1:d) {
        if(!is.null(z_design)) {
          B[l, j] = psi_vec[l]*prod(z_vec)
        } else {
          B[l, j] = psi_vec[l]
        }
      }
    }
  }
  return(B)
}
