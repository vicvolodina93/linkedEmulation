#' Function to compute \eqn{J} matrix of Linked emulator with
#' squared exponential correlation function.
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
#' @return an \eqn{m} by \eqn{m} matrix \eqn{J}.
#' @export
#'
#' @examples
Jmatrix <- function(w, z = NULL, mu, Sigma, z_design = NULL, params,
                    ind = "TRUE") {

  m <- dim(w)[1] # number of training (design) points
  d <- dim(w)[2] # number of dimensions
  d1 <- dim(z_design)[2]

  J <- matrix(nrow = m, ncol = m)

  for(i in 1:m) {
    for(j in i:m) {
      zeta_vec <- zeta_fun(w[i, ], w[j, ], mu, Sigma, params, ind)
      if(!is.null(z_design)) {
        b1 <- sapply(1:d1, function(x)
          exp(-(z_design[i, x] - z[x])^2/(params$delta_par[d+x])^2))
        b2 <- sapply(1:d1, function(x)
          exp(-(z_design[j, x] - z[x])^2/(params$delta_par[d+x])^2))
        J[i, j] = prod(b1)*prod(b2)*prod(zeta_vec)
      } else {
        J[i, j] = prod(zeta_vec)
      }
      J[j, i] <- J[i, j]
    }
  }
  return(J)
}
