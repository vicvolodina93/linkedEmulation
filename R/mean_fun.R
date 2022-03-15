#' Function to compute the analytical mean \eqn{mu_L} of
#' the Linked Emulator under squared exponential correlation function
#'
#' @param w an \eqn{m} (number of design points) by \eqn{d} (dimensions)
#' matrix of design points
#' @param z a vector of length d1 of exogenous inputs.
#' @param y a vector length \eqn{m} (number of design points) responses.
#' @param mu a vector of length \eqn{d} of predictive mean.
#' @param Sigma a \eqn{d} by \eqn{d} variance-covariance matrix.
#' @param z_design a design matrix of exogenous inputs.
#' @param params a list with parameter values.
#' @param FastParts a list of precomputed objects.
#' @param ind logical argument indicating if we adopt independence
#' assumption for \eqn{W}.
#'
#' @return a value of the analytical mean.
#' @export
#'
#' @examples
mean_fun <- function(w, z = NULL, y, mu, Sigma, z_design = NULL,
                     params, FastParts, ind = "TRUE") {


  m <- dim(w)[1] # number of training points
  d <- dim(w)[2] # number of dimensions
  d_t <- length(params$beta)

  if(ind == "TRUE") {
    b <- bvec(w, z, mu, Sigma, z_design, params, ind)
  } else {
    Lambda.mat <- diag(params$delta_par^2/2, nrow = d, ncol = d)
    LS.inv <- solve(Lambda.mat + Sigma)
    det_comp <- 1/sqrt(det((Lambda.mat+Sigma)%*%solve(Lambda.mat)))

    b <- bvec(w, z, mu, Sigma, z_design, params,
              ind, det_comp, LS.inv)
  }
  txtA <- backsolve(FastParts$QA, b, transpose = TRUE)
  if(is.null(z_design)) {
    res <- t(mu)%*%params$beta +crossprod(txtA, FastParts$Ldiff)
  } else {
    res <- t(mu)%*%params$beta[1:d] +
      t(as.matrix(c(1, z), nrow=1))%*%matrix(params$beta[(d+1):d_t],
                                      ncol = 1) +
      crossprod(txtA, FastParts$Ldiff)
  }
  return(res)
}
