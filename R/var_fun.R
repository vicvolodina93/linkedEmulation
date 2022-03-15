#' Function to compute analytical variance \eqn{\sigma_L^2} of Linked Emulator
#' under squared exponential correlation function.
#'
#' @param w w an \eqn{m} (number of design points) by \eqn{d} (dimensions)
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
#' @return a value of the analytical variance.
#'
#' @examples
var_fun <- function(w, z = NULL, y, mu, Sigma, z_design = NULL, params,
                    FastParts, ind = "TRUE") {

  m <- dim(w)[1]
  d <- dim(w)[2]

  b <-  bvec(w, z, mu, Sigma, z_design, params, ind)
  B <- Bmatrix(w, z, mu, Sigma, z_design, params, ind)
  J <- Jmatrix(w, z, mu, Sigma, z_design, params, ind)

  AA <- FastParts$R%*%FastParts$diff
  part1 <- t(AA)%*%(J-b%*%t(b))%*%AA
  part2 <- 2*params$beta[1:d]%*%(B-mu%*%t(b))%*%AA
  part3 <- sum(diag(params$beta[1:d]%*%t(params$beta[1:d])%*%Sigma))


  if(is.null(z)) {
    tilde_H <- w
    G <- matrix(mu, nrow = d)
    K <- t(B)
    P <- Sigma
  } else {
    tilde_H <- cbind(w, rep(1, m), z_design)
    G <- matrix(c(mu, 1, z), nrow = length(c(mu, 1, z)))
    K <- cbind(t(B), b%*%c(1, z))
    #P <- diag(c(diag(Sigma), rep(0, 1+length(z))),
    #          nrow = dim(Sigma)[1]+1+length(z), ncol = dim(Sigma)[1]+1+length(z))
    P <- diag(0, nrow = dim(Sigma)[1]+1+length(z),
              ncol = dim(Sigma)[1]+1+length(z))
    P[1:dim(Sigma)[1], 1:dim(Sigma)[2]] = Sigma

  }


  Q <- FastParts$R%*%tilde_H%*%solve(t(tilde_H)%*%FastParts$R%*%tilde_H)%*%t(tilde_H)%*%FastParts$R -FastParts$R

  C <- solve(t(tilde_H)%*%FastParts$R%*%tilde_H)

  part4 <- params$sigma^2*(1+params$nugget+sum(diag(Q%*%J)) + t(G)%*%C%*%G +sum(diag(C%*%P - 2*C%*%t(tilde_H)%*%FastParts$R%*%K)))
  res <- sum(c(part1, part2, part3, part4))
  if(res < 0) { warning( "Negative variance" ) }
  return(res)
}
