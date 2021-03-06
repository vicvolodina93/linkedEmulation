% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bmatrix.R
\name{Bmatrix}
\alias{Bmatrix}
\title{Function to compute \eqn{B} matrix of Linked Emulator with squared
exponential correlation function.}
\usage{
Bmatrix(w, z = NULL, mu, Sigma, z_design = NULL, params, ind = "TRUE")
}
\arguments{
\item{w}{an \eqn{m} (number of design points) by \eqn{d} (dimensions)
matrix of design points}

\item{z}{a vector of length d1 of exogenous inputs.}

\item{mu}{a vector of length \eqn{d} of predictive mean.}

\item{Sigma}{a \eqn{d} by \eqn{d} variance-covariance matrix.}

\item{z_design}{a design matrix of exogenous inputs.}

\item{params}{a list with parameter values.}

\item{ind}{logical argument indicating if we adopt independence
assumption for \eqn{W}.}
}
\value{
an \eqn{d} by \eqn{m} matrix \eqn{B}.
}
\description{
Function to compute \eqn{B} matrix of Linked Emulator with squared
exponential correlation function.
}
