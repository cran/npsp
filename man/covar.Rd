\name{covar}
\alias{covar}
\title{Covariance values}
\usage{
  covar(x, h, sill = x$sill)
}
\arguments{
  \item{x}{variogram model (\code{\link{svarmod}} object).}

  \item{h}{vector (isotropic case) or matrix of lags
  values.}

  \item{sill}{variance \eqn{C(0) = \sigma^2} or pseudo-sill
  (unbounded variograms).}
}
\value{
  A vector of (pseudo) covariance values \eqn{C(h_i) =
  \sigma^2 - \gamma(h_i)}.
}
\description{
  Computes covariance values given a variogram model (or
  pseudo-covariances for unbounded variograms).
}
\seealso{
  \code{\link{sv}}.
}

