% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov.R
\name{varcov}
\alias{varcov}
\alias{varcov.isotropic}
\alias{varcov.np.svar}
\title{Covariance matrix}
\usage{
varcov(x, coords, ...)

\method{varcov}{isotropic}(
  x,
  coords,
  sill = x$sill,
  range.taper,
  discretize = nrow(coords) > 256,
  ...
)

\method{varcov}{np.svar}(x, coords, sill = max(x$est), range.taper = x$grid$max, ...)
}
\arguments{
\item{x}{variogram model (\code{\link{svarmod}} object) or semivariogram estimate.}

\item{coords}{matrix of coordinates (columns correspond with dimensions and rows with data).}

\item{...}{further arguments passed to or from other methods.}

\item{sill}{(theoretical or estimated) variance \eqn{C(0) = \sigma^2} or pseudo-sill (unbounded variograms).}

\item{range.taper}{(optional) if provided, covariances corresponding to 
distances larger than this value are set to 0.}

\item{discretize}{logical. If \code{TRUE} (default), the variogram is (previously) discretized.}
}
\value{
The covariance matrix of the data.
}
\description{
Computes the covariance matrix a corresponding to a set of spatial locations 
given a variogram model or a semivariogram estimate.
}
\seealso{
\code{\link{sv}}, \code{\link{covar}}.
}
