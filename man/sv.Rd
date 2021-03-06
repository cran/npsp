% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svarmod.R
\name{sv}
\alias{sv}
\alias{sv.default}
\alias{sv.svarmod}
\alias{sv.svar.grid}
\alias{sv.sb.iso}
\title{Evaluate a semivariogram model}
\usage{
sv(x, h, ...)

\method{sv}{default}(x, h, ...)

\method{sv}{svarmod}(x, h, ...)

\method{sv}{svar.grid}(x, h, ...)

\method{sv}{sb.iso}(x, h, discretize = FALSE, ...)
}
\arguments{
\item{x}{variogram model (\code{\link{svarmod}} object).}

\item{h}{vector (isotropic case) or matrix of lag values.}

\item{...}{further arguments passed to or from other methods.}

\item{discretize}{logical. If \code{TRUE} the variogram is previously discretized.}
}
\value{
A vector of semivariance values \eqn{\gamma(h_i)}.
}
\description{
Evaluates an \code{svarmod} object \code{x} at lags \code{h} (S3 generic function).
}
\seealso{
\code{\link{covar}}
}
