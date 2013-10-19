\name{svar.bin}
\alias{iso.svar}
\alias{svar.bin}
\alias{svar.bin-class}
\alias{svar.bin.default}
\alias{svariso}
\title{Linear binning of semivariances}
\usage{
  svar.bin(x, ...)

  \method{svar.bin}{default} (x, y, maxlag = NULL,
    nlags = NULL, minlag = maxlag/nlags,
    estimator = c("classical", "modulus"), ...)

  svariso(x, y, maxlag = NULL, nlags = NULL,
    minlag = maxlag/nlags,
    estimator = c("classical", "modulus"), ...)
}
\arguments{
  \item{x}{a (data) object used to select a method.}

  \item{...}{further arguments passed to or from other
  methods.}

  \item{estimator}{character, estimator name (e.g.
  "classical"). See "Details" below.}

  \item{y}{vector of data (response variable).}

  \item{maxlag}{maximum lag. Defaults to 55\% of largest
  lag.}

  \item{nlags}{number of lags. Defaults to 101.}

  \item{minlag}{minimun lag.}
}
\value{
  Returns an S3 object of class \code{svar.bin} (extends
  \code{\link{bin.data}}), a \code{\link{data.grid}} object
  with the following 4 components: \item{biny}{array
  (dimension \code{nlags}) with the binned semivariances. }
  \item{binw}{array (dimension \code{nlags}) with the bin
  counts (weights).} \item{grid}{a
  \code{\link{grid.par}}-\code{\link{class}} object with
  the grid parameters.} \item{data}{a list with 3
  components: \itemize{ \item{\code{x} argument \code{x}.}
  \item{\code{y} argument \code{y}.} \item{\code{med}
  (weighted) mean of the (binned) semivariances.} }}
  \item{svar}{a list of 2 components: \itemize{
  \item{\code{type} character, type of estimation (e.g.
  "isotropic").} \item{\code{estimator} character,
  estimator name (e.g. "classical").} }}
}
\description{
  Creates a \code{svar.bin} (binned semivar. + grid
  parameters) object with linearly binned semivariances.
}
\details{
  Currently, only isotropic semivariogram estimation is
  supported.

  If parameter \code{nlags} is not specified is set to
  \code{101}.
}
\seealso{
  \code{\link{svarisonp}}, \code{\link{np.svar}},
  \code{\link{data.grid}}, \code{\link{binning}},
  \code{\link{locpol}}.
}

