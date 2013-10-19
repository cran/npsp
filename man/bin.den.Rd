\name{bin.den}
\alias{as.bin.den}
\alias{as.bin.den.bin.data}
\alias{bin.den}
\alias{bin.den-class}
\title{Linear binning for density estimation}
\usage{
  bin.den(x, nbin = NULL)

  as.bin.den(object, ...)

  \method{as.bin.den}{bin.data} (object, ...)
}
\arguments{
  \item{x}{vector or matrix of covariates (e.g. spatial
  coordinates).  Columns correspond with dimensions and
  rows with observations.}

  \item{nbin}{vector with the number of bins
  (\eqn{intervals = nbin + 1}) on each dimension.}

  \item{object}{(gridded data) used to select a method.}

  \item{...}{further arguments passed to or from other
  methods.}
}
\value{
  Returns an S3 object of \code{\link{class}}
  \code{bin.den} (extends \code{\link{data.grid}}).  A list
  with the following 3 components: \item{binw}{vector or
  array (dimension \code{nbin}) with the bin counts
  (weights).} \item{grid}{a
  \code{\link{grid.par}}-\code{\link{class}} object with
  the grid parameters.} \item{data}{a list with a component
  \code{$x} with argument \code{x}.}
}
\description{
  Creates a \code{bin.den}-\code{\link{class}} (gridded
  binned density) object with linear binning counts.
}
\details{
  If parameter \code{nbin} is not specified is set to
  \code{rep(25, ncol(x))}.
}
\seealso{
  \code{\link{data.grid}}, \code{\link{bin.data}},
  \code{\link{locpol}}.
}

