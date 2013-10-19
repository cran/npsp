\name{binning}
\alias{bin.data}
\alias{bin.data-class}
\alias{binning}
\title{Linear binning}
\usage{
  binning(x, y = NULL, nbin = NULL, set.NA = FALSE)
}
\arguments{
  \item{x}{vector or matrix of covariates (e.g. spatial
  coordinates). Columns correspond with covariates
  (coordinate dimension) and rows with data.}

  \item{y}{vector of data (response variable).}

  \item{nbin}{vector with the number of bins on each
  dimension.}

  \item{set.NA}{logical. If \code{TRUE}, sets binning cells
  without data to missing.}
}
\value{
  If \code{y != NULL}, an S3 object of \code{\link{class}}
  \code{bin.data} (gridded binned data; extends
  \code{\link{bin.den}}) is returned. A
  \code{\link{data.grid}} object with the following 4
  components: \item{biny}{vector or array (dimension
  \code{nbin}) with the bin averages. } \item{binw}{vector
  or array (dimension \code{nbin}) with the bin counts
  (weights).} \item{grid}{a
  \code{\link{grid.par}}-\code{\link{class}} object with
  the grid parameters.} \item{data}{a list with 3
  components: \itemize{ \item{\code{x} argument \code{x}.}
  \item{\code{y} argument \code{y}.} \item{\code{med}
  (weighted) mean of the (binned) data.} }}

  If \code{y == NULL}, \code{\link{bin.den}} is called and
  a \code{\link{bin.den}}-\code{\link{class}} object is
  returned.
}
\description{
  Discretizes the data into a regular grid (computes a
  binned approximation) using the multivariate linear
  binning technique described in Wand (1994).
}
\details{
  If parameter \code{nbin} is not specified is set to
  \code{rep(25, ncol(x))}.

  Setting \code{set.NA = TRUE} (equivalent to
  \code{biny[binw == 0] <- NA}) may be useful for plotting
  the binned averages \code{$biny} (the hat matrix should
  be handle with care when using \code{\link{locpol}}).
}
\examples{
bin <- binning(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(30,30), set.NA = TRUE)

coorvs <- coordvalues(bin)
ns <- names(coorvs)                           # dimnames(bin$grid)
image( coorvs[[1]], coorvs[[2]], bin$biny, main = 'Binning averages', xlab = ns[1], ylab = ns[2])
with(earthquakes, points(lon, lat, pch = 20))
}
\references{
  Wand M.P. (1994) Fast Computation of Multivariate Kernel
  Estimators.  \emph{Journal of Computational and Graphical
  Statistics}, \bold{3}, 433-445.
}
\seealso{
  \code{\link{data.grid}}, \code{\link{locpol}},
  \code{\link{bin.den}}.
}

