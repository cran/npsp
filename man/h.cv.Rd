\name{h.cv}
\alias{h.cv}
\alias{h.cv.bin.data}
\alias{h.cv.bin.den}
\alias{hcv.data}
\title{Cross-validation methods for bandwidth selection}
\usage{
  h.cv(bin, ...)

  \method{h.cv}{bin.data} (bin,
    objective = c("CV", "GCV", "MASE"), h.start = NULL,
    h.lower = NULL, h.upper = NULL, degree = 1,
    ncv = ifelse(objective == "GCV", 0, 1), cov.bin = NULL,
    DEalgorithm = FALSE, warn = FALSE, ...)

  \method{h.cv}{bin.den} (bin, h.start = NULL,
    h.lower = NULL, h.upper = NULL, degree = 1, ncv = 1,
    DEalgorithm = FALSE, warn = FALSE, ...)

  hcv.data(bin, objective = c("CV", "GCV", "MASE"),
    h.start = NULL, h.lower = NULL, h.upper = NULL,
    degree = 1, ncv = ifelse(objective == "GCV", 0, 1),
    cov = NULL, DEalgorithm = FALSE, warn = FALSE, ...)
}
\arguments{
  \item{bin}{object used to select a method (binned data,
  binned density or binned semivariogram).}

  \item{...}{further arguments passed to or from other
  methods (e.g. parameters of the optimization routine).}

  \item{objective}{character; optimal criterion to be used
  ("CV", "GCV" or "MASE").}

  \item{h.start}{vector; initial values for the parameters
  (diagonal elements) to be optimized over. If
  \code{DEalgorithm == FALSE} (otherwise not used),
  defaults to \code{(3 + ncv) * lag}, where \code{lag =
  bin$grid$lag}.}

  \item{h.lower}{vector; lower bounds on each parameter
  (diagonal elements) to be optimized. Defaults to
  \code{(1.5 + ncv) * bin$grid$lag}.}

  \item{h.upper}{vector; upper bounds on each parameter
  (diagonal elements) to be optimized.  Defaults to
  \code{1.5 * dim(bin) * bin$grid$lag}.}

  \item{DEalgorithm}{logical; if \code{TRUE}, the
  differential evolution optimization algorithm in package
  \pkg{DEoptim} is used.}

  \item{ncv}{integer; determines the number of cells leaved
  out in each dimension. (0 to GCV considering all the
  data, \eqn{>0} to traditional or modified
  cross-validation). See "Details" bellow.}

  \item{cov.bin}{covariance matrix of the binned data.
  Defaults to identity.}

  \item{warn}{logical; sets the handling of warning
  messages (normally due to the lack of data in some
  neighborhoods). If \code{FALSE} (the default) all
  warnings are ignored.}

  \item{cov}{covariance matrix of the data. Defaults to
  identity (uncorrelated data).}

  \item{degree}{degree of the local polynomial used.
  Defaults to 1 (local linear estimation).}
}
\value{
  Returns a list containing the following 3 components:
  \item{h}{the best (diagonal) bandwidth matrix found.}
  \item{value}{the value of the objective function
  corresponding to \code{h}.} \item{objective}{the
  criterion used.}
}
\description{
  Selects the bandwidth of a local polynomial kernel
  (regression, density or variogram) estimator using
  (standart or modified) CV, GCV or MASE criteria.
}
\details{
  Currently, only diagonal windows are supported.

  \code{h.cv} methods use binning approximations to the
  objective function values. If \code{ncv > 0}, estimates
  are computed by leaving out binning cells with indexes
  within the intervals \eqn{[x_i - ncv + 1, x_i + ncv -
  1]}, at each dimension i, where \eqn{x} denotes the index
  of the estimation position. \eqn{ncv = 1} corresponds
  with traditional cross-validation and \eqn{ncv > 1} with
  modified CV (see e.g. Chu and Marron, 1991, for the one
  dimensional case). For standard GCV, set \code{ncv = 0}
  (the full data is used). For theoretical MASE, set
  \code{y = trend.teor}, \code{cov = cov.teor} and
  \code{ncv = 0}.

  If \code{DEalgorithm == FALSE}, the \code{"L-BFGS-B"}
  method in \code{\link{optim}} is used.

  \code{hcv.data} evaluates the objective functions at the
  original data (combining a binning approximation to the
  nonparametric estimates with a linear interpolation). If
  \code{ncv > 1} (modified CV), a similar algorithm to that
  in \code{h.cv.bin.data} is used, estimates are computed
  by leaving out binning cells with indexes within the
  intervals \eqn{[x_i - ncv + 1, x_i + ncv - 1]}.
}
\examples{
bin <- binning(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(30,30))
hcv <- h.cv(bin, ncv = 2)
lp <- locpol(bin, h = hcv$h)
## Alternatively:
## lp <- locpolhcv(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(30,30), ncv = 2)

simage(lp, main = 'Smoothed magnitude')
contour(lp, add = TRUE)
with(earthquakes, points(lon, lat, pch = 20))

## Density estimation
hden <- h.cv(as.bin.den(bin))
den <- np.den(bin, h = hden$h)

plot(den, main = 'Estimated log(density)')
}
\references{
  Chu, C.K. and Marron, J.S. (1991) Comparison of Two
  Bandwidth Selectors with Dependent Errors. \emph{The
  Annals of Statistics}, \bold{19}, 1906-1918.

  Francisco-Fernandez M. and Opsomer J.D. (2005) Smoothing
  parameter selection methods for nonparametric regression
  with spatially correlated errors.  \emph{Canadian Journal
  of Statistics}, \bold{33}, 539-558.
}
\seealso{
  \code{\link{locpol}}, \code{\link{locpolhcv}},
  \code{\link{binning}}, \code{\link{np.svar}}.
}

