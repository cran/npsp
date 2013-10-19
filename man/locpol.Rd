\name{locpol}
\alias{locpol}
\alias{locpol.bin}
\alias{locpol.bin-class}
\alias{locpol.bin.data}
\alias{locpol.bin.den}
\alias{locpol.default}
\alias{locpol.svar.bin}
\alias{locpolhcv}
\title{Local polynomial estimation}
\usage{
  locpol(x, ...)

  \method{locpol}{default} (x, y, h = NULL, nbin = NULL,
    degree = 1 + as.numeric(drv), drv = FALSE,
    hat.bin = FALSE, ncv = 0, set.NA = FALSE, ...)

  \method{locpol}{bin.data} (x, h = NULL,
    degree = 1 + as.numeric(drv), drv = FALSE,
    hat.bin = FALSE, ncv = 0, ...)

  \method{locpol}{svar.bin} (x, h = NULL, degree = 1,
    drv = FALSE, hat.bin = FALSE, ncv = 0, ...)

  \method{locpol}{bin.den} (x, h = NULL,
    degree = 1 + as.numeric(drv), drv = FALSE, ncv = 0,
    ...)

  locpolhcv(x, y, nbin = NULL,
    objective = c("CV", "GCV", "MASE"),
    degree = 1 + as.numeric(drv), drv = FALSE,
    hat.bin = FALSE, set.NA = FALSE,
    ncv = ifelse(objective == "GCV", 0, 1), cov = NULL,
    ...)
}
\arguments{
  \item{x}{a (data) object used to select a method.}

  \item{...}{further arguments passed to or from other
  methods (e.g. to \code{\link{hcv.data}}).}

  \item{y}{vector of data (response variable).}

  \item{h}{(full) bandwidth matrix (controls the degree of
  smoothing).}

  \item{nbin}{vector with the number of bins on each
  dimension.}

  \item{degree}{degree of the local polynomial used.
  Defaults to 1 (local linear estimation).}

  \item{drv}{logical; if TRUE, the matrix of estimated
  first derivatives is returned.}

  \item{hat.bin}{logical; if TRUE, the hat matrix of the
  binned data is returned.}

  \item{ncv}{integer; determines the number of cells leaved
  out in each dimension. Defaults to 0 (the full data is
  used) and it is not normally changed by the user in this
  setting. See "Details" below.}

  \item{set.NA}{logical. If \code{TRUE}, sets binning cells
  without data to missing.}

  \item{objective}{character; optimal criterion to be used
  ("CV", "GCV" or "MASE").}

  \item{cov}{covariance matrix of the data. Defaults to
  identity (uncorrelated data).}
}
\value{
  Returns an S3 object of class \code{locpol.bin} (locpol +
  bin data + grid par.). A \code{\link{bin.data}} object
  with the additional (some optional) 3 components:
  \item{est}{vector or array (dimension \code{nbin}) with
  the local polynomial estimates. } \item{locpol}{a list
  with 7 components: \itemize{ \item{\code{degree} degree
  of the polinomial.} \item{\code{h} bandwidth matrix.}
  \item{\code{rm} residual mean.} \item{\code{rss} sum of
  squared residuals.} \item{\code{ncv} number of cells
  ignored in each direction.} \item{\code{hat} (if
  requested) hat matrix of the binned data.}
  \item{\code{nrl0} (if appropriate) number of cells with
  data (\code{binw > 0}) and missing estimate (\code{est ==
  NA}).} }} \item{deriv}{(if requested) matrix of first
  derivatives.}

  \code{locpol.svar.bin} returns an S3 object of class
  \code{\link{np.svar}} (locpol semivar + bin semivar +
  grid par.).

  \code{locpol.bin.den} returns an S3 object of class
  \code{\link{np.den}} (locpol den + bin den + grid par.).
}
\description{
  Estimates a multidimensional regression function (and its
  first derivatives) using local polynomial kernel
  smoothing of linearly binned data.
}
\details{
  Standard generic function with a default method
  (interface to the fortran routine \code{lp_raw}), in
  which argument \code{x} is a vector or matrix of
  covariates (e.g. spatial coordinates).

  If parameter \code{nbim} is not specified is set to
  \code{rep(25, ncol(x))}.

  A multiplicative triweight kernel is used to compute the
  weights.

  If \code{ncv > 0}, estimates are computed by leaving out
  cells with indexes within the intervals \eqn{[x_i - ncv +
  1, x_i + ncv - 1]}, at each dimension i, where \eqn{x}
  denotes the index of the estimation position. \eqn{ncv =
  1} corresponds with traditional cross-validation and
  \eqn{ncv > 1} with modified CV (see e.g. Chu and Marron,
  1991, for the one dimensional case).

  Setting \code{set.NA = TRUE} (equivalent to
  \code{biny[binw == 0] <- NA}) may be useful for plotting
  the binned averages \code{$biny} (the hat matrix should
  be handle with care).

  \code{locpolhcv} calls \code{\link{hcv.data}} to obtain
  an "optimal" bandwith (additional arguments \code{...}
  are passed to this function). Argument \code{ncv} is only
  used here at bandwith selection stage (estimation is done
  with all the data).
}
\examples{
bin <- binning(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(30,30))
hcv <- h.cv(bin)
lp <- locpol(bin, h = hcv$h)
## Equivalent to:
## lp <- locpolhcv(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(30,30))

coorvs <- coordvalues(lp)
ns <- names(coorvs)                           # dimnames(lp$grid)
image( coorvs[[1]], coorvs[[2]], lp$est, main = 'Smoothed magnitude',
            xlab = ns[1], ylab = ns[2])
contour(coorvs[[1]], coorvs[[2]], log(lp$est), add = TRUE)
with(earthquakes, points(lon, lat, pch = 20))
}
\references{
  Chu, C.K. and Marron, J.S. (1991) Comparison of Two
  Bandwidth Selectors with Dependent Errors. \emph{The
  Annals of Statistics}, \bold{19}, 1906-1918.

  Rupert D. and Wand M.P. (1994) Multivariate locally
  weighted least squares regression.  \emph{The Annals of
  Statistics}, \bold{22}, 1346-1370.
}
\seealso{
  \code{\link{binning}}, \code{\link{data.grid}},
  \code{\link{svarisonp}}, \code{\link{svar.bin}},
  \code{\link{np.den}}, \code{\link{bin.den}}.
}

