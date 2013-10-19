\name{np.svar}
\alias{iso.np.svar}
\alias{np.svar}
\alias{np.svar-class}
\alias{np.svar.default}
\alias{np.svar.svar.bin}
\alias{svarisohcv}
\alias{svarisonp}
\title{Local polynomial estimation of the semivariogram}
\usage{
  np.svar(x, ...)

  \method{np.svar}{default} (x, y, h = NULL, maxlag = NULL,
    nlags = NULL, minlag = maxlag/nlags, degree = 1,
    drv = FALSE, hat.bin = FALSE, ncv = 0, ...)

  \method{np.svar}{svar.bin} (x, h = NULL, degree = 1,
    drv = FALSE, hat.bin = FALSE, ncv = 0, ...)

  svarisonp(x, y, h = NULL, maxlag = NULL, nlags = NULL,
    minlag = maxlag/nlags, degree = 1, drv = FALSE,
    hat.bin = FALSE, ncv = 0, ...)

  svarisohcv(x, y, maxlag = NULL, nlags = NULL,
    minlag = maxlag/nlags,
    objective = c("CV", "GCV", "MASE"),
    ncv = ifelse(objective == "GCV", 0, 1), cov.bin = NULL,
    ...)
}
\arguments{
  \item{x}{a (data) object used to select a method.}

  \item{...}{further arguments passed to or from other
  methods.}

  \item{y}{vector of data (response variable).}

  \item{maxlag}{maximum lag. Defaults to 55\% of largest
  lag.}

  \item{nlags}{number of lags. Defaults to 101.}

  \item{minlag}{minimun lag.}

  \item{hat.bin}{logical; if TRUE, the hat matrix of the
  binned semivariances is returned.}

  \item{cov.bin}{covariance matrix of the binned
  semivariances. Defaults to identity.}

  \item{h}{(full) bandwidth matrix (controls the degree of
  smoothing).}

  \item{degree}{degree of the local polynomial used.
  Defaults to 1 (local linear estimation).}

  \item{drv}{logical; if TRUE, the matrix of estimated
  first derivatives is returned.}

  \item{ncv}{integer; determines the number of cells leaved
  out in each dimension. Defaults to 0 (the full data is
  used) and it is not normally changed by the user in this
  setting. See "Details" below.}

  \item{objective}{character; optimal criterion to be used
  ("CV", "GCV" or "MASE").}
}
\value{
  Returns an S3 object of class \code{np.svar} (locpol svar
  + binned svar + grid par.), extends
  \code{\link{svar.bin}}, with the additional (some
  optional) 3 components: \item{est}{vector or array with
  the local polynomial semivariogram estimates. }
  \item{locpol}{a list of 6 components: \itemize{
  \item{\code{degree} degree of the local polinomial used.}
  \item{\code{h} smoothing matrix.} \item{\code{rm} mean of
  residual semivariances.} \item{\code{rss} sum of squared
  residual semivariances.} \item{\code{ncv} number of cells
  ignored in each direction.} \item{\code{hat} (if
  requested) hat matrix of the binned semivariances.}
  \item{\code{nrl0} (if appropriate) number of cells with
  \code{binw > 0} and \code{est == NA}.} }}
  \item{deriv}{(if requested) matrix of estimated first
  semivariogram derivatives.}
}
\description{
  Estimates a multidimensional semivariogram (and its first
  derivatives) using local polynomial kernel smoothing of
  linearly binned semivariances.
}
\details{
  Currently, only isotropic semivariogram estimation is
  supported.

  If parameter \code{nlags} is not specified is set to
  \code{101}.

  A multiplicative triweight kernel is used to compute the
  weights.
}
\references{
  Fernandez Casal R., Gonzalez Manteiga W. and Febrero
  Bande M. (2003) Space-time dependency modeling using
  general classes of flexible stationary variogram models,
  \emph{J. Geophys. Res.}, \bold{108}, 8779,
  doi:10.1029/2002JD002909.

  Garcia-Soidan P.H., Gonzalez-Manteiga W. and
  Febrero-Bande M. (2003) Local linear regression
  estimation of the variogram, \emph{Stat. Prob. Lett.},
  \bold{64}, 169-179.
}
\seealso{
  \code{\link{svar.bin}}, \code{\link{data.grid}},
  \code{\link{locpol}}.
}

