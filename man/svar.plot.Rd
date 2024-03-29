% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svar.plot.R
\name{svar.plot}
\alias{svar.plot}
\alias{plot.fitsvar}
\alias{plot.svar.bin}
\alias{plot.np.svar}
\title{Plot a semivariogram object}
\usage{
\method{plot}{fitsvar}(
  x,
  y = NULL,
  legend = TRUE,
  xlab = "distance",
  ylab = "semivariance",
  xlim = NULL,
  ylim = c(0, 1.25 * max(x$fit$sv, na.rm = TRUE)),
  lwd = c(1, 2),
  add = FALSE,
  ...
)

\method{plot}{svar.bin}(
  x,
  y = NULL,
  xlab = "distance",
  ylab = "semivariance",
  xlim = NULL,
  ylim = c(0, max(x$biny, na.rm = TRUE)),
  add = FALSE,
  ...
)

\method{plot}{np.svar}(
  x,
  y = NULL,
  xlab = "distance",
  ylab = "semivariance",
  xlim = NULL,
  ylim = c(0, max(x$biny, na.rm = TRUE)),
  add = FALSE,
  ...
)
}
\arguments{
\item{x}{a variogram object. Typically an output of functions
\code{\link{np.svariso}} or \code{\link{fitsvar.sb.iso}}.}

\item{y}{ignored argument.}

\item{legend}{logical; if \code{TRUE} (default), a legend is added to the plot.}

\item{xlab}{label for the x axis (defaults to "distance").}

\item{ylab}{label for the y axis (defaults to "semivariance").}

\item{xlim}{x-limits.}

\item{ylim}{y-limits.}

\item{lwd}{line widths for points (estimates) and lines (fitted model) respectively.}

\item{add}{logical; if \code{TRUE} the semivariogram plot is just added 
to the existing plot.}

\item{...}{additional graphical parameters (see \code{\link{par}}).}
}
\value{
No return value, called for side effects (generate the plot).
}
\description{
Utilities for plotting pilot semivariograms or fitted models.

\code{plot.fitsvar} plots a fitted variogram model.

\code{plot.svar.bin} plots the binned semivariances.

\code{plot.np.svar} plots a local polynomial estimate of the semivariogram.
}
\seealso{
\code{\link{svariso}}, \code{\link{np.svariso}}, \code{\link{fitsvar.sb.iso}}.
}
