% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mask.R
\name{npsp.tolerance}
\alias{npsp.tolerance}
\title{npsp Tolerances}
\usage{
npsp.tolerance(level = 2, warn = TRUE)
}
\arguments{
\item{level}{numerical,}

\item{warn}{logical; If \code{TRUE} (the default) a warning message is issued
when \code{level < 1}.}
}
\value{
Returns \code{.Machine$double.eps^(1/level)} if \code{level >= 1},
in other case \code{1 - .Machine$double.eps}.
}
\description{
Returns a (convergence, taper, approximation,...) tolerance.
Defaults to \code{.Machine$double.eps^(1/level)}, typically about \code{1e-8}.
}
\examples{
curve(npsp.tolerance, 1, 1000)
abline(h = npsp.tolerance(0, FALSE), lty = 2)
}
\seealso{
\code{\link{.Machine}}
}
