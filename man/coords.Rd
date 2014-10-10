\name{coords}
\alias{coords}
\alias{coords.data.grid}
\alias{coords.grid.par}
\title{(spatial) coordinates}
\usage{
  coords(x, ...)

  \method{coords}{grid.par} (x, ...)

  \method{coords}{data.grid} (x, ...)
}
\arguments{
  \item{x}{a (spatial) object used to select a method.}

  \item{...}{further arguments passed to or from other
  methods.}
}
\value{
  A matrix of coordinates (columns correspond with
  dimensions and rows with data).
}
\description{
  Retrieves the (spatial) coordinates of the object.
}
\seealso{
  \code{\link{coordvalues}}.
}
\keyword{internal}

