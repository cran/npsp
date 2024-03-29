% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interp.R
\name{interp}
\alias{interp}
\alias{interp.grid.par}
\alias{interp.data.grid}
\alias{predict.locpol.bin}
\alias{predict.np.den}
\title{Fast linear interpolation of a regular grid}
\usage{
interp(object, ...)

\method{interp}{grid.par}(object, data, newx, ...)

\method{interp}{data.grid}(object, data.ind = 1, newx, ...)

\method{predict}{locpol.bin}(object, newx = NULL, hat.data = FALSE, ...)

\method{predict}{np.den}(object, newx = NULL, ...)
}
\arguments{
\item{object}{(gridded data) object used to select a method.}

\item{...}{further arguments passed to or from other methods.}

\item{data}{vector or array of data values.}

\item{newx}{vector or matrix with the (irregular) locations to interpolate. 
Columns correspond with dimensions and rows with data.}

\item{data.ind}{integer (or character) with the index (or name) of the data component.}

\item{hat.data}{logical; if \code{TRUE} (and possible), the hat matrix corresponding 
to the (original) data is returned.}
}
\value{
A list with two components:
\item{x}{interpolation locations.}
\item{y}{interpolated values.}

If \code{newx == NULL}, \code{predict.locpol.bin} returns the estimates 
(and optionally the hat matrix) corresponding to the data
(otherwise \code{interp.data.grid} is called).
}
\description{
Computes a linear interpolation of multidimensional regularly gridded data.
}
\details{
\code{interp} methods are interfaces to the fortran routine \code{interp_data_grid}  
   (in \code{grid_module.f90}).

\code{predict.locpol.bin} is an interface to the fortran routine 
\code{predict_lp} (in \code{lp_module.f90}).
}
\note{
Linear extrapolation is performed from the end nodes of the grid.

WARNING: May fail with missing values (especially if \code{object$locpol$ncv > 0}).
}
\seealso{
\code{\link[fields]{interp.surface}}.
}
