#--------------------------------------------------------------------
#   binning.R (npsp package demo)
#--------------------------------------------------------------------
# EXAMPLES:
#   binning   S3 generic
#   bin.data  S3 class and methods
#   locpol    S3 generic  (locpol.bin.data)
#
#   (c) R. Fernandez-Casal         Last revision: Oct 2013
#--------------------------------------------------------------------
library(npsp)
library(fields)                               # required for drape.plot()


# 2D data example (regularly spaced)
set.seed(1)
nx <- c(40, 40)                               # ndata =  prod(nx)
x1 <- seq(-1, 1, length.out=nx[1])
x2 <- seq(-1, 1, length.out=nx[2])
x <- as.matrix(expand.grid(x1 = x1, x2 = x2)) # regularly spaced two-dimensional grid

f2d <- function(x) x[1]^2 - x[2]^2
y <- apply(x, 1, f2d) + rnorm(prod(nx), 0, 0.1)
## Equivalently: f2d <- function(x,y) x^2 - y^2 ; y <- outer(x1, x2, f2d)

drape.plot( x1, x2, matrix(y,nrow=nx[1]), main = 'Data values',
            xlab = 'x1', ylab = 'x2', zlab = 'y')


# Binning
bin <- binning(x, y)
dim(bin)
str(bin)
str(coords(bin))

coorvs <- coordvalues(bin)
ns <- names(coorvs)                           # dimnames(bin$grid)
drape.plot( coorvs[[1]], coorvs[[2]], bin$biny, main = 'Binning surface',
            xlab = ns[1], ylab = ns[2], zlab = 'response means')


# Local polynomial kernel regression
# (currently local linear...)
lp <- locpol(bin, h = diag(0.3, ncol=2, nrow=2), hat.bin = TRUE)
## Equivalent to:  lp <- locpol(x, y, h = diag(0.3, ncol=2, nrow=2)) # not equal...
str(lp)

drape.plot( coorvs[[1]], coorvs[[2]], lp$est, main = 'locpol surface',
            xlab = ns[1], ylab = ns[2], zlab = 'trend estimates')


# Estimation with (binning) hat matrix
if (!is.null(lp$locpol$hat)) {
    y2 <- apply(x, 1, f2d) + rnorm(prod(nx), 0, 0.1)  # new y-data
    bin <- binning(x, y2)

  est2 <- lp$locpol$hat %*% as.numeric(bin$biny)
  drape.plot( coorvs[[1]], coorvs[[2]], matrix(est2, nrow=dim(lp)[1]), 
              main = 'locpol surface (binning hat matrix)', 
              xlab = ns[1], ylab = ns[2], zlab = 'trend estimates')
} else cat("'locpol.bin' object must be created with parameter 'hat.bin = TRUE'")
