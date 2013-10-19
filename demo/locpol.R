#--------------------------------------------------------------------
#   locpol.R (npsp package demo)
#--------------------------------------------------------------------
# EXAMPLES:
#   locpol.bin  S3 class and methods
#   locpol()    S3 generic
#
#   (c) R. Fernandez-Casal         Last revision: Oct 2013
#--------------------------------------------------------------------
library(npsp)
library(fields)                               # required for drape.plot()

#--------------------------------------------------------------------
# 2D data example (regularly spaced)
#--------------------------------------------------------------------

nx <- c(40, 40)                               # ndata =  prod(nx)
x1 <- seq(-1, 1, length.out=nx[1])
x2 <- seq(-1, 1, length.out=nx[2])
x <- as.matrix(expand.grid(x1 = x1, x2 = x2)) # regularly spaced two-dimensional grid

f2d <- function(x) x[1]^2 - x[2]^2
y <- apply(x, 1, f2d) + rnorm(prod(nx), 0, 0.1)
## Equivalent to: f2d <- function(x,y) x^2 - y^2 ; y <- outer(x1, x2, f2d) + rnorm(prod(nx), 0, 0.1)

drape.plot(x1, x2, matrix(y,nrow=nx[1]), main = 'Data values', xlab = 'x1', ylab = 'x2', zlab = 'y')

#  # Binning
#  bin <- binning(x, y)
#  dim(bin)
#  str(bin)
#  coords(bin)
#  
#  coorvs <- coordvalues(bin)
#  ns <- names(coorvs)
#  drape.plot( coorvs[[1]], coorvs[[2]], bin$biny, main = 'Binning surface',
#              xlab = ns[1], ylab = ns[2], zlab = 'response means')


# Local polynomial kernel regression 
lp <- locpol(x, y, h = diag(0.3, ncol=2, nrow=2))
str(lp)
dim(lp)
str(coords(lp))

coorvs <- coordvalues(lp)
ns <- names(coorvs)                           # dimnames(lp$grid)
drape.plot( coorvs[[1]], coorvs[[2]], lp$est, main = 'locpol surface', 
            xlab = ns[1], ylab = ns[2], zlab = 'trend estimates')


# lopol from binned data
lp2 <- locpol(lp, h = diag(0.1, ncol=2, nrow=2))  # avoids binning
drape.plot( coorvs[[1]], coorvs[[2]], lp2$est, main = 'locpol surface (h=0.1)', 
            xlab = ns[1], ylab = ns[2], zlab = 'trend estimates')


#--------------------------------------------------------------------
# Regularly spaced 1D data
#--------------------------------------------------------------------
# one-dimensional data grid
nx <- 1000
x <- seq(0, 1, length.out = nx)
f1d <- function(x) sin( 2*pi*x )
y <- f1d(x) + rnorm(nx, 0, 0.5)

plot(x, y, type='l', col='darkgray')
curve(f1d, 0, 1, add = TRUE)

lp <- locpol(x, y, h = 0.05, nbin = 100)
lines(coords(lp), lp$biny, lwd = 2, lty = 2)
lines(coords(lp), lp$est, lwd = 2)

str(lp)

# cuadratic fit + derivatives

lp <- locpol(x, y, h = 0.2, nbin = 100, drv = TRUE, hat.bin = FALSE)
lines(coords(lp), lp$est, lwd = 2, col = 2)

plot(coords(lp), lp$deriv, type = "l", main = "Estimated first derivative")