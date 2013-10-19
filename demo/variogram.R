#--------------------------------------------------------------------
#   variogram.R (npsp package demo)
#--------------------------------------------------------------------
# EXAMPLES:
#   svar.bin        S3 class and methods
#   np.svar         S3 class and methods
#   svariso()
#   svarisonp()     S3 generic
#   as.variogram()  S3 generic
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2012
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# svarisonp
#--------------------------------------------------------------------

# Cargar librerias
library(npsp)
library(geoR)

# Stationary data (data(s100) geoR)
summary(s100)
plot(s100)

# Empirical variogram
vario.geor <- variog(s100, max.dist=0.6) # geoR variog()
str(vario.geor)

# Local linear variogram 
svarnp <- svarisonp(s100$coords, s100$data, h = 0.2, maxlag = 0.6)
str(svarnp)

# Graphical comparison
oldpar <- par(mfrow=c(1,2))
plot(vario.geor, main="Empirical semivariogram")
plot(coords(svarnp), svarnp$biny, ylim = c(0,max(svarnp$biny)), main = "Nonparametric semivariogram",
            xlab = "distance", ylab = "semivariance")
lines(coords(svarnp), svarnp$est)
par(oldpar) 


#--------------------------------------------------------------------
# as.variogram
#--------------------------------------------------------------------

svar.geor <- as.variogram(svarnp)

vario.ols <- variofit(svar.geor, ini = c(1, 0.5), weights = "equal")  # OLS
vario.ols
vario.wls <- variofit(svar.geor, ini = c(1, 0.5), weights = "cressie")  # WLS
vario.wls

# Representacion grafica
plot(svar.geor, main = "Nonparametric estimates and fitted models")
lines(vario.ols, max.dist = 0.6)
lines(vario.wls, lty = 2, max.dist = 0.6)
legend(0.3, 0.3, legend = c("OLS","WLS"), lty = c(1, 2))


#--------------------------------------------------------------------
# svariso
#--------------------------------------------------------------------

# Empirical variogram
vario.geor <- variog(s100, max.dist=0.6)
varior.geor <- variog(s100, estimator.type = "modulus", max.dist=0.6)
str(vario.geor)

# Linearly binned variogram
svar <- svariso(s100$coords, s100$data, maxlag = 0.6, nlags = 13)
svarr <- svariso(s100$coords, s100$data, maxlag = 0.6, nlags = 13, estimator = "modulus")
str(svar)

# Representacion grafica
oldpar <- par(mfrow=c(2,2))
plot(vario.geor, main="Empirical semivariogram")
plot(varior.geor, main="Robust semivariogram")
plot(coords(svar), svar$biny, main = "Binned semivariogram",
            xlab = "distance", ylab = "semivariance",
            xlim = c(0,max(coords(svar))), ylim = c(0,max(svar$biny)))
plot(coords(svarr), svarr$biny, main = "Binned robust semivariogram",
            xlab = "distance", ylab = "semivariance",
            xlim = c(0,max(coords(svarr))), ylim = c(0,max(svarr$biny)))
par(oldpar) #Restaurar opciones de graficos



#--------------------------------------------------------------------
# locpol.svar.bin

svarnp2 <- locpol(svar, h = 0.2)
str(svarnp2)
plot(coords(svar), svar$biny, main = "Binned and nonparametric semivariogram",
            xlab = "distance", ylab = "semivariance")
lines(coords(svarnp2), svarnp2$est, lty = 2)