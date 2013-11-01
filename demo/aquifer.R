#--------------------------------------------------------------------
# Examples package "npsp", release = 0.2-3
# Wolfcamp aquifer data
#--------------------------------------------------------------------

library(npsp)
library(fields)

# ?aquifer
str(aquifer)
summary(aquifer)
with(aquifer, plot(lon, lat))

#------------------------------------------
# Trend estimation

lp <- locpol(aquifer[,1:2], aquifer$head, nbin = c(50,50), h = diag(100, 2))
# Warning: missing values were generated

# Delete grid nodes far from data
coorx <- coords(lp)
lp$est[(coorx[,1] + 150)/100 < (coorx[,2] - 100)/80] <- NA

# Plot
coorvs <- coordvalues(lp)
ns <- names(coorvs)
drape.plot( coorvs[[1]], coorvs[[2]], lp$est, main = 'Trend estimates',
            xlab = ns[1], ylab = ns[2], zlab = 'piezometric-head levels', 
            ticktype = 'detailed', theta = 120)

#------------------------------------------
# Variogram estimation

lp.resid <- lp$data$y - predict(lp)
# maxlag <- 0.55*sqrt(sum(diff(apply(aquifer[,1:2], 2, range))^2))

esvar <- svarisonp(aquifer[,1:2], lp.resid, maxlag = 150, nlags = 60, h = 60)
svm <- fitsvar.sb.iso(esvar)  # dk = 2

with(svm$fit, plot(u, sv, main = "Nonparametric semivariogram and fitted model", 
    xlab = "distance", ylab = "semivariance", ylim = c(0, 25000)))
with(svm$fit, lines(u, fitted.sv))

# Note that the direct use of the residuals introduces a bias in the estimation 
# of the variogram. This bias is usually negative and higher at large lags 
# (e.g. Cressie, 1993, section 3.4.3).
# A correction for this bias is proposed in:
# Fernandez-Casal R. and Francisco-Fernandez M. (2013) 
# Nonparametric bias-corrected variogram estimation under non-constant trend. 
# Stoch. Environ. Res. Ris. Assess. (Accepted for publication).


#------------------------------------------
# Kriging

library(gstat)

spdf <- SpatialPointsDataFrame(aquifer[,1:2], data.frame(y = aquifer$head, r = lp.resid))
newdata <- SpatialPoints(coords(lp))
krig <- krige(r ~ 1, locations = spdf, newdata = newdata, model = as.vgm(svm))

pred <- lp$est + krig@data$var1.pred

drape.plot( coorvs[[1]], coorvs[[2]], pred, main = 'Kriging predictions',
            xlab = ns[1], ylab = ns[2], zlab = 'piezometric-head', 
            ticktype = 'detailed', theta = 120)

# kvar <- matrix(krig@data$var1.var, nrow=dim(lp)[1])
# image.plot(coorvs[[1]], coorvs[[2]], kvar, main = 'Kriging variances',
#            xlab = ns[1], ylab = ns[2])

