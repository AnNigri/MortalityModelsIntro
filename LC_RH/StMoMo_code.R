##############################################################
## Il presente codice prende spunto direttamente            ##
## dalle lezioni tenute dalla prof.ssa Susanna Levantesi    ##
##                                                          ##
##############################################################


library(StMoMo)

# -----------------  Dati ITALIA  --------------------#
age=c(60:90)
ex_italia =read.csv("Ex_ITALIA.csv",header=TRUE,sep=";")[,-1][,-43][age,] 
qx_italia=read.csv("qx.csv",header=TRUE,sep=";")[age,]
mx_italia <- -log(1-qx_italia)
dx_italia <- mx_italia*ex_italia
yv=c(1974:2015)
wxt <- genWeightMat(ages = age, years = yv, clip = 3)

#--------------------- mortality reduction factors --------------#
qx_RF <- matrix(0, nrow=length(age), ncol=(length(yv)-1))  # reduction factors
for(i in 1:(length(yv)-1)){
  qx_RF[,i] <- qx_italia[,i+1]/qx_italia[,i]
}
#----------------	Heatmap RF	---------------------#
require(grDevices) # for colours
require(graphics)
library(lattice)
Y <- as.numeric(yv[2:length(yv)])
A <- as.numeric(age)
par(mar=c(5,5,5,5))  
filled.contour(A, Y, qx_RF, color.palette = colorRampPalette(c("white", "yellow", "orange", "red", "blue", "darkblue", "black"), space = "Lab"),
               plot.title = title(main = "Mortality reduction factors", xlab = "età", ylab = "anni"),
               plot.axes = {	
                 axis(1, seq(A[1], A[length(A)], by = 5))
                 axis(2, seq(Y[1], Y[length(Y)], by = 5)) 
               },
               key.title = title(main=""),
               key.axes = axis(4, seq(min(qx_RF), max(qx_RF), by = ((max(qx_RF)-min(qx_RF))/20))))

# -----------------  Modello Lee-Carter  --------------------#
LCfit <- fit(lc(), Dxt = dx_italia, Ext = ex_italia, ages = age, years = yv)
plot(LCfit)
names(LCfit)
BIC(LCfit); AIC(LCfit)
LCres <- residuals(LCfit)
plot(LCres, type="scatter")
plot(LCres, type = "signplot")
plot(LCres, type = "colourmap", reslim = c(-3.5, 3.5))

# The period indexes are forecasted using a Multivariate Random Walk with Drift (MRWD). 
# The cohort index is forecasted using an ARIMA(p; d; q) . By default: ARIMA(1; 1; 0).
# h= number of years ahead to forecast
# level=  confidence level for prediction intervals of the period and cohort indices.
# jumpchoice option to select the jump-off rates, i.e. the rates from the final year of observation, to use in projections of mortality rates.  
# "fit"(default) uses the fitted rates and  "actual" uses the actual rates from the final year.
# kt.lookback optional argument to specify the look-back window to use in the estimation of the MRWD for period indexes. 
# By default all the estimated values are used in estimating the MRWD. If  kt.lookback is provided then the last  kt.lookback years of kt are used.
LCfor <- forecast(LCfit, h=50, level=95)	# proiezione 
names(LCfor)
LCfor$rates
plot(LCfor)


# Simulate future sample paths from a Bootstrapped Stochastic Mortality Model. The period indexes kt are modelled using a Multivariate RandomWalk with Drift. 
# The cohort index is modelled using an ARIMA(p; d; q) . By default:  ARIMA(1; 1; 0).
LCsim <- simulate(LCfit, nsim = 1000, h=50) # simulazione di scenari di moratlità
names(LCsim)
LCsim$rates

par(mfrow=c(1, 3))
plot(LCfit$years, (LCfit$Dxt / LCfit$Ext)["65", ], xlim = range(LCfit$years, LCsim$years), ylim = range((LCfit$Dxt / LCfit$Ext)["65", ], LCsim$rates["65", , ]), type = "l", xlab = "year t", ylab = "q(x,t)", main = "Lee-Carter: tasso di mortalità simulato (x=65)")
matlines(LCsim$years, LCsim$rates["65", , ], type = "l", lty = 1)
plot(LCfit$years, (LCfit$Dxt / LCfit$Ext)["75", ], xlim = range(LCfit$years, LCsim$years), ylim = range((LCfit$Dxt / LCfit$Ext)["75", ], LCsim$rates["75", , ]), type = "l", xlab = "year t", ylab = "q(x,t)", main = "Lee-Carter: tasso di mortalità simulato (x=75)")
matlines(LCsim$years, LCsim$rates["75", , ], type = "l", lty = 1)
plot(LCfit$years, (LCfit$Dxt / LCfit$Ext)["85", ], xlim = range(LCfit$years, LCsim$years), ylim = range((LCfit$Dxt / LCfit$Ext)["85", ], LCsim$rates["85", , ]), type = "l", xlab = "year t", ylab = "q(x,t)", main = "Lee-Carter: tasso di mortalità simulato (x=85)")
matlines(LCsim$years, LCsim$rates["85", , ], type = "l", lty = 1)


# -----------------  Modello Renshaw-Haberman  --------------------#
RHfit <- fit(rh(), Dxt = dx_italia, Ext = ex_italia, ages = age, years = yv)
RHres <- residuals(RHfit)
plot(RHres, type = "colourmap", reslim = c(-3.5, 3.5))
RHfor <- forecast(RHfit)
RHsim <- simulate(RHfit, nsim = 1000, h=50)

plot(LCfor, only.kt = TRUE)
plot(RHfor, only.kt = TRUE)
plot(RHfor, only.gc = TRUE)

LCrates995 <- apply(LCsim$rates,c(1,2),quantile,prob=0.995)
RHrates995 <- apply(RHsim$rates,c(1,2),quantile,prob=0.995)
par(mfrow=c(1,2))
matplot(LCrates995,type="l", main="LC")
matplot(RHrates995,type="l", main="RH")

