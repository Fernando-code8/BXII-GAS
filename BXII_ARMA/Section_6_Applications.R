################################################################################
# PAPER: Quantile-based dynamic modeling of asymmetric data: a novel Burr XII
#        approach for positive continuous random variables
# SECTION: 6. Applications
# GOAL: Application of the proposed model to real data
# AUTHORS: Fernando Jose Monteiro de Araujo, Renata Rojas Guerra and 
#          Fernando Arturo Pena-Ramirez
# LAST UPDATE: May 25, 2024
################################################################################

rm(list=ls()) 

# Packages
library(readr)
library(forecast)
library(tseries)
library(lmtest)
library(extraDistr)
library(SuppDists)
library(Rcpp)
library(xtable)
library(actuar)
library(timetk)
library(ggplot2)
library(numDeriv)
library(readxl)
library(rlang)
library(PTSR)
source("bxiiarma.fit.r")
source("bxiiarmaCOV.fit.r")
source("predict_bxii.R")

################################################################################
#########           Banco Bradesco S.A. (BBD)- trading volume         ##########
################################################################################

# Data set
dados<-read_csv("BBD3-Time series model for finance data.csv")
attach(dados)

# View(dados)
dados1<-dados$Volume

h1<-30 ## Steps to be predicted
h2<-h1-1

dados2<-dados1[1:(length(dados1)-h1)]
dadosh<-dados1[(length(dados1)-h2):length(dados1)]


div<-100000000
dadosh<-dadosh/div
yy<-dados2/div

inds <- seq(as.Date("2022-07-06"), as.Date("2023-02-10"), by = "day")

## Create a time series object
dados.ts <- ts(yy,     # random data
              start = c(2022, as.numeric(format(inds[1], "%j"))),
              frequency = 365)

y <- dados.ts

# Table of descriptive measures (Table 2)
a<-c(summary(dados.ts,na.rm=T),
     var(dados.ts))
round(a,4)
a<-as.matrix(t(a))
xtable(a,digits = 4)

# Initial graphics
w1<-6
h11<-4

# Serie Plot (Figure 3)
postscript(file = "plotserie.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
plot(dados.ts,ylab="Volume")
dev.off()  

w1<-5
h11<-4

# ACF Plot (Figure 4.a)
postscript(file = "ACF.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
acf(dados.ts, main="")
dev.off()  

# PACF Plot (Figure 4.b)
postscript(file = "PACF.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
pacf(dados.ts, main="")
dev.off()  
  
# Stationarity of the series
y <- dados.ts
pp.test(y, alternative ="stationary") 
kpss.test(dados.ts)

#-------------------------------------------------------------------
# BXII-ARMA(1,4) model                               
#-------------------------------------------------------------------

# Autoregressive and moving averages

q1<-c(1:3,6)
p1<-c(7)

BXIIARMA.mod <- bxiiarma.fit(y,ma=q1,ar=p1,h=30,mod=3,resid = 3) 

### choosing the best model based on AIC BIC and HQ criterion
BXIIARMA.mod2 <- bxiiarma.fit(y,ma=c(1),ar=c(1),h=30,mod=3,resid = 3) 
BXIIARMA.mod3 <- bxiiarma.fit(y,ma=c(1),ar=c(2),h=30,mod=3,resid = 3)
BXIIARMA.mod4 <- bxiiarma.fit(y,ma=c(1),ar=c(3),h=30,mod=3,resid = 3)
BXIIARMA.mod5 <- bxiiarma.fit(y,ma=c(1),ar=c(4),h=30,mod=3,resid = 3)
BXIIARMA.mod6 <- bxiiarma.fit(y,ma=c(2),ar=c(1),h=30,mod=3,resid = 3)
BXIIARMA.mod7 <- bxiiarma.fit(y,ma=c(2),ar=c(2),h=30,mod=3,resid = 3)
BXIIARMA.mod8 <- bxiiarma.fit(y,ma=c(2),ar=c(3),h=30,mod=3,resid = 3)
BXIIARMA.mod9 <- bxiiarma.fit(y,ma=c(2),ar=c(4),h=30,mod=3,resid = 3)
BXIIARMA.mod10 <- bxiiarma.fit(y,ma=c(3),ar=c(1),h=30,mod=3,resid = 3)
BXIIARMA.mod11 <- bxiiarma.fit(y,ma=c(3),ar=c(2),h=30,mod=3,resid = 3)
BXIIARMA.mod12 <- bxiiarma.fit(y,ma=c(3),ar=c(3),h=30,mod=3,resid = 3)
BXIIARMA.mod13 <- bxiiarma.fit(y,ma=c(3),ar=c(4),h=30,mod=3,resid = 3)
BXIIARMA.mod14 <- bxiiarma.fit(y,ma=c(4),ar=c(1),h=30,mod=3,resid = 3)
BXIIARMA.mod15 <- bxiiarma.fit(y,ma=c(4),ar=c(2),h=30,mod=3,resid = 3)
BXIIARMA.mod16 <- bxiiarma.fit(y,ma=c(4),ar=c(3),h=30,mod=3,resid = 3)
BXIIARMA.mod17 <- bxiiarma.fit(y,ma=c(4),ar=c(4),h=30,mod=3,resid = 3)
BXIIARMA.mod18 <- bxiiarma.fit(y,ma=c(1),h=30,mod=2,resid = 3)
BXIIARMA.mod19 <- bxiiarma.fit(y,ma=c(2),h=30,mod=2,resid = 3)
BXIIARMA.mod20 <- bxiiarma.fit(y,ma=c(3),h=30,mod=2,resid = 3)
BXIIARMA.mod21 <- bxiiarma.fit(y,ma=c(4),h=30,mod=2,resid = 3)
BXIIARMA.mod22 <- bxiiarma.fit(y,ar=c(1),h=30,mod=1,resid = 3)
BXIIARMA.mod23 <- bxiiarma.fit(y,ar=c(1),h=30,mod=1,resid = 3)
BXIIARMA.mod24 <- bxiiarma.fit(y,ar=c(1),h=30,mod=1,resid = 3)
BXIIARMA.mod25 <- bxiiarma.fit(y,ar=c(1),h=30,mod=1,resid = 3)

AIC<-c(BXIIARMA.mod$aic,BXIIARMA.mod2$aic,BXIIARMA.mod3$aic,BXIIARMA.mod4$aic,BXIIARMA.mod5$aic,
       BXIIARMA.mod6$aic,BXIIARMA.mod7$aic,BXIIARMA.mod8$aic,BXIIARMA.mod9$aic,BXIIARMA.mod10$aic,
       BXIIARMA.mod11$aic,BXIIARMA.mod12$aic,BXIIARMA.mod13$aic,BXIIARMA.mod14$aic,BXIIARMA.mod15$aic,
       BXIIARMA.mod16$aic,BXIIARMA.mod17$aic,BXIIARMA.mod18$aic,BXIIARMA.mod19$aic,BXIIARMA.mod20$aic,
       BXIIARMA.mod21$aic,BXIIARMA.mod22$aic,BXIIARMA.mod23$aic,BXIIARMA.mod24$aic,BXIIARMA.mod25$aic)

BIC<-c(BXIIARMA.mod$bic,BXIIARMA.mod2$bic,BXIIARMA.mod3$bic,BXIIARMA.mod4$bic,BXIIARMA.mod5$bic,
       BXIIARMA.mod6$bic,BXIIARMA.mod7$bic,BXIIARMA.mod8$bic,BXIIARMA.mod9$bic,BXIIARMA.mod10$bic,
       BXIIARMA.mod11$bic,BXIIARMA.mod12$bic,BXIIARMA.mod13$bic,BXIIARMA.mod14$bic,BXIIARMA.mod15$bic,
       BXIIARMA.mod16$bic,BXIIARMA.mod17$bic,BXIIARMA.mod18$bic,BXIIARMA.mod19$bic,BXIIARMA.mod20$bic,
       BXIIARMA.mod21$bic,BXIIARMA.mod22$bic,BXIIARMA.mod23$bic,BXIIARMA.mod24$bic,BXIIARMA.mod25$bic)

HQ<-c(BXIIARMA.mod$hq,BXIIARMA.mod2$hq,BXIIARMA.mod3$hq,BXIIARMA.mod4$hq,BXIIARMA.mod5$hq,
      BXIIARMA.mod6$hq,BXIIARMA.mod7$hq,BXIIARMA.mod8$hq,BXIIARMA.mod9$hq,BXIIARMA.mod10$hq,
      BXIIARMA.mod11$hq,BXIIARMA.mod12$hq,BXIIARMA.mod13$hq,BXIIARMA.mod14$hq,BXIIARMA.mod15$hq,
      BXIIARMA.mod16$hq,BXIIARMA.mod17$hq,BXIIARMA.mod18$hq,BXIIARMA.mod19$hq,BXIIARMA.mod20$hq,
      BXIIARMA.mod21$hq,BXIIARMA.mod22$hq,BXIIARMA.mod23$hq,BXIIARMA.mod24$hq,BXIIARMA.mod25$hq)

# best model
BXIIARMA.mod

bb<-BXIIARMA.mod$fitted
a<-BXIIARMA.mod$forecast
datatest<-dadosh
n<-length(dados2)

# forecast estimation using a one-step-ahead approach
bxii_forecast <- predict_bxii(BXIIARMA.mod)
pred<-bxii_forecast

# Residuals
residq<-BXIIARMA.mod$resid3

# Ljung-box test
Box.test(residq,lag = 10, type =  "Ljung-Box", fitdf = 0)

# accuracy measures BXII_ARMA
BXII_prev =(dadosh-pred)
BXII_prev_1=(dadosh[2:h1]-dadosh[1:h2])

mse_BXII_ARMA = (sum(BXII_prev^2))/length(BXII_prev)
mape_BXII_ARMA = sum( abs(BXII_prev)/abs(dadosh) )/ length(BXII_prev)
mase_BXII_ARMA = sum( abs(BXII_prev)/(sum(abs(BXII_prev_1))*(1/(length(BXII_prev)-1))))/ length(BXII_prev)

c<-c(BXIIARMA.mod$aic,BXIIARMA.mod$bic,
     BXIIARMA.mod$hq,
     mse_BXII_ARMA,mape_BXII_ARMA,mase_BXII_ARMA)

#-------------------------------------------------------------------
# ARMA(3,1) model
#-------------------------------------------------------------------
ARMA.mod<-auto.arima(y)
result<-summary(ARMA.mod)
coeftest(ARMA.mod)

# forecast
ARMA.mod.prev<-forecast(ARMA.mod,h=h1)
ARMA.mod.prev<-ARMA.mod.prev$mean

# accuracy measures ARMA
ARMA_prev =(dadosh-ARMA.mod.prev)
ARMA_prev_1=(dadosh[2:h1]-dadosh[1:h2])

mse_ARMA = (sum(ARMA_prev^2))/length(ARMA_prev)
mape_ARMA = sum( abs(ARMA_prev)/abs(dadosh) )/ length(ARMA_prev)
mase_ARMA = sum( abs(ARMA_prev)/(sum(abs(ARMA_prev_1))*(1/(length(ARMA_prev)-1))))/ length(ARMA_prev)

HQC.ARMA<--2*ARMA.mod$loglik+log(log(length(y)))*length(ARMA.mod$coef)

a<-c(ARMA.mod$aic,BIC(ARMA.mod),HQC.ARMA
     ,mse_ARMA,mape_ARMA,mase_ARMA)

#-------------------------------------------------------------------
# Gamma-ARMA(1,4) model with no regressors
#-------------------------------------------------------------------
n<-length(y)

# ### choosing the best model based on AIC BIC and HQ criterion
# b<-matrix(NA,25,3)

# i=0
# i=i+1
# 
p<-4
q<-4
gARMA.mod = ptsr.fit(start = c(0,rep(0,p+q),10), yt =yy ,
                     fit.alpha = TRUE, p = p, q = q,
                     ddist = d.gamma, link1 = "log",
                     link2 = "log", method = "BFGS")
summary(gARMA.mod)
n<-length(y)
AIC.gARMA<--2*gARMA.mod$sll+2*length(gARMA.mod$coefficients)
BIC.gARMA<--2*gARMA.mod$sll+log(length(y))*length(gARMA.mod$coefficients)
HQC.gARMA<--2*gARMA.mod$sll+log(log(length(y)))*length(gARMA.mod$coefficients)

# b[i,]<-c(AIC.gARMA,BIC.gARMA,HQC.gARMA)
# b
# i=i+1


# gARMA[1,1] -365.1871 -351.6126 -366.4463
# gARMA[1,2] -364.0854 -347.1172 -365.6593
# gARMA[1,3] -381.9657 -361.6040 -383.8544
# gARMA[1,4] -379.9677 -356.2123 -382.1712
# gARMA[2,1] -372.8606 -355.8925 -374.4346
# gARMA[2,2] -386.8273 -366.4656 -388.7160
# gARMA[2,3] -385.2035 -361.4481 -387.4070
# gARMA[2,4] -385.9835 -358.8344 -388.5017
# gARMA[3,1] -382.2368 -361.8751 -384.1255
# gARMA[3,2] -385.2358 -361.4804 -387.4393
# gARMA[3,3] -383.3163 -356.1673 -385.8346
# gARMA[3,4] -385.5206 -354.9779 -388.3536
# gARMA[4,1] -380.2966 -356.5412 -382.5000
# gARMA[4,2] -383.7522 -356.6032 -386.2704
# gARMA[4,3] -381.7745 -351.2318 -384.6075
# gARMA[4,4] -388.0148# -354.0785 -391.1626 # best model
# gARMA[0,1] -343.8477 -333.6668 -344.7921
# gARMA[0,2] -348.6908 -335.1163 -349.9499
# gARMA[0,3] -383.8791 -366.9109# -385.4530
# gARMA[0,4] -381.8851 -361.5233 -383.7738
# gARMA[1,0] -365.5218 -355.3410 -366.4662
# gARMA[2,0] -372.5263 -358.9518 -373.7854
# gARMA[3,0] -373.3119 -356.3438 -374.8859
# gARMA[4,0] -371.6408 -351.2790 -373.5295

# gARMA.mod forecast
gARMA.mod.prev<-predict(gARMA.mod,nnew = h1)

# accuracy measures ARMA
gARMA_prev =(dadosh-gARMA.mod.prev$forecast)
gARMA_prev_1=(dadosh[2:h1]-dadosh[1:h2])

mse_gARMA = (sum(gARMA_prev^2))/length(gARMA_prev)
mape_gARMA = sum( abs(gARMA_prev)/abs(dadosh) )/ length(gARMA_prev)
mase_gARMA = sum( abs(gARMA_prev)/(sum(abs(gARMA_prev_1))*(1/(length(gARMA_prev)-1))))/ length(gARMA_prev)

b<-c(AIC.gARMA,BIC.gARMA,HQC.gARMA
     ,mse_gARMA,mape_gARMA,mase_gARMA)

#-------------------------------------------------------------------
# RARMA(0,3) model with no regressors
#-------------------------------------------------------------------
n<-length(y)

RARMA.mod = ptsr.fit(start = c(0,rep(0,5)), yt =yy ,
                     fit.alpha = TRUE, p = 1, q = 4,
                     ddist =d.ray , link1 = "log",
                     link2 = "identity", method = "BFGS")
summary(RARMA.mod)

# ### choosing the best model based on AIC BIC and HQ criterion

# d<-matrix(NA,25,3)

# i=1

p<-0
q<-3
RARMA.mod = ptsr.fit(start = c(0,rep(0,p+q)), yt =yy ,
                     fit.alpha = TRUE, p = p, q = q,
                     ddist =d.ray , link1 = "log",
                     link2 = "identity", method = "BFGS")
summary(RARMA.mod)
n<-length(y)
AIC.RARMA<--2*RARMA.mod$sll+2*length(RARMA.mod$coefficients)
BIC.RARMA<--2*RARMA.mod$sll+log(length(y))*length(RARMA.mod$coefficients)
HQC.RARMA<--2*RARMA.mod$sll+log(log(length(y)))*length(RARMA.mod$coefficients)



# d[i,]<-c(AIC.RARMA,BIC.RARMA,HQC.RARMA)
# d
# i=i+1

# RARMA[1,1] -253.6320 -243.4511 -254.5763
# RARMA[1,2] -254.2006 -240.6261 -255.4597
# RARMA[1,3] -261.7150 -244.7468 -263.2889
# RARMA[1,4] -259.7182 -239.3564 -261.6069
# RARMA[2,1] -253.7387 -240.1642 -254.9979
# RARMA[2,2] -256.7484 -239.7803 -258.3223
# RARMA[2,3] -259.7470 -239.3852 -261.6357
# RARMA[2,4] -257.7502 -233.9948 -259.9537
# RARMA[3,1] -259.2371 -242.2689 -260.8110
# RARMA[3,2] -257.3646 -237.0028 -259.2533
# RARMA[3,3] -257.7415 -233.9861 -259.9450
# RARMA[3,4] -256.3034 -229.1544 -258.8217
# RARMA[4,1] -257.2978 -236.9360 -259.1865
# RARMA[4,2] -255.3603 -231.6049 -257.5637
# RARMA[4,3] -257.3013 -230.1523 -259.8195
# RARMA[4,4] -260.8339 -230.2913 -263.6670
# RARMA[0,1] -246.2790 -239.4918 -246.9086
# RARMA[0,2] -247.4768 -237.2960 -248.4212
# RARMA[0,3] -263.3870# -249.8125# -264.6461 # best model
# RARMA[0,4] -261.3878 -244.4196 -262.9617
# RARMA[1,0] -254.1342 -247.3470 -254.7638
# RARMA[2,0] -255.6068 -245.4259 -256.5511
# RARMA[3,0] -255.4490 -241.8744 -256.7081
# RARMA[4,0] -253.5836 -236.6155 -255.1575

# RARMA - best model
print(RARMA.mod, digits = max(3L, getOption("digits") - 3L),
      signif.stars = getOption("show.signif.stars"))

n<-length(y)

# RARMA forecast
RARMA.mod.prev<-predict(RARMA.mod,nnew = h1)
RARMA_fitted<-RARMA.mod$fitted.values

# accuracy measures ARMA
RARMA_prev =(dadosh-RARMA.mod.prev$forecast)
RARMA_prev_1=(dadosh[2:h1]-dadosh[1:h2])

mse_RARMA = (sum(RARMA_prev^2))/length(RARMA_prev)
mape_RARMA = sum( abs(RARMA_prev)/abs(dadosh) )/ length(RARMA_prev)
mase_RARMA = sum( abs(RARMA_prev)/(sum(abs(RARMA_prev_1))*(1/(length(RARMA_prev)-1))))/ length(RARMA_prev)

d<-c(AIC.RARMA,BIC.RARMA,HQC.RARMA
     ,mse_RARMA,mape_RARMA,mase_RARMA)

###### Models (Table 3)
mod1<-round(coeftest(ARMA.mod)[,-3],4)
mod2<-round(summary(gARMA.mod)$coefficients[,-3],4)
mod3<-round(summary(RARMA.mod)$coefficients[,-3],4)
mod4<-round(BXIIARMA.mod$model[,-3],4)
print("ARMA")
mod1
print("gARMA")
mod2
print("RARMA")
mod3
print("BXIIARMA")
mod4

###### Accuracy measures (Table 4 e 5)
final<-rbind(a,b,d,c)
row.names(final)<-c("ARMA","gARMA","RARMA","BXIIARMA")
colnames(final)<-c("AIC","BIC","HQ","MSE","MAPE","MASE")
round(final,4)



# Graph: BXII-ARMA and RARMA forecast (Figure 6)
postscript(file = "res_prev_wind2.eps",width = 6, height = 4,family = "Times")
par(mfrow=c(1,1))
par(mar=c(2.8, 2.7, 1, 1))
par(mgp=c(1.7, 0.45, 0))
plot(c(BXIIARMA.mod$fitted,pred),type="l",lty=2,col="red", ylim=c(-0.2,1.6),ylab="Volume",xlab="Times")
lines(c(y,dadosh))
lines(c(RARMA_fitted,RARMA_prev),col = "blue",type="l",lty=2)
abline(v = 220,lty = 3)
legend("topleft",c("Observed data","Fitted and forecast values - BXII-ARMA(1,4)","Fitted and forecast values - RARMA(0,3)"),#pch=vpch,
       pt.bg="white", lty=c(1,2,3), bty="n",col=c(1,"red","blue"),cex=0.75
)
dev.off()


################################################################################
####                            WIND SPEED YENAGOA                        ######
################################################################################

# Data set
dados<-read.csv("Wind_Time series model for meteorological data.csv",header=T,sep=",")
attach(dados)

# View(dados)
dados1<-dados$YENAGOA
h1<-12 ## Steps to be predicted
h2<-h1-1
dados2<-dados1[1:(length(dados1)-h1)]
dadosh<-dados1[(length(dados1)-h2):length(dados1)]

div<-10
dadosh<-dadosh/div
yy<-dados2/div

## Create a time series object
dados.ts<- ts(yy,start = c(2017,1), frequency = 12)
y <- dados.ts

# Table of descriptive measures (Table 6)
a<-c(summary(dados.ts,na.rm=T),
     var(dados.ts))
round(a,4)
a<-as.matrix(a)
xtable(t(a),digits = 4)

# Initial graphics
w1<-6
h11<-4

# Serie Plot (Figure 7)
postscript(file = "plotseriewind.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
plot(dados.ts,ylab="Wind Speed")
dev.off()

w1<-5
h11<-4

# ACF Plot (Figure 8.a)
postscript(file = "ACF_wind.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
acf(dados.ts,main="")
dev.off()

# PACF Plot (Figure 8.b)
postscript(file = "PACF_wind.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
pacf(dados.ts,main="")
dev.off()

# Stationarity of the series
y <- dados.ts
adf.test(y, alternative ="stationary") # não é estacionária
pp.test(y, alternative ="stationary") # é estacionária
kpss.test(dados.ts) # é estacionária

# # tendency
n<-length(dados$YENAGOA)
# t <- 1:(n-12) # in sample
# t_hat <- ((n-12)+1):((n-12)+h1) 
# 
# # deterministic seazonality
# C<-cos(2*pi*t/12)  
# C_hat<-cos(2*pi*t_hat/12) 

# Model adjustments

#-------------------------------------------------------------------
# BXII-ARMA(1,2) model
#-------------------------------------------------------------------

# Autoregressive and moving averages 
p1<-3
q1<-1:2

# With covariates
# source("bxiiarmaCOV.fit.r")
# BXIIARMA.mod <- bxiiarmaCOV.fit(y,ar=p1,
#                               ma=q1,X = as.matrix(C),X_hat=as.matrix(C_hat),
#                               h=1,diag=0,tau=0.5)


BXIIARMA.mod <- bxiiarma.fit(y,ar=p1,
                                ma=q1,h=1,mod=3,resid = 3)

### choosing the best model based on AIC BIC and HQ criterion
BXIIARMA.mod2 <- bxiiarma.fit(y,ma=c(1),ar=c(1),h=30,mod=3,resid = 3) 
BXIIARMA.mod3 <- bxiiarma.fit(y,ma=c(1),ar=c(2),h=30,mod=3,resid = 3)
BXIIARMA.mod4 <- bxiiarma.fit(y,ma=c(1),ar=c(3),h=30,mod=3,resid = 3)
BXIIARMA.mod5 <- bxiiarma.fit(y,ma=c(1),ar=c(4),h=30,mod=3,resid = 3)
BXIIARMA.mod6 <- bxiiarma.fit(y,ma=c(2),ar=c(1),h=30,mod=3,resid = 3)
BXIIARMA.mod7 <- bxiiarma.fit(y,ma=c(2),ar=c(2),h=30,mod=3,resid = 3)
BXIIARMA.mod8 <- bxiiarma.fit(y,ma=c(2),ar=c(3),h=30,mod=3,resid = 3)
BXIIARMA.mod9 <- bxiiarma.fit(y,ma=c(2),ar=c(4),h=30,mod=3,resid = 3)
BXIIARMA.mod10 <- bxiiarma.fit(y,ma=c(3),ar=c(1),h=30,mod=3,resid = 3)
BXIIARMA.mod11 <- bxiiarma.fit(y,ma=c(3),ar=c(2),h=30,mod=3,resid = 3)
BXIIARMA.mod12 <- bxiiarma.fit(y,ma=c(3),ar=c(3),h=30,mod=3,resid = 3)
BXIIARMA.mod13 <- bxiiarma.fit(y,ma=c(3),ar=c(4),h=30,mod=3,resid = 3)
BXIIARMA.mod14 <- bxiiarma.fit(y,ma=c(4),ar=c(1),h=30,mod=3,resid = 3)
BXIIARMA.mod15 <- bxiiarma.fit(y,ma=c(4),ar=c(2),h=30,mod=3,resid = 3)
BXIIARMA.mod16 <- bxiiarma.fit(y,ma=c(4),ar=c(3),h=30,mod=3,resid = 3)
BXIIARMA.mod17 <- bxiiarma.fit(y,ma=c(4),ar=c(4),h=30,mod=3,resid = 3)
BXIIARMA.mod18 <- bxiiarma.fit(y,ma=c(1),h=30,mod=2,resid = 3)
BXIIARMA.mod19 <- bxiiarma.fit(y,ma=c(2),h=30,mod=2,resid = 3)
BXIIARMA.mod20 <- bxiiarma.fit(y,ma=c(3),h=30,mod=2,resid = 3)
BXIIARMA.mod21 <- bxiiarma.fit(y,ma=c(4),h=30,mod=2,resid = 3)
BXIIARMA.mod22 <- bxiiarma.fit(y,ar=c(1),h=30,mod=1,resid = 3)
BXIIARMA.mod23 <- bxiiarma.fit(y,ar=c(1),h=30,mod=1,resid = 3)
BXIIARMA.mod24 <- bxiiarma.fit(y,ar=c(1),h=30,mod=1,resid = 3)
BXIIARMA.mod25 <- bxiiarma.fit(y,ar=c(1),h=30,mod=1,resid = 3)

AIC<-c(BXIIARMA.mod$aic,BXIIARMA.mod2$aic,BXIIARMA.mod3$aic,BXIIARMA.mod4$aic,BXIIARMA.mod5$aic,
       BXIIARMA.mod6$aic,BXIIARMA.mod7$aic,BXIIARMA.mod8$aic,BXIIARMA.mod9$aic,BXIIARMA.mod10$aic,
       BXIIARMA.mod11$aic,BXIIARMA.mod12$aic,BXIIARMA.mod13$aic,BXIIARMA.mod14$aic,BXIIARMA.mod15$aic,
       BXIIARMA.mod16$aic,BXIIARMA.mod17$aic,BXIIARMA.mod18$aic,BXIIARMA.mod19$aic,BXIIARMA.mod20$aic,
       BXIIARMA.mod21$aic,BXIIARMA.mod22$aic,BXIIARMA.mod23$aic,BXIIARMA.mod24$aic,BXIIARMA.mod25$aic)

BIC<-c(BXIIARMA.mod$bic,BXIIARMA.mod2$bic,BXIIARMA.mod3$bic,BXIIARMA.mod4$bic,BXIIARMA.mod5$bic,
       BXIIARMA.mod6$bic,BXIIARMA.mod7$bic,BXIIARMA.mod8$bic,BXIIARMA.mod9$bic,BXIIARMA.mod10$bic,
       BXIIARMA.mod11$bic,BXIIARMA.mod12$bic,BXIIARMA.mod13$bic,BXIIARMA.mod14$bic,BXIIARMA.mod15$bic,
       BXIIARMA.mod16$bic,BXIIARMA.mod17$bic,BXIIARMA.mod18$bic,BXIIARMA.mod19$bic,BXIIARMA.mod20$bic,
       BXIIARMA.mod21$bic,BXIIARMA.mod22$bic,BXIIARMA.mod23$bic,BXIIARMA.mod24$bic,BXIIARMA.mod25$bic)

HQ<-c(BXIIARMA.mod$hq,BXIIARMA.mod2$hq,BXIIARMA.mod3$hq,BXIIARMA.mod4$hq,BXIIARMA.mod5$hq,
      BXIIARMA.mod6$hq,BXIIARMA.mod7$hq,BXIIARMA.mod8$hq,BXIIARMA.mod9$hq,BXIIARMA.mod10$hq,
      BXIIARMA.mod11$hq,BXIIARMA.mod12$hq,BXIIARMA.mod13$hq,BXIIARMA.mod14$hq,BXIIARMA.mod15$hq,
      BXIIARMA.mod16$hq,BXIIARMA.mod17$hq,BXIIARMA.mod18$hq,BXIIARMA.mod19$hq,BXIIARMA.mod20$hq,
      BXIIARMA.mod21$hq,BXIIARMA.mod22$hq,BXIIARMA.mod23$hq,BXIIARMA.mod24$hq,BXIIARMA.mod25$hq)

# best model - BXII-ARMA
BXIIARMA.mod$model

residc<-BXIIARMA.mod$resid3

# max_r<- max(res,na.rm=T)
# min_r<- min(res,na.rm=T)
# residc<-res

datatest<-dadosh
n<-length(dados2)


bxii_forecast <- predict_bxii_wind(BXIIARMA.mod)

previsões<-bxii_forecast

# Residuals
residq<-BXIIARMA.mod$resid3
Box.test(residq,lag = 10, type =  "Ljung-Box", fitdf = 0)

# accuracy measures BXII_ARMA
BXII_prev =(dadosh-previsões) 
BXII_prev_1=(dadosh[2:h1]-dadosh[1:h2])

mse_BXII_ARMA = (sum(BXII_prev^2))/length(BXII_prev)
mape_BXII_ARMA = sum( abs(BXII_prev)/abs(dadosh) )/ length(BXII_prev)
mase_BXII_ARMA = sum( abs(BXII_prev)/(sum(abs(BXII_prev_1))*(1/(length(BXII_prev)-1))))/ length(BXII_prev)

c<-c(BXIIARMA.mod$aic,BXIIARMA.mod$bic,
     BXIIARMA.mod$hq,
     mse_BXII_ARMA,mape_BXII_ARMA,mase_BXII_ARMA)

#-------------------------------------------------------------------
# ARMA(1,0) model
#-------------------------------------------------------------------
ARMA.mod<-arima(y,order=c(1,0,2))
ARMA.mod<-auto.arima(y)
ARMA.mod<-arima(y,order=c(1,0,0))

# Best model-ARMA
result<-summary(ARMA.mod)
coeftest(ARMA.mod)
ARMA.mod.prev<-forecast(ARMA.mod,h=h1)

# accuracy measures ARMA
ARMA_prev =(dadosh-ARMA.mod.prev$mean)
ARMA_prev_1=(dadosh[2:h1]-dadosh[1:h2])

mse_ARMA = (sum(ARMA_prev^2))/length(ARMA_prev)
mape_ARMA = sum( abs(ARMA_prev)/abs(dadosh) )/ length(ARMA_prev)
mase_ARMA = sum( abs(ARMA_prev)/(sum(abs(ARMA_prev_1))*(1/(length(ARMA_prev)-1))))/ length(ARMA_prev)

HQC.ARMA<--2*ARMA.mod$loglik+log(log(length(y)))*length(ARMA.mod$coef)

a<-c(ARMA.mod$aic,BIC(ARMA.mod),HQC.ARMA
     ,mse_ARMA,mape_ARMA,mase_ARMA)

################## Fitted
# accuracy measures ARMA

accuracy(result)

#-------------------------------------------------------------------
# Gamma-ARMA(3,3) model with no regressors
#-------------------------------------------------------------------
n<-length(y)

# b<-matrix(NA,24,3)

# i=1

p<-3
q<-3

gARMA.mod = ptsr.fit(start = c(0,rep(0,p+q),10), yt =yy ,
                     fit.alpha = TRUE, p = p, q = q,
                     ddist = d.gamma, link1 = "log",
                     link2 = "log", method = "BFGS")
summary(gARMA.mod)

n<-length(y)
AIC.gARMA<--2*gARMA.mod$sll+2*length(gARMA.mod$coefficients)
BIC.gARMA<--2*gARMA.mod$sll+log(length(y))*length(gARMA.mod$coefficients)
HQC.gARMA<--2*gARMA.mod$sll+log(log(length(y)))*length(gARMA.mod$coefficients)


# b[i,]<-c(AIC.gARMA,BIC.gARMA,HQC.gARMA)
# b
# i=i+1

# best model - gARMA
gARMA.mod
gARMA.mod.prev<-predict(gARMA.mod,nnew = h1)

# accuracy measures gARMA
gARMA_prev =(dadosh-gARMA.mod.prev$forecast)
gARMA_prev_1=(dadosh[2:h1]-dadosh[1:h2])

mse_gARMA = (sum(gARMA_prev^2))/length(gARMA_prev)
mape_gARMA = sum( abs(gARMA_prev)/abs(dadosh) )/ length(gARMA_prev)
mase_gARMA = sum( abs(gARMA_prev)/(sum(abs(gARMA_prev_1))*(1/(length(gARMA_prev)-1))))/ length(gARMA_prev)

b<-c(AIC.gARMA,BIC.gARMA,HQC.gARMA
     ,mse_gARMA,mape_gARMA,mase_gARMA)

#-------------------------------------------------------------------
# RARMA(2,2) model with no regressors
#-------------------------------------------------------------------
n<-length(y)

# d<-matrix(NA,24,3)

# i=1

p<-2
q<-2


RARMA.mod = ptsr.fit(start = c(0,rep(1,p+q)), yt =yy ,
                     fit.alpha = TRUE, p = p ,q = q,
                     ddist =d.ray , link1 = "log",
                     link2 = "identity", method = "BFGS")
summary(RARMA.mod)

n<-length(y)
AIC.RARMA<--2*RARMA.mod$sll+2*length(RARMA.mod$coefficients)
BIC.RARMA<--2*RARMA.mod$sll+log(length(y))*length(RARMA.mod$coefficients)
HQC.RARMA<--2*RARMA.mod$sll+log(log(length(y)))*length(RARMA.mod$coefficients)

# d[i,]<-c(AIC.RARMA,BIC.RARMA,HQC.RARMA)
# d
# i=i+1

# best model
RARMA.mod
RARMA.mod.prev<-predict(RARMA.mod,nnew = h1)
RARMA_fitted<-RARMA.mod$fitted.values

# accuracy measures ARMA
RARMA_prev =(dadosh-RARMA.mod.prev$forecast)
RARMA_prev_1=(dadosh[2:h1]-dadosh[1:h2])

mse_RARMA = (sum(RARMA_prev^2))/length(RARMA_prev)
mape_RARMA = sum( abs(RARMA_prev)/abs(dadosh) )/ length(RARMA_prev)
mase_RARMA = sum( abs(RARMA_prev)/(sum(abs(RARMA_prev_1))*(1/(length(RARMA_prev)-1))))/ length(RARMA_prev)

d<-c(AIC.RARMA,BIC.RARMA,HQC.RARMA
     ,mse_RARMA,mape_RARMA,mase_RARMA)

###### Models (Table 7)
mod1<-round(coeftest(ARMA.mod)[,-3],4)
mod2<-round(summary(gARMA.mod)$coefficients[,-3],4)
mod3<-round(summary(RARMA.mod)$coefficients[,-3],4)
mod4<-round(BXIIARMA.mod$model[,-3],4)
print("ARMA")
mod1
print("gARMA")
mod2
print("RARMA")
mod3
print("BXIIARMA")
mod4

###### Accuracy measures (Table 8 e 9)
final<-rbind(a,b,d,c)
row.names(final)<-c("ARMA","gARMA","RARMA","BXIIARMA")
colnames(final)<-c("AIC","BIC","HQ","MSE","MAPE","MASE")
round(final,4)

# Graph: BXII-ARMA and RARMA forecast (Figure 10)
postscript(file = "res_prev_wind2.eps",width = 6, height = 4,family = "Times")
par(mar=c(5,6,4,1)+.1)
plot(c(BXIIARMA.mod$fitted,previsões),type="l",lty=2,col="red", ylim=c(-0.15,0.4),ylab="Volume",xlab="Times")
lines(c(y,dadosh))
lines(c(RARMA_fitted,RARMA_prev),col = "blue",type="l",lty=2)
abline(v = 49,lty = 3)
legend("bottomleft",c("Observed data","Fitted and forecast values - BXII-ARMA(1,2)","Fitted and forecast values - RARMA(2,2)"),#pch=vpch,
       pt.bg="white", lty=c(1,2,3), bty="n",col=c(1,"red","blue"),cex=0.75
)
dev.off()

