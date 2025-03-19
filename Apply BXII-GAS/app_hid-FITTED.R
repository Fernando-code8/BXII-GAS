# Reference: BXII-GAS
# Created by Fernando José Monteiro de Araújo (nandomonteiro418@gmail.com), September/2024

library(tidyverse)
library(dplyr)
library(zoo)
library(forecast)
library(psych)
library(readr)
library(xtable)
setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/BXII-GAS")
source("GASBXIIFit.R")
setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/Ray-GAS/Ray-GAS-master/Ray-GAS")
source("GASRayFit.R")
setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/gamma-GAS")
source("GASGFit.R")

############################
#### initial quantities ####
############################
w1<-8 # width for plots 
h11<-6 # height for plots

## Data
Dados <- read.delim("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/Apply BXII-GAS/Data/Dados.txt")
dados<-Dados
dados <- dados %>% filter(X == "Selecione Tipo DH")
dados<-as.data.frame(dados[,-2])
dados<-t(dados)
colnames(dados) <- as.character(dados[1, ])
dados <- dados[-c(1,2), ]
dados<-data.frame(dados)
str(dados)
dados <- data.frame(lapply(dados, as.numeric))
dim(dados)
# View(dados)
sum(is.na(dados))

#Descriptive statistics
Desc<-describe(dados)
Desc

#Fitted
dados<-dados[,-c(11,21,30:33)]
dados<-dados[,-c(7,9:11,14:15,21:22,24:25,27)]

datas <- seq(as.Date("2010-01-01"), as.Date("2024-08-31"), by = "day")
R<-1:dim(dados)[2]
mat<-matrix(c(1,1,1,1,1,1,2,2,1,1,1,2,2,2,1,1,2,2,2,2,2,2,1,2,1,2,1,1,1,2,2,2,
              1,2,1,2,1,1,3,3,1,1,1,3,2,2,3,3,2,2,1,3,1,2,3,3,1,2,1,3,3,3,1,1,
              3,3,2,2,3,3,1,2,1,3,1,1,1,3,2,2,1,3,1,2,1,3,3,3,3,3,3,3,3,3,1,3),ncol=4,byrow=T)

acuracia1<-array(NA,c(3,length(R),5))
pv_LB_BXII<-pv_LB_Ray<-pv_LB_GAMMA<-
  nome<-modelo<-modelo_Ray<-modelo_GAMMA<-c()
final<-final_Ray<-final_GAMMA<-list()
start_time <- Sys.time()
j<-0

setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/Apply BXII-GAS")
for(i in R){
  ################################
  #### defining the variables ####
  ################################
  y<-na.omit(dados[,i])
  s<-12
  h=12
  y <- zoo(y, order.by = datas)
  y <- aggregate(y, as.yearmon, mean)
  y<-y/1000
  y<-ts(y[1:(length(y)-h)],frequency = s,start=c(2010,1))
  y[y==0] <- mean(y, na.rm = TRUE)

  t<-1:length(y)
  t_hat<-length(y)+1
  Xsin<-sin(2*pi*t/s)
  Xsin_hat<-sin(2*pi*t_hat/s)
  Xcos<-cos(2*pi*t/s)
  Xcos_hat<-cos(2*pi*t_hat/s)
  X<-cbind(sin(2*pi*t/s),cos(2*pi*t/s))
  X_hat<-cbind(sin(2*pi*t_hat/s),cos(2*pi*t_hat/s))
  
  ##########################
  #### fitting BXII-GAS ####
  ##########################
  
  hidsin_l<-hidcos_l<-hidsincos_l<-hid_l<-list()
  minBXIIsin<-minBXIIcos<-minBXIIsincos<-minBXIIs<-c()
  
  ## Sin
  for(l in 1:nrow(mat)){
    print(c(l,"BXIIsin"))
    pq<-mat[l,]  
    hidsin<-try(BXIIGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=Xsin,X_hat=Xsin_hat,h1=1),silent=T)
    if(class(hidsin)=="try-error" || sum(sum(hidsin$model[-1,4]<0.11,na.rm = FALSE)< length(hidsin$model[-1,4])) > 0 || hidsin$conv != 0 || sum(is.nan(hidsin$model[,4]))>0 || sum(is.infinite(hidsin$residuals))>0){hidsin$aic<-.Machine$double.xmax}
    hidsin_l[[l]]<-hidsin
    minBXIIsin[l]<-hidsin$aic
  }
  
  minBXIIsin1<-min(minBXIIsin)
  id_sin<-which(minBXIIsin==minBXIIsin1)
  if(length(id_sin)>1){hidsin_final<-hidsin_l[[1]]}else{hidsin_final<-hidsin_l[[id_sin]]}
  
  ## Cos
  for(l in 1:nrow(mat)){
    pq<-mat[l,]
    print(c(l,"BXIIcos"))
    hidcos<-try(BXIIGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=Xcos,X_hat=Xcos_hat,h1=1),silent=T)
    if(class(hidcos)=="try-error" || sum(sum(hidcos$model[-1,4]<0.11,na.rm = FALSE)< length(hidcos$model[-1,4])) > 0 || hidcos$conv != 0 || sum(is.nan(hidcos$model[,4]))>0 || sum(is.infinite(hidcos$residuals))>0){hidcos$aic<-.Machine$double.xmax}
    hidcos_l[[l]]<-hidcos
    minBXIIcos[l]<-hidcos$aic
  }
  minBXIIcos1<-min(minBXIIcos)
  id_cos<-which(minBXIIcos==minBXIIcos1)
  if(length(id_cos)>1){hidcos_final<-hidcos_l[[1]]}else{hidcos_final<-hidcos_l[[id_cos]]}
  
  ## Sin&Cos
  for(l in 1:nrow(mat)){
    print(c(l,"BXIIsincos"))
    pq<-mat[l,]  
    hidsincos<-try(BXIIGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=X,X_hat=X_hat,h1=1),silent=T)
    if(class(hidsincos)=="try-error" || sum(sum(hidsincos$model[-1,4]<0.11,na.rm = FALSE)< length(hidsincos$model[-1,4])) > 0 || hidsincos$conv != 0 || sum(is.nan(hidsincos$model[,4]))>0 || sum(is.infinite(hidsincos$residuals))>0){hidsincos$aic<-.Machine$double.xmax}
    hidsincos_l[[l]]<-hidsincos
    minBXIIsincos[l]<-hidsincos$aic
  }
  minBXIIsincos1<-min(minBXIIsincos)
  id_sincos<-which(minBXIIsincos==minBXIIsincos1)
  if(length(id_sincos)>1){hidsincos_final<-hidsincos_l[[1]]}else{hidsincos_final<-hidsincos_l[[id_sincos]]}
  
  ## sem cov
  for(l in 1:nrow(mat)){
    print(c(l,"BXII"))
    pq<-mat[l,]  
    hid<-try(BXIIGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],h1=1),silent=T)
    if(class(hid)=="try-error" || sum(sum(hid$model[-1,4]<0.11,na.rm = FALSE)< length(hid$model[-1,4])) > 0 || hid$conv != 0 || sum(is.nan(hid$model[,4]))>0 || sum(is.infinite(hid$residuals))>0){hid$aic<-.Machine$double.xmax}
    hid_l[[l]]<-hid
    minBXIIs[l]<-hid$aic
  }
  minBXIIs1<-min(minBXIIs)
  id_s<-which(minBXIIs==minBXIIs1)
  if(length(id_s)>1){hid_final<-hid_l[[1]]}else{hid_final<-hid_l[[id_s]]}
  
  minBXII<-min(hidsin_final$aic,hidcos_final$aic,hidsincos_final$aic,hid_final$aic)
  if(minBXII==hidsin_final$aic){final1<-hidsin_final;modelo1<-"BXII with sin covariate"} else{
    if(minBXII==hidcos_final$aic){final1<-hidcos_final; modelo1<-"BXII with cos covariate"} else{
      if(minBXII==hidsincos_final$aic){final1<-hidsincos_final; modelo1<-"BXII with sin and cos covariate"} else{
        if(minBXII==hid_final$aic){final1<-hid_final; modelo1<-"BXII without covariate"} 
      }}}
  
  #############################
  ###  fitting Ray-GAS   ######
  #############################
  
  hidRaysin_l<-hidRaycos_l<-hidRaysincos_l<-hidRay_l<-list()
  minRaysin<-minRaycos<-minRaysincos<-minRays<-c()
  
  ## Sin
  for(l in 1:nrow(mat)){
    print(c(l,"Raysin"))
    pq<-mat[l,]
    hidsin<-try(RayGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=Xsin,X_hat=Xsin_hat,h1=1),silent=T)
    if(class(hidsin)=="try-error" || sum(sum(hidsin$model[-1,4]<0.11,na.rm = FALSE)< length(hidsin$model[-1,4])) > 0 || hidsin$conv != 0 || sum(is.nan(hidsin$model[,4]))>0 || sum(is.infinite(hidsin$residuals))>0){hidsin$aic<-.Machine$double.xmax}
    hidRaysin_l[[l]]<-hidsin
    minRaysin[l]<-hidsin$aic
  }
  minRaysin1<-min(minRaysin)
  id_sin<-which(minRaysin==minRaysin1)
  if(length(id_sin)>1){hidRaysin_final<-hidRaysin_l[[1]]}else{hidRaysin_final<-hidRaysin_l[[id_sin]]}
  
  ## Cos
  for(l in 1:nrow(mat)){
    print(c(l,"Raycos"))
    pq<-mat[l,]  
    hidcos<-try(RayGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=Xcos,X_hat=Xcos_hat,h1=1),silent=T)
    if(class(hidcos)=="try-error" || sum(sum(hidcos$model[-1,4]<0.11,na.rm = FALSE)< length(hidcos$model[-1,4])) > 0 || hidcos$conv != 0 || sum(is.nan(hidcos$model[,4]))>0 || sum(is.infinite(hidcos$residuals))>0){hidcos$aic<-.Machine$double.xmax}
    hidRaycos_l[[l]]<-hidcos
    minRaycos[l]<-hidcos$aic
  }
  minRaycos1<-min(minRaycos)
  id_cos<-which(minRaycos==minRaycos1)
  if(length(id_cos)>1){hidRaycos_final<-hidRaycos_l[[1]]}else{hidRaycos_final<-hidRaycos_l[[id_cos]]}
  
  ## Sin&Cos
  for(l in 1:nrow(mat)){
    print(c(l,"Raysincos"))
    pq<-mat[l,]  
    hidsincos<-try(RayGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=X,X_hat=X_hat,h1=1),silent=T)
    if(class(hidsincos)=="try-error" || sum(sum(hidsincos$model[-1,4]<0.11,na.rm = FALSE)< length(hidsincos$model[-1,4])) > 0 || hidsincos$conv != 0 || sum(is.nan(hidsincos$model[,4]))>0 || sum(is.infinite(hidsincos$residuals))>0){hidsincos$aic<-.Machine$double.xmax}
    hidRaysincos_l[[l]]<-hidsincos
    minRaysincos[l]<-hidsincos$aic
  }
  minRaysincos1<-min(minRaysincos)
  id_sincos<-which(minRaysincos==minRaysincos1)
  if(length(id_sincos)>1){hidRaysincos_final<-hidRaysincos_l[[1]]}else{hidRaysincos_final<-hidRaysincos_l[[id_sincos]]}
  
  ## sem cov
  for(l in 1:nrow(mat)){
    print(c(l,"Ray"))
    pq<-mat[l,]  
    hid<-try(RayGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],h1=1),silent=T)
    if(class(hid)=="try-error" || sum(sum(hid$model[-1,4]<0.11,na.rm = FALSE)< length(hid$model[-1,4])) > 0 || hid$conv != 0 || sum(is.nan(hid$model[,4]))>0 || sum(is.infinite(hid$residuals))>0){hid$aic<-.Machine$double.xmax}
    hidRay_l[[l]]<-hid
    minRays[l]<-hid$aic
  }
  minRays1<-min(minRays)
  id_s<-which(minRays==minRays1)
  if(length(id_s)>1){hidRay_final<-hidRay_l[[1]]}else{hidRay_final<-hidRay_l[[id_s]]}
  
  minRay<-min(hidRaysin_final$aic,hidRaycos_final$aic,hidRaysincos_final$aic,hidRay_final$aic)
  if(minRay==hidRaysin_final$aic){final_Ray1<-hidRaysin_final; modelo_Ray1<-"Ray with sin covariate"} else{
    if(minRay==hidRaycos_final$aic){final_Ray1<-hidRaycos_final;modelo_Ray1<-"Ray with cos covariate"} else{
      if(minRay==hidRaysincos_final$aic){final_Ray1<-hidRaysincos_final; modelo_Ray1<-"Ray with sin and cos covariate"} else{
        if(minRay==hidRay_final$aic){final_Ray1<-hidRay_final; modelo_Ray1<-"Ray without covariate"} 
      }}}
  
  #############################
  ###  fitting gamma-GAS   ####
  #############################
  
  hidgammasin_l<-hidgammacos_l<-hidgammasincos_l<-hidgamma_l<-list()
  mingammasin<-mingammacos<-mingammasincos<-mingammas<-c()
  
  ## Sin
  for(l in 1:nrow(mat)){
    print(c(l,"gammasin"))
    pq<-mat[l,]
    hidsin<-try(GGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=Xsin,X_hat=Xsin_hat,h1=1),silent=T)
    if(class(hidsin)=="try-error" || sum(sum(hidsin$model[-1,4]<0.11,na.rm = FALSE)< length(hidsin$model[-1,4])) > 0 || hidsin$conv != 0 || sum(is.nan(hidsin$model[,4]))>0 || sum(is.infinite(hidsin$residuals))>0){hidsin$aic<-.Machine$double.xmax}
    hidgammasin_l[[l]]<-hidsin
    mingammasin[l]<-hidsin$aic
  }
  mingammasin1<-min(mingammasin)
  id_sin<-which(mingammasin==mingammasin1)
  if(length(id_sin)>1){hidgammasin_final<-hidgammasin_l[[1]]}else{hidgammasin_final<-hidgammasin_l[[id_sin]]}
  
  ## Cos
  for(l in 1:nrow(mat)){
    print(c(l,"gammacos"))
    pq<-mat[l,]  
    hidcos<-try(GGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=Xcos,X_hat=Xcos_hat,h1=1),silent=T)
    if(class(hidcos)=="try-error" || sum(sum(hidcos$model[-1,4]<0.11,na.rm = FALSE)< length(hidcos$model[-1,4])) > 0 || hidcos$conv != 0 || sum(is.nan(hidcos$model[,4]))>0 || sum(is.infinite(hidcos$residuals))>0){hidcos$aic<-.Machine$double.xmax}
    hidgammacos_l[[l]]<-hidcos
    mingammacos[l]<-hidcos$aic
  }
  mingammacos1<-min(mingammacos)
  id_cos<-which(mingammacos==mingammacos1)
  if(length(id_cos)>1){hidgammacos_final<-hidgammacos_l[[1]]}else{hidgammacos_final<-hidgammacos_l[[id_cos]]}
  
  ## Sin&Cos
  for(l in 1:nrow(mat)){
    print(c(l,"gammasincos"))
    pq<-mat[l,]  
    hidsincos<-try(GGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=X,X_hat=X_hat,h1=1),silent=T)
    if(class(hidsincos)=="try-error" || sum(sum(hidsincos$model[-1,4]<0.11,na.rm = FALSE)< length(hidsincos$model[-1,4])) > 0 || hidsincos$conv != 0 || sum(is.nan(hidsincos$model[,4]))>0 || sum(is.infinite(hidsincos$residuals))>0){hidsincos$aic<-.Machine$double.xmax}
    hidgammasincos_l[[l]]<-hidsincos
    mingammasincos[l]<-hidsincos$aic
  }
  mingammasincos1<-min(mingammasincos)
  id_sincos<-which(mingammasincos==mingammasincos1)
  if(length(id_sincos)>1){hidgammasincos_final<-hidgammasincos_l[[1]]}else{hidgammasincos_final<-hidgammasincos_l[[id_sincos]]}
  
  ## sem cov
  for(l in 1:nrow(mat)){
    print(c(l,"gamma"))
    pq<-mat[l,]  
    hid<-try(GGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],h1=1),silent=T)
    if(class(hid)=="try-error" || sum(sum(hid$model[-1,4]<0.11,na.rm = FALSE)< length(hid$model[-1,4])) > 0 || hid$conv != 0 || sum(is.nan(hid$model[,4]))>0 || sum(is.infinite(hid$residuals))>0){hid$aic<-.Machine$double.xmax}
    hidgamma_l[[l]]<-hid
    mingammas[l]<-hid$aic
  }
  mingammas1<-min(mingammas)
  id_s<-which(mingammas==mingammas1)
  if(length(id_s)>1){hidgamma_final<-hidgamma_l[[1]]}else{hidgamma_final<-hidgamma_l[[id_s]]}
  
  mingamma<-min(hidgammasin_final$aic,hidgammacos_final$aic,hidgammasincos_final$aic,hidgamma_final$aic)
  if(mingamma==hidgammasin_final$aic){final_GAMMA1<-hidgammasin_final; modelo_GAMMA1<-"gamma with sin covariate"} else{
    if(mingamma==hidgammacos_final$aic){final_GAMMA1<-hidgammacos_final;modelo_GAMMA1<-"gamma with cos covariate"} else{
      if(mingamma==hidgammasincos_final$aic){final_GAMMA1<-hidgammasincos_final; modelo_GAMMA1<-"gamma with sin and cos covariate"} else{
        if(mingamma==hidgamma_final$aic){final_GAMMA1<-hidgamma_final; modelo_GAMMA1<-"gamma without covariate"} 
      }}}
  
  ###########################
  #### residual analysis ####  
  ###########################
  
  pv_LB_BXII1<-Box.test(final1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_Ray1<-Box.test(final_Ray1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_GAMMA1<-Box.test(final_GAMMA1$residuals, lag = 20, type = "Ljung")$p.value
  
  #############################
  #### selecting resevoirs ####  
  #############################
  
  # if(((pv_LB_BXII1>0.05) #&
      # (pv_LB_Ray1>0.05)&
      # (pv_LB_GAMMA1>0.05)
      # )>0){ 
    j<-j+1
    final[[j]]<-final1
    final_Ray[[j]]<-final_Ray1
    final_GAMMA[[j]]<-final_GAMMA1
    modelo[j]<-modelo1
    modelo_Ray[j]<-modelo_Ray1
    modelo_GAMMA[j]<-modelo_GAMMA1
    nome[j]<-substring(colnames(dados)[i],1,nchar(colnames(dados)[i])#-10
                       )
    print(nome[j])
    pv_LB_BXII[j]<-Box.test(final1$residuals, lag = 20, type = "Ljung")$p.value
    pv_LB_Ray[j]<-Box.test(final_Ray1$residuals, lag = 20, type = "Ljung")$p.value
    pv_LB_GAMMA[j]<-Box.test(final_GAMMA1$residuals, lag = 20, type = "Ljung")$p.value
    
    acuracia1[1,j,]<-c(accuracy(final[[j]]$fitted, y)[,c(2,3,5)],final[[j]]$aic,final[[j]]$bic)
    acuracia1[2,j,]<-c(accuracy(final_Ray[[j]]$fitted, y)[,c(2,3,5)],final_Ray[[j]]$aic,final_Ray[[j]]$bic)
    acuracia1[3,j,]<-c(accuracy(final_GAMMA[[j]]$fitted, y)[,c(2,3,5)],final_GAMMA[[j]]$aic,final_GAMMA[[j]]$bic)
    
    ###############
    #### plots ####  
    ###############
    #---------------------------------
    # time series
    
    time_series<-paste0("Plots/time_series",i,".pdf")
    pdf(time_series,width = w1, height = h11)
    par(mfrow=c(1,1))
    plot(y,main=nome[j],ylim=c(min(y),max(y)+0.15))
    # plot(y,main="")
    lines(final[[j]]$fitted, col=2,lty=2,lwd=2)
    lines(final_Ray[[j]]$fitted, col=5,lty=3,lwd=1.7)
    lines(final_GAMMA[[j]]$fitted, col=6,lty=4,lwd=1.6)
    legend("topright", 
           c("Original","BXII-GAS","Ray-GAS","gamma-ARMA"),
           col = c(1,2,5,6),
           lty= c(1,2,3,4),
           lwd = c(1,2,1.7,1.6), bty="n", cex = 1)
    dev.off()
    # ---------------------------------
    # seasonality
    months<-paste0("Plots/months",i,".pdf")
    pdf(months,width = w1, height = h11)
    par(mfrow=c(1,1))
    monthplot(y,main=nome[j])
    dev.off()
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plot",i,".pdf")
    pdf(resid_plot,width = w1*2, height = h11)
    # acf BXII residuals 
    par(mfrow=c(1,3))
    acf(final[[j]]$residuals,main="BXII-GAS")
    # acf Ray residuals 
    acf(final_Ray[[j]]$residuals,main="Ray-GAS")
    # acf GAMMA residuals 
    acf(final_GAMMA[[j]]$residuals,main="gamma-GAS")
    dev.off()
    # }
}

########################
#### final measures ####  
########################
acuracia1<-acuracia1[,1:j,]
rankRMSE<-apply((apply(acuracia1[,,1],2,rank)==1),1,sum)
rankMAE<-apply((apply(acuracia1[,,2],2,rank)==1),1,sum)
rankMAPE<-apply((apply(acuracia1[,,3],2,rank)==1),1,sum)
rankAIC<-apply((apply(acuracia1[,,4],2,rank)==1),1,sum)
rankBIC<-apply((apply(acuracia1[,,5],2,rank)==1),1,sum)

result<-rbind(rankAIC,rankBIC,
              rankMAPE,
              rankRMSE, rankMAE
)
prop.table(result,1)

type_model_BXII<-table(modelo)
type_model_Ray<-table(modelo_Ray)
type_model_GAMMA<-table(modelo_GAMMA)

count_LB_BXII<-sum(round(pv_LB_BXII,4)>0.05)
count_LB_Ray<-sum(round(na.omit(pv_LB_Ray),4)>0.05)
count_LB_GAMMA<-sum(round(pv_LB_GAMMA,4)>0.05)

results<-data.frame(
  AIC_BXII=acuracia1[1,,4],
  AIC_Ray=acuracia1[2,,4],
  AIC_GAMMA=acuracia1[3,,4], 
  BIC_BXII=acuracia1[1,,5],
  BIC_Ray=acuracia1[2,,5],
  BIC_GAMMA=acuracia1[3,,5], 
  MAE_BXII=acuracia1[1,,2],
  MAE_Ray=acuracia1[2,,2],
  MAE_GAMMA=acuracia1[3,,2],
  MAPE_BXII=acuracia1[1,,3],
  MAPE_RARMA=acuracia1[2,,3],
  MAPE_GAMMA=acuracia1[3,,3],
  RMSE_BXII=acuracia1[1,,1],
  RMSE_RARMA=acuracia1[2,,1],
  RMSE_GAMMA=acuracia1[3,,1],
  pv_LB_BXII=round(pv_LB_BXII,4),
  pv_LB_RARMA=round(pv_LB_Ray,4),
  pv_LB_GAMMA=round(pv_LB_GAMMA,4)
  )
rownames(results)<-nome

setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/Apply BXII-GAS")
saving<-paste0("app_hid_new.RData")
save.image(saving)
end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)

resultados<-results
resultados[7,2]<-0
resultados<-round(resultados,4)
resultados[7,2]<-"_"
# library(gt)

# Criar uma tabela estilizada
resultados[c(1,4,5,8,11,12,13,14,16),-c(1:6,16:18)] %>%
  gt()
