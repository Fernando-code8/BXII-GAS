# Reference: BXII-GAS
# Created by Fernando José Monteiro de Araújo (nandomonteiro418@gmail.com), September/2024

library(tidyverse)
library(dplyr)
library(zoo)
library(forecast)
library(psych)
library(readr)
library(xtable)
library(VGAM)
library(stats)
setwd("~/GitHub/BXII-GAS/BXII-GAS")
source("GASBXIIFit.R")
setwd("~/GitHub/BXII-GAS/Ray-GAS/Ray-GAS-master/Ray-GAS")
source("GASRayFit.R")
setwd("~/GitHub/BXII-GAS/gamma-GAS")
source("GASGFit.R")
setwd("~/GitHub/BXII-GAS/BXII_ARMA")
source("bxiiarmaCOV.fit.r")
source("bxiiarma.fit.r")


############################
#### initial quantities ####
############################
w1<-8 # width for plots 
h11<-6 # height for plots

## Data
Dados <- read.delim("~/GitHub/BXII-GAS/Apply BXII-GAS/Data/Dados.txt")
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
dados<-dados[,-c(3,5:9,13,15)]
dados<-dados[,-c(3,4,7)]


datas <- seq(as.Date("2010-01-01"), as.Date("2024-08-31"), by = "day")
R<-1:dim(dados)[2]
mat<-matrix(c(rep(1,1*4)),ncol=4,byrow=T)

acuracia1<-array(NA,c(7,length(R),5))
pv_LB_BXII<-pv_LB_Ray<-pv_LB_GAMMA<-pv_LB_BXIIARMA<-pv_LB_ETS<-pv_LB_ARMA<-pv_LB_SARMA<-
  nome<-modelo<-modelo_Ray<-modelo_GAMMA<-modeloBXIIARMA<-c()
final<-final_Ray<-final_GAMMA<-finalBXIIARMA<-final_ETS<-final_ARMA<-final_SARMA<-list()
start_time <- Sys.time()
j<-0

setwd("~/GitHub/BXII-GAS/Apply BXII-GAS")
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
    if(class(hidsin)=="try-error" || Box.test(hidsin$residuals, lag = 20, type = "Ljung")$p.value<0.05 || sum(sum(hidsin$model[-1,4]<0.11,na.rm = FALSE)< length(hidsin$model[-1,4])) > 0 || hidsin$conv != 0 || sum(is.nan(hidsin$model[,4]))>0 || sum(is.infinite(hidsin$residuals))>0){hidsin$aic<-.Machine$double.xmax}
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
    if(class(hidcos)=="try-error" || Box.test(hidcos$residuals, lag = 20, type = "Ljung")$p.value<0.05 || sum(sum(hidcos$model[-1,4]<0.11,na.rm = FALSE)< length(hidcos$model[-1,4])) > 0 || hidcos$conv != 0 || sum(is.nan(hidcos$model[,4]))>0 || sum(is.infinite(hidcos$residuals))>0){hidcos$aic<-.Machine$double.xmax}
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
    if(class(hidsincos)=="try-error" || Box.test(hidsincos$residuals, lag = 20, type = "Ljung")$p.value<0.05 || sum(sum(hidsincos$model[-1,4]<0.11,na.rm = FALSE)< length(hidsincos$model[-1,4])) > 0 || hidsincos$conv != 0 || sum(is.nan(hidsincos$model[,4]))>0 || sum(is.infinite(hidsincos$residuals))>0){hidsincos$aic<-.Machine$double.xmax}
    hidsincos_l[[l]]<-hidsincos
    minBXIIsincos[l]<-hidsincos$aic
  }
  minBXIIsincos1<-min(minBXIIsincos)
  id_sincos<-which(minBXIIsincos==minBXIIsincos1)
  if(length(id_sincos)>1){hidsincos_final<-hidsincos_l[[1]]}else{hidsincos_final<-hidsincos_l[[id_sincos]]}
  
  ## sem cov
  # for(l in 1:nrow(mat)){
  #   print(c(l,"BXII"))
  #   pq<-mat[l,]  
  #   hid<-try(BXIIGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],h1=1),silent=T)
  #   if(class(hid)=="try-error" || Box.test(hid$residuals, lag = 20, type = "Ljung")$p.value<0.05 || sum(sum(hid$model[-1,4]<0.11,na.rm = FALSE)< length(hid$model[-1,4])) > 0 || hid$conv != 0 || sum(is.nan(hid$model[,4]))>0 || sum(is.infinite(hid$residuals))>0){hid$aic<-.Machine$double.xmax}
  #   hid_l[[l]]<-hid
  #   minBXIIs[l]<-hid$aic
  # }
  # minBXIIs1<-min(minBXIIs)
  # id_s<-which(minBXIIs==minBXIIs1)
  # if(length(id_s)>1){hid_final<-hid_l[[1]]}else{hid_final<-hid_l[[id_s]]}
  
  hid_final<-c()
  hid_final$aic<-.Machine$double.xmax
  
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
    if(class(hidsin)=="try-error" #|| Box.test(hidsin$residuals, lag = 20, type = "Ljung")$p.value<0.05 
       || sum(sum(hidsin$model[-1,4]<0.11,na.rm = FALSE)< length(hidsin$model[-1,4])) > 0 || hidsin$conv != 0 || sum(is.nan(hidsin$model[,4]))>0 || sum(is.infinite(hidsin$residuals))>0){hidsin$aic<-.Machine$double.xmax}
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
    if(class(hidcos)=="try-error" #|| Box.test(hidcos$residuals, lag = 20, type = "Ljung")$p.value<0.05 
       || sum(sum(hidcos$model[-1,4]<0.11,na.rm = FALSE)< length(hidcos$model[-1,4])) > 0 || hidcos$conv != 0 || sum(is.nan(hidcos$model[,4]))>0 || sum(is.infinite(hidcos$residuals))>0){hidcos$aic<-.Machine$double.xmax}
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
    if(class(hidsincos)=="try-error" #|| Box.test(hidsincos$residuals, lag = 20, type = "Ljung")$p.value<0.05 
       || sum(sum(hidsincos$model[-1,4]<0.11,na.rm = FALSE)< length(hidsincos$model[-1,4])) > 0 || hidsincos$conv != 0 || sum(is.nan(hidsincos$model[,4]))>0 || sum(is.infinite(hidsincos$residuals))>0){hidsincos$aic<-.Machine$double.xmax}
    hidRaysincos_l[[l]]<-hidsincos
    minRaysincos[l]<-hidsincos$aic
  }
  minRaysincos1<-min(minRaysincos)
  id_sincos<-which(minRaysincos==minRaysincos1)
  if(length(id_sincos)>1){hidRaysincos_final<-hidRaysincos_l[[1]]}else{hidRaysincos_final<-hidRaysincos_l[[id_sincos]]}
  
  # ## sem cov
  # for(l in 1:nrow(mat)){
  #   print(c(l,"Ray"))
  #   pq<-mat[l,]  
  #   hid<-try(RayGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],h1=1),silent=T)
  #   if(class(hid)=="try-error" || Box.test(hid$residuals, lag = 20, type = "Ljung")$p.value<0.05 || sum(sum(hid$model[-1,4]<0.11,na.rm = FALSE)< length(hid$model[-1,4])) > 0 || hid$conv != 0 || sum(is.nan(hid$model[,4]))>0 || sum(is.infinite(hid$residuals))>0){hid$aic<-.Machine$double.xmax}
  #   hidRay_l[[l]]<-hid
  #   minRays[l]<-hid$aic
  # }
  # minRays1<-min(minRays)
  # id_s<-which(minRays==minRays1)
  # if(length(id_s)>1){hidRay_final<-hidRay_l[[1]]}else{hidRay_final<-hidRay_l[[id_s]]}
  
  hidRay_final<-c()
  hidRay_final$aic<-.Machine$double.xmax
  
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
  mingammasin2<-mingammacos2<-mingammasincos2<-mingammas2<-c()
  
  ## Sin
  for(l in 1:nrow(mat)
  ){
    print(c(l,"gammasin"))
    pq<-mat[l,]  
    hidsin<-try(GGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=Xsin,X_hat=Xsin_hat,h1=1),silent=T)
    if(class(hidsin)=="try-error" || hidsin$conv != 0 || sum(is.nan(hidsin$model[,4]))>0 || sum(is.infinite(hidsin$residuals))>0|| hidsin$conv != 0 || sum(is.nan(hidsin$model[,4]))>0 || sum(is.infinite(hidsin$residuals))>0){hidsin$aic<-.Machine$double.xmax}
    mingammasin2[l]<-hidsin$aic
    if(sum(hidsin$model[-1,4] < 0.11, na.rm = TRUE) < length(hidsin$model[-1,4])){hidsin$aic <- .Machine$double.xmax}
    hidgammasin_l[[l]]<-hidsin
    mingammasin[l]<-hidsin$aic
  }
  mingammasin1<-min(mingammasin)
  id_sin<-which(mingammasin==mingammasin1)
  if(length(id_sin)>1){
    mingammasin3<-min(mingammasin2)
    id_sin2<-which(mingammasin2==mingammasin3)
    hidgammasin_final<-hidgammasin_l[[id_sin2]]
    hidgammasin_final$aic<-mingammasin3
  }else{hidgammasin_final<-hidgammasin_l[[id_sin]]}  
  
  ## Cos
  for(l in 1:nrow(mat)
  ){
    print(c(l,"gammacos"))
    pq<-mat[l,]  
    hidcos<-try(GGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=Xcos,X_hat=Xcos_hat,h1=1),silent=T)
    if(class(hidcos)=="try-error" || hidcos$conv != 0 || sum(is.nan(hidcos$model[,4]))>0 || sum(is.infinite(hidcos$residuals))>0|| hidcos$conv != 0 || sum(is.nan(hidcos$model[,4]))>0 || sum(is.infinite(hidcos$residuals))>0){hidcos$aic<-.Machine$double.xmax}
    mingammacos2[l]<-hidcos$aic
    if(sum(hidcos$model[-1,4] < 0.11, na.rm = TRUE) < length(hidcos$model[-1,4])){hidcos$aic <- .Machine$double.xmax}
    hidgammacos_l[[l]]<-hidcos
    mingammacos[l]<-hidcos$aic
  }
  mingammacos1<-min(mingammacos)
  id_cos<-which(mingammacos==mingammacos1)
  if(length(id_cos)>1){
    mingammacos3<-min(mingammacos2)
    id_cos2<-which(mingammacos2==mingammacos3)
    hidgammacos_final<-hidgammacos_l[[id_cos2]]
    hidgammacos_final$aic<-mingammacos3
  }else{hidgammacos_final<-hidgammacos_l[[id_cos]]}  
  
  ## Sin&Cos
  for(l in 1:nrow(mat)
  ){
    print(c(l,"gammasincos"))
    pq<-mat[l,]  
    hidsincos<-try(GGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],X=X,X_hat=X_hat,h1=1),silent=T)
    if(class(hidsincos)=="try-error" || hidsincos$conv != 0 || sum(is.nan(hidsincos$model[,4]))>0 || sum(is.infinite(hidsincos$residuals))>0|| hidsincos$conv != 0 || sum(is.nan(hidsincos$model[,4]))>0 || sum(is.infinite(hidsincos$residuals))>0){hidsincos$aic<-.Machine$double.xmax}
    mingammasincos2[l]<-hidsincos$aic
    if(sum(hidsincos$model[-1,4] < 0.11, na.rm = TRUE) < length(hidsincos$model[-1,4])){hidsincos$aic <- .Machine$double.xmax}
    hidgammasincos_l[[l]]<-hidsincos
    mingammasincos[l]<-hidsincos$aic
  }
  mingammasincos1<-min(mingammasincos)
  id_sincos<-which(mingammasincos==mingammasincos1)
  if(length(id_sincos)>1){
    mingammasincos3<-min(mingammasincos2)
    id_sincos2<-which(mingammasincos2==mingammasincos3)
    hidgammasincos_final<-hidgammasincos_l[[id_sincos2]]
    hidgammasincos_final$aic<-mingammasincos3
  }else{hidgammasincos_final<-hidgammasincos_l[[id_sincos]]}  
  
  # ## sem cov
  # for(l in 1:nrow(mat)){
  #   print(c(l,"gamma"))
  #   pq<-mat[l,]  
  #   hid<-try(GGAS.fit(y,pq[1]:pq[2],pq[3]:pq[4],h1=1),silent=T)
  #   if(class(hid)=="try-error" #|| Box.test(hid$residuals, lag = 20, type = "Ljung")$p.value<0.05 
  # || sum(sum(hid$model[-1,4]<0.11,na.rm = FALSE)< length(hid$model[-1,4])) > 0 || hid$conv != 0 || sum(is.nan(hid$model[,4]))>0 || sum(is.infinite(hid$residuals))>0){hid$aic<-.Machine$double.xmax}
  #   hidgamma_l[[l]]<-hid
  #   mingammas[l]<-hid$aic
  # }
  # mingammas1<-min(mingammas)
  # id_s<-which(mingammas==mingammas1)
  # if(length(id_s)>1){hidgamma_final<-hidgamma_l[[1]]}else{hidgamma_final<-hidgamma_l[[id_s]]}
  
  hidgamma_final<-c()
  hidgamma_final$aic<-.Machine$double.xmax
  
  mingamma<-min(hidgammasin_final$aic,hidgammacos_final$aic,hidgammasincos_final$aic,hidgamma_final$aic)
  if(mingamma==hidgammasin_final$aic){final_GAMMA1<-hidgammasin_final; modelo_GAMMA1<-"gamma with sin covariate"} else{
    if(mingamma==hidgammacos_final$aic){final_GAMMA1<-hidgammacos_final;modelo_GAMMA1<-"gamma with cos covariate"} else{
      if(mingamma==hidgammasincos_final$aic){final_GAMMA1<-hidgammasincos_final; modelo_GAMMA1<-"gamma with sin and cos covariate"} else{
        if(mingamma==hidgamma_final$aic){final_GAMMA1<-hidgamma_final; modelo_GAMMA1<-"gamma without covariate"} 
      }}}
  
  ##########################
  ### fitting BXII-ARMA ####
  ##########################
  
  hidBXIIARMAsin_l<-hidBXIIARMAcos_l<-hidBXIIARMAsincos_l<-hidBXIIARMA_l<-list()
  minBXIIARMAsin<-minBXIIARMAcos<-minBXIIARMAsincos<-minBXIIARMAs<-c()
  matBXII<-matrix(c(1,1,1,1,1,1,2,2,1,1,1,2,2,2,1,1,2,2,2,2,2,2,1,2,1,2,1,1,1,2,2,2,
                    1,2,1,2,1,1,3,3,1,1,1,3,2,2,3,3,2,2,1,3,1,2,3,3,1,2,1,3,3,3,1,1,
                    3,3,2,2,3,3,1,2,1,3,1,1,1,3,2,2,1,3,1,2,1,3,3,3,3,3,3,3,3,3,1,3
  ),ncol=4,byrow=T)
  
  ## Sin
  for(l in 1:nrow(matBXII)){
    print(c(l,"BXIIARMAsin"))
    pqbxii<-matBXII[l,]
    hidsin<-try(bxiiarmaCOV.fit(y,ar=pqbxii[1]:pqbxii[2],ma=pqbxii[3]:pqbxii[4],X = as.matrix(Xsin),X_hat=as.matrix(Xsin_hat),h=1,diag=0,tau=0.5),silent=T)
    if(class(hidsin)=="try-error" || sum(sum(hidsin$model[-1,4]<0.11,na.rm = FALSE)< length(hidsin$model[-1,4])) > 0 || hidsin$conv != 0 || sum(is.nan(hidsin$model[,4]))>0 || sum(is.infinite(hidsin$residuals))>0){hidsin$aic<-.Machine$double.xmax}
    hidBXIIARMAsin_l[[l]]<-hidsin
    minBXIIARMAsin[l]<-hidsin$aic
  }

  minBXIIARMAsin1<-min(minBXIIARMAsin)
  id_sin<-which(minBXIIARMAsin==minBXIIARMAsin1)
  if(length(id_sin)>1){hidBXIIARMAsin_final<-hidBXIIARMAsin_l[[1]]}else{hidBXIIARMAsin_final<-hidBXIIARMAsin_l[[id_sin]]}
  
  #hidBXIIARMAsin_final<-c()
  #hidBXIIARMAsin_final$aic<-.Machine$double.xmax
  
  ## Cos
  for(l in 1:nrow(matBXII)){
    pqbxii<-matBXII[l,]
    print(c(l,"BXIIARMAcos"))
    hidcos<-try(bxiiarmaCOV.fit(y,ar=pqbxii[1]:pqbxii[2],ma=pqbxii[3]:pqbxii[4],X = as.matrix(Xcos),X_hat=as.matrix(Xcos_hat),h=1,diag=0,tau=0.5),silent=T)
    if(class(hidcos)=="try-error" || sum(sum(hidcos$model[-1,4]<0.11,na.rm = FALSE)< length(hidcos$model[-1,4])) > 0 || hidcos$conv != 0 || sum(is.nan(hidcos$model[,4]))>0 || sum(is.infinite(hidcos$residuals))>0){hidcos$aic<-.Machine$double.xmax}
    hidBXIIARMAcos_l[[l]]<-hidcos
    minBXIIARMAcos[l]<-hidcos$aic
  }
  minBXIIARMAcos1<-min(minBXIIARMAcos)
  id_cos<-which(minBXIIARMAcos==minBXIIARMAcos1)
  if(length(id_cos)>1){hidBXIIARMAcos_final<-hidBXIIARMAcos_l[[1]]}else{hidBXIIARMAcos_final<-hidBXIIARMAcos_l[[id_cos]]}
  
  #hidBXIIARMAcos_final<-c()
  #hidBXIIARMAcos_final$aic<-.Machine$double.xmax
  
  ## Sin&Cos
  for(l in 1:nrow(matBXII)){
    print(c(l,"BXIIARMAsincos"))
    pqbxii<-matBXII[l,]
    hidsincos<-try(bxiiarmaCOV.fit(y,ar=pqbxii[1]:pqbxii[2],ma=pqbxii[3]:pqbxii[4],X = as.matrix(X),X_hat=as.matrix(X_hat),h=1,diag=0,tau=0.5),silent=T)
    if(class(hidsincos)=="try-error" || sum(sum(hidsincos$model[-1,4]<0.11,na.rm = FALSE)< length(hidsincos$model[-1,4])) > 0 || hidsincos$conv != 0 || sum(is.nan(hidsincos$model[,4]))>0 || sum(is.infinite(hidsincos$residuals))>0){hidsincos$aic<-.Machine$double.xmax}
    hidBXIIARMAsincos_l[[l]]<-hidsincos
    minBXIIARMAsincos[l]<-hidsincos$aic
  }
  minBXIIARMAsincos1<-min(minBXIIARMAsincos)
  id_sincos<-which(minBXIIARMAsincos==minBXIIARMAsincos1)
  if(length(id_sincos)>1){hidBXIIARMAsincos_final<-hidBXIIARMAsincos_l[[1]]}else{hidBXIIARMAsincos_final<-hidBXIIARMAsincos_l[[id_sincos]]}
  
  #hidBXIIARMAsincos_final<-c()
  #hidBXIIARMAsincos_final$aic<-.Machine$double.xmax
  
  # sem cov

  # for(l in 1:nrow(matBXII)){
  #  print(c(l,"BXIIARMA"))
  #  pqbxii<-matBXII[l,]
  #  hid<-try(bxiiarma.fit(y,ar=pqbxii[1]:pqbxii[2],ma=pqbxii[3]:pqbxii[4],h1=1,resid = 3,diag1=1),silent=T)
  #  if(class(hid)=="try-error" || sum(sum(hid$model[-1,4]<0.11,na.rm = FALSE)< length(hid$model[-1,4])) > 0 || hid$conv != 0 || sum(is.nan(hid$model[,4]))>0 || sum(is.infinite(hid$resid1))>0){hid$aic<-.Machine$double.xmax}
  #  hidBXIIARMA_l[[l]]<-hid
  #  minBXIIARMAs[l]<-hid$aic
  # }
  # minBXIIARMAs1<-min(minBXIIARMAs)
  # id_s<-which(minBXIIARMAs==minBXIIARMAs1)
  # if(length(id_s)>1){hidBXIIARMA_final<-hidBXIIARMA_l[[1]]}else{hidBXIIARMA_final<-hidBXIIARMA_l[[id_s]]}
  
  hidBXIIARMA_final<-c()
  hidBXIIARMA_final$aic<-.Machine$double.xmax
  
  minBXIIARMA<-min(hidBXIIARMAsin_final$aic,hidBXIIARMAcos_final$aic,hidBXIIARMAsincos_final$aic,hidBXIIARMA_final$aic)
  if(minBXIIARMA==hidBXIIARMAsin_final$aic){finalBXIIARMA1<-hidBXIIARMAsin_final;modeloBXIIARMA1<-"BXIIARMA with sin covariate"} else{
    if(minBXIIARMA==hidBXIIARMAcos_final$aic){finalBXIIARMA1<-hidBXIIARMAcos_final; modeloBXIIARMA1<-"BXIIARMA with cos covariate"} else{
      if(minBXIIARMA==hidBXIIARMAsincos_final$aic){finalBXIIARMA1<-hidBXIIARMAsincos_final; modeloBXIIARMA1<-"BXIIARMA with sin and cos covariate"} else{
        if(minBXIIARMA==hidBXIIARMA_final$aic){finalBXIIARMA1<-hidBXIIARMA_final; modeloBXIIARMA1<-"BXIIARMA without covariate"} 
      }}}
  
  #############################
  ###     classic models    ###
  #############################
  
  # Holt-Winters (HW)
  
  final_ETS1 <- try(ets(y, model = "AAA"), silent = TRUE)
  
  aic_ETS <- if(inherits(final_ETS1, "try-error")){
    Inf
    } else {
      final_ETS1$aic
    }
  
  # ARMA
  arma_l <- list()
  minARMA <- numeric(nrow(mat))
  
  for (l in 1:nrow(mat)) {
    pq <- mat[l, ]
    print(c(l, "ARMA"))
    
    fit <- try(
      Arima(y, order = c(pq[1], 0, pq[3])),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) {
      minARMA[l] <- Inf
      arma_l[[l]] <- NULL
    } else {
      minARMA[l] <- fit$aic
      arma_l[[l]] <- fit
    }
  }
  
  id_arma <- which.min(minARMA)
  final_ARMA1 <- arma_l[[id_arma]]
  
  # SARMA
  sarma_l <- list()
  minSARMA <- c()
  matsarma<-matrix(c(1,1,1,1,1,1,0,1,1,1,1,0),ncol=4,byrow=T)
  for (l in 1:nrow(matsarma)) {
    pqs <- matsarma[l,]
    print(c(l, "SARMA"))
    
    fit <- try(
      Arima(
        y,
        order = c(pqs[1], 0, pqs[2]),
        seasonal = list(order = c(pqs[3], 0, pqs[4]), period = s)
      ),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) {
      aic <- Inf
    } else {
      aic <- fit$aic
    }
    
    sarma_l[[l]] <- fit
    minSARMA[l] <- aic
  }
  
  id_sarma <- which.min(minSARMA)
  final_SARMA1 <- sarma_l[[id_sarma]]
  
  
  ###########################
  #### residual analysis ####  
  ###########################
  
  pv_LB_BXII1<-Box.test(final1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_Ray1<-Box.test(final_Ray1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_GAMMA1<-Box.test(final_GAMMA1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_BXIIARMA1<-Box.test(finalBXIIARMA1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_ETS1<-Box.test(final_ETS1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_ARMA1<-Box.test(final_ARMA1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_SARMA1<-Box.test(final_SARMA1$residuals, lag = 20, type = "Ljung")$p.value
  
  
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
    finalBXIIARMA[[j]]<-finalBXIIARMA1
    final_ETS[[j]]<-final_ETS1
    final_ARMA[[j]]<-final_ARMA1
    final_SARMA[[j]]<-final_SARMA1
    modelo[j]<-modelo1
    modelo_Ray[j]<-modelo_Ray1
    modelo_GAMMA[j]<-modelo_GAMMA1
    modeloBXIIARMA[j]<-modeloBXIIARMA1
    nome[j]<-substring(colnames(dados)[i],1,nchar(colnames(dados)[i])#-10
                       )
    print(nome[j])
    pv_LB_BXII[j]<-Box.test(final1$residuals, lag = 20, type = "Ljung")$p.value
    pv_LB_Ray[j]<-Box.test(final_Ray1$residuals, lag = 20, type = "Ljung")$p.value
    pv_LB_GAMMA[j]<-Box.test(final_GAMMA1$residuals, lag = 20, type = "Ljung")$p.value
    pv_LB_BXIIARMA[j]<-Box.test(finalBXIIARMA1$residuals, lag = 20, type = "Ljung")$p.value
    pv_LB_ETS[j]<-Box.test(final_ETS1$residuals, lag = 20, type = "Ljung")$p.value
    pv_LB_ARMA[j]<-Box.test(final_ARMA1$residuals, lag = 20, type = "Ljung")$p.value
    pv_LB_SARMA[j]<-Box.test(final_SARMA1$residuals, lag = 20, type = "Ljung")$p.value
    
    acuracia1[1,j,]<-c(accuracy(final[[j]]$fitted, y)[,c(2,3,5)],final[[j]]$aic,final[[j]]$bic)
    acuracia1[2,j,]<-c(accuracy(final_Ray[[j]]$fitted, y)[,c(2,3,5)],final_Ray[[j]]$aic,final_Ray[[j]]$bic)
    acuracia1[3,j,]<-c(accuracy(final_GAMMA[[j]]$fitted, y)[,c(2,3,5)],final_GAMMA[[j]]$aic,final_GAMMA[[j]]$bic)
    acuracia1[4,j,]<-c(accuracy(finalBXIIARMA[[j]]$fitted, y)[,c(2,3,5)],finalBXIIARMA[[j]]$aic,finalBXIIARMA[[j]]$bic)
    acuracia1[5,j,]<-c(accuracy(final_ETS[[j]]$fitted, y)[,c(2,3,5)],final_ETS[[j]]$aic,final_ETS[[j]]$bic)
    acuracia1[6,j,]<-c(accuracy(final_ARMA[[j]]$fitted, y)[,c(2,3,5)],final_ARMA[[j]]$aic,final_ARMA[[j]]$bic)
    acuracia1[7,j,]<-c(accuracy(final_SARMA[[j]]$fitted, y)[,c(2,3,5)],final_SARMA[[j]]$aic,final_SARMA[[j]]$bic)
    
    
    ###############
    #### plots ####  
    ###############
    setwd("~/GitHub/BXII-GAS/Apply BXII-GAS")
    
    # time series plot
    time_series_plot<-paste0("Plots/time_series_plot",i,".pdf")
    pdf(time_series_plot,width = w1, height = h11)
    par(mfrow=c(1,1))
    plot(y,main=nome[j],ylim=c(min(y),max(y)+0.15))
    #plot(y,main="")
    dev.off()
    
    # Acf plot
    Acf_plot<-paste0("Plots/Acf_plot",i,".pdf")
    pdf(Acf_plot,width =4, height = 4)
    # acf  
    par(mfrow=c(1,1))
    acf(y,main=nome[j]#,main=""
        )
    dev.off()
    #---------------------------------
    # time series plot - fitted models
    time_series<-paste0("Plots/time_series",i,".pdf")
    pdf(time_series,width = w1, height = h11)
    par(mfrow=c(1,1))
    plot(y,main=nome[j],ylim=c(min(y),max(y)+0.15))
    # plot(y,main="")
    lines(final[[j]]$fitted, col=2,lty=2,lwd=2)
    lines(final_Ray[[j]]$fitted, col=5,lty=3,lwd=1.7)
    lines(finalBXIIARMA[[j]]$fitted, col=4,lty=7,lwd=1.5)
    lines(final_ETS[[j]]$fitted, col=7,lty=5,lwd=1.8)
    lines(final_ARMA[[j]]$fitted, col=6,lty=4,lwd=1.6)
    lines(final_SARMA[[j]]$fitted, col=3,lty=6,lwd=1.9)
    legend("topright", 
           c("Original","BXII-GAS","Ray-GAS","BXII-ARMA","ETS","ARMA","SARMA"),
           col = c(1,2,5,4,7,6,3),
           lty= c(1,2,3,7,5,4,6),
           lwd = c(1,2,1.7,1.5,1.8,1.6,1.9), bty="n", cex = 1)
    dev.off()
    # ---------------------------------
    # seasonality
    months<-paste0("Plots/months",i,".pdf")
    pdf(months,width = 10, height = 6)
    par(mfrow=c(1,1))
    monthplot(y,main=nome[j],ylab = "RF",base = "median")
    dev.off()
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotBXII",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf BXII residuals 
    acf(final[[j]]$residuals,main="")
    dev.off()
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotRay",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf Ray residuals 
    acf(final_Ray[[j]]$residuals,main="")
    dev.off()
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotgamma",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf GAMMA residuals 
    acf(final_GAMMA[[j]]$residuals,main="")
    dev.off()
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotBXIIARMA",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf BXIIARMA residuals 
    acf(finalBXIIARMA[[j]]$residuals,main="")
    dev.off()
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotETS",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf ETS residuals 
    acf(final_ETS[[j]]$residuals,main="ETS")
    dev.off()
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotARMA",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf ARMA residuals 
    acf(final_ARMA[[j]]$residuals,main="ARMA")
    dev.off()
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotSARMA",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf SARMA residuals 
    acf(final_SARMA[[j]]$residuals,main="SARMA")
    dev.off()
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
type_model_BXIIARMA<-table(modeloBXIIARMA)

count_LB_BXII<-sum(round(pv_LB_BXII,4)>0.05)
count_LB_Ray<-sum(round(na.omit(pv_LB_Ray),4)>0.05)
count_LB_GAMMA<-sum(round(pv_LB_GAMMA,4)>0.05)
count_LB_BXIIARMA<-sum(round(pv_LB_BXIIARMA,4)>0.05)
count_LB_ETS<-sum(round(na.omit(pv_LB_ETS),4)>0.05)
count_LB_ARMA<-sum(round(pv_LB_ARMA,4)>0.05)
count_LB_SARMA<-sum(round(pv_LB_SARMA,4)>0.05)

results<-data.frame(
  AIC_BXII=acuracia1[1,,4],
  AIC_Ray=acuracia1[2,,4],
  AIC_GAMMA=acuracia1[3,,4],
  AIC_BXIIARMA=acuracia1[4,,4],
  AIC_ETS=acuracia1[5,,4],
  AIC_ARMA=acuracia1[6,,4],
  AIC_SARMA=acuracia1[7,,4],
  BIC_BXII=acuracia1[1,,5],
  BIC_Ray=acuracia1[2,,5],
  BIC_GAMMA=acuracia1[3,,5],
  BIC_BXIIARMA=acuracia1[4,,5], 
  BIC_ETS=acuracia1[5,,5], 
  BIC_ARMA=acuracia1[6,,5], 
  BIC_SARMA=acuracia1[7,,5], 
  MAE_BXII=acuracia1[1,,2],
  MAE_Ray=acuracia1[2,,2],
  MAE_GAMMA=acuracia1[3,,2],
  MAE_BXIIARMA=acuracia1[4,,2],
  MAE_ETS=acuracia1[5,,2],
  MAE_ARMA=acuracia1[6,,2],
  MAE_SARMA=acuracia1[7,,2],
  MAPE_BXII=acuracia1[1,,3],
  MAPE_Ray=acuracia1[2,,3],
  MAPE_GAMMA=acuracia1[3,,3],
  MAPE_BXIIARMA=acuracia1[4,,3],
  MAPE_ETS=acuracia1[5,,3],
  MAPE_ARMA=acuracia1[6,,3],
  MAPE_SARMA=acuracia1[7,,3],
  RMSE_BXII=acuracia1[1,,1],
  RMSE_Ray=acuracia1[2,,1],
  RMSE_GAMMA=acuracia1[3,,1],
  RMSE_BXIIARMA=acuracia1[4,,1],
  RMSE_ETS=acuracia1[5,,1],
  RMSE_ARMA=acuracia1[6,,1],
  RMSE_SARMA=acuracia1[7,,1],
  pv_LB_BXII=round(pv_LB_BXII,4),
  pv_LB_Ray=round(pv_LB_Ray,4),
  pv_LB_GAMMA=round(pv_LB_GAMMA,4),
  pv_LB_BXIIARMA=round(pv_LB_BXIIARMA,4),
  pv_LB_ETS=round(pv_LB_ETS,4),
  pv_LB_ARMA=round(pv_LB_ARMA,4),
  pv_LB_SARMA=round(pv_LB_SARMA,4)
  )
rownames(results)<-nome

setwd("~/GitHub/BXII-GAS/Apply BXII-GAS")
saving<-paste0("app_hid_new.RData")
save.image(saving)
end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)

results[,c(4:7,11:14,18:21,25:28,32:35,39:42)]


resultados<-results
resultados[2,4]<-0
resultados<-round(resultados,4)
resultados[2,4]<-"_"
library(gt)

# Criar uma tabela estilizada
resultados[,] %>%
  gt()

#library(gt)

#results[results > 1.797693e+307] <- 0
#resultados<-round(results,4)

#resultados[c(12,14,15,16),c(7:15)] %>%
#  gt()


#results_out[c(12,14,15,16),] %>%
#  gt()
 