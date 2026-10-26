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
# GAS
library("gasmodel")
library(forecast)
library(lmtest)
library(xgboost)

# Main project folder
if(Sys.info()["sysname"]=="Windows"){
  path <- "C:/Users/Fernando/Documents/GitHub/BXII-GAS" # Windowns
}else{
  path <- "~/Fernando/BXII-GAS" # Ubuntu
}

setwd(file.path(path, "BXII-GAS"))
source("GASBXIIFit.R")
setwd(file.path(path, "Ray-GAS/Ray-GAS-master/Ray-GAS"))
source("GASRayFit.R")
setwd(file.path(path, "gamma-GAS"))
source("GASGFit.R")
setwd(file.path(path, "BXII_ARMA"))
source("bxiiarmaCOV.fit.R")
source("bxiiarma.fit.R")


############################
#### initial quantities ####
############################
w1<-8 # width for plots 
h11<-6 # height for plots

## Data
setwd(file.path(path, "Apply BXII-GAS/Data"))
Dados <- read.delim("Dados.txt")
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
#dados<-dados[,-c(11,21,30:33)]
#dados<-dados[,-c(7,9:11,14:15,21:22,24:25,27)]
#dados<-dados[,-c(3,5:9,13,15)]
#dados<-dados[,-c(3,4,7)]
#dados<-dados[,c(1,2,18,19,28)]
View(dados)

datas <- seq(as.Date("2010-01-01"), as.Date("2024-08-31"), by = "day")
R<-1:dim(dados)[2]

mat<-matrix(c(rep(1,1*4)),ncol=4,byrow=T)
mat<-matrix(c(rep(1,3),rep(2,3),rep(3,3),rep(1:3,3)),ncol=2,byrow=F)

# IMPORTANT: Candidate grid for the gamma
gamma_orders <- mat

acuracia1<-array(NA,c(7,length(R),5))

pv_LB_BXII<-pv_LB_Ray<-pv_LB_GAMMA<-pv_LB_BXIIARMA<-pv_LB_ETS<-pv_LB_SARMA<-pv_LB_XGBoost<-c()
nome<-modelo<-modelo_Ray<-modelo_GAMMA<-modeloBXIIARMA<-modelo_ETS<-modelo_XGBoost<-c()
final<-final_Ray<-final_GAMMA<-finalBXIIARMA<-final_ETS<-final_SARMA<-final_XGBoost<-list()
res_XGBoost<-pred_XGBoost<-list()

start_time <- Sys.time()
j<-0
pvalor<-0.10
tau<-0.25
setwd(file.path(path, "Apply BXII-GAS"))
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
  
  hidsin_l<-hidcos_l<-hidsincos_l<-list()
  
  aic_BXII_sin<-pv_LB_BXIIsin<-int_sig_BXIIsin<-model_sig_BXIIsin<-RMSE_BXIIsin<-MAPE_BXIIsin<-MAE_BXIIsin<-c()
  aic_BXII_cos<-pv_LB_BXIIcos<-int_sig_BXIIcos<-model_sig_BXIIcos<-RMSE_BXIIcos<-MAPE_BXIIcos<-MAE_BXIIcos<-c()
  aic_BXII_sincos<-pv_LB_BXIIsincos<-int_sig_BXIIsincos<-model_sig_BXIIsincos<-RMSE_BXIIsincos<-MAPE_BXIIsincos<-MAE_BXIIsincos<-c()
  
  for(l in 1:nrow(mat)){
    
    pq<-mat[l,]
    
    ## Sin
    print(c(l,"BXIIsin",1))
    hidsin<-try(BXIIGAS.fit(y,1:pq[1],1:pq[2],X=Xsin,X_hat=Xsin_hat,h1=1,tau = tau),silent=T)
    hidsin_l[[l]]<-hidsin
    
    ## Cos
    print(c(l,"BXIIcos",2))
    hidcos<-try(BXIIGAS.fit(y,1:pq[1],1:pq[2],X=Xcos,X_hat=Xcos_hat,h1=1,tau = tau),silent=T)
    hidcos_l[[l]]<-hidcos
    
    ## Sin&Cos
    print(c(l,"BXIIsincos",3))
    hidsincos<-try(BXIIGAS.fit(y,1:pq[1],1:pq[2],X=X,X_hat=X_hat,h1=1,tau = tau),silent=T)
    hidsincos_l[[l]]<-hidsincos
    
    ## Calculate model statistics
    for(tipo in c("sin","cos","sincos")){
      
      if(tipo=="sin"){
        fit<-hidsin
      }else if(tipo=="cos"){
        fit<-hidcos
      }else{
        fit<-hidsincos
      }
      
      if(inherits(fit,"try-error")){
        aic<-.Machine$double.xmax
        pv<-NA
        int_sig<-0
        model_sig<-0
        RMSE<-MAE<-MAPE<-NA
      }else{
        aic<-round(fit$aic,4)
        pv<-round(Box.test(fit$residuals,lag=20,type="Ljung")$p.value,4)
        int_sig<-as.numeric(fit$model[1,4]<pvalor)
        model_sig<-as.numeric(all(fit$model[2:nrow(fit$model),4]<pvalor))
        RMSE<-accuracy(y,fit$fitted)[,2]
        MAE<-accuracy(y,fit$fitted)[,3]
        MAPE<-accuracy(y,fit$fitted)[,5]
      }
      
      if(tipo=="sin"){
        aic_BXII_sin[l]<-aic
        pv_LB_BXIIsin[l]<-pv
        int_sig_BXIIsin[l]<-int_sig
        model_sig_BXIIsin[l]<-model_sig
        RMSE_BXIIsin[l]<-RMSE
        MAE_BXIIsin[l]<-MAE
        MAPE_BXIIsin[l]<-MAPE
      }else if(tipo=="cos"){
        aic_BXII_cos[l]<-aic
        pv_LB_BXIIcos[l]<-pv
        int_sig_BXIIcos[l]<-int_sig
        model_sig_BXIIcos[l]<-model_sig
        RMSE_BXIIcos[l]<-RMSE
        MAE_BXIIcos[l]<-MAE
        MAPE_BXIIcos[l]<-MAPE
      }else{
        aic_BXII_sincos[l]<-aic
        pv_LB_BXIIsincos[l]<-pv
        int_sig_BXIIsincos[l]<-int_sig
        model_sig_BXIIsincos[l]<-model_sig
        RMSE_BXIIsincos[l]<-RMSE
        MAE_BXIIsincos[l]<-MAE
        MAPE_BXIIsincos[l]<-MAPE
      }
    }
  }
  
  ## Results table
  Top_BXIIsin<-data.frame(
    Index=1:nrow(mat),AIC=aic_BXII_sin,pvalue_LB=pv_LB_BXIIsin,
    int_sig=int_sig_BXIIsin,model_sig=model_sig_BXIIsin,
    RMSE=RMSE_BXIIsin,MAPE=MAPE_BXIIsin,MAE=MAE_BXIIsin
  )
  
  Top_BXIIcos<-data.frame(
    Index=1:nrow(mat),AIC=aic_BXII_cos,pvalue_LB=pv_LB_BXIIcos,
    int_sig=int_sig_BXIIcos,model_sig=model_sig_BXIIcos,
    RMSE=RMSE_BXIIcos,MAPE=MAPE_BXIIcos,MAE=MAE_BXIIcos
  )
  
  Top_BXIIsincos<-data.frame(
    Index=1:nrow(mat),AIC=aic_BXII_sincos,pvalue_LB=pv_LB_BXIIsincos,
    int_sig=int_sig_BXIIsincos,model_sig=model_sig_BXIIsincos,
    RMSE=RMSE_BXIIsincos,MAPE=MAPE_BXIIsincos,MAE=MAE_BXIIsincos
  )
  
  Top_BXII<-rbind(
    data.frame(Type="Sin",Top_BXIIsin),
    data.frame(Type="Cos",Top_BXIIcos),
    data.frame(Type="SinCos",Top_BXIIsincos)
  )
  
  rownames(Top_BXII)<-NULL
  Top_BXII<-Top_BXII[order(Top_BXII$AIC),]
  
  ## Select the best model
  # 1. Selecionar os 10 melhores modelos pelo menor AIC
  Top_BXII <- Top_BXII %>%
    arrange(AIC) %>%
    slice_head(n = 10)

  # 2. Selecionar modelos com resíduos adequados e modelo significativo
  Select_BXII <- Top_BXII %>%
    filter(pvalue_LB > 0.05, model_sig == 1)

  # 3. Se houver modelos que atendam aos critérios,
  # selecionar aquele com menor RMSE
  if (nrow(Select_BXII) > 0) {
    Select_BXII <- Select_BXII %>%
      filter(RMSE == min(RMSE, na.rm = TRUE))
  } else {
    # 4. Se nenhum modelo atender aos critérios,
    # selecionar o modelo com o menor AIC
    Select_BXII <- Top_BXII %>%
      filter(AIC == min(AIC, na.rm = TRUE))
  }
  
  
  if(nrow(Select_BXII)>0){
    if(Select_BXII$Type[1]=="Sin"){
      final1<-hidsin_l[[Select_BXII$Index[1]]]
      modelo1<-"BXII Sin"
    }else if(Select_BXII$Type[1]=="Cos"){
      final1<-hidcos_l[[Select_BXII$Index[1]]]
      modelo1<-"BXII Cos"
    }else{
      final1<-hidsincos_l[[Select_BXII$Index[1]]]
      modelo1<-"BXII SinCos"
    }
  }
  
  ##########################
  #### fitting Ray-GAS  ####
  ##########################
  
  hidRaysin_l<-hidRaycos_l<-hidRaysincos_l<-list()
  
  aic_Ray_sin<-pv_LB_Raysin<-int_sig_Raysin<-model_sig_Raysin<-RMSE_Raysin<-MAPE_Raysin<-MAE_Raysin<-c()
  aic_Ray_cos<-pv_LB_Raycos<-int_sig_Raycos<-model_sig_Raycos<-RMSE_Raycos<-MAPE_Raycos<-MAE_Raycos<-c()
  aic_Ray_sincos<-pv_LB_Raysincos<-int_sig_Raysincos<-model_sig_Raysincos<-RMSE_Raysincos<-MAPE_Raysincos<-MAE_Raysincos<-c()
  
  for(l in 1:nrow(mat)){
    
    pq<-mat[l,]
    
    ## Sin
    print(c(l,"Raysin",1))
    hidsin<-try(RayGAS.fit(y,1:pq[1],1:pq[2],X=Xsin,X_hat=Xsin_hat,h1=1),silent=T)
    hidRaysin_l[[l]]<-hidsin
    
    ## Cos
    print(c(l,"Raycos",2))
    hidcos<-try(RayGAS.fit(y,1:pq[1],1:pq[2],X=Xcos,X_hat=Xcos_hat,h1=1),silent=T)
    hidRaycos_l[[l]]<-hidcos
    
    ## Sin&Cos
    print(c(l,"Raysincos",3))
    hidsincos<-try(RayGAS.fit(y,1:pq[1],1:pq[2],X=X,X_hat=X_hat,h1=1),silent=T)
    hidRaysincos_l[[l]]<-hidsincos
    
    ## Calculate model statistics
    for(tipo in c("sin","cos","sincos")){
      
      if(tipo=="sin"){
        fit<-hidsin
      }else if(tipo=="cos"){
        fit<-hidcos
      }else{
        fit<-hidsincos
      }
      
      if(inherits(fit,"try-error")){
        aic<-.Machine$double.xmax
        pv<-NA
        int_sig<-0
        model_sig<-0
        RMSE<-MAE<-MAPE<-NA
      }else{
        aic<-round(fit$aic,4)
        pv<-round(Box.test(fit$residuals,lag=20,type="Ljung")$p.value,4)
        int_sig<-as.numeric(fit$model[1,4]<pvalor)
        model_sig<-as.numeric(all(fit$model[2:nrow(fit$model),4]<pvalor))
        RMSE<-accuracy(y,fit$fitted)[,2]
        MAE<-accuracy(y,fit$fitted)[,3]
        MAPE<-accuracy(y,fit$fitted)[,5]
      }
      
      if(tipo=="sin"){
        aic_Ray_sin[l]<-aic
        pv_LB_Raysin[l]<-pv
        int_sig_Raysin[l]<-int_sig
        model_sig_Raysin[l]<-model_sig
        RMSE_Raysin[l]<-RMSE
        MAE_Raysin[l]<-MAE
        MAPE_Raysin[l]<-MAPE
      }else if(tipo=="cos"){
        aic_Ray_cos[l]<-aic
        pv_LB_Raycos[l]<-pv
        int_sig_Raycos[l]<-int_sig
        model_sig_Raycos[l]<-model_sig
        RMSE_Raycos[l]<-RMSE
        MAE_Raycos[l]<-MAE
        MAPE_Raycos[l]<-MAPE
      }else{
        aic_Ray_sincos[l]<-aic
        pv_LB_Raysincos[l]<-pv
        int_sig_Raysincos[l]<-int_sig
        model_sig_Raysincos[l]<-model_sig
        RMSE_Raysincos[l]<-RMSE
        MAE_Raysincos[l]<-MAE
        MAPE_Raysincos[l]<-MAPE
      }
    }
  }
  
  ## Results table
  Top_Raysin<-data.frame(
    Index=1:nrow(mat),AIC=aic_Ray_sin,pvalue_LB=pv_LB_Raysin,
    int_sig=int_sig_Raysin,model_sig=model_sig_Raysin,
    RMSE=RMSE_Raysin,MAPE=MAPE_Raysin,MAE=MAE_Raysin
  )
  
  Top_Raycos<-data.frame(
    Index=1:nrow(mat),AIC=aic_Ray_cos,pvalue_LB=pv_LB_Raycos,
    int_sig=int_sig_Raycos,model_sig=model_sig_Raycos,
    RMSE=RMSE_Raycos,MAPE=MAPE_Raycos,MAE=MAE_Raycos
  )
  
  Top_Raysincos<-data.frame(
    Index=1:nrow(mat),AIC=aic_Ray_sincos,pvalue_LB=pv_LB_Raysincos,
    int_sig=int_sig_Raysincos,model_sig=model_sig_Raysincos,
    RMSE=RMSE_Raysincos,MAPE=MAPE_Raysincos,MAE=MAE_Raysincos
  )
  
  Top_Ray<-rbind(
    data.frame(Type="Sin",Top_Raysin),
    data.frame(Type="Cos",Top_Raycos),
    data.frame(Type="SinCos",Top_Raysincos)
  )
  
  rownames(Top_Ray)<-NULL
  Top_Ray<-Top_Ray[order(Top_Ray$AIC),]
  
  ## Select the best model
  # 1. Selecionar os 10 melhores modelos pelo menor AIC
  Top_Ray <- Top_Ray %>%
    arrange(AIC) %>%
    slice_head(n = 10)

  # 2. Selecionar modelos com resíduos adequados e modelo significativo
  Select_Ray <- Top_Ray %>%
    filter(pvalue_LB > 0.05, model_sig == 1)

  # 3. Se houver modelos que atendam aos critérios,
  # selecionar aquele com menor RMSE
  if (nrow(Select_Ray) > 0) {
    Select_Ray <- Select_Ray %>%
      filter(RMSE == min(RMSE, na.rm = TRUE))
  } else {
    # 4. Se nenhum modelo atender aos critérios,
    # selecionar o modelo com o menor AIC
    Select_Ray <- Top_Ray %>%
      filter(AIC == min(AIC, na.rm = TRUE))
  }
  
  if(nrow(Select_Ray)>0){
    if(Select_Ray$Type[1]=="Sin"){
      final_Ray1<-hidRaysin_l[[Select_Ray$Index[1]]]
      modelo_Ray1<-"Ray Sin"
    }else if(Select_Ray$Type[1]=="Cos"){
      final_Ray1<-hidRaycos_l[[Select_Ray$Index[1]]]
      modelo_Ray1<-"Ray Cos"
    }else{
      final_Ray1<-hidRaysincos_l[[Select_Ray$Index[1]]]
      modelo_Ray1<-"Ray SinCos"
    }
  }
  
  ###########################
  #### fitting gamma-GAS ####
  ###########################
  
  hidgammasin_l<-hidgammacos_l<-hidgammasincos_l<-list()
  
  aic_gamma_sin<-pv_LB_gammasin<-int_sig_gammasin<-model_sig_gammasin<-RMSE_gammasin<-MAPE_gammasin<-MAE_gammasin<-c()
  aic_gamma_cos<-pv_LB_gammacos<-int_sig_gammacos<-model_sig_gammacos<-RMSE_gammacos<-MAPE_gammacos<-MAE_gammacos<-c()
  aic_gamma_sincos<-pv_LB_gammasincos<-int_sig_gammasincos<-model_sig_gammasincos<-RMSE_gammasincos<-MAPE_gammasincos<-MAE_gammasincos<-c()
  
  for(l in 1:nrow(mat)){
    
    pq<-mat[l,]
    
    ## Sin
    print(c(l,"gammasin",1))
    hidsin<-try(GGAS.fit(y,pq[1],pq[2],X=Xsin,X_hat=Xsin_hat,h1=1),silent=T)
    hidgammasin_l[[l]]<-hidsin
    
    ## Cos
    print(c(l,"gammacos",2))
    hidcos<-try(GGAS.fit(y,pq[1],pq[2],X=Xcos,X_hat=Xcos_hat,h1=1),silent=T)
    hidgammacos_l[[l]]<-hidcos
    
    ## Sin&Cos
    print(c(l,"gammasincos",3))
    hidsincos<-try(GGAS.fit(y,pq[1],pq[2],X=X,X_hat=X_hat,h1=1),silent=T)
    hidgammasincos_l[[l]]<-hidsincos
    
    ## Calculate model statistics
    for(tipo in c("sin","cos","sincos")){
      
      if(tipo=="sin"){
        fit<-hidsin
      }else if(tipo=="cos"){
        fit<-hidcos
      }else{
        fit<-hidsincos
      }
      
      if(inherits(fit,"try-error")){
        aic<-.Machine$double.xmax
        pv<-NA
        int_sig<-0
        model_sig<-0
        RMSE<-MAE<-MAPE<-NA
      }else{
        aic<-round(fit$aic,4)
        pv<-round(Box.test(fit$residuals,lag=20,type="Ljung")$p.value,4)
        int_sig<-as.numeric(fit$model[1,4]<pvalor)
        model_sig<-as.numeric(all(fit$model[2:nrow(fit$model),4]<pvalor))
        RMSE<-accuracy(y,fit$fitted)[,2]
        MAE<-accuracy(y,fit$fitted)[,3]
        MAPE<-accuracy(y,fit$fitted)[,5]
      }
      
      if(tipo=="sin"){
        aic_gamma_sin[l]<-aic
        pv_LB_gammasin[l]<-pv
        int_sig_gammasin[l]<-int_sig
        model_sig_gammasin[l]<-model_sig
        RMSE_gammasin[l]<-RMSE
        MAE_gammasin[l]<-MAE
        MAPE_gammasin[l]<-MAPE
      }else if(tipo=="cos"){
        aic_gamma_cos[l]<-aic
        pv_LB_gammacos[l]<-pv
        int_sig_gammacos[l]<-int_sig
        model_sig_gammacos[l]<-model_sig
        RMSE_gammacos[l]<-RMSE
        MAE_gammacos[l]<-MAE
        MAPE_gammacos[l]<-MAPE
      }else{
        aic_gamma_sincos[l]<-aic
        pv_LB_gammasincos[l]<-pv
        int_sig_gammasincos[l]<-int_sig
        model_sig_gammasincos[l]<-model_sig
        RMSE_gammasincos[l]<-RMSE
        MAE_gammasincos[l]<-MAE
        MAPE_gammasincos[l]<-MAPE
      }
    }
  }
  
  ## Results table
  Top_gammasin<-data.frame(
    Index=1:nrow(mat),AIC=aic_gamma_sin,pvalue_LB=pv_LB_gammasin,
    int_sig=int_sig_gammasin,model_sig=model_sig_gammasin,
    RMSE=RMSE_gammasin,MAPE=MAPE_gammasin,MAE=MAE_gammasin
  )
  
  Top_gammacos<-data.frame(
    Index=1:nrow(mat),AIC=aic_gamma_cos,pvalue_LB=pv_LB_gammacos,
    int_sig=int_sig_gammacos,model_sig=model_sig_gammacos,
    RMSE=RMSE_gammacos,MAPE=MAPE_gammacos,MAE=MAE_gammacos
  )
  
  Top_gammasincos<-data.frame(
    Index=1:nrow(mat),AIC=aic_gamma_sincos,pvalue_LB=pv_LB_gammasincos,
    int_sig=int_sig_gammasincos,model_sig=model_sig_gammasincos,
    RMSE=RMSE_gammasincos,MAPE=MAPE_gammasincos,MAE=MAE_gammasincos
  )
  
  Top_gamma<-rbind(
    data.frame(Type="Sin",Top_gammasin),
    data.frame(Type="Cos",Top_gammacos),
    data.frame(Type="SinCos",Top_gammasincos)
  )
  
  rownames(Top_gamma)<-NULL
  Top_gamma<-Top_gamma[order(Top_gamma$AIC),]
  
  ## Select the best model
  
  # 1. Selecionar os 10 melhores modelos pelo menor AIC
  Top_gamma <- Top_gamma %>%
    arrange(AIC) %>%
    slice_head(n = 10)

  # 2. Selecionar modelos com resíduos adequados e modelo significativo
  Select_gamma <- Top_gamma %>%
    filter(pvalue_LB > 0.05, model_sig == 1)

  # 3. Se houver modelos que atendam aos critérios,
  # selecionar aquele com menor RMSE
  if (nrow(Select_gamma) > 0) {
    Select_gamma <- Select_gamma %>%
      filter(RMSE == min(RMSE, na.rm = TRUE))
  } else {
    # 4. Se nenhum modelo atender aos critérios,
    # selecionar o modelo com o menor AIC
    Select_gamma <- Top_gamma %>%
      filter(AIC == min(AIC, na.rm = TRUE))
  }
  
  
  if(nrow(Select_gamma)>0){
    if(Select_gamma$Type[1]=="Sin"){
      final_GAMMA1<-hidgammasin_l[[Select_gamma$Index[1]]]
      modelo_GAMMA1<-"gamma Sin"
    }else if(Select_gamma$Type[1]=="Cos"){
      final_GAMMA1<-hidgammacos_l[[Select_gamma$Index[1]]]
      modelo_GAMMA1<-"gamma Cos"
    }else{
      final_GAMMA1<-hidgammasincos_l[[Select_gamma$Index[1]]]
      modelo_GAMMA1<-"gamma SinCos"
    }
  }
  
  ###########################
  #### fitting BXII-ARMA ####
  ###########################
  
  hidBXIIARMAsin_l<-hidBXIIARMAcos_l<-hidBXIIARMAsincos_l<-list()
  
  aic_BXIIARMA_sin<-pv_LB_BXIIARMAsin<-int_sig_BXIIARMAsin<-model_sig_BXIIARMAsin<-RMSE_BXIIARMAsin<-MAPE_BXIIARMAsin<-MAE_BXIIARMAsin<-c()
  aic_BXIIARMA_cos<-pv_LB_BXIIARMAcos<-int_sig_BXIIARMAcos<-model_sig_BXIIARMAcos<-RMSE_BXIIARMAcos<-MAPE_BXIIARMAcos<-MAE_BXIIARMAcos<-c()
  aic_BXIIARMA_sincos<-pv_LB_BXIIARMAsincos<-int_sig_BXIIARMAsincos<-model_sig_BXIIARMAsincos<-RMSE_BXIIARMAsincos<-MAPE_BXIIARMAsincos<-MAE_BXIIARMAsincos<-c()
  
  for(l in 1:nrow(mat)){
    
    pq<-mat[l,]
    
    ## Sin
    print(c(l,"BXIIARMAsin",1))
    hidsin<-try(bxiiarmaCOV.fit(y,ar=1:pq[1],ma=1:pq[2],X = as.matrix(Xsin),X_hat=as.matrix(Xsin_hat),h=1,diag=0,tau=tau),silent=T)
    hidBXIIARMAsin_l[[l]]<-hidsin
    
    ## Cos
    print(c(l,"BXIIARMAcos",2))
    hidcos<-try(bxiiarmaCOV.fit(y,ar=1:pq[1],ma=1:pq[2],X = as.matrix(Xcos),X_hat=as.matrix(Xcos_hat),h=1,diag=0,tau=tau),silent=T)
    hidBXIIARMAcos_l[[l]]<-hidcos
    
    ## Sin&Cos
    print(c(l,"BXIIARMAsincos",3))
    hidsincos<-try(bxiiarmaCOV.fit(y,ar=1:pq[1],ma=1:pq[2],X = as.matrix(X),X_hat=as.matrix(X_hat),h=1,diag=0,tau=tau),silent=T)
    hidBXIIARMAsincos_l[[l]]<-hidsincos
    
    ## Calculate model statistics
    for(tipo in c("sin","cos","sincos")){
      
      if(tipo=="sin"){
        fit<-hidsin
      }else if(tipo=="cos"){
        fit<-hidcos
      }else{
        fit<-hidsincos
      }
      
      if(inherits(fit,"try-error")){
        aic<-.Machine$double.xmax
        pv<-NA
        int_sig<-0
        model_sig<-0
        RMSE<-MAE<-MAPE<-NA
      }else{
        aic<-round(fit$aic,4)
        pv<-round(Box.test(fit$residuals,lag=20,type="Ljung")$p.value,4)
        int_sig<-as.numeric(fit$model[1,4]<pvalor)
        model_sig<-as.numeric(all(fit$model[2:nrow(fit$model),4]<pvalor))
        RMSE<-accuracy(y,fit$fitted)[,2]
        MAE<-accuracy(y,fit$fitted)[,3]
        MAPE<-accuracy(y,fit$fitted)[,5]
      }
      
      if(tipo=="sin"){
        aic_BXIIARMA_sin[l]<-aic
        pv_LB_BXIIARMAsin[l]<-pv
        int_sig_BXIIARMAsin[l]<-int_sig
        model_sig_BXIIARMAsin[l]<-model_sig
        RMSE_BXIIARMAsin[l]<-RMSE
        MAE_BXIIARMAsin[l]<-MAE
        MAPE_BXIIARMAsin[l]<-MAPE
      }else if(tipo=="cos"){
        aic_BXIIARMA_cos[l]<-aic
        pv_LB_BXIIARMAcos[l]<-pv
        int_sig_BXIIARMAcos[l]<-int_sig
        model_sig_BXIIARMAcos[l]<-model_sig
        RMSE_BXIIARMAcos[l]<-RMSE
        MAE_BXIIARMAcos[l]<-MAE
        MAPE_BXIIARMAcos[l]<-MAPE
      }else{
        aic_BXIIARMA_sincos[l]<-aic
        pv_LB_BXIIARMAsincos[l]<-pv
        int_sig_BXIIARMAsincos[l]<-int_sig
        model_sig_BXIIARMAsincos[l]<-model_sig
        RMSE_BXIIARMAsincos[l]<-RMSE
        MAE_BXIIARMAsincos[l]<-MAE
        MAPE_BXIIARMAsincos[l]<-MAPE
      }
    }
  }
  
  ## Results table
  Top_BXIIARMAsin<-data.frame(
    Index=1:nrow(mat),AIC=aic_BXIIARMA_sin,pvalue_LB=pv_LB_BXIIARMAsin,
    int_sig=int_sig_BXIIARMAsin,model_sig=model_sig_BXIIARMAsin,
    RMSE=RMSE_BXIIARMAsin,MAPE=MAPE_BXIIARMAsin,MAE=MAE_BXIIARMAsin
  )
  
  Top_BXIIARMAcos<-data.frame(
    Index=1:nrow(mat),AIC=aic_BXIIARMA_cos,pvalue_LB=pv_LB_BXIIARMAcos,
    int_sig=int_sig_BXIIARMAcos,model_sig=model_sig_BXIIARMAcos,
    RMSE=RMSE_BXIIARMAcos,MAPE=MAPE_BXIIARMAcos,MAE=MAE_BXIIARMAcos
  )
  
  Top_BXIIARMAsincos<-data.frame(
    Index=1:nrow(mat),AIC=aic_BXIIARMA_sincos,pvalue_LB=pv_LB_BXIIARMAsincos,
    int_sig=int_sig_BXIIARMAsincos,model_sig=model_sig_BXIIARMAsincos,
    RMSE=RMSE_BXIIARMAsincos,MAPE=MAPE_BXIIARMAsincos,MAE=MAE_BXIIARMAsincos
  )
  
  Top_BXIIARMA<-rbind(
    data.frame(Type="Sin",Top_BXIIARMAsin),
    data.frame(Type="Cos",Top_BXIIARMAcos),
    data.frame(Type="SinCos",Top_BXIIARMAsincos)
  )
  
  rownames(Top_BXIIARMA)<-NULL
  Top_BXIIARMA<-Top_BXIIARMA[order(Top_BXIIARMA$AIC),]
  
  ## Select the best model
    
  # 1. Selecionar os 10 melhores modelos pelo menor AIC
  Top_BXIIARMA <- Top_BXIIARMA %>%
    arrange(AIC) %>%
    slice_head(n = 10)

  # 2. Selecionar modelos com resíduos adequados e modelo significativo
  Select_BXIIARMA <- Top_BXIIARMA %>%
    filter(pvalue_LB > 0.05, model_sig == 1)

  # 3. Se houver modelos que atendam aos critérios,
  # selecionar aquele com menor RMSE
  if (nrow(Select_BXIIARMA) > 0) {
    Select_BXIIARMA <- Select_BXIIARMA %>%
      filter(RMSE == min(RMSE, na.rm = TRUE))
  } else {
    # 4. Se nenhum modelo atender aos critérios,
    # selecionar o modelo com o menor AIC
    Select_BXIIARMA <- Top_BXIIARMA %>%
      filter(AIC == min(AIC, na.rm = TRUE))
  }
  
  if(nrow(Select_BXIIARMA)>0){
    if(Select_BXIIARMA$Type[1]=="Sin"){
      finalBXIIARMA1<-hidBXIIARMAsin_l[[Select_BXIIARMA$Index[1]]]
      modeloBXIIARMA1<-"BXIIARMA Sin"
    }else if(Select_BXIIARMA$Type[1]=="Cos"){
      finalBXIIARMA1<-hidBXIIARMAcos_l[[Select_BXIIARMA$Index[1]]]
      modeloBXIIARMA1<-"BXIIARMA Cos"
    }else{
      finalBXIIARMA1<-hidBXIIARMAsincos_l[[Select_BXIIARMA$Index[1]]]
      modeloBXIIARMA1<-"BXIIARMA SinCos"
    }
  }
  
  #############################
  ###     classic models    ###
  #############################
  
  ##########################
  #### fitting ETS #########
  ##########################
  
  ets_models <- c("ANN","AAN","AAdN","AAA","AAdA",
                  "MNN","MAN","MAdN","MNA","MAM","MAdA","ZZZ")
  
  ets_l <- list()
  aic_ETS <- RMSE_ETS <- MAE_ETS <- MAPE_ETS <- c()
  
  for(l in 1:length(ets_models)){
    
    print(c(l,ets_models[l]))
    
    fit <- try(
      ets(y, model=ets_models[l]),
      silent=T
    )
    
    ets_l[[l]] <- fit
    
    if(inherits(fit,"try-error")){
      aic_ETS[l] <- Inf
      RMSE_ETS[l] <- MAE_ETS[l] <- MAPE_ETS[l] <- NA
    }else{
      aic_ETS[l] <- round(fit$aic,4)
      
      acc <- accuracy(fit)
      
      RMSE_ETS[l] <- acc[1,"RMSE"]
      MAE_ETS[l] <- acc[1,"MAE"]
      MAPE_ETS[l] <- acc[1,"MAPE"]
    }
  }
  
  ## Results table
  Top_ETS <- data.frame(
    Index=1:length(ets_models),
    Model=ets_models,
    AIC=aic_ETS,
    RMSE=RMSE_ETS,
    MAE=MAE_ETS,
    MAPE=MAPE_ETS
  )
  
  Top_ETS <- Top_ETS[order(Top_ETS$AIC),]
  
  ## Select the model with the lowest AIC
  
  if(all(is.infinite(aic_ETS))){
    final_ETS1 <- NULL
    modelo_ETS1 <- NA
    warning("All ETS models failed.")
    }else{
      id_ETS <- which.min(aic_ETS)
      final_ETS1 <- ets_l[[id_ETS]]
      modelo_ETS1 <- ets_models[id_ETS]
      }
  
  ##########################
  #### fitting SARMA #######
  ##########################
  
  sarma_l<-list()
  
  aic_SARMA<-pv_LB_SARMAmodel<-int_sig_SARMA<-model_sig_SARMA<-
    RMSE_SARMA<-MAPE_SARMA<-MAE_SARMA<-c()
  
  matsarma<-matrix(c(
    1,1,0,1,
    1,1,1,0,
    1,1,1,1,
    1,2,0,1,
    1,2,1,0,
    1,2,1,1,
    1,3,0,1,
    1,3,1,0,
    1,3,1,1,
    2,1,0,1,
    2,1,1,0,
    2,1,1,1,
    2,2,0,1,
    2,2,1,0,
    2,2,1,1,
    2,3,0,1,
    2,3,1,0,
    2,3,1,1,
    3,1,0,1,
    3,1,1,0,
    3,1,1,1,
    3,2,0,1,
    3,2,1,0,
    3,2,1,1,
    3,3,0,1,
    3,3,1,0,
    3,3,1,1
  ),ncol=4,byrow=T)
  
  for(l in 1:nrow(matsarma)){
    
    pqs<-matsarma[l,]
    print(c(l,"SARMA"))
    
    fit<-try(
      Arima(
        y,
        order=c(pqs[1],0,pqs[2]),
        seasonal=list(order=c(pqs[3],0,pqs[4]),period=s)
      ),
      silent=T
    )
    
    sarma_l[[l]]<-fit
    
    if(inherits(fit,"try-error")){
      aic_SARMA[l]<-Inf
      pv_LB_SARMAmodel[l]<-NA
      int_sig_SARMA[l]<-0
      model_sig_SARMA[l]<-0
      RMSE_SARMA[l]<-MAE_SARMA[l]<-MAPE_SARMA[l]<-NA
      
    }else{
      
      aic_SARMA[l]<-round(fit$aic,4)
      pv_LB_SARMAmodel[l]<-round(Box.test(residuals(fit),lag=20,type="Ljung")$p.value,4)
      int_sig_SARMA[l]<-as.numeric(coeftest(fit)[1,4]<pvalor)
      model_sig_SARMA[l]<-as.numeric(all(coeftest(fit)[2:nrow(coeftest(fit)),4]<pvalor))
      acc<-accuracy(fit)
      RMSE_SARMA[l]<-acc[1,"RMSE"]
      MAE_SARMA[l]<-acc[1,"MAE"]
      MAPE_SARMA[l]<-acc[1,"MAPE"]
    }
  }
  
  ## Results table
  Top_SARMA<-data.frame(
    Index=1:nrow(matsarma),
    AIC=aic_SARMA,
    pvalue_LB=pv_LB_SARMAmodel,
    int_sig=int_sig_SARMA,
    model_sig=model_sig_SARMA,
    RMSE=RMSE_SARMA,
    MAPE=MAPE_SARMA,
    MAE=MAE_SARMA
  )
  
  Top_SARMA<-Top_SARMA[order(Top_SARMA$AIC),]

  ## Select the best model
    
  # 1. Selecionar os 10 melhores modelos pelo menor AIC
  Top_SARMA <- Top_SARMA %>%
    arrange(AIC) %>%
    slice_head(n = 10)

  # 2. Selecionar modelos com resíduos adequados e modelo significativo
  Select_SARMA <- Top_SARMA %>%
    filter(pvalue_LB > 0.05, model_sig == 1)

  # 3. Se houver modelos que atendam aos critérios,
  # selecionar aquele com menor RMSE
  if (nrow(Select_SARMA) > 0) {
    Select_SARMA <- Select_SARMA %>%
      filter(RMSE == min(RMSE, na.rm = TRUE))
  } else {
    # 4. Se nenhum modelo atender aos critérios,
    # selecionar o modelo com o menor AIC
    Select_SARMA <- Top_SARMA %>%
      filter(AIC == min(AIC, na.rm = TRUE))
  }
  
  if(nrow(Select_SARMA)>0){
    final_SARMA1<-sarma_l[[Select_SARMA$Index[1]]]
    }else{
      final_SARMA1<-NULL
      warning("No SARMA model passed the Ljung-Box criterion.")
      }

  #########################
  #### Fitting XGBoost ####
  #########################

  # 1. Construct temporal covariates
  n <- length(y)
  Index <- seq_len(n)

  X_XGBoost <- data.frame(
    Index = Index,

    # Lags
    lag1  = c(NA, y[-n]),
    lag2  = c(NA, NA, y[seq_len(n - 2)]),
    lag3 = c(rep(NA, 3), y[seq_len(n - 3)]),

    # Moving averages, using only past observations
    MA2 = c(
      rep(NA, 2),
      vapply(3:n, function(t)
        mean(y[(t - 2):(t - 1)]), numeric(1))
    ),

    MA3 = c(
      rep(NA, 3),
      vapply(4:n, function(t)
        mean(y[(t - 3):(t - 1)]), numeric(1))
    ),

    # Annual seasonality
    Xsin = Xsin,
    Xcos = Xcos
  )

  dados_XGBoost <- data.frame(
    y = as.numeric(y),
    X_XGBoost
  )

  # Keep the original time index for each observation
  dados_XGBoost$Index_original <- dados_XGBoost$Index

  dados_XGBoost <- na.omit(dados_XGBoost)

  y_XGBoost <- dados_XGBoost$y

  X_XGBoost <- as.matrix(
    dados_XGBoost[, c(
      "Index", "lag1", "lag2", "lag3",
      "MA2", "MA3", "Xsin", "Xcos"
    )]
  )

  ############################
  #### Candidate models ####
  ############################

  # Hyperparameter combinations
  mat_XGBoost <- expand.grid(
    max_depth = c(2, 4, 6),
    eta = c(0.01, 0.05, 0.10),
    nrounds = c(100, 300)
  )

  hids_XGBoost_l <- list()

  RMSE_XGBoost <- rep(NA_real_, nrow(mat_XGBoost))
  MAE_XGBoost  <- rep(NA_real_, nrow(mat_XGBoost))
  MAPE_XGBoost <- rep(NA_real_, nrow(mat_XGBoost))
  pv_LB_XGBoost <- rep(NA_real_, nrow(mat_XGBoost))

  # Expanding-window validation
  n_obs <- length(y_XGBoost)
  initial <- floor(0.60 * n_obs)
  n_folds <- 5

  if (initial < 20 || n_obs - initial < n_folds) {
    stop("Not enough observations for the specified temporal validation.")
  }

  folds <- split(
    seq.int(initial + 1, n_obs),
    cut(
      seq.int(initial + 1, n_obs),
      breaks = n_folds,
      labels = FALSE
    )
  )

  for (l in seq_len(nrow(mat_XGBoost))) {

    pq1 <- mat_XGBoost[l, ]

    print(c(l, "XGBoost"))

    pred_cv <- rep(NA_real_, n_obs)

    # Rolling-origin cross-validation
    for (k in seq_along(folds)) {

      id_valid <- folds[[k]]
      id_train <- seq_len(min(id_valid) - 1)

      fit <- try(
        xgboost(
          data = X_XGBoost[id_train, , drop = FALSE],
          label = y_XGBoost[id_train],
          objective = "reg:squarederror",
          max_depth = pq1$max_depth,
          eta = pq1$eta,
          nrounds = pq1$nrounds,
          subsample = 0.8,
          colsample_bytree = 0.8,
          nthread = 2,
          verbosity = 0
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        pred_cv[id_valid] <- predict(
          fit,
          X_XGBoost[id_valid, , drop = FALSE]
        )
      }
    }

    # Calculate statistics from out-of-fold predictions
    valid <- which(!is.na(pred_cv))

    if (length(valid) > 0) {

      res <- y_XGBoost[valid] - pred_cv[valid]

      RMSE_XGBoost[l] <- sqrt(mean(res^2))
      MAE_XGBoost[l] <- mean(abs(res))

      if (all(y_XGBoost[valid] != 0)) {
        MAPE_XGBoost[l] <- mean(
          abs(res / y_XGBoost[valid])
        ) * 100
      }

      if (length(res) > 20) {
        pv_LB_XGBoost[l] <- Box.test(
          res,
          lag = 20,
          type = "Ljung-Box"
        )$p.value
      }
    }

    # Fit candidate using the complete modelling sample
    hids_XGBoost_l[[l]] <- try(
      xgboost(
        data = X_XGBoost,
        label = y_XGBoost,
        objective = "reg:squarederror",
        max_depth = pq1$max_depth,
        eta = pq1$eta,
        nrounds = pq1$nrounds,
        subsample = 0.8,
        colsample_bytree = 0.8,
        nthread = 2,
        verbosity = 0
      ),
      silent = TRUE
    )
  }

  ############################
  #### Results table ####
  ############################

  Top_XGBoost <- data.frame(
    Index = seq_len(nrow(mat_XGBoost)),
    max_depth = mat_XGBoost$max_depth,
    eta = mat_XGBoost$eta,
    nrounds = mat_XGBoost$nrounds,
    pvalue_LB = pv_LB_XGBoost,
    RMSE = RMSE_XGBoost,
    MAPE = MAPE_XGBoost,
    MAE = MAE_XGBoost
  )

  # Keep models that fitted successfully and have a valid RMSE
  valid_models <- vapply(
    hids_XGBoost_l,
    function(fit) !inherits(fit, "try-error"),
    logical(1)
  )

  Top_XGBoost <- Top_XGBoost %>%
    filter(valid_models, is.finite(RMSE)) %>%
    arrange(RMSE) %>%
    slice_head(n = 10)

  ############################
  #### Select best model ####
  ############################

  # Prefer models with uncorrelated forecast errors
  Select_XGBoost <- Top_XGBoost %>%
    filter(pvalue_LB > 0.05)

  if (nrow(Select_XGBoost) > 0) {

    Select_XGBoost <- Select_XGBoost %>%
      filter(RMSE == min(RMSE, na.rm = TRUE))

  } else {

    # Fallback: lowest cross-validation RMSE
    Select_XGBoost <- Top_XGBoost %>%
      filter(RMSE == min(RMSE, na.rm = TRUE))
  }

  ############################
  #### Recover selected model ####
  ############################

  if (nrow(Select_XGBoost) > 0) {

    id_XGBoost <- Select_XGBoost$Index[1]

    final_XGBoost1 <- hids_XGBoost_l[[id_XGBoost]]

    modelo_XGBoost1 <- "XGBoost"
  }

  ############################
  #### Covariate importance ####
  ############################

  if (exists("final_XGBoost")) {

    imp_XGBoost <- xgb.importance(
      model = final_XGBoost1
    )

    print(imp_XGBoost)

    xgb.plot.importance(
      importance_matrix = imp_XGBoost,
      measure = "Gain",
      top_n = 15
    )
  }

  
  
  ###########################
  #### residual analysis ####  
  ###########################
  
  pv_LB_BXII1<-Box.test(final1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_Ray1<-Box.test(final_Ray1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_GAMMA1 <- try(
    Box.test(final_GAMMA1$residuals, lag = 20, type = "Ljung")$p.value,
    silent = TRUE
  )
  
  if (inherits(pv_LB_GAMMA1, "try-error")) {
    pv_LB_GAMMA1 <- NA
  }
  pv_LB_BXIIARMA1<-Box.test(finalBXIIARMA1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_ETS1<-Box.test(final_ETS1$residuals, lag = 20, type = "Ljung")$p.value
  pv_LB_SARMA1<-Box.test(final_SARMA1$residuals, lag = 20, type = "Ljung")$p.value
  # Predicted values
  pred_XGBoost1 <- predict(
    final_XGBoost1,
    X_XGBoost
  )
  # Residuals
  res_XGBoost1 <- y_XGBoost - pred_XGBoost1
  pv_LB_XGBoost1<-Box.test(res_XGBoost1, lag = 20, type = "Ljung")$p.value

  
  #############################
  #### selecting resevoirs ####  
  #############################

    j<-j+1
    final[[j]]<-final1
    final_Ray[[j]]<-final_Ray1
    final_GAMMA[[j]]<-final_GAMMA1
    finalBXIIARMA[[j]]<-finalBXIIARMA1
    final_ETS[[j]]<-final_ETS1
    final_SARMA[[j]]<-final_SARMA1
    final_XGBoost[[j]]<-final_XGBoost1
    modelo[j]<-modelo1
    modelo_Ray[j]<-modelo_Ray1
    modelo_GAMMA[j]<-modelo_GAMMA1
    modeloBXIIARMA[j]<-modeloBXIIARMA1
    modelo_ETS[j]<-modelo_ETS1
    modelo_XGBoost[j]<-modelo_XGBoost1
    nome[j]<-substring(colnames(dados)[i],1,nchar(colnames(dados)[i])#-10
                       )
    print(nome[j])
    pv_LB_BXII[j]<- pv_LB_BXII1
    pv_LB_Ray[j]<- pv_LB_Ray1
    pv_LB_GAMMA[j] <- try(
      pv_LB_GAMMA1,
      silent = TRUE
    )
    
    if (inherits(pv_LB_GAMMA[j], "try-error")) {
      pv_LB_GAMMA[j] <- NA
    }
    pv_LB_BXIIARMA[j]<-pv_LB_BXIIARMA1
    pv_LB_ETS[j]<-pv_LB_ETS1
    pv_LB_SARMA[j]<-pv_LB_SARMA1
    pv_LB_XGBoost[j]<-pv_LB_XGBoost1


    acuracia1[1,j,]<-c(accuracy(final[[j]]$fitted, y)[,c(2,3,5)],final[[j]]$aic,final[[j]]$bic)
    acuracia1[2,j,]<-c(accuracy(final_Ray[[j]]$fitted, y)[,c(2,3,5)],final_Ray[[j]]$aic,final_Ray[[j]]$bic)
    if (inherits(final_GAMMA[[j]], "try-error")){
      acuracia1[3, j, ] <- NA
    }else{
      acuracia1[3, j, ] <- c(accuracy(final_GAMMA[[j]]$fitted, y)[, c(2,3,5)],
      final_GAMMA[[j]]$aic,
      final_GAMMA[[j]]$bic)
    }
    acuracia1[4,j,]<-c(accuracy(finalBXIIARMA[[j]]$fitted, y)[,c(2,3,5)],finalBXIIARMA[[j]]$aic,finalBXIIARMA[[j]]$bic)
    acuracia1[5,j,]<-c(accuracy(final_ETS[[j]]$fitted, y)[,c(2,3,5)],final_ETS[[j]]$aic,final_ETS[[j]]$bic)
    acuracia1[6,j,]<-c(accuracy(final_SARMA[[j]]$fitted, y)[,c(2,3,5)],final_SARMA[[j]]$aic,final_SARMA[[j]]$bic)
    res_XGBoost[[j]]<-res_XGBoost1
    pred_XGBoost[[j]]<-pred_XGBoost1
    acuracia1[7,j,]<-c(accuracy(pred_XGBoost[[j]], y_XGBoost)[,c(2,3,5)],NA,NA)

    
    ###############
    #### plots ####  
    ###############
    setwd(file.path(path, "Apply BXII-GAS"))
    
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
    lines(final_Ray[[j]]$fitted, col=5,lty=3,lwd=1.7)
    # lines(ts(final_GAMMA[[j]]$fitted,end = c(2023, 8),frequency = 12), col=8,lty=3.5,lwd=1.5)
    lines(finalBXIIARMA[[j]]$fitted, col=6,lty=4,lwd=1.6)
    lines(final_SARMA[[j]]$fitted, col=7,lty=5,lwd=1.8)
    lines(ts(pred_XGBoost[[j]],end = c(2023, 8),frequency = 12), col=4,lty=7,lwd=1.4)
    lines(final_ETS[[j]]$fitted, col=3,lty=6,lwd=1.9)
    lines(final[[j]]$fitted, col=2,lty=2,lwd=2)
    legend("topleft", 
           c("Original","BXII-GAS","Ray-GAS",# "gamma-GAS",
           "BXII-ARMA","SARMA",
             "ETS","XGBoost" 
             ),
           col = c(1,2,5,# 8,
            6,7,3,4),
           lty= c(1,2,3,# 3.5,
            4,5,6,7),
           lwd = c(1,2,1.7,# 1.5,
            1.6,1.8,1.9,1.4), bty="n", cex = 1)
    dev.off()
    # ---------------------------------
    # seasonality
    months<-paste0("Plots/months",i,".pdf")
    pdf(months,width = 10, height = 6)
    par(mfrow=c(1,1))
    monthplot(y,main=nome[j],ylab = "RF",base = "median")
    dev.off()
    
    # ---------------------------------
    if(sum(final[[j]]$residuals)!=(-Inf)){
    resid_plot<-paste0("Plots/resid_plotBXII",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf BXII residuals 
    acf(final[[j]]$residuals,main="")
    dev.off()
    }else{"BXII-GAS: It is not possible to create the chart."}
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotRay",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf Ray residuals 
    acf(final_Ray[[j]]$residuals,main="")
    dev.off()
    
    # ---------------------------------
    
    if (inherits(final_GAMMA[[j]], "try-error")) {
      print("gamma-GAS: It is not possible to create the chart.")
    }else{pdf(resid_plot, width = 4, height = 4)
      resid_plot <- paste0("Plots/resid_plotgamma", i, ".pdf")
      acf(final_GAMMA[[j]]$residuals, main = "")
      dev.off()
    }
    
    # ---------------------------------
    if(sum(finalBXIIARMA[[j]]$residuals)!=(-Inf)){
    resid_plot<-paste0("Plots/resid_plotBXIIARMA",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf BXIIARMA residuals 
    acf(finalBXIIARMA[[j]]$residuals,main="")
    dev.off()
    }else{"BXII-ARMA: It is not possible to create the chart."}
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotETS",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf ETS residuals 
    acf(final_ETS[[j]]$residuals,main="")
    dev.off()
    
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotSARMA",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf SARMA residuals 
    acf(final_SARMA[[j]]$residuals,main="")
    dev.off()
  
    # ---------------------------------
    resid_plot<-paste0("Plots/resid_plotXGBoost",i,".pdf")
    pdf(resid_plot,width = 4, height = 4)
    # acf SARMA residuals 
    acf(ts(res_XGBoost[[j]],end = c(2023, 8),frequency = 12),main="")
    dev.off()
    }

########################
#### final measures ####  
########################

results<-data.frame(
  AIC_BXII=acuracia1[1,,4],
  AIC_Ray=acuracia1[2,,4],
  AIC_GAMMA=acuracia1[3,,4],
  AIC_BXIIARMA=acuracia1[4,,4],
  AIC_ETS=acuracia1[5,,4],
  AIC_SARMA=acuracia1[6,,4],
  AIC_XGBoost=c(rep(NA,33)),
  BIC_BXII=acuracia1[1,,5],
  BIC_Ray=acuracia1[2,,5],
  BIC_GAMMA=acuracia1[3,,5],
  BIC_BXIIARMA=acuracia1[4,,5], 
  BIC_ETS=acuracia1[5,,5], 
  BIC_SARMA=acuracia1[6,,5], 
  BIC_XGBoost=c(rep(NA,33)),
  MAE_BXII=acuracia1[1,,2],
  MAE_Ray=acuracia1[2,,2],
  MAE_GAMMA=acuracia1[3,,2],
  MAE_BXIIARMA=acuracia1[4,,2],
  MAE_ETS=acuracia1[5,,2],
  MAE_SARMA=acuracia1[6,,2],
  MAE_XGBoost=acuracia1[7,,2],
  MAPE_BXII=acuracia1[1,,3],
  MAPE_Ray=acuracia1[2,,3],
  MAPE_GAMMA=acuracia1[3,,3],
  MAPE_BXIIARMA=acuracia1[4,,3],
  MAPE_ETS=acuracia1[5,,3],
  MAPE_SARMA=acuracia1[6,,3],
  MAPE_XGBoost=acuracia1[7,,3],
  RMSE_BXII=acuracia1[1,,1],
  RMSE_Ray=acuracia1[2,,1],
  RMSE_GAMMA=acuracia1[3,,1],
  RMSE_BXIIARMA=acuracia1[4,,1],
  RMSE_ETS=acuracia1[5,,1],
  RMSE_SARMA=acuracia1[6,,1],
  RMSE_XGBoost=acuracia1[7,,1],
  pv_LB_BXII=round(pv_LB_BXII,4),
  pv_LB_Ray=round(pv_LB_Ray,4),
  pv_LB_GAMMA=round(pv_LB_GAMMA,4),
  pv_LB_BXIIARMA=round(pv_LB_BXIIARMA,4),
  pv_LB_ETS=round(pv_LB_ETS,4),
  pv_LB_SARMA=round(pv_LB_SARMA,4),
  pv_LB_XGBoost=round(pv_LB_XGBoost,4)
  )
rownames(results)<-nome

setwd(file.path(path, "Apply BXII-GAS"))
saving<-paste0("app_hid_new","tau",tau,".RData")
save.image(saving)
end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)

resultados<-results
resultados<-round(resultados,4)
library(gt)

# Criar uma tabela estilizada
resultados[,] %>%
  gt()
