# Reference: BXII-GAS
# Created by Fernando José Monteiro de Araújo (nandomonteiro418@gmail.com), november/2024

####################
#### R packages ####
####################
library(tidyverse)
library(dplyr)
library(zoo)
library(forecast)
library(psych)
library(readr)
library(PTSR)
library(xtable)
setwd("~/GitHub/BXII-GAS/BXII-GAS")
source("GASBXIIFit.R")
source("predict_BXII_GAS.R")
setwd("~/GitHub/BXII-GAS/Ray-GAS/Ray-GAS-master/Ray-GAS")
source("GASRayFit.R")
source("predict_Ray_GAS.R")
setwd("~/GitHub/BXII-GAS/gamma-GAS")
source("GASGFit.R")
source("predict_G_GAS.R")
setwd("~/GitHub/BXII-GAS/BXII_ARMA")
source("bxiiarmaCOV.fit.r")
source("bxiiarma.fit.r")
source("predict_BXII_ARMA.R")
setwd("~/GitHub/BXII-GAS/classic_models")
source("predict_ARMA.R")
source("predict_ETS.R")
source("predict_SARMA.R")

############################
#### initial quantities ####
############################

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
R<-dim(dados)[2]
h<-12
s<-12
out_forecast<-array(NA,c(h,7,R))
ac<-array(NA,c(7,R,3))
nome<-c()
model<-c("BXII with sin and cos covariate",        
         "BXII with sin and cos covariate",
         "BXII with sin and cos covariate",        
         "BXII with sin and cos covariate",
         "BXII with sin and cos covariate")

matBXII<-matrix(c(rep(1,1*20)),ncol=4,byrow=T)


model_Ray<-c("Ray with sin and cos covariate",        
             "Ray with sin and cos covariate",
             "Ray with sin and cos covariate",        
             "Ray with sin and cos covariate",        
             "Ray with sin and cos covariate")

matRay<-matrix(c(rep(1,1*20)),ncol=4,byrow=T)

model_GAMMA<-c("gamma with sin covariate", 
               "gamma with sin covariate",
               "gamma with sin covariate", 
               "gamma with sin covariate", 
               "gamma with cos covariate")

matgamma<-matrix(c(rep(1,1*20)),ncol=4,byrow=T)

modeloBXIIARMA<-c("BXIIARMA with sin covariate",
                  "BXIIARMA with sin covariate",
                  "BXIIARMA with sin covariate",
                  "BXIIARMA with sin covariate",
                  "BXIIARMA with sin covariate")

matBXIIARMA<-matrix(c(2,2,1,2,
                      1,1,1,1,
                      2,2,1,2,
                      3,3,1,1,
                      2,2,1,1),ncol=4,byrow=T)

matARMA<-matrix(c(rep(1,1*20)),ncol=4,byrow=T)

matSARMA<-matrix(c(1,1,1,1,
                   1,1,0,1,
                   1,1,1,1,
                   1,1,1,1,
                   1,1,1,0),ncol=4,byrow=T)


start_time <- Sys.time()


for(i in 1:R){
  nome[i]<-substring(colnames(dados)[i],1,nchar(colnames(dados)[i])#-10
                     )
  y1<-na.omit(dados[,i])
  y1 <- zoo(y1, order.by = datas)
  y1 <- aggregate(y1, as.yearmon, mean)
  y1<-y1/1000
  y1[y1==0] <- mean(y1, na.rm = TRUE)
  
  y_prev<-y1[(length(y1)[1]-h+1):length(y1)[1]]
  print(i)
  #for(j in 1:h){
  #  print(c(i,j))
    y<-ts(y1[1:(length(y1)[1]-h)],frequency = s,start = c(2010,1))
    t<-1:length(y)
    t_hat<-c((length(y)+1):(length(y)+h))
    if(model[i]=="BXII with sin covariate"){
      X<-sin(2*pi*t/s)
      X_hat<-sin(2*pi*t_hat/s)
    }else{
      if(model[i]=="BXII with sin and cos covariate"){
        X<-cbind(sin(2*pi*t/s),cos(2*pi*t/s))
        X_hat<-cbind(sin(2*pi*t_hat/s),cos(2*pi*t_hat/s))
      }else{
        if(model[i]=="BXII without covariate"){
          X<-X_hat<-NA
        }else{
          if(model[i]=="BXII with cos covariate")
            X<-cos(2*pi*t/s)
            X_hat<-cos(2*pi*t_hat/s)
        }
      }
    }
    
    if(model_Ray[i]=="Ray with sin covariate"){
      XRay<-sin(2*pi*t/s)
      XRay_hat<-sin(2*pi*t_hat/s)
      r=1
    }else{
      if(model_Ray[i]=="Ray with sin and cos covariate"){
        XRay<-cbind(sin(2*pi*t/s),cos(2*pi*t/s))
        XRay_hat<-cbind(sin(2*pi*t_hat/s),cos(2*pi*t_hat/s))
        r=2
      }else{
        if(model_Ray[i]=="Ray without covariate"){
          XRay<-XRay_hat<-NA
          r=0
        }else{
          if(model_Ray[i]=="Ray with cos covariate")
            XRay<-cos(2*pi*t/s)
            XRay_hat<-cos(2*pi*t_hat/s)
            r=1
        }
      }
    }
    
    if(model_GAMMA[i]=="gamma with sin covariate"){
      XGAMMA<-sin(2*pi*t/s)
      XGAMMA_hat<-sin(2*pi*t_hat/s)
    }else{
      if(model_GAMMA[i]=="gamma with sin and cos covariate"){
        XGAMMA<-cbind(sin(2*pi*t/s),cos(2*pi*t/s))
        XGAMMA_hat<-cbind(sin(2*pi*t_hat/s),cos(2*pi*t_hat/s))
      }else{
        if(model_GAMMA[i]=="gamma without covariate"){
          XGAMMA<-XGAMMA_hat<-NA
        }else{
          if(model_GAMMA[i]=="gamma with cos covariate")
            XGAMMA<-cos(2*pi*t/s)
          XGAMMA_hat<-cos(2*pi*t_hat/s)
        }
      }
    }
    
    if(modeloBXIIARMA[i]=="BXIIARMA with sin covariate"){
      XBXIIARMA<-sin(2*pi*t/s)
      XBXIIARMA_hat<-sin(2*pi*t_hat/s)
    }else{
      if(modeloBXIIARMA[i]=="BXIIARMA with sin and cos covariate"){
        XBXIIARMA<-cbind(sin(2*pi*t/s),cos(2*pi*t/s))
        XBXIIARMA_hat<-cbind(sin(2*pi*t_hat/s),cos(2*pi*t_hat/s))
      }else{
        if(modeloBXIIARMA[i]=="BXIIARMA without covariate"){
          XBXIIARMA<-XBXIIARMA_hat<-NA
        }else{
          if(modeloBXIIARMA[i]=="BXIIARMA with cos covariate")
            XBXIIARMA<-cos(2*pi*t/s)
          XBXIIARMA_hat<-cos(2*pi*t_hat/s)
        }
      }
    }
    
    # train
    hid<-BXIIGAS.fit(y,matBXII[i,1]:matBXII[i,2],matBXII[i,3]:matBXII[i,4],X=X,X_hat=X_hat,h1=1)
    hidRay<-RayGAS.fit(y,matRay[i,1]:matRay[i,2],matRay[i,3]:matRay[i,4],X=XRay,X_hat=XRay_hat,h1=1)
    hidGAMMA<-GGAS.fit(y,matgamma[i,1]:matgamma[i,2],matgamma[i,3]:matgamma[i,4],X=XGAMMA,X_hat=XGAMMA_hat,h1=1)
    if(modeloBXIIARMA[i]=="BXIIARMA with sin covariate" || modeloBXIIARMA[i]=="BXIIARMA with sin and cos covariate" || modeloBXIIARMA[i]=="BXIIARMA with cos covariate"){
      hidBXIIARMA<-bxiiarmaCOV.fit(y,ar=matBXIIARMA[i,1]:matBXIIARMA[i,2],ma=matBXIIARMA[i,3]:matBXIIARMA[i,4],X = as.matrix(XBXIIARMA),X_hat=as.matrix(XBXIIARMA),h=1,diag=0,tau=0.5)
    } else {
      hidBXIIARMA<-bxiiarma.fit(y,ar=matBXIIARMA[i,2],ma=matBXIIARMA[i,4],h1=1,resid = 3,diag1=1)
    }

    # test
    prev_BXIIGAS<-predict_BXIIGAS(hid,y_prev,X=X,X_hat=X_hat)
    prev_RayGAS<-predict_RayGAS(hidRay,y_prev,X=XRay,X_hat=XRay_hat)
    prev_GGAS<-predict_GGAS(hidGAMMA,y_prev,X=XGAMMA,X_hat=XGAMMA_hat)
    prev_BXIIARMA<-predict_BXIIARMA(hidBXIIARMA,y_prev,X = as.matrix(XBXIIARMA),X_hat=as.matrix(XBXIIARMA_hat),serie=y)
    prev_ETS <- predict_ETS_1step(
      y_train = y,
      y_test  = y_prev,
      model   = "AAA"
    )
    
    prev_ARMA <- predict_ARMA_1step(
      y_train = y,
      y_test  = y_prev,
      order   = c(matARMA[i,1], 0, matARMA[i,3])
    )
    
    prev_SARMA <- predict_SARMA_1step(
      y_train = y,
      y_test  = y_prev,
      order   = c(matSARMA[i,1], 0, matSARMA[i,2]),
      seasonal = c(matSARMA[i,3],0,matSARMA[i,4]),
      period  = s
    )
    
    out_forecast[,1,i] <- prev_BXIIGAS
    out_forecast[,2,i] <- prev_RayGAS
    out_forecast[,3,i] <- prev_GGAS
    out_forecast[,4,i] <- prev_BXIIARMA
    out_forecast[,5,i] <- prev_ETS
    out_forecast[,6,i] <- prev_ARMA
    out_forecast[,7,i] <- prev_SARMA
    #}
  # Accuracy
  assign(paste0("y_prev",i),y_prev)
  ac[1,i,]<- accuracy(out_forecast[,1,i], y_prev)[,c(2,3,5)]
  ac[2,i,]<- accuracy(out_forecast[,2,i], y_prev)[,c(2,3,5)]
  ac[3,i,]<- accuracy(out_forecast[,3,i], y_prev)[,c(2,3,5)]
  ac[4,i,]<- accuracy(out_forecast[,4,i], y_prev)[,c(2,3,5)]
  ac[5,i,]<- accuracy(out_forecast[,5,i], y_prev)[,c(2,3,5)]
  ac[6,i,]<- accuracy(out_forecast[,6,i], y_prev)[,c(2,3,5)]
  ac[7,i,]<- accuracy(out_forecast[,7,i], y_prev, na.rm = TRUE)[,c(2,3,5)]
  
  setwd("~/GitHub/BXII-GAS/Apply BXII-GAS")
  forecast<-paste0("Plots/forecast",i,".pdf")
  pdf(forecast,width = 8, height = 4)
  plot(y_prev,ylab="RF",ylim=c(min(y_prev,out_forecast[,1,i]),0.25+max(y_prev,out_forecast[,1,i])),type = "l")
  lines(ts(out_forecast[,1,i],frequency = s,start=c(2023,9)),col=2,lty=2,lwd=2)
  legend("topright", 
         c("BXII-GAS"),
         col = c(2),
         lty= c(2),
         lwd = c(2), bty="n", cex = 1)
  dev.off()
}

rankRMSE_out<-apply((apply(ac[,,1],2,rank)==1),1,sum)
rankMAE_out<-apply((apply(ac[,,2],2,rank)==1),1,sum)
rankMAPE_out<-apply((apply(ac[,,3],2,rank)==1),1,sum)

result_out<-rbind(rankRMSE_out,
                  rankMAE_out,
                  rankMAPE_out)
print(result_out)
print(prop.table(result_out,1))
xtable(t(matrix(paste0(round(result_out,3)," (",
                       round(prop.table(result_out,1),4)*100,")"),3,7)))

results_out<-data.frame(
  MAE_BXII_out=ac[1,,2],
  MAE_Ray_out=ac[2,,2],
  MAE_GAMMA_out=ac[3,,2],
  MAE_BXIIARMA_out=ac[4,,2],
  MAE_ETS_out=ac[5,,2],
  MAE_ARMA_out=ac[6,,2],
  MAE_SARMA_out=ac[7,,2],
  RMSE_BXII_out=ac[1,,1],
  RMSE_Ray_out=ac[2,,1],
  RMSE_GAMMA_out=ac[3,,1],
  RMSE_BXIIARMA_out=ac[4,,1],
  RMSE_ETS_out=ac[5,,1],
  RMSE_ARMA_out=ac[6,,1],
  RMSE_SARMA_out=ac[7,,1],
  MAPE_BXII_out=ac[1,,3],
  MAPE_Ray_out=ac[2,,3],
  MAPE_GAMMA_out=ac[3,,3],
  MAPE_BXIIARMA_out=ac[4,,3],
  MAPE_ETS_out=ac[5,,3],
  MAPE_ARMA_out=ac[6,,3],
  MAPE_SARMA_out=ac[7,,3]
)

rownames(results_out)<-nome
print(results_out)

setwd("~/GitHub/BXII-GAS/Apply BXII-GAS")
saving<-paste0("app_hid-FORECAST1stepnew.RData")
save.image(saving)

end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)

# Criar uma tabela estilizada
results_out %>%
  gt()

results_out[,c(4:7,11:14,18:21)]
