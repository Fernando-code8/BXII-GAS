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
setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/BXII-GAS")
source("GASBXIIFit.R")
setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/Ray-GAS/Ray-GAS-master/Ray-GAS")
source("GASRayFit.R")
setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/gamma-GAS")
source("GASGFit.R")

############################
#### initial quantities ####
############################

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
R<-dim(dados)[2]
h<-12
s<-12
out_forecast<-array(NA,c(h,3,R))
ac<-array(NA,c(3,R,3))
nome<-c()
model<-c(  "BXII with cos covariate",        
           "BXII with cos covariate" ,       
           "BXII with sin covariate" ,       
           "BXII with cos covariate" ,       
           "BXII with sin and cos covariate",
           "BXII with sin and cos covariate",
           "BXII with sin covariate"    ,    
           "BXII with sin covariate"  ,      
           "BXII with sin covariate"    ,    
           "BXII with sin and cos covariate",
           "BXII with cos covariate"        ,
           "BXII with sin and cos covariate",
           "BXII with sin covariate"   ,     
           "BXII with sin covariate"   ,     
           "BXII with cos covariate"   ,     
           "BXII with sin and cos covariate")

matBXII<-matrix(c(1,3,1,2,1,2,1,1,1,2,3,3,1,3,1,2,1,1,1,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,1,1,1,2,1,1,1,3,1,2,1,1,1,3,1,3,1,1,1,2,3,3,2,2,1,2,1,1,1,1),ncol=4,byrow=T)


model_Ray<-c(    "Ray with cos covariate",        
                 "Ray with sin and cos covariate",
                 "Ray without covariate"     ,    
                 "Ray with sin and cos covariate",
                 "Ray with sin covariate"    ,    
                 "Ray with sin and cos covariate",
                 "Ray with sin covariate"   ,     
                 "Ray with sin covariate"  ,      
                 "Ray with sin covariate"    ,    
                 "Ray with sin and cos covariate",
                 "Ray with sin and cos covariate",
                 "Ray with sin and cos covariate",
                 "Ray with sin covariate"  ,      
                 "Ray with cos covariate"  ,      
                 "Ray with sin covariate"   ,     
                 "Ray with sin and cos covariate")

matRay<-matrix(c(1,3,1,2,1,2,1,3,1,3,1,1,1,2,2,2,1,1,1,3,1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,2,3,3,1,2,1,3,1,1,1,3,1,2,1,3,3,3,1,3,1,2,1,3,1,1,1,1),ncol=4,byrow=T)

model_GAMMA<-c(   "gamma without covariate" ,
                  "gamma with cos covariate",
                  "gamma with cos covariate",
                  "gamma without covariate" ,
                  "gamma without covariate" ,
                  "gamma without covariate" ,
                  "gamma without covariate" ,
                  "gamma without covariate" ,
                  "gamma with cos covariate",
                  "gamma with cos covariate",
                  "gamma without covariate" ,
                  "gamma without covariate" ,
                  "gamma without covariate", 
                  "gamma without covariate" ,
                  "gamma without covariate" ,
                  "gamma with sin covariate")

matgamma<-matrix(c(2,2,1,2,2,2,2,2,2,2,3,3,1,1,1,2,2,2,3,3,1,1,1,1,1,1,1,3,2,2,3,3,2,2,3,3,1,2,1,3,1,1,2,2,1,1,1,2,1,1,1,2,1,1,1,3,1,1,1,2,1,2,1,3),ncol=4,byrow=T)
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
  for(j in 1:h){
    print(c(i,j))
    y<-ts(y1[1:(length(y1)[1]-h+j-1)],frequency = s,start = c(2010,1))
    t<-1:length(y)
    t_hat<-length(y)+1
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
    
    hid<-BXIIGAS.fit(y,matBXII[i,1]:matBXII[i,2],matBXII[i,3]:matBXII[i,4],X=X,X_hat=X_hat,h1=1)
    hidRay<-RayGAS.fit(y,matRay[i,1]:matRay[i,2],matRay[i,3]:matRay[i,4],X=XRay,X_hat=XRay_hat,h1=1)
    hidGAMMA<-GGAS.fit(y,matgamma[i,1]:matgamma[i,2],matgamma[i,3]:matgamma[i,4],X=XGAMMA,X_hat=XGAMMA_hat,h1=1)
    
    out_forecast[j,1,i]<-hid$forecast
    out_forecast[j,2,i]<-hidRay$forecast
    out_forecast[j,3,i]<-hidGAMMA$forecast
  }
  assign(paste0("y_prev",i),y_prev)
  ac[1,i,]<- accuracy(out_forecast[,1,i], y_prev)[,c(2,3,5)]
  ac[2,i,]<- accuracy(out_forecast[,2,i], y_prev)[,c(2,3,5)]
  ac[3,i,]<- accuracy(out_forecast[,3,i], y_prev)[,c(2,3,5)]
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
                       round(prop.table(result_out,1),4)*100,")"),3,3)))

results_out<-data.frame(
  MAE_BXII_out=ac[1,,2],
  MAE_Ray_out=ac[2,,2],
  MAE_GAMMA_out=ac[3,,2],
  RMSE_BXII_out=ac[1,,1],
  RMSE_Ray_out=ac[2,,1],
  RMSE_GAMMA_out=ac[3,,1],
  MAPE_BXII_out=ac[1,,3],
  MAPE_Ray_out=ac[2,,3],
  MAPE_GAMMA_out=ac[3,,3]
)

rownames(results_out)<-nome
print(results_out)

setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/Apply BXII-GAS")
saving<-paste0("app_hid-FORECAST1stepnew.RData")
save.image(saving)

end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)

results_out<-res

# Criar uma tabela estilizada
results_out %>%
  gt()
