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

# Main project folder
if(Sys.info()["sysname"]=="Windows"){
  path <- "C:/Users/nando/OneDrive - BTH/Fernando/BXII-GAS-1" # Windowns
}else{
  path <- "~/Fernando/BXII-GAS" # Ubuntu
}

setwd(file.path(path, "BXII-GAS"))
source("GASBXIIFit.R")
source("predict_BXII_GAS.R")
setwd(file.path(path, "Ray-GAS/Ray-GAS-master/Ray-GAS"))
source("GASRayFit.R")
source("predict_Ray_GAS.R")
setwd(file.path(path, "gamma-GAS"))
source("GASGFit.R")
source("predict_G_GAS.R")
setwd(file.path(path, "BXII_ARMA"))
source("bxiiarmaCOV.fit.r")
source("bxiiarma.fit.r")
source("predict_BXII_ARMA.R")
setwd(file.path(path, "classic_models"))
source("predict_ARMA.R")
source("predict_ETS.R")
source("predict_SARMA.R")

############################
#### initial quantities ####
############################

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

datas <- seq(as.Date("2010-01-01"), as.Date("2024-08-31"), by = "day")
R<-dim(dados)[2]
h<-12
s<-12
out_forecast<-array(NA,c(h,6,R))
ac<-array(NA,c(6,R,3))
nome<-c()

start_time <- Sys.time()

for(i in 1:R){
  nome[i]<-substring(colnames(dados)[i],1,nchar(colnames(dados)[i]))
  y1<-na.omit(dados[,i])
  y1 <- zoo(y1, order.by = datas)
  y1 <- aggregate(y1, as.yearmon, sum)
  y1<-y1/100000
  y1[y1==0] <- mean(y1, na.rm = TRUE)
  y<-y1
  y_prev<-y[(length(y)[1]-h+1):length(y)[1]]
  print(i)

  prev_BXIIGAS<-predict_BXIIGAS(hid,y_prev,X=X,X_hat=X_hat)
  prev_RayGAS<-predict_RayGAS(hidRay,y_prev,X=XRay,X_hat=XRay_hat)
  prev_BXIIARMA<-try(predict_BXIIARMA(hidBXIIARMA,y_prev,X = as.matrix(XBXIIARMA),X_hat=as.matrix(XBXIIARMA_hat),serie=y),silent = TRUE)
  prev_ETS <- predict_ETS_1step(y_train = y, y_test  = y_prev, model   = "AAA")
  prev_ARMA <- predict_ARMA_1step( y_train = y, y_test  = y_prev, order   = c(matARMA[i,1], 0, matARMA[i,3]))
  prev_SARMA <- predict_SARMA_1step(y_train = y, y_test  = y_prev, order   = c(matSARMA[i,1], 0, matSARMA[i,2]), seasonal = c(matSARMA[i,3],0,matSARMA[i,4]), period  = s)

  X_XGBoost <- data.frame(
    # Lags
    lag1  = c(NA, y[-n]),
    # Moving averages, using only past observations
    MA2 = c(
      rep(NA, 2), vapply(3:n, function(t) mean(y[(t - 2):(t - 1)]), numeric(1))),
    # Annual seasonality
    Xsin = Xsin_hat,
    Xcos = Xcos_hat)

  X_XGBoost <- as.matrix( dados_XGBoost[, c("lag1", "MA2", "Xsin", "Xcos")])
  prev_XGBoost <- predict(final_XGBoost[[1]], X_XGBoost)

  out_forecast[,1,i] <- prev_BXIIGAS
  out_forecast[,2,i] <- prev_RayGAS
  if (inherits(prev_BXIIARMA, "try-error")) { out_forecast[,3,i] <- NA} else {out_forecast[,3,i] <- prev_BXIIARMA}
  out_forecast[,4,i] <- prev_ETS
  out_forecast[,5,i] <- prev_ARMA
  out_forecast[,6,i] <- prev_SARMA
  out_forecast[,7,i] <- prev_XGBoost

  # Accuracy
  assign(paste0("y_prev",i),y_prev)
  ac[1,i,]<- accuracy(out_forecast[,1,i], y_prev)[,c(2,3,5)]
  ac[2,i,]<- accuracy(out_forecast[,2,i], y_prev)[,c(2,3,5)]
  ac[3,i,]<-accuracy(out_forecast[,3,i], y_prev, na.rm = TRUE)[,c(2,3,5)]
  ac[4,i,]<- accuracy(out_forecast[,4,i], y_prev)[,c(2,3,5)]
  ac[5,i,]<- accuracy(out_forecast[,5,i], y_prev)[,c(2,3,5)]
  ac[6,i,]<- accuracy(out_forecast[,6,i], y_prev, na.rm = TRUE)[,c(2,3,5)]
  ac[7,i,]<- accuracy(out_forecast[,7,i], y_prev, na.rm = TRUE)[,c(2,3,5)]
  
  setwd("~/GitHub/BXII-GAS/Apply BXII-GAS")

  forecast <- paste0("Plots/forecast", i, ".pdf")
  if (sum(is.na(out_forecast[,1,i]))==0 &&
      sum(is.infinite(out_forecast[,1,i]))==0) {
    pdf(forecast,width = 8, height = 4)
    plot(y_prev,ylab="RF",ylim=c(min(y_prev,out_forecast[,1,i]),0.25+max(y_prev,out_forecast[,1,i])),type = "l")
    lines(ts(out_forecast[,1,i],frequency = s,start=c(2023,9)),col=2,lty=2,lwd=2)
    legend("topright",c("BXII-GAS"),col = c(2),lty= c(2),lwd = c(2), bty="n", cex = 1)
    dev.off()
  }
}

rankRMSE_out<-apply((apply(ac[,,1],2,rank)==1),1,sum)
rankMAE_out<-apply((apply(ac[,,2],2,rank)==1),1,sum)
rankMAPE_out<-apply((apply(ac[,,3],2,rank)==1),1,sum)

result_out<-rbind(rankRMSE_out, rankMAE_out, rankMAPE_out)
print(result_out)

results_out<-data.frame(
  MAE_BXII_out=ac[1,,2],
  MAE_Ray_out=ac[2,,2],
  MAE_BXIIARMA_out=ac[3,,2],
  MAE_ETS_out=ac[4,,2],
  MAE_ARMA_out=ac[5,,2],
  MAE_SARMA_out=ac[6,,2],
  RMSE_BXII_out=ac[1,,1],
  RMSE_Ray_out=ac[2,,1],
  RMSE_BXIIARMA_out=ac[3,,1],
  RMSE_ETS_out=ac[4,,1],
  RMSE_ARMA_out=ac[5,,1],
  RMSE_SARMA_out=ac[6,,1],
  MAPE_BXII_out=ac[1,,3],
  MAPE_Ray_out=ac[2,,3],
  MAPE_BXIIARMA_out=ac[3,,3],
  MAPE_ETS_out=ac[4,,3],
  MAPE_ARMA_out=ac[5,,3],
  MAPE_SARMA_out=ac[6,,3]
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

install.packages("openxlsx")
library(openxlsx)
write.xlsx(results, "my_file2.xlsx")
