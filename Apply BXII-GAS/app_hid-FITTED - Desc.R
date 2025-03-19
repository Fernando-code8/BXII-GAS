# Reference: BXII-GAS
# Created by Fernando José Monteiro de Araújo (nandomonteiro418@gmail.com), September/2024

library(tidyverse)
library(dplyr)
library(zoo)
library(forecast)
library(psych)
library(readr)
library(xtable)

############################
#### initial quantities ####
############################
w1<-10 # width for plots 
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

#Fitted
dados<-dados[,-c(11,21,30:33)]
dados<-dados[,-c(7,9:11,14:15,21:22,24:25,27)]
datas <- seq(as.Date("2010-01-01"), as.Date("2024-08-31"), by = "day")
dados<-dados[,-c(3,5:9,13,15)]
data<-dados
data<-data[,-c(3,4,7)]

# Table of descriptive measures
i=4
dat <-data[,i]
dat <- zoo(dat, order.by = datas)
dat <- aggregate(dat, as.yearmon, mean)
dat<-ts(dat,start = c(2010,1),frequency = 12)

a<-c(summary(dat,na.rm=T),
     var(dat,na.rm=T))
round(a,4)
a<-as.matrix(t(a))
xtable(a,digits = 4)

min1<-a[1]
ids<-which(dat==min1)

Sul<-rowMeans(data,na.rm=T)
Sul.ts<-ts(Sul,start = c(2010,1), frequency = 12)
y <- zoo(Sul.ts, order.by = datas)
y <- aggregate(y, as.yearmon, mean)
y<-y#/1000
# ---------------------------------
# seasonality
setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/Apply BXII-GAS")
months<-paste0("Plots/monthsS",".pdf")
pdf(months,width = w1, height = h11)
par(mfrow=c(1,1))
monthplot(ts(y,start = c(2010,1), frequency = 12),main="",ylab = "RF")
dev.off()
 
?monthplot()

