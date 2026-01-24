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

#Fitted
dados<-dados[,-c(11,21,30:33)]
dados<-dados[,-c(7,9:11,14:15,21:22,24:25,27)]
datas <- seq(as.Date("2010-01-01"), as.Date("2024-08-31"), by = "day")
dados<-dados[,-c(3,5:9,13,15)]
data<-dados
data<-data[,-c(3,4,7)]

# Table of descriptive measures
cols <- 1:ncol(data)
res_list <- list()

for (i in cols) {
  
  dat <- data[, i] / 1000
  dat <- zoo(dat, order.by = datas)
  dat <- aggregate(dat, as.yearmon, mean)
  dat <- ts(dat, start = c(2010, 1), frequency = 12)
  
  a <- c(
    summary(dat, na.rm = TRUE),
    sd(dat, na.rm = TRUE),
    (sd(dat, na.rm = TRUE) / summary(dat, na.rm = TRUE)[4]) * 100
  )
  
  a <- round(a[c(1, 3, 4, 6, 7, 8)], 4)
  
  res_list[[colnames(data)[i]]] <- a
}

res_mat <- do.call(rbind, res_list)
colnames(res_mat) <- c("Min", "1st Qu.", "Median", "Mean", "SD", "CV (%)")
xtable(res_mat, digits = 4)
