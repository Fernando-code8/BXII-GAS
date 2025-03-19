# Reference: UG-GAS
# Created by Fernando J M de Araújo (nandomonteiro418@gmail.com), April/2024
rm(list=ls())

setwd("C:/Users/nando/OneDrive/Desktop/mestrado UFRGS/Dissertação/Scripts Atual BXII-GAS/BXII-GAS")
source("simu_BXII-GAS.R")
source("GASBXIIFit.R")
library(beepr)

############################
#### initial quantities ####
############################
start_time <- Sys.time()
vn<-c(
  100,200,300,500,1000
  )
R<-1000
bug<-0

# w<-0.1;A=.2;B=.1;c=0.9;beta=c(0.1,0.5) 

w<-0.1;A=0.2;B=0.4;c=5;beta=c(-0.1#,.5
                              ) 

# X<-sin(2*pi*t/12) # in sample
# X_hat<-sin(2*pi*t_hat/12) # out of sample
theta=c(w,A,B,beta,
        c)
# cov<-"sin&cos"
cov<-"sin"
# cov<-"sin&cos"
# cov<-NA
# X=NA
# X_hat=NA
# theta=c(w,A,B,psi)
estim<-err<-ICi<-ICs<-array(NA,c(R,length(theta),length(vn)))
comega<-calpha1<-cbeta1<-cc<-cdelta1<-cdelta2<-0
RMSE<-mean_MLEs<-RB<-TC<-matrix(NA,length(vn),length(theta))
# For the confidence intervals
alpha_erro <-0.05
quantil <-1-alpha_erro/2
z<-qnorm(quantil)
##########################
#### simulation study ####
##########################

set.seed(2024)
for(j in 1:length(vn)){
  n<-vn[j]
  t=1:n
  t_hat=n+1
  i<-1
  # X<-cbind(sin(2*pi*t/12),cos(2*pi*t/12))
  # X_hat<-cbind(sin(2*pi*t_hat/12),cos(2*pi*t_hat/12))
  X<-sin(2*pi*t/12) # in sample
  X_hat<-sin(2*pi*t_hat/12) # out of sample
  # X=NA
  # X_hat=NA
  set.seed(10)
  comega<-calpha1<-cbeta1<-cc<-cdelta1<-cdelta2<-0
  bug<-0
  while(i<=R){
    print(c(i,j))
    y <- simu.BXIIGAS(n,A=A,B=B,w=w,c=c,beta=beta,tau=0.5,X=cov)
    # y <- simu.BXIIGAS(n,A=A,B=B,w=w,c=c,tau=0.5)
    result <- try(BXIIGAS.fit(y,ar=1,ma=1,X=X,X_hat=X_hat,
                              h=length(t_hat)),silent = T)
    if(class(result)=="try-error" || 
       result$conv != 0 ||
       sum(is.na(result$stderror))!=0){
      bug<-bug+1
      print("Erro de Estimativa")
    }else{
      estim[i,,j]<-result$model[,1]
      err[i,,j] <- result$model[,2]
      ICi[i,,j]<-estim[i,,j]-(z*err[i,,j])
      ICs[i,,j]<-estim[i,,j]+(z*err[i,,j])
      
      if (ICi[i,1,j]<=theta[1] && ICs[i,1,j]>=theta[1])
      {
        comega<-comega+1
      }
      if (ICi[i,2,j]<= theta[2] && ICs[i,2,j]>=theta[2])
      {
        calpha1<-calpha1+1
      }
      if (ICi[i,3,j]<= theta[3] && ICs[i,3,j]>=theta[3])
      {
        cbeta1<-cbeta1+1
      }
      if (ICi[i,4,j]<= theta[4] && ICs[i,4,j]>=theta[4])
      {
         cdelta1<-cdelta1+1
      }
      # if (ICi[i,5,j]<= theta[5] && ICs[i,5,j]>=theta[5])
      # {
      #   cdelta2<-cdelta2+1
      # }
      if (ICi[i,5,j]<= theta[5] && ICs[i,5,j]>=theta[5])
      {
        cc<-cc+1
      }
      i<-i+1
    }
  }
  mean_MLEs[j,] <- apply(estim[,,j], 2, mean)
  RB[j,] <- (mean_MLEs[j,]-theta)/theta*100 
  RMSE[j,] <- sqrt(apply(estim[,,j],2,var)+(mean_MLEs[j,]-theta)^2)
  TC[j,]<- c(comega,calpha1,cbeta1,cdelta1,#cdelta2,
             cc
             )/R

  }
scen=2

setwd("C:/Users/nando/OneDrive/Desktop/mestrado UFRGS/Dissertação/Scripts Atual BXII-GAS/Simulation")
saving<-paste0("results_simu/simu_scensin",scen,"R",R,".RData")
save.image(saving)
end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)
beep(8)
