# Reference: Ray-GAS
# Created by Miguel Peña-Ramirez (miguepe@gmail.com), March/2023
rm(list = ls())
source("Ray-Gas.R")
source("Ray-Gas_fit.R")
source("simu_rarma.R")
source("rarma.r")
library(forecast)
############################
#### initial quantities ####
############################
start_time <- Sys.time()
vn<-c(100, 200, 500) 
R<-10
bug<-0
true<-"rarma"
#true<-"Ray"
scenarios<-c("Forest","Lake")
scenario=scenarios[2]
if(scenario==scenarios[1]){
  scen=1
  w=-1.24; A=.76; B=.39 #scenario Forest 
  theta=c(w,A,B)
}else{
  if(scenario==scenarios[2]){
    scen=2
    w=-1.73; A=.64; B=.24 #scenario Lake 
    theta=c(w,A,B)
  }
}
estimrarma<-estim<-array(NA,c(R,length(theta),length(vn)))
AIC<-BIC<-RMSER<-MAPE<-array(NA,c(R,2,length(vn)))
mean_MLEs<-RMSE<-RB<-matrix(NA,length(vn),length(theta))
contAIC<-contBIC<-contRMSER<-contMAPE<-matrix(NA,length(vn),2)
##########################
#### simulation study ####
##########################
for(j in 1:length(vn)){
  n<-vn[j]
  print(paste0("Sample size: ",n))
  t=1:n
  t_hat=n+1
  i<-1
  X=NA
  X_hat=NA
  set.seed(10)
  while(i<=R){
    if(true=="Ray"){
      y <- simu.GASRay(n,A=A,B=B,w=w)
    }else{
      y <- simu.Rarma(n,phi=A,theta=B,alpha=w)
    }
    result <- try(GASRAY.fit(y,1,1,h=25),silent = T)
    resultrarma <- try(rarma(y,1,1,h=25,diag=0,resid=1),silent = T)
    if(class(result)=="try-error" ||  result$conv != 0 ||
       class(resultrarma)=="try-error" ||  resultrarma$conv != 0){
      bug<-bug+1
    }else{
      modelray <- cbind(round(result$coeff,4),round(result$stderror,4),
                        round(result$zstat,4),round(result$pvalues,4))
      colnames(modelray)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
      estim[i,,j]<-modelray[,1]
      modelrarma <- cbind(round(resultrarma$coeff,4),round(resultrarma$stderror,4),
                          round(resultrarma$zstat,4),round(resultrarma$pvalues,4))
      colnames(modelrarma)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
      estimrarma[i,,j]<-modelrarma[,1]
      AIC[i,,j]<-c(result$aic,resultrarma$aic)
      BIC[i,,j]<-c(result$bic,resultrarma$bic)
      RMSER[i,,j]<-c(accuracy(result$fitted, y)[2],
                     accuracy(resultrarma$fitted, y)[2])
      MAPE[i,,j]<-c(accuracy(result$fitted, y)[5],
                    accuracy(resultrarma$fitted, y)[5])
      i<-i+1
    }
  }
  mean_MLEs[j,] <- apply(estim[,,j], 2, mean)
  RB[j,] <- (mean_MLEs[j,]-theta)/theta*100 
  RMSE[j,] <- (apply(estim[,,j],2,var)+(mean_MLEs[j,]-theta)^2)
  contAIC[j,]<-apply(apply(AIC[,,j],1,rank)==1,1,sum)
  contBIC[j,]<-apply(apply(BIC[,,j],1,rank)==1,1,sum)
  contRMSER[j,]<-apply(apply(RMSER[,,j],1,rank)==1,1,sum)
  contMAPE[j,]<-apply(apply(MAPE[,,j],1,rank)==1,1,sum)
}
dimnames(contAIC) <- list(vn,c("Ray-GAS", "RARMA"))
dimnames(contBIC) <- list(vn,c("Ray-GAS", "RARMA"))
resultados<-rbind(contAIC*100/R,contBIC*100/R)
print(resultados)
if(!file.exists("simulation")){
  dir.create("simulation") 
}
saving<-paste0("simulation/simuTable2_",true,"_scen",scen,"_vn",min(vn),"-",max(vn),"_R",R,".RData")
save.image(saving)
end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)