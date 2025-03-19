# Reference: Ray-GAS
# Created by Miguel Pena-Ramirez (miguepe@gmail.com), March/2023 
rm(list = ls())
source("Ray-Gas.R")
source("Ray-Gas_fit.R")

############################
#### initial quantities ####
############################
time <- Sys.time()
set.seed(6)
bug<-0
scenarios<-c("Forest","Lake")
R <- 10000
n_samp=500 #100,200,500,1000
scenario=scenarios[1]
ns=0.05
COEF<-matrix(NA,nrow=R,ncol = 3) 
RESQ<-matrix(NA,nrow=R,ncol = 3)
BIAS<-matrix(NA,nrow=R,ncol = 3)
LI  <-matrix(NA,nrow=R,ncol = 3)
LS  <-matrix(NA,nrow=R,ncol = 3)

if (scenario==scenarios[1]){
  #Scenario Forest	
  w=-1.24; A=.76; B=.39   
}else{
  #Scenario Lake	
  w=-1.73; A=.64; B=.24  
}
##########################
#### simulation study ####
##########################
i<-1
while(i<=R){
  y<-simu.GAS(n_samp,A=A,B=B,w=w)  
  res<-try(GASRAY.fit(y,1,1,h1=1),silent=T)
  if(class(res)=="try-error" ||  res$conv != 0 || anyNA(res$stderror)){
    bug<-bug+1
  }else{
   COEF[i,]<-res$coef
   BIAS[i,]<-res$coef-c(w,A,B) 
   RESQ[i,]<-(res$coef-c(w,A,B))^2 
   #IC
   LI[i,]<-res$coeff+qnorm(ns/2)*sqrt(res$stderror)
   LS[i,]<-res$coeff+qnorm(1-ns/2)*sqrt(res$stderror)
   print(i)
   i<-i+1
  }
}
param <- c(w,A,B)
m_MLEs <- apply(COEF, 2, mean) #mean_MLEs
RB <- (apply(COEF, 2, mean)-param)/param*100 #RB
MSE <- (apply(COEF,2,var)+(m_MLEs-param)^2) #MSE
round(RB,3)
round(MSE,3)

######## Coverage rates#############
CR<-c()
CR[1]<-sum(w>LI[,1] & w<LS[,1])/R*100
CR[2]<-sum(A>LI[,2] & A<LS[,2])/R*100
CR[3]<-sum(B>LI[,3] & B<LS[,3])/R*100
####################################
#anyNA(res$stderror)
#apply(is.na(LS),2,sum)
#which(is.na(LS[,3]))
#apply(is.na(COEF),2,sum)
#which(is.na(LS[,3]))
#LI[21637,]
resultados<-rbind(c("w","A","B"),round(RB,3),round(MSE,3),round(CR,3))
print(resultados)
timef <- Sys.time() - time
if(!file.exists("simulation")){
  dir.create("simulation") 
}
saving<-paste0("simulation/simuTable1_",scenario,"_R",R,"_",n_samp,".RData")
print(timef)
save.image(saving)
