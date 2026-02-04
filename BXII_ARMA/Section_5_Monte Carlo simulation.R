################################################################################
# PAPER: Quantile-based dynamic modeling of asymmetric data: a novel Burr XII
#        approach for positive continuous random variables
# SUBSECTION: 5. Monte Carlo simulation 
# GOAL: Conducting a simulation study for the BXII-ARMA model
# AUTHORS: Fernando Jose Monteiro de Araujo, Renata Rojas Guerra and 
#          Fernando Arturo Pena-Ramirez
# LAST UPDATE: May 25, 2024
################################################################################

rm(list = ls())
source("bxiiarma.fit.r")
source("simu_bxiiarma.R")
library(beepr)

# Regression Simulation (Table 1)

R <- 10000
vn <- c(70,150,300,500,1000)

alpha <- 1
phi<-c(.5)
theta <- c(0.2)
c<-.5
model <- c(alpha,
             phi,
           theta,
           c)

p1<- 1:length(phi)
q1<- 1:length(theta)

set.seed(2024)
tempo.inicio = Sys.time()
for(n in vn){

# For the confidence intervals
alpha_erro <-0.05
quantil <-1-alpha_erro/2
z<-qnorm(quantil)
t<-qt(quantil, n-1)
  
# To save the results
estim <-ICi<-ICs<- err <- matrix(NA, nrow = R, ncol = length(model))
calpha<-cphi1<-cphi2<-ctheta1<-ctheta2<-cc<-0
i<-0


### simulation
while(i<R) {
  y <- simu.bxiiarma(n,phi=phi,theta=theta,alpha=alpha,c=c,tau = 0.5)
  result <- try(bxiiarma.fit(y,ar=p1,ma=q1,mod=3,tau = 0.5,diag1=1), silent = T)
  if(class(result) == "try-error" )
  {
    result$conv<-1
    convergencia<-1
  }
  
  if(result$conv == 0)
  {
    convergencia<-  result$conv
  }else{
    
    warning("FUNCTION DID NOT CONVERGE!")
  }
  
  if(convergencia==0 && sum(is.na(result$stderror))==0)
  {
    i<-i+1
    print(c("i=",i))
    estim[i,] <- result$model[,1]
    err[i,] <- result$model[,2]
    ICi[i,]<-estim[i,]-(z*err[i,])
    ICs[i,]<-estim[i,]+ (z*err[i,])
    
    if (ICi[i,1]<=alpha && ICs[i,1]>=alpha)
    {
      calpha<-calpha+1
    }
    
    if (ICi[i,2]<= phi[1] && ICs[i,2]>=phi[1])
    {
      cphi1<-cphi1+1
    }
    # if (ICi[i,3]<= phi[2] && ICs[i,3]>=phi[2])
    # {
    #   cphi2<-cphi2+1
    # }
    if (ICi[i,3]<= theta[1] && ICs[i,3]>=theta[1])
    {
      ctheta1<-ctheta1+1
    }
    # if (ICi[i,5]<= theta[2] && ICs[i,5]>=theta[2])
    # {
    #   ctheta2<-ctheta2+1
    # }
    if (ICi[i,4]<= c && ICs[i,4]>=c)
    {
      cc<-cc+1
    }
  }
} 

### mean
m <- apply(estim, 2, mean)
### bias
bias <- (model-m)
### relative percentage bias
biasP <- bias/model *100
### SD
erro <- apply(estim, 2, sd)
### MSE
MSE <- apply(estim, 2, var)+bias^2
### 
TC<-c(calpha,
      cphi1,#cphi2,
      ctheta1,
      # ctheta2,
      cc)/R
## final results
results <- rbind(m, bias, biasP, erro, MSE,TC)
rownames(results) <- c("Mean", "Bias","RB%", "SE", "MSE","TC")
print(c("Tamanho da Amostra:",n))
print(round(results,4))
}
mnjtempo.fim = Sys.time()
tempo.exec = mnjtempo.fim-tempo.inicio
print(tempo.exec)
beep(8)



# (Figure 2) Average coverage rates of confidence intervals for the several 

postscript(file = "ACR",width = 6, height = 4,family = "Times")
par(mfrow=c(1,1))
par(mar=c(2.8, 2.5, 0.5, 0.5))
par(mgp=c(1.3, 0.45, 0))

ARMA11<-c(0.9356,0.9417,0.9420,0.9430,0.9420)
ARMA10<-c(0.9478,0.9491,0.9481,0.9480,0.9489)
ARMA01<-c(0.9392,0.9455,0.9488,0.9482,0.9488)

plot(ARMA11,xaxs="r",type="l",
     ylab=expression("ACR"),
     xlab=expression(plain(italic(n))),
     axes=FALSE,
     ylim=c(0.92,0.97))

lines(ARMA11, lty=1,col=1)
lines(ARMA10, lty=1,col=7)
lines(ARMA01, lty=1,col=4)

points(ARMA11,pch=1,col=1)
points(ARMA10,pch=2,col=7)
points(ARMA01,pch=3,col=4)

abline(h=0.95,col=2,lty=3)

legend("top",c(expression(plain(BXII-ARMA(1,1))),
                    expression(plain(BXII-ARMA(1,0))), 
                    expression(BXII-ARMA(0,1))),
       lty=c(1,1,1), pch=c(1,2,4),col=c(1,2,4),bty="n",cex=1)
axis(1,1:5,c(70,150,300,500,1000))
box();axis(2)
dev.off() 