# Reference: Ray-GAS
# Created by Miguel Peña Ramirez (miguepe@ufsm.br), february/2023
# set.seed(5)
# n=100
# A=.4
# b=.6
# w=1
# link="log"
# X=NA
# beta=0
source("rarma.r")
simu.Rarma <- function(n,phi=NA,theta=NA,alpha=0.0,freq=0,link="log",beta=NA,X=NA){
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("log", "probit", "cloglog"))){
    stats <- make.link(linktemp)
  }  else {
    stop(paste(linktemp, "link not available, available links are \"log\", ","\"probit\" and \"cloglog\""))
  }
  
  link = linktemp 
  linkfun = stats$linkfun
  linkinv = stats$linkinv 
  mu.eta = stats$mu.eta 
  diflink = mu.eta
  ar<-NA
  ma<-NA
  beta<-as.matrix(beta,1,2)
  
  ##### X definitions
  if(is.na(X)){
    X<-matrix(0, c(n,1))
  }else{
    if(X=="cos"){
      X=as.matrix(cos(2*pi*(1:n)/12))
      if(beta==0) stop("Inform the value of beta")
    }else
      if(X=="sin"){
        X=as.matrix(sin(2*pi*(1:n)/12))
        if(beta==0) stop("Inform the value of beta")
      }else
        if(X=="sin&cos"){
          X=cbind(sin(2*pi*(1:n)/12),cos(2*pi*(1:n)/12))
          if(dim(beta)[1]!=2) stop("Inform the value of beta 2")
        }
  }
  if (any(is.na(phi) ==F ))
  {
    ar <- 1:length(phi)
  }
  # ##### GASq model
  if(any(is.na(theta) ==F ))
  {
    ma <-1 : length(theta)
  }
  # ##### GASpq model
  if(any(is.na(phi) ==F ) && any (is.na(theta) ==F ))
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)
    
    error<-rep(0,n+m) # E(error)=0 
    eta<- y <- rep(NA,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rr(mu[i])
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]   
    }
   }  
    return(ts(y[(m+1):(n+m)]))
  } # RARMA model
# res<-simu.Rarma(200,0.64,0.24,-1.73)
# plot(res)
