# Reference: BXII-GAS
# Created by Fernando J M de Araújo (nandomonteiro418@gmail.com), June/2024

#############################################
# calculating score vector for the BXII-GAS #
#############################################
 
## function for st with fixed mu
st.funcBXII <- function(mu0,y1,c,tau,link="log"){
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (linktemp == "log"){stats<-VGAM::loglink
  } else if (linktemp == "sqrt"){stats<-VGAM::sqrtlink
  } else {
    stop(paste(linktemp, "link not available, available links are \"log\" and \"sqrt\""))
  }
  ###########################################################
  link = linktemp 
  diflink<-function(x)stats(x,inverse = F,deriv = 1)
  ########################
  ut=(1+mu0^c)
  ht= log(1-tau)/log(ut)
  at<-(c^2*mu0^(2*(c-1)))/((ut^2)*(log(ut)^2)) ## loglik/dmu2
  
  nablat<-((-c*mu0^(c-1))/(ut*log(ut)))*(1+(ht*log(1+y1^c)))
  St<-1/at
  nablat*St*diflink(mu0)
  
}

### functions for the analytic derivative
nablat<- function(mu0,y,c,tau,link="log"){
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (linktemp == "log"){stats<-VGAM::loglink
  } else if (linktemp == "sqrt"){stats<-VGAM::sqrtlink
  } else {
    stop(paste(linktemp, "link not available, available links are \"log\", ","\"sqrt\" and \"cloglog\""))
  }
  ###########################################################
  link = linktemp 
  diflink<-function(x)stats(x,inverse = F,deriv = 1)
  ###########################################################
  ut<-(1+mu0^c)
  ht<- (log(1-tau)/log(ut))
  result<-(((-c*mu0^(c-1))/(ut*log(ut)))*(1+(ht*log(1+y^c))))/diflink(mu0)
  return(result)
} 

st_q<-function(mu0,y,c,tau,link=link){
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (linktemp == "log"){stats<-VGAM::loglink
  } else if (linktemp == "sqrt"){stats<-VGAM::sqrtlink
  } else {
    stop(paste(linktemp, "link not available, available links are \"log\", ","\"sqrt\" and \"cloglog\""))
  }
  ###########################################################
  link = linktemp 
  diflink<-function(x)stats(x,inverse = F,deriv = 1)
  diflink2<-function(x)stats(x,inverse = F,deriv = 2)
  ###########################################################
  ut<-(1+mu0^c)
  ht<- (log(1-tau)/log(ut))
  result<-((-((1+(log(ut)*(ut-c))/(c*mu0^c))*(1+ht*log(1+y^c)))+(log(1-tau)*log(1+y^c)))*diflink(mu0)-(((ut*log(ut)*diflink2(mu0))/(c*mu0^(c-1)))*(1+ht*log(1+y^c))))*(1/diflink(mu0))#diflink para completar o s_t^mu
  return(result)
}

dll.dc<- function(mu0,y1,c,tau){
  ut<-(1+mu0^c)
  ht<-(log(1-tau)/log(ut))
  result<- 1/c+log(y1)-(((1-ht)*y1^c*log(y1))/(1+y1^c))-((mu0^c*log(mu0))/(ut*log(ut)))*(1+ht*log(1+y1^c))
  return(result)
}

dst.c <- function(mu0,y,c,tau){
  ut<-(1+mu0^c)
  ht<- (log(1-tau)/log(ut))
  if(mu0==0 || y==0) return(0) else
    result<-(-((mu0*log(ut))/(c^2*mu0^c))*(1+ht*log(1+y^c))*(((c*mu0^c*log(mu0))/(log(ut)))-mu0^c-c*log(mu0)-1)+((ht)/(c*mu0^(c-1)*(y^c+1)))*(-ut*log(ut)*y^c*log(y)+mu0^c*log(mu0)*(1+y^c)*log(1+y^c)))
  return(result)
}

BXIIGAS.score <- function(w,A,B,beta=0,c,y,
                        tau=0.5, ar=NA, ma=NA, X=NA, link = "log")
{
  ##### link function definitions
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (linktemp == "log"){stats<-VGAM::loglink
  } else if (linktemp == "sqrt"){stats<-VGAM::sqrtlink
  } else {
    stop(paste(linktemp, "link not available, available links are \"log\", ","\"sqrt\" and \"cloglog\""))
  }
  ###########################################################
  link = linktemp 
  linkinv = function(x)stats(x,inverse = T)
  diflink<-function(x)stats(x,inverse = F,deriv = 1)
  diflink2<-function(x)stats(x,inverse = F,deriv = 2)
  ###########################################################  
  ##### X definitions
  if(any(is.na(X))==FALSE){
    X<-as.matrix(X)
    if(any(beta==0)) stop("You need to inform beta")
    if(length(beta)!=dim(X)[2]) stop("The length of beta must be equal to the number of columns of X")
    # k = ncol(X)
  }else{
    if(beta!=0) stop("You need to inform X")
    X <- matrix(0, c(n,1))
  }
  ##### defining the lags for st
  if(any(is.na(ar)==T)){
    if (any(is.na(A)==F))
    {
      stop("You must inform the order of A")
    } else {A <- ar <- 0}
  } else {
    if (any(is.na(A)==T)) stop("You must inform A")
    if (!isTRUE(all(ar == floor(ar)))) 
      stop("'ar' must only contain integer values")
    if(length(A)!=length(ar))
      stop("The lengths of A and ar must be equal")
  }
  ##### defining the lags for ft
  if(any(is.na(ma)==T)){
    if(any(is.na(B) ==F ))
    {
      stop("You must inform the order of B")
    } else {B <- ma <- 0}
  } else {
    if (any(is.na(A)==T)) stop("You must inform B")
    if (!isTRUE(all(ma == floor(ma)))) 
      stop("'ma' must only contain integer values")
    if(length(B)!=length(ma))
      stop("The lengths of B and ma must be equal")
  }
  ###########################################################
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  ## initializations
  st<- mu<- f<- ll <- dst.dmu <-nabla<-
    dst<- df.dw <- df.dB <- 
    dc <- df.dc <- rep(0,n)
  
  df.dA <- matrix(0,n,length(A))
  df.dB <- matrix(0,n,length(B))
  df.dbeta <- matrix(0,n,length(beta))
  
  for(i in (m+1):n)
  {
    f[i]  <- w + as.numeric(A%*%st[i-ar]) + 
      as.numeric(B%*%f[i-ma]) + X[i,]%*%beta
    mu[i]   <- linkinv(f[i])
    st[i] <- st.funcBXII(mu[i],y[i],c,tau,link=link)
    dst[i] <-st_q(mu[i],y[i],c,tau,link=link)
    if(eval(mu[i])==0) dfi <- 0 else dfi <- diflink(eval(mu[i]))
    df.dw[i] <- 1+(A*dst[i-ar])%*%df.dw[i-ar]+B%*%df.dw[i-ma]
    df.dA[i,] <- st[i-ar]+t(A*dst[i-ar])%*%df.dA[i-ar,]+B%*%df.dA[i-ma,]
    df.dB[i,]<- f[i-ma]+(A*dst[i-ar])%*%df.dB[i-ar,]+B%*%df.dB[i-ma,]
    df.dbeta[i,] <- X[i,]+(A*dst[i-ar])%*%df.dbeta[i-ar,]+B%*%df.dbeta[i-ma,]
    df.dc[i]<- A%*%dc[i-ar]+B%*%df.dc[i-ma]
    dc[i] <- dst.c(mu[i],y[i],c,tau=tau)%*%dfi+dst[i]*df.dc[i]
  }
  
  Uw <- sum(nablat(mu[-c(1:m)],y[-c(1:m)],c,tau)*df.dw[-c(1:m)])
  UA <- apply(
    as.matrix(nablat(mu[-c(1:m)],y[-c(1:m)],c,tau)*df.dA[-c(1:m),]),
    2,sum)
  UB <- apply(
    as.matrix(nablat(mu[-c(1:m)],y[-c(1:m)],c,tau)*df.dB[-c(1:m),]),
    2,sum)
  Ubeta <- apply(
    as.matrix(nablat(mu[-c(1:m)],y[-c(1:m)],c,tau)*df.dbeta[-c(1:m),]),
    2,sum)
  Uc <- sum(dll.dc(mu[-c(1:m)],y[-c(1:m)],c,tau)+
                   nablat(mu[-c(1:m)],y[-c(1:m)],c,tau)*df.dc[-c(1:m)])
  
  if(any(beta==0)) {rval <- c(Uw,UA,UB,Uc)
  names(rval) <- c("Uw","UA","UB","Uc")
  }
  else {rval <- rval <- c(Uw,UA,UB,Ubeta,Uc)
  names(rval) <- c("Uw",c(paste0("UA", 1:length(A))),
                   c(paste0("UB", 1:length(B))),
                   c(paste0("Ubeta", 1:length(beta))),
                   "Uc")
  }
  
  
  return(rval)
}