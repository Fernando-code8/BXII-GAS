rm(list=ls())
gc()
set.seed(12)
n=7;c=10; tau=.5; y=rgamma(n,.2,.3)
w=-0.1; A1= 0.5; A2=0.2; A=c(A1,A2)
B1= 0.4; B2=0.2; B=c(B1,B2)
beta1=-0.22; X1<-cos(2*pi*1:n/12)
beta2=0.3; X2<-sin(2*pi*1:n/12)
beta=c(beta1,beta2); X <- cbind(X1,X2)
# beta=0; X <- rep(0,n)
## logit link
link<-VGAM::loglink
linkinv = function(x)link(x,inverse = T)
diflink<-function(x)link(x,inverse = F,deriv = 1)
diflink2<-function(x)link(x,inverse = F,deriv = 2)
###
ar<-1:length(A)
ma<-1:length(B)
p <- max(ar)
q <- max(ma)
n <- length(y)
m <- max(p,q,na.rm=T)
p1 <- length(ar)
q1 <- length(ma)
### functions R derivative


st.func <- function(mu0,y1){
  # ut<-(1+mu0^c)
  # ht<- (log(1-tau)/log(ut))
  -(((1+mu0^c)*log((1+mu0^c)))/(c*mu0^(c-1)))*(1+(log(1-tau)/log((1+mu0^c)))*log(1+y1^c))*(1/mu0)
}
ft.func <- function(w,A,B,beta=0,c,y,X=NA,tau=.5){
  ###########################################################
  if(any(is.na(X))==FALSE){
    X<-as.matrix(X)
    if(any(beta==0)) stop("You need to inform beta")
    # k = ncol(X)
  }else{
    if(beta!=0) stop("You need to inform X")
    X <- matrix(0, c(n,1))
  }
  ###########################################################  
  z <- c()
  p1 <- length(A)
  q1 <- length(B)
  p <- 1:p1#max(ar)
  q <- 1:q1#max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  st<-f<-rep(0,n) 
  mu<- rep(0,n)
  
  for(i in (m+1):n)
  {
    f[i]  <- w + as.numeric(A%*%st[i-p]) + 
      as.numeric(B%*%f[i-q]) + X[i,]%*%beta
    mu[i]   <- linkinv(f[i])
    st[i] <- st.func(mu[i],y[i])
  }
  z$f<-f[-c(1:m)]
  z$mu<-mu[-c(1:m)]
  z$st<-st[-c(1:m)]
  
  return(z)
}
s=1
f.expr <- expression(w + A1*st1 + A2*st0 + B1*f1 + B2*f0  + beta1*x1 + beta2*x2)
mu.expr<- expression(exp(ft1))
ll.expr<- expression(
  log((log(1/(1-tau))*c)/(log(1+(mu0^c))))+(c-1)*(log(y1))
  +(log(1-tau)/log(1+mu0^c)-1)*((log(1+y1^c)))
)
### functions for the analytic derivative
nablat<- function(mu0,y,c){
  ut<-(1+mu0^c)
  ht<- (log(1-tau)/log(ut))
  result<-(((-c*mu0^(c-1))/(ut*log(ut)))*(1+(ht*log(1+y^c))))/diflink(mu0)
  return(result)
} 

# ##Completa
# ut<-(1+mu0^c)
# ht<- (log(1-tau)/log(ut))
# stqc<-(((1+(log(ut)*(ut-c))/(c*mu0^c))*(1+ht*log(1+y^c)))+((-log(1-tau)*(log(1+y^c)*c*mu0^(c-1))/(ut*log(ut)))*((ut*log(ut))/(c*mu0^(c-1)))))+(((ut*log(ut))/(c*mu0^(c-1)))*(1+ht*log(1+y^c)))
# ##Simplificada
# stqs<-(((1+(log(ut)*(ut-c))/(c*mu0^c))*(1+ht*log(1+y^c)))-(log(1-tau)*log(1+y^c)))+(((ut*log(ut))/(c*mu0^(c-1)))*(1+ht*log(1+y^c)))

st_q<-function(mu0,y,c){
  ut<-(1+mu0^c)
  ht<- (log(1-tau)/log(ut))
  result<-((-((1+(log(ut)*(ut-c))/(c*mu0^c))*(1+ht*log(1+y^c)))+(log(1-tau)*log(1+y^c)))*diflink(mu0)-(((ut*log(ut)*diflink2(mu0))/(c*mu0^(c-1)))*(1+ht*log(1+y^c))))*(1/diflink(mu0))#diflink para completar o s_t^mu
  return(result)
}

dll.dc<- function(mu0,y1,c){
  ut<-(1+mu0^c)
  ht<-(log(1-tau)/log(ut))
  result<- 1/c+log(y1)-(((1-ht)*y1^c*log(y1))/(1+y1^c))-((mu0^c*log(mu0))/(ut*log(ut)))*(1+ht*log(1+y1^c))
  return(result)
}


##Completa
# ut<-(1+mu0^c)
# ht<- (log(1-tau)/log(ut))
# stdc<- -(((mu0^(1-c))/(c^2))*(mu0^c*log(ut)-c*mu0^c*log(mu0)+c*log(mu0)*log(ut)+log(ut)))*(1+ht*log(1+y^c))+((ht)/(ut*log(ut)*(y^c+1)))*((ut*log(ut))/(c*mu0^(c-1)))*(ut*log(ut)*y^c*log(y)-mu0^c*log(mu0)*(1+y^c)*log(1+y^c))
# ##Simplificada
# stds<-((mu0*log(ut))/(c^2*mu0^c))*(1+ht*log(1+y^c))*(((c*mu0^c*log(mu0))/(log(ut)))-mu0^c-c*log(mu0)-1)+((ht)/(c*mu0^(c-1)*(y^c+1)))*(ut*log(ut)*y^c*log(y)-mu0^c*log(mu0)*(1+y^c)*log(1+y^c))

dst.c <- function(mu0,y,c){
  ut<-(1+mu0^c)
  ht<- (log(1-tau)/log(ut))
  if(mu0==0 || y==0) return(0) else
  result<-(-((mu0*log(ut))/(c^2*mu0^c))*(1+ht*log(1+y^c))*(((c*mu0^c*log(mu0))/(log(ut)))-mu0^c-c*log(mu0)-1)+((ht)/(c*mu0^(c-1)*(y^c+1)))*(-ut*log(ut)*y^c*log(y)+mu0^c*log(mu0)*(1+y^c)*log(1+y^c)))
  return(result)
  }

## initializations
st<- mu<- f<- ll <- dst.dmu <-nabla<-
  dst<- df.dw <- df.dA1 <- df.dA2 <- 
  df.dB1 <- df.dB2 <- df.dbeta1 <-
  df.dbeta2 <- dc <- df.dc <-
  Uw<- UA1<- UA2<- UB1 <- UB2 <- Ubeta1 <- Ubeta2 <-
  Uc<- rep(0,n)

df.dA <- matrix(0,n,length(A))

for(i in (m+1):n)
{
  y1<-y[i]
  L1<-list(st0=st[i-2][[1]],st1=st[i-1][[1]],
           f0=f[i-2][[1]], f1=f[i-1][[1]],
           x1=X[i,1],x2=X[i,2])
  f[i] <-as.expression(do.call(substitute, list(f.expr[[1]], L1)))
  # f[i]  <- w + A*st[i-p][[1]] + B*f[i-q] 
  L2<-list(ft1=f[i][[1]])
  mu[i]   <- as.expression(do.call(substitute, list(mu.expr[[1]], L2)))
  L3<-list(mu0=mu[i][[1]],y1=y[i])
  st[i] <- as.expression(do.call(substitute, list(body(st.func)[[2]], L3)))
  ll[i] <- as.expression(do.call(substitute, list(ll.expr[[1]], L3)))
  Uw[i] <-eval(as.expression(D(ll[i],"w")))
  UA1[i] <-eval(as.expression(D(ll[i],"A1")))
  UA2[i] <-eval(as.expression(D(ll[i],"A2")))
  UB1[i] <-eval(as.expression(D(ll[i],"B1")))
  UB2[i] <-eval(as.expression(D(ll[i],"B2")))
  Ubeta1[i]<-eval(as.expression(D(ll[i],"beta1")))
  Ubeta2[i]<-eval(as.expression(D(ll[i],"beta2")))
  Uc[i] <-eval(as.expression(D(ll[i],"c")))
  L4 <-list(mu0=mu[i-1][[1]],y1=y[i-1])
  dst[i] <-st_q(eval(mu[i]),y1,c)
  if(eval(mu[i])==0) dfi <- 0 else dfi <- diflink(eval(mu[i]))
  df.dw[i] <- 1+(A*dst[i-ar])%*%df.dw[i-ar]+B%*%df.dw[i-ma]
  df.dA1[i] <-  eval(st[i-1])+(A1*dst[i-1])*df.dA1[i-1]+
    (A2*dst[i-2])*df.dA1[i-2]+B%*%df.dA1[i-ma]
  df.dA2[i] <- eval(st[i-2])+(A1*dst[i-1])*df.dA2[i-1]+
    (A2*dst[i-2])*df.dA2[i-2]+B1*df.dA2[i-1]+B2*df.dA2[i-2]
  df.dA[i,] <- c(eval(st[i-1]),eval(st[i-2]))+
    t(A*dst[i-ar])%*%df.dA[i-ar,]+B%*%df.dA[i-ma,]
  df.dB1[i]<- eval(f[i-1])+(A*dst[i-ar])%*%df.dB1[i-ar]+
    B1*df.dB1[i-1]+B2*df.dB1[i-2]
  df.dB2[i]<- eval(f[i-2])+(A*dst[i-ar])%*%df.dB2[i-ar]+
    B%*%df.dB2[i-ma]
  df.dbeta1[i] <- X[i,1]+(A*dst[i-ar])%*%df.dbeta1[i-ar]+B%*%df.dbeta1[i-ma]
  df.dbeta2[i] <- X[i,2]+(A*dst[i-ar])%*%df.dbeta2[i-ar]+B%*%df.dbeta2[i-ma]
  df.dc[i]<- A%*%dc[i-ar]+B%*%df.dc[i-ma]
  dc[i] <- dst.c(eval(mu[i]),y[i],c)%*%dfi+dst[i]*df.dc[i]
}
values<-ft.func(w,A,B,beta=beta,c=c,y,X)

round(nablat(values$mu,y[-(1:m)],c)*df.dw[-(1:m)],6)==round(Uw[-(1:m)],6)
round(nablat(values$mu,y[-(1:m)],c)*df.dA1[-(1:m)],6)==
  round(UA1[-(1:m)],6)
round(nablat(values$mu,y[-(1:m)],c)*df.dA2[-(1:m)],6)==
  round(UA2[-(1:m)],6)
round(nablat(values$mu,y[-(1:m)],c)*df.dB1[-(1:m)],6)==
  round(UB1[-(1:m)],6)
round(nablat(values$mu,y[-(1:m)],c)*df.dB2[-(1:m)],6)==
  round(UB2[-(1:m)],6)
round(nablat(values$mu,y[-(1:m)],c)*df.dbeta1[-(1:m)],6)==
  round(Ubeta1[-(1:m)],6)
round(nablat(values$mu,y[-(1:m)],c)*df.dbeta2[-(1:m)],6)==
  round(Ubeta2[-(1:m)],6)
round(dll.dc(values$mu,y[-(1:m)],c)+
        nablat(values$mu,y[-(1:m)],c)*df.dc[-(1:m)],6)==
  round(Uc[-(1:m)],6)

setwd("C:/Users/nando/OneDrive/Desktop/mestrado UFRGS/Dissertação/All scripts/UG-GAS e BXII-GAS new version/Scripts Atual BXII-GAS/Git/BXII-GAS/BXII-GAS")
source("BXIIGASfuncs.R")
round(BXIIGAS.score(w,A=A,B=B,beta=beta,c,y,tau=.5,
                  ar=ar,ma=ma,X=X,link = "log"),6)==
  round(c(sum(Uw),sum(UA1),sum(UA2),sum(UB1),sum(UB2),sum(Ubeta1),sum(Ubeta2),sum(Uc)),6)

