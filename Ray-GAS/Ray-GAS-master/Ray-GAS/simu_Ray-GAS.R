library(VGAM)
source("Ray-mu.R")
source("RayGASfuncs.R")
simu.RayGAS <- function(n, w=1 , A=NA, B=NA, beta=0,
                        ar=NA, ma=NA, X=NA, link = "log")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (linktemp == "log"){stats<-loglink
  } else if (linktemp == "sqrt"){stats<-sqrtlink
  } else {
    stop(paste(linktemp, "link not available, available links are \"log\" and \"sqrt\""))
  }
  
  link = linktemp 
  linkinv = function(x)stats(x,inverse = T)
  diflink<-function(x)stats(x,inverse = F,deriv = 1)
  diflink2<-function(x)stats(x,inverse = F,deriv = 2)
  beta<-as.matrix(beta,1,2)
  
  ##### X definitions
  if(is.na(X)){
    X<-matrix(0, c(n,1))
    if(beta!=0) stop("Inform the value of X")
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
  ##### defining the lags for st
  if(any(is.na(ar)==T)){
    if (any(is.na(A) ==F ))
    {
      ar <- 1:length(A)
    } else {A <- ar <- 0}
  } else {
    if (!isTRUE(all(ar == floor(ar)))) 
      stop("'ar' must only contain integer values")
  }
  
  ##### defining the lags for ft
  if(any(is.na(ma)==T)){
    if(any(is.na(B) ==F ))
    {
      ma <-1 : length(B)
    } else {B <- ma <- 0}
  } else {
    if (!isTRUE(all(ma == floor(ma)))) 
      stop("'ma' must only contain integer values")
  }
  
  
  # ##### GASpq model
  {
    p <- max ( ar )
    q <- max ( ma )
    m <- max ( p, q, na.rm=T )
    p1 <- length( ar )
    q1 <- length( ma )
    # inicializa os vetores dos dados
    st <- f <- rep(0, n+m) # E(erro) =0
    mu <- y <- rep(NA, n+m)
    
    for ( i in (m+1) : (n+m)){
      f[i] <- w + as.numeric(A%*%st[i-ar]) + as.numeric(B%*%f[i-ma]) + X[i-m,]%*%beta
      mu[i] <- linkinv(f[i])
      y[i] <- rr(1,mu[i])
      
      ########################
      # ut[i]=(1+mu[i]^c)
      # ht[i]= log(1-tau)/log(ut[i])
      # vh[i] <- ((-c*mu[i]^(c-1))/(ut[i]*log(ut[i])))*(1+(ht[i]*log(1+y[i]^c))) ## loglik/dmu
      # at[i]<-(c^2*mu[i]^(2*(c-1)))/((ut[i]^2)*(log(ut[i])^2))
      # ########################
      # 
      # st[i] <- at[i]*vh[i]*diflink(mu[i])
        st[i] <- st.funcRay(mu[i],y[i],link=link)
      
      
    }
  }
  return(ts(y[-(1:m)],frequency = 12))
  
}



# set.seed(123)
# y<-simu.RayGAS(1000,.1,.1)
# y
# plot(y)
