# Reference: gamma-GAS
# Created by Fernando J M de Araújo (nandomonteiro418@gmail.com), June/2024

source("GGAS-alpha-lambda.R")
source("GGASfuncs.R")

GGAS.fit <- function (y, ar, ma,X=NA, X_hat=NA, tau=0.5 ,link = "log", h1=1)
{
  if (min(y) < 0)
    stop("OUT OF SUPPORT!")
  
  ##### defining the lags for st
  if(any(is.na(ar)==T)){
    ar <- 0
  } else {
    if (!isTRUE(all(ar == floor(ar)))) 
      stop("'ar' must only contain integer values")
  }
  ##### defining the lags for ft
  if(any(is.na(ma)==T)){
    ma <- 0
  } else {
    if (!isTRUE(all(ma == floor(ma)))) 
      stop("'ma' must only contain integer values")
  }
  
  maxit1<-10000
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
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
    stop(paste(linktemp, "link not available, available links are \"log\" and \"cloglog\""))
  }
  
  link = linktemp 
  linkinv = function(x)stats(x,inverse = T)
  diflink<-function(x)stats(x,inverse = F,deriv = 1)
  diflink2<-function(x)stats(x,inverse = F,deriv = 2)
  
  
  ###########################################################
  names_A <- c(paste("A", ar, sep = ""))
  names_B <- c(paste("B", ma, sep = ""))
  
  if(any(is.na(X))==FALSE){
    if(any(is.na(X_hat))==TRUE) 
      stop("You need to inform X_hat")
    X<-as.matrix(X)
    X_hat<-as.matrix(X_hat)
    k = ncol(X)
    names_beta <- c(paste("beta", 1:k, sep = ""))
  }else{
    X <- matrix(0, c(n,1))
    X_hat<- as.matrix(rep(0,h1+1))
    k=0
    names_beta <- NULL
  }
  ###########################################################  
  ######### GAS pq model
  loglik <- function(z) 
  {
    w <- z[1]
    A <- z[2:(p1+1)]
    B <- z[(p1+2):(p1+q1+1)]
    if(k==0)  beta <- as.matrix(0) else beta <- as.matrix(z[(p1+q1+2):(p1+q1+1+k)])
    lambda <- z[p1+q1+k+2]
    st<-f<-rep(0,n) 
    alpha<- rep(0,n)
    
    for(i in (m+1):n)
    {
      f[i]  <- w + as.numeric(A%*%st[i-ar]) + as.numeric(B%*%f[i-ma]) + X[i,]%*%beta
      alpha[i]   <- linkinv(f[i])
      st[i] <- st.funcGGAS(alpha[i],y[i],lambda,tau,link=link)
    }
    alpha   <- linkinv(f[(m+1):n])
    y1<-y[(m+1):n]
    
    ll <- suppressWarnings(log(dgammaGAS(y1, alpha, lambda)))
    
    sum(ll)
  } 

  
  opt <- try(
    optim(c(rep(0,(p1+q1+k+2)))
          , loglik, #score.func, 
          method = "BFGS", hessian = T,
          control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    ,silent = T)
  
  
  
  if(class(opt)=="try-error"){
    opt <- suppressWarnings(
      optim(c(rep(.1,(p1+q1+k+2))), 
            loglik, #score.func, 
            method = "BFGS", hessian = T,
            control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    )
    
    
  }
  
  if (opt$conv != 0)
  {
    warning("FUNCTION DID NOT CONVERGE WITH ANALYTICAL GRADIENT!")
    opt <- optim(c(rep(.001,(p1+q1+k+2))), loglik, method = "BFGS", hessian = T,
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    
    
    if (opt$conv != 0)
    {
      warning("THE FUNCTION DID NOT CONVERGE WITH EITHER NUMERICAL OR ANALYTICAL GRADIENTS!")
    }else{
      warning("IT WORKS WITH NUMERICAL GRADIENT!")
    }
  }
  ########################################################### 
  #### SAÍDAS DO Z
  z <- c()
  z$conv <- opt$conv
  coef <- opt$par
  w <- coef[1]
  A <- coef[2:(p1+1)]
  B <- coef[(p1+2):(p1+q1+1)]
  if(k==0)  beta <- as.matrix(0) else beta <- coef[(p1+q1+2):(p1+q1+1+k)]
  lambda <- coef[p1+q1+k+2]
  
  names(coef)<-c("omega",names_A,names_B,names_beta,"lambda")
  z$coeff <- coef
  J_inv <- solve(-(opt$hessian))
  
  z$w <- w
  z$A <- A
  z$B <- B
  z$beta <- beta
  z$lambda <- lambda
  
  z$stderror<-sqrt(diag(J_inv))
  z$zstat <- abs(z$coef/z$stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat) )
  
  z$loglik<-opt$value
  z$k<- (p1+q1+k+2)
  z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
  z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
  z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  
  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  z$model <- model_presentation
  
  ########################################################### 
  #### FORECAST
  fhat <- sthat <-rep(0,n) 
  alphahat<-rep(NA,n)
  
  for(i in (m+1):n)
  {
    fhat[i] <- w + as.numeric(A%*%sthat[i-ar]) + as.numeric(B%*%fhat[i-ma]) + X[i,]%*%beta
    alphahat[i]   <- linkinv(fhat[i])
    sthat[i] <- st.funcGGAS(alphahat[i],y[i],lambda,tau,link=link)
  }
  y1<-y[(m+1):n]
 
  
  z$fitted <- ts(alphahat,start=start(y),frequency=frequency(y))
  z$fhat <- fhat
  z$sthat <- sthat
  z$serie <- y
  
  ################################################
  # ynew_prev <- c(y[n],rep(NA,h1))
  # 
  # fhatf <- c(fhat[n],rep(0,h1))
  # sthatf <-c(sthat[n],rep(0,h1))
  # alphahatf<- c(alphahat[n],rep(0,h1))
  # 
  # for(i in 2:(h1+1))
  # {
  #   fhatf[i] <-  w + as.numeric(A%*%sthatf[i-ar]) +
  #     as.numeric(B%*%fhatf[i-q]) +  X_hat[i-1,]%*%beta
  #   alphahatf[i]   <- linkinv(fhatf[i])
  #   ynew_prev[i] <- alphahatf[i]
  #   sthatf[i] <-  st.funcGGAS(alphahatf[i],y[i],lambda,tau,link=link)
  # }
  # z$forecast <- ynew_prev[2:(h1+1)]
  
  
  ###############################
  y_prev1 <- c(rep(NA,(n+h1)))
  fhatf1 <- c(rep(NA,(n+h1)))
  sthatf1 <-c(rep(NA,(n+h1)))
  alphahat1<- c(rep(NA,(n+h1)))
  
  fhatf1[1:n] <- fhat
  sthatf1[1:n] <-sthat
  alphahat1[1:n]<- alphahat
  
  #### Forecasting
  ynew_prev1 <- c(y,rep(NA,h1))
  ynew_fhatf1 <- c(fhat,rep(NA,h1))
  ynew_sthatf1 <- c(sthat,rep(NA,h1))
  ynew_alphahat1 <- c(alphahat,rep(NA,h1))
  y_prev1[1:n] <- z$fitted
  
  X_prev<- rbind(X,X_hat)
  
  for(i in 1:h1)
  {
    ynew_fhatf1[n+i] <-  w + as.numeric(A%*%ynew_sthatf1[n-i-ar]) + as.numeric(B%*%ynew_fhatf1[n-i-ma]) +  X_prev[n+i,]%*%beta
    ynew_alphahat1[n+i]   <- linkinv(ynew_fhatf1[n+i])
    ynew_prev1[n+i] <- ynew_alphahat1[n+i]
    ynew_sthatf1[n+i] <- st.funcGGAS(ynew_alphahat1[n+i],ynew_prev1[n+i],lambda,tau,link=link)
  }
  
  z$forecast <- ynew_prev1[(n+1):(n+h1)] 
  #######################################################
  
  # residuals
  z$residuals<-rep(0,n) 
  for(i in (m+1):n){
    z$residuals[i] <- qnorm(pgammaGAS(y[i],z$fitted[i],lambda))
  }
  return(z)
}