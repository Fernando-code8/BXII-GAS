# Reference: Ray-GAS
# Created by Miguel Peña-Ramírez (miguepe@gmail.com), march/2023

########################################
#### R packages and source archives ####
########################################
source("rayleigh.R")
library("rootSolve")

##################
#### Function ####
##################
GASRAY.fit <- function (y, ar, ma,X=NA,names_A=paste0("A",length(ar)),names_B=paste0("B",length(ma)), link = "log", h1=1)#, names_phi,names_theta,names_beta,diag,h1,X,X_hat,resid)
{
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("log", "sqrt", "identity"))){
    stats <- make.link(linktemp)
  }  else {
    stop(paste(linktemp, "link not available, available links are \"log\", ","\"sqrt\" and \"identity\""))
  }
  
  link = linktemp 
  linkfun = stats$linkfun
  linkinv = stats$linkinv 
  mu.eta = stats$mu.eta 
  diflink = mu.eta
  ynew = linkfun(y)
  
  maxit1<-10000
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)

  ###########################################################
  if(any(is.na(X))==FALSE){
    if(any(is.na(X_hat))==TRUE) 
      stop("You need to inform X_hat")
    X<-as.matrix(X)
    X_hat<-as.matrix(X_hat)
    k = ncol(X)
    names_theta <- c(paste("theta", 1:k, sep = ""))
  }else{
    X <- matrix(0, c(n,1))
    X_hat<- as.matrix(rep(0,h1+1))
    k=0
    names_theta <- NULL
  }
  ###########################################################  
  
  # ##### GASp model
  if (any(is.na(p) ==F ))
  {
    arp <- 1:p
  }
  # ##### GASq model
  if(any(is.na(q) ==F ))
  {
    maq <-1 :q
  }
  ######### GAS pq model
  { 
    loglik <- function(z){
      w <- z[1]
      A <- z[2:(p1+1)]
      B <- z[(p1+2):(p1+q1+1)]
      
      st<-f<-rep(0,n+m) 
      mu<- rep(0,n+m) 
      
      for(i in (m+1):n){
        f[i]  <- w + as.numeric(A%*%st[i-ar]) + as.numeric(B%*%f[i-ma])
        mu[i]   <- linkinv(f[i])
        st[i] <- 1/4*(pi*y[i]^2/(2*mu[i]^2)-2)
      }
      mu <- linkinv(f[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- sum (
        log(pi/2)+log(y1)-log(mu^2)-(pi*y1^2)/(4*(mu^2))
      )
    }
    
    opt <- try(
      opt <- optim(c(log(mean(y)),0,0), loglik, method = "BFGS", hessian = T,
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
      ,silent = T)
    if(class(opt)=="try-error"){
      opt <- suppressWarnings(
        opt <- optim(c(log(mean(y)),.1,.1), loglik,method = "BFGS", hessian = T,
                     control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
      )
    }
    
    #### SAÍDAS DO Z
    z <- c()
    z$fgrad <- gradient(loglik, opt$par)
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p+q+1)]
    z$coeff <- coef
    z$linkfun <- linkfun #tmp
    z$linkinv <-linkinv #tmp
    w <- coef[1]
    A <- coef[2:(p+1)]
    B <- coef[(p+2):(p+q+1)]
    z$w <- w
    z$A <- A
    z$B <- B
    J_inv <- solve(-(opt$hessian))
    z$stderror<-sqrt(diag(J_inv))
    z$zstat <-abs(z$coef/z$stderror)
    z$pvalues <- 2*(1 -pnorm(z$zstat))
    
    z$loglik<-opt$value
    z$k<- (p1+q1+1+k)
    
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    
    model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),
                                round(z$zstat,4),round(z$pvalues,4))
    colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
    z$model <- model_presentation
    ####### FORECAST
    fhat <- sthat <-rep(0,n) 
    muhat<-rep(NA,n+m)

    for(i in (m+1):n)
    {
      fhat[i] <- w + as.numeric(A%*%sthat[i-p]) + as.numeric(B%*%fhat[i-q])
      muhat[i]   <- linkinv(fhat[i])
      sthat[i] <- 1/4*(pi*y[i]^2/(2*muhat[i]^2)-2)
    }
    y1<-y[(m+1):n]
    
    z$fitted <- ts(muhat,start=start(y),frequency=frequency(y))
    z$fhat <- fhat
    z$sthat <- sthat
    z$muhat <- muhat #tmp
    
    names_par <- c("w",names_A,names_B)
    z$serie <- y
    z$GAS <- names_par
    
    ynew_prev <- c(y[n],rep(NA,h1))
    
    fhatf <- c(fhat[n],rep(0,h1))
    sthatf <-c(sthat[n],rep(0,h1))
    muhatf<- c(muhat[n],rep(0,h1))
    
    for(i in 1:(h1))
    {
      fhatf[n+i] <-  w + as.numeric(A%*%sthatf[n+i-arp]) + as.numeric(B%*%fhatf[n+i-maq])
      muhatf[n+i]   <- linkinv(fhatf[n+i])
      ynew_prev[n+i] <- muhatf[n+i]
      nabla <-  (pi*ynew_prev[n+i]^2)/(2*muhatf[n+i]^2)-(2)
      St <- 1/4
      sthatf[n+i] <-nabla*St
      
    }
    z$forecast <- ynew_prev[n+1:(n+h1)]
    
    ##############################################
    # residuals
    z$residuals<-rep(0,n) 
    for(i in (m+1):n){
      z$residuals[i] <- qnorm(pr(y[i],z$fitted[i])) # quantile residuals
    }
  }
  return(z)
#}
}


