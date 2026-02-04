################################################################################
# BXII-ARMA - Forecast one-step-ahead approach
################################################################################

predict_BXIIARMA<-function(fit_BXIIARMA,y_test,X=X,X_hat=X_hat,tau=0.5,link="log",serie=y){
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
  
  link = linktemp 
  linkinv = function(x)stats(x,inverse = T)
  
  alpha<-fit_BXIIARMA$coeff[1]
  phi<-fit_BXIIARMA$phi
  theta<-fit_BXIIARMA$theta
  c<-fit_BXIIARMA$c
  y_train<-serie
  n<-length(y_train)
  h1<-length(y_test)
  #### Forecasting
  ar<-1:length(phi)
  ma<-1:length(theta)
  namesphi <- names(fit_BXIIARMA$coeff[2+ar]);orderphi <- as.numeric(gsub("phi", "", namesphi))
  namestheta <- names(fit_BXIIARMA$coeff[2+length(phi)+ma]);ordertheta <- as.numeric(gsub("theta", "", namestheta))
  m<-max(c(orderphi,ordertheta))
  
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
  if(k==0)  beta <- as.matrix(0) else beta<-fit_BXIIARMA$beta
  
  
  ###############################
  
  y_t<-c(y_train[(n-m):(n-1)],y_test) # taking the last value in the train test
  
  #### Forecasting
  errorhat<- c(fit_BXIIARMA$errorhat[(n-m):(n-1)],rep(NA,h1))
  ynew_prev <-y_prev <- c()
  
  if(m==1)  X_prev<- rbind(rbind(X[(n-m):(n-1),]),cbind(X_hat)) else X_prev<- rbind(cbind(X[(n-m):(n-1),]),cbind(X_hat))
  
  for(i in 1:h1){
    ynew_prev[i] <- alpha + X_prev[m+i,]%*%as.matrix(beta) +
      (phi%*%(y_t[m+i-ar]-X_prev[m+i-ar,]%*%as.matrix(beta))) +
      (theta%*%errorhat[m+i-ma])
    errorhat[m+i]<- y_t[m+i]-ynew_prev[i] # predictor scale
    y_prev[i] <- linkinv(ynew_prev[i])
  }
  return(y_prev)
}
