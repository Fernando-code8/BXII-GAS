################################################################################
# gamma-GAS - Forecast one-step-ahead approach
################################################################################

predict_GGAS<-function(fit_GGAS,y_test,X=XGAMMA,X_hat=XGAMMA_hat,tau=0.5,link="log"){
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

  omega<-fit_GGAS$coeff[1]
  A<-fit_GGAS$A
  B<-fit_GGAS$B
  lambda<-fit_GGAS$lambda
  y_train<-fit_GGAS$serie
  n<-length(y_train)
  h1<-length(y_test)
  #### Forecasting
  ar<-1:length(A)
  ma<-1:length(B)
  namesA <- names(fit_GGAS$coeff[1+ar]);orderA <- as.numeric(gsub("A", "", namesA))
  namesB <- names(fit_GGAS$coeff[1+length(A)+ma]);orderB <- as.numeric(gsub("B", "", namesB))
  m<-max(c(orderA,orderB))

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
  if(k==0)  beta <- as.matrix(0) else beta<-fit_GGAS$beta


  ###############################

  y_t<-c(y_train[(n-m):(n-1)],y_test) # taking the last value in the train test

  #### Forecasting

  muhat<- c(fit_GGAS$fitted[(n-m):(n-1)],rep(NA,h1))
  sthat<- c(fit_GGAS$sthat[(n-m):(n-1)],rep(NA,h1))
  fthat<- c(fit_GGAS$fhat[(n-m):(n-1)],rep(NA,h1))
  ynew_prev <- c()
  
  if(m==1)  X_prev<- rbind(rbind(X[(n-m):(n-1),]),cbind(X_hat)) else X_prev<- rbind(cbind(X[(n-m):(n-1),]),cbind(X_hat))
  
  for(i in 1:h1){
    ynew_prev[i] <- omega + as.numeric(A%*%sthat[m+i-orderA]) + as.numeric(B%*%fthat[m+i-orderB]) +  X_prev[m+i,]%*%beta
    fthat[m+i] <- ynew_prev[i]
    muhat[m+i] <- linkinv(ynew_prev[i])
    sthat[m+i] <- st.funcGGAS(muhat[m+i],y_t[m+i],lambda=lambda,tau=tau,link=link)
  }
  
  return(muhat[-c(1:m)])
}