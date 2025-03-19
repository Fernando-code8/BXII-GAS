# Reference: Ray-GAS
# Created by Miguel Peña-Ramírez (miguepe@gmail.com), march/2023

predict_rayGAS<-function(fit_rayGAS,y_test){
  w<-fit_rayGAS$w
  A<-fit_rayGAS$A
  B<-fit_rayGAS$B
  y_train<-ray_GAS$serie
  n<-length(y_train)
  h1<-length(y_test)
  #### Forecasting
  ar<-1:length(A)
  ma<-1:length(B)
  m<-max(length(ar),length(ma))
  y_t<-c(y_train[(n-m):(n-1)],y_test) # taking the last value in the train test
  muhat<- c(fit_rayGAS$muhat[(n-m):(n-1)],rep(NA,h1)) 
  sthat<- c(fit_rayGAS$sthat[(n-m):(n-1)],rep(NA,h1)) 
  fthat<- c(fit_rayGAS$fhat[(n-m):(n-1)],rep(NA,h1)) 
  ynew_prev <- c()
  for(i in 1:h1){
    ynew_prev[i] <- w + as.numeric(A%*%sthat[m+i-ar]) + 
      as.numeric(B%*%fthat[m+i-ma])
    fthat[m+i] <- ynew_prev[i]
    muhat[m+i] <- fit_rayGAS$linkinv(ynew_prev[i])
    sthat[m+i] <- 1/4*(pi*y_t[m+i]^2/(2*muhat[m+i]^2)-2)
                  
  }
  return(muhat[-1])
}





