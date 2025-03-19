# Reference: Ray-GAS
# Created by Miguel Peña-Ramírez (miguepe@gmail.com), march/2023

predict_rarma<-function(fit_rarma,y_test){
  alpha<-fit_rarma$alpha
  phi<-fit_rarma$phi
  theta<-fit_rarma$theta
  y_train<-fit_rarma$serie
  n=length(y_train)
  n<-length(y_train)
  h1<-length(y_test)
  #### Forecasting
  ar<-1:length(phi)
  ma<-1:length(theta)
  m<-max(length(ar),length(ma))
  y_t<-fit_rarma$link$linkfun(c(y_train[(n-m):(n-1)],y_test)) # taking the last value in the train test
  errorhat<- c(fit_rarma$errorhat[(n-m):(n-1)],rep(NA,h1))
  ynew_prev <-y_prev <- c()
  
  for(i in 1:h1){
    ynew_prev[i] <- alpha + (phi%*%y_t[m+i-ar]) + (theta%*%errorhat[m+i-ma])
    errorhat[m+i]<- y_t[m+i]-ynew_prev[i] # predictor scale
    y_prev[i] <- fit_rarma$link$linkinv(ynew_prev[i])
  }
  return(y_prev)
}



