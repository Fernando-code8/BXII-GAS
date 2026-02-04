################################################################################
# PAPER: Quantile-based dynamic modeling of asymmetric data: a novel Burr XII
#        approach for positive continuous random variables
# SUBSECTION: 6. Applications
# GOAL: Forecast estimated using a one-step-ahead approach
# AUTHORS: Fernando Jose Monteiro de Araujo, Renata Rojas Guerra and 
#          Fernando Arturo Pena-Ramirez
# LAST UPDATE: May 25, 2024
################################################################################


# 6.1 Time seriesmodel for finance data

predict_bxii<-function(fit_bxii){
  alpha<-fit_bxii$coeff[1]
  phi<-fit_bxii$phi
  theta<-fit_bxii$theta
  #### Forecasting
  h1<-length(datatest)
  ar<-1:length(phi)
  ma<-1:length(theta)
  m<-max(length(ma),length(ar))
  y_t<-fit_bxii$stats$linkfun(c(yy[(n-m):(n-2)],dadosh))
  errorhat<- c(fit_bxii$errorhat[(n-m):(n-1)],rep(NA,h1))
  ynew_prev <-y_prev <- c()
  
  for(i in 1:h1){
    ynew_prev[i] <- alpha  + (phi%*%y_t[m+i-ar])+(theta%*%errorhat[m+i-ma])
    errorhat[m+i]<- y_t[m+i]-ynew_prev[i] # predictor scale
    y_prev[i] <- fit_bxii$stats$linkinv(ynew_prev[i])
  }
  return(y_prev)
}

# 6.2 Time seriesmodel formeteorological data

predict_bxii_wind<-function(fit_bxii){
  alpha<-fit_bxii$coeff[1]
  phi<-fit_bxii$phi
  theta<-fit_bxii$theta
  #### Forecasting
  h1<-length(datatest)
  ar<-1:length(phi)
  ma<-1:length(theta)
  m<-max(length(ma))
  y_t<-fit_bxii$stats$linkfun(c(yy[(n-m):(n-2)],dadosh))
  errorhat<- c(fit_bxii$errorhat[(n-m):(n-1)],rep(NA,h1))
  ynew_prev <-y_prev <- c()
  
  for(i in 1:h1){
    ynew_prev[i] <- alpha +(phi%*%y_t[m+i-ar])+(theta%*%errorhat[m+i-ma])
    errorhat[m+i]<- y_t[m+i]-ynew_prev[i] # predictor scale
    y_prev[i] <- fit_bxii$stats$linkinv(ynew_prev[i])
  }
  return(y_prev)
}