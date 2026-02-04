################################################################################
# ARMA - Forecast one-step-ahead approach
################################################################################

predict_ARMA_1step <- function(y_train, y_test, order = c(1,0,1)) {
  
  h1 <- length(y_test)
  y_full <- c(y_train, y_test)
  
  prev <- numeric(h1)
  
  for (i in 1:h1) {
    
    y_est <- ts(
      y_full[1:(length(y_train) + i - 1)],
      frequency = frequency(y_train),
      start = start(y_train)
    )
    
    fit <- Arima(y_est, order = order)
    prev[i] <- forecast(fit, h = 1)$mean
  }
  
  return(prev)
}