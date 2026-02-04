################################################################################
# ETS - Forecast one-step-ahead approach
################################################################################

predict_ETS_1step <- function(y_train, y_test, model = "AAA") {
  
  h1 <- length(y_test)
  y_full <- c(y_train, y_test)
  
  prev <- numeric(h1)
  
  for (i in 1:h1) {
    
    y_est <- ts(
      y_full[1:(length(y_train) + i - 1)],
      frequency = frequency(y_train),
      start = start(y_train)
    )
    
    fit <- ets(y_est, model = model)
    prev[i] <- forecast(fit, h = 1)$mean
  }
  
  return(prev)
}
