################################################################################
# SARMA - Forecast one-step-ahead approach
################################################################################

predict_SARMA_1step <- function(
    y_train,
    y_test,
    order = c(1,0,1),
    seasonal = c(1,0,1),
    period
) {
  
  h1 <- length(y_test)
  y_full <- c(y_train, y_test)
  
  prev <- rep(NA, h1)
  
  for (i in 1:h1) {
    
    y_est <- ts(
      y_full[1:(length(y_train) + i - 1)],
      frequency = frequency(y_train),
      start = start(y_train)
    )
    
    fit <- try(
      Arima(
        y_est,
        order = order,
        seasonal = list(order = seasonal, period = period),
        method = "ML"
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      prev[i] <- forecast(fit, h = 1)$mean
    }
  }
  
  return(prev)
}
