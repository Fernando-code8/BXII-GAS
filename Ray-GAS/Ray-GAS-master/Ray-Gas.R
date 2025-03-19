# Reference: Ray-GAS
# Created by Miguel Peña-Ramírez (miguepe@gmail.com), march/2023

###########################
#### R source archives ####
###########################
source("rayleigh.R")

##################
#### Function ####
##################
simu.GAS <- function(n, A=.1, B=.1, w=1 , link = "log"){
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
  ar<-NA
  ma<-NA
  
  # ##### GASp model
  if (any(is.na(A) ==F ))
  {
    ar <- 1:length(A)
  }
  # ##### GASq model
  if(any(is.na(B) ==F ))
  {
    ma <-1 : length(B)
  }
  # ##### GAS model
  if(any(is.na(A) ==F ) && any (is.na(B) ==F ))
  {
    p <- max ( ar )
    q <- max ( ma )
    m <- max ( p, q, na.rm=T )
    p1 <- length( ar )
    q1 <- length( ma )
    st <- f <- rep(0, n+m) 
    mu <- y <- rep(NA, n+m)
    
    for ( i in (m+1) : (n+m)){
      f[i] <- w + as.numeric(A%*%st[i-ar]) + as.numeric(B%*%f[i-ma])
      mu[i] <- linkinv(f[i])
      y[i] <- rr(1,mu[i])
      st[i] <- 1/4*(pi*y[i]^2/(2*mu[i]^2)-2)
    }
  }
  return(ts(y[-(1:m)],frequency = 12))
}

#y<-simu.GAS(100,.1,.5,.1)
#plot(y)
#plot(ts(y[-(1:m)],frequency = 12))
#par(new=T)
#plot(z,col="blue")