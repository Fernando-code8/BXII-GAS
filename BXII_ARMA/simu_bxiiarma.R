################################################################################
# PAPER: Quantile-based dynamic modeling of asymmetric data: a novel Burr XII
#        approach for positive continuous random variables
# SECTION: 5. Monte Carlo simulation
# GOAL: Function to generate random samples of the BXII-ARMA model.
# AUTHORS: Fernando Jose Monteiro de Araujo, Renata Rojas Guerra and 
#          Fernando Arturo Pena-Ramirez
# LAST UPDATE: May 25, 2024
################################################################################

simu.bxiiarma <- function(n,phi=NA,theta=NA, alpha=1,c=1,tau=0.5,freq=12,
                          link="log")
{
  source("bxii-funcs.r")
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }else{
    ar<-0
    phi<-0
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }else{
    ma<-0
    theta<-0
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("log", "sqrt")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\" 
               and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv ))
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(NA,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) 
                 + as.numeric(theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rbxii(1,mu[i],c)
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]   
      
    }
    
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } 
}
