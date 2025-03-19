################################################################
#
# Rayleigh autorregressive moveing average model (RARMA)
#
# Implemented by Fabio M Bayer (bayer@ufsm.br) 
#
# February, 2019
#
################################################################

rr<-function(mu) # metodo da inversao
{
  n<-length(mu)
  
  u<- runif(n)
  y<- 2*mu*sqrt(-log(1-u)/pi) # metodo da inversao
  
  y
}

qr<-function(alpha,mu=1) # quantile function
{
  q<- 2*mu*sqrt((-log(1-alpha))/pi)
  q
}

dr<-function(x,mu=1) # density function
{
  d<- pi*x/(2*mu^2)*exp(-(pi*x^2)/(4*mu^2))
  d
}

pr<-function(x,mu=1) # cumulative function
{
  p<- 1- exp((-pi*x^2)/(4*mu^2))
  p
}

rarma<- function (y, ar=NA, ma=NA, link = "log",diag=1,h=6,X=NA,X_hat=NA,resid=1,pw=0,PI=NA)
{  
  
  source("rarma.fit.r")
  
  if (min(y) < 0)
    stop("OUT OF RANGE R+")
  
  # if(is.ts(y)==T)
  # {
  #   freq<-frequency(y)
  # }else stop("data can be a time-series object")
  
  
  if(any(is.na(ar))==F) names_phi<-c(paste("phi",ar,sep=""))
  
  if(any(is.na(ma))==F) names_theta<-c(paste("theta",ma,sep=""))
  
  if(any(is.na(X))==F)
  {
    names_beta<-c(paste("beta",1:ncol( as.matrix(X) ),sep=""))
  }
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  
  m <- max(p,q,na.rm=T)
  
  p1 <- length(ar)
  q1 <- length(ma)
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("log", "sqrt", "identity")))
    stats <- make.link(linktemp)
  else stop(paste(linktemp, "link not available, available links are \"log\", ",
                  "\"sqrt\" and \"identity\""))
  
  link1 <- structure(list(link = linktemp, 
                          linkfun = stats$linkfun,
                          linkinv = stats$linkinv, 
                          mu.eta = stats$mu.eta
                          # diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t))
                          #                          )
  )
  )
  
  fit1 <- rarma.fit(y, ar, ma, link1, names_phi, names_theta, names_beta, diag, h, X, X_hat,resid=resid,pw,PI) # model estimation
  
  return(fit1)
}

