# Reference: gamma-GAS
# Created by Fernando J M de Araújo (nandomonteiro418@gmail.com), november/2024

# density function
dgammaGAS<-function(y,alpha,lambda)
{
  d<-(1/(gamma(alpha)*(alpha^(-1)*lambda)^alpha))*y^(alpha-1)*exp(-y*(alpha/lambda))
  d
}

# cumulative distribution function
pgammaGAS<-function(y,alpha,lambda)
{
  p<- pgamma(y,alpha,alpha/lambda)
  p
}

# quantile function
qgammaGAS<-function(u,alpha,lambda)
{
  q<- qgamma(u,alpha,alpha/lambda)
  q 
}

# inversion method for randon generation
rgammaGAS<-function(n,alpha,lambda)
{
  u<- runif(n)
  y<- rgamma(u,alpha,alpha/lambda)
  y
}
