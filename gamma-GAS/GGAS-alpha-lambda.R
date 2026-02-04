# Reference: gamma-GAS
# Created by Fernando J M de Araújo (nandomonteiro418@gmail.com), november/2024

# density function
dgammaGAS<-function(y,alpha,lambda)
{
  d<-(1/(gamma(alpha)*(alpha^(-1)*lambda)^alpha))*y^(alpha-1)*exp(-y*(alpha/lambda))
  d
}

# # density function
# alpha=1
# lambda=2
# dgammaGAS<-function(y){
#   d<-(1/(gamma(alpha)*(alpha^(-1)*lambda)^alpha))*y^(alpha-1)*exp(-y*(alpha/lambda))
#   d
# }
# 
# integrate(dgammaGAS,0,2)

# cumulative distribution function
pgammaGAS<-function(y,alpha,lambda)
{
  p<- pgamma(y,alpha,alpha/lambda)
  p
}

# pgammaGAS(2,alpha=1,
#           lambda=2)

# quantile function
qgammaGAS<-function(u,alpha,lambda)
{
  q<- qgamma(u,alpha,alpha/lambda)
  q 
}

# qgammaGAS(0.6321206,alpha=1,
#                     lambda=2)

# inversion method for randon generation
rgammaGAS<-function(n,alpha,lambda)
{
  u<- runif(n)
  y<- rgamma(u,alpha,alpha/lambda)
  y
}
