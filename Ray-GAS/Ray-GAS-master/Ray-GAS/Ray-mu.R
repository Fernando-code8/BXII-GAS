# Reference: BXII-GAS
# Created by Fernando J M de Araújo (nandomonteiro418@gmail.com), March/2024

# density function
dr<-function(y,mu)
{
  d<- pi*y/(2*mu^2)*exp(-(pi*y^2)/(4*mu^2))
  d
}

# cumulative distribution function
pr<-function(y,mu)
{
  p<- 1- exp((-pi*y^2)/(4*mu^2))
  p
}

# quantile function
qr<-function(u,mu)
{
  q<- 2*mu*sqrt((-log(1-u))/pi)
  q
}

# inversion method for randon generation
rr<-function(n,mu)
{
  u<- runif(n)
  y<- 2*mu*sqrt(-log(1-u)/pi)
  y
}