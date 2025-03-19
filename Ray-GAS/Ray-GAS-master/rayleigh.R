
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