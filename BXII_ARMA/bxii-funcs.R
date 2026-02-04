################################################################################
# PAPER: Quantile-based dynamic modeling of asymmetric data: a novel Burr XII
#        approach for positive continuous random variables
# SECTION: 6. Applications
# GOAL: It provides some useful functions for the BXII-ARMA model, 
#       such as: cumulative distribution function (fda), probability density 
#       function (pdf), quantile function (qf) and function to generate 
#       occurrences by inversion method.
# AUTHORS: Fernando Jose Monteiro de Araujo, Renata Rojas Guerra and 
#          Fernando Arturo Pena-Ramirez
# LAST UPDATE: May 25, 2024
################################################################################

################################################################################
#####################  BXII-ARMA model - USEFUL FUNCTIONS  #####################
################################################################################

# density function
dbxii<-function(y,mu,c,tau=.5)
{
  s=1
  d<-log(1/(1-tau))*(c*y^(c-1))/(s^(c)*log(1+(mu/s)^c))*
    (1+(y/s)^c)^(log(1-tau)/(log(1+(mu/s)^c))-1)
  d
}

# cumulative distribution function
pbxii<-function(y,mu,c,tau=.5)
{
  s=1
  p<- 1-(1+(y/s)^c)^(log(1-tau)/(log(1+(mu/s)^c)))
  p
}

# quantile function
qbxii<-function(u,mu,c,tau=.5)
{
  s=1
  q<- s*((1-u)^(log(1+(mu/s)^c)/log(1-tau))-1)^(1/c)
  q
}

# inversion method for randon generation
rbxii<-function(n,mu,c,tau=.5)
{
  s=1
  u<- runif(n)
  y<- s*((1-u)^(log(1+(mu/s)^c)/log(1-tau))-1)^(1/c)
  y
}

