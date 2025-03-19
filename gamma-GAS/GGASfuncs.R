# Reference: gamma-GAS
# Created by Fernando J M de Araújo (nandomonteiro418@gmail.com), november/2024

##############################################
# calculating score vector for the gamma-GAS #
##############################################
 
## function for st with fixed alpha
st.funcGGAS <- function(alpha0,y1,lambda,tau,link="log"){
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (linktemp == "log"){stats<-VGAM::loglink
  } else if (linktemp == "sqrt"){stats<-VGAM::sqrtlink
  } else {
    stop(paste(linktemp, "link not available, available links are \"log\" and \"sqrt\""))
  }
  ###########################################################
  link = linktemp 
  diflink<-function(x)stats(x,inverse = F,deriv = 1)
  ########################
  d=1
  St<-alpha0^d*((y1/lambda)-1)
  St*diflink(alpha0)
  # St<-alpha0*(log(y1)-log(lambda)+psigamma(alpha0, deriv = 1))*(alpha0^2)*psigamma(alpha0, deriv = 2)
  # St
  
}
