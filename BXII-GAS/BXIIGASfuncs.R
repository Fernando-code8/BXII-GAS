# Reference: BXII-GAS
# Created by Fernando J M de Araújo (nandomonteiro418@gmail.com), June/2024

#############################################
# calculating score vector for the BXII-GAS #
#############################################
 
## function for st with fixed mu
st.funcBXII <- function(mu0,y1,c,tau,link="log"){
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
  ut=(1+mu0^c)
  ht= log(1-tau)/log(ut)
  at<-(c^2*mu0^(2*(c-1)))/((ut^2)*(log(ut)^2)) ## loglik/dmu2
  
  nablat<-((-c*mu0^(c-1))/(ut*log(ut)))*(1+(ht*log(1+y1^c)))
  St<-1/at
  nablat*St*diflink(mu0)
  
}