###########################################
# calculating score vector for the Kw-GAS #
###########################################
 
## function for st with fixed qt
st.funcRay <- function(mu0,y1,link="log"){
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
  
  1/4*(pi*y1^2/(2*mu0^2)-2)
  
}


### CONTINUAR AQUI

# ### functions for the analytic derivative
# nablat<- function(mu0,y,c,tau,link="log"){
#   linktemp <- substitute(link)
#   if (!is.character(linktemp))
#   {
#     linktemp <- deparse(linktemp)
#     if (linktemp == "link")
#       linktemp <- eval(link)
#   }
#   if (linktemp == "log"){stats<-VGAM::loglink
#   } else if (linktemp == "sqrt"){stats<-VGAM::sqrtlink
#   } else {
#     stop(paste(linktemp, "link not available, available links are \"log\" and \"sqrt\""))
#   }
#   ###########################################################
#   link = linktemp 
#   diflink<-function(x)stats(x,inverse = F,deriv = 1)
#   ###########################################################
#   delta<-(log(1-tau)/log(1-qt0^varphi))
#   ct<-(qt0^(varphi-1))/
#     ((1-qt0^varphi)*log(1-qt0^varphi))*
#     (delta*log(1-y^varphi)+1)
#   ###########################################################
#   result<-varphi*ct/diflink(qt0)
#   return(result)
# } 
# 
# st_q<-function(qt0,y,varphi,tau,link=link){
#   linktemp <- substitute(link)
#   if (!is.character(linktemp))
#   {
#     linktemp <- deparse(linktemp)
#     if (linktemp == "link")
#       linktemp <- eval(link)
#   }
#   if (linktemp == "logit"){stats<-VGAM::logitlink
#   } else if (linktemp == "probit"){stats<-VGAM::probitlink
#   } else if (linktemp == "cloglog"){stats<-VGAM::clogloglink
#   } else {
#     stop(paste(linktemp, "link not available, available links are \"logit\", ","\"probit\" and \"cloglog\""))
#   }
#   ###########################################################
#   link = linktemp 
#   diflink<-function(x)stats(x,inverse = F,deriv = 1)
#   diflink2<-function(x)stats(x,inverse = F,deriv = 2)
#   ###########################################################
#   delta<-(log(1-tau)/log(1-qt0^varphi))
#   ct<-(qt0^(varphi-1))/
#     ((1-qt0^varphi)*log(1-qt0^varphi))*
#     (delta*log(1-y^varphi)+1)
#   ###########################################################
#   result<-  (1+delta*log(1-y^varphi))^2*(
#     (1-varphi)/
#       (varphi*ct*qt0)-log(1-qt0^varphi)/
#       (1+delta*log(1-y^varphi))
#     +diflink2(qt0)/(varphi*ct*diflink(qt0))
#   )-1
#   return(result)
# }
# 
# dll.dvarphi<- function(qt0,y1,varphi,tau){
#   delta<-(log(1-tau)/log(1-qt0^varphi))
#   ct<-(qt0^(varphi-1))/
#     ((1-qt0^varphi)*log(1-qt0^varphi))*
#     (delta*log(1-y1^varphi)+1)
#   ###########################################################
#   result<- (1/varphi+log(y1)+ct*qt0*log(qt0)-
#               (delta-1)*y1^varphi*log(y1)/(1-y1^varphi))
#   return(result)
# }
# 
# dst.varphi <- function(qt0,y,varphi,tau){
#   delta<-(log(1-tau)/log(1-qt0^varphi))
#   ct<-(qt0^(varphi-1))/
#     ((1-qt0^varphi)*log(1-qt0^varphi))*
#     (delta*log(1-y^varphi)+1)
#   p1<- delta*qt0/varphi*(
#     log(qt0)*log(1-y^varphi)-
#       y^varphi*
#       log(y)*(1-qt0^varphi)*
#       log(1-qt0^varphi)/(qt0^(varphi)*(1-y^varphi))
#   )
#   ###########################################################
#   result <- (p1-ct*(1-qt0^varphi)*log(1-qt0^varphi)/
#                (varphi*qt0^(2*varphi-2))*(
#                  log(qt0)*(qt0^varphi+log(1-qt0^varphi))+
#                    1/varphi*(1-qt0^varphi)*log(1-qt0^varphi)
#                ))
#   ###########################################################
#   if(qt0==0 || y==0) return(0) else
#     return(result)
# }
# 
# KwGAS.score <- function(w,A,B,beta=0,varphi,y,
#                         tau=0.5, ar=NA, ma=NA, X=NA, link = "logit")
# {
#   ##### link function definitions
#   linktemp <- substitute(link)
#   if (!is.character(linktemp))
#   {
#     linktemp <- deparse(linktemp)
#     if (linktemp == "link")
#       linktemp <- eval(link)
#   }
#   if (linktemp == "logit"){stats<-VGAM::logitlink
#   } else if (linktemp == "probit"){stats<-VGAM::probitlink
#   } else if (linktemp == "cloglog"){stats<-VGAM::clogloglink
#   } else {
#     stop(paste(linktemp, "link not available, available links are \"logit\", ","\"probit\" and \"cloglog\""))
#   }
#   ###########################################################
#   link = linktemp 
#   linkinv = function(x)stats(x,inverse = T)
#   diflink<-function(x)stats(x,inverse = F,deriv = 1)
#   diflink2<-function(x)stats(x,inverse = F,deriv = 2)
#   ###########################################################  
#   ##### X definitions
#   if(any(is.na(X))==FALSE){
#     X<-as.matrix(X)
#     if(any(beta==0)) stop("You need to inform beta")
#     if(length(beta)!=dim(X)[2]) stop("The length of beta must be equal to the number of columns of X")
#     # k = ncol(X)
#   }else{
#     if(beta!=0) stop("You need to inform X")
#     X <- matrix(0, c(n,1))
#   }
#   ##### defining the lags for st
#   if(any(is.na(ar)==T)){
#     if (any(is.na(A)==F))
#     {
#       stop("You must inform the order of A")
#     } else {A <- ar <- 0}
#   } else {
#     if (any(is.na(A)==T)) stop("You must inform A")
#     if (!isTRUE(all(ar == floor(ar)))) 
#       stop("'ar' must only contain integer values")
#     if(length(A)!=length(ar))
#       stop("The lengths of A and ar must be equal")
#   }
#   ##### defining the lags for ft
#   if(any(is.na(ma)==T)){
#     if(any(is.na(B) ==F ))
#     {
#       stop("You must inform the order of B")
#     } else {B <- ma <- 0}
#   } else {
#     if (any(is.na(A)==T)) stop("You must inform B")
#     if (!isTRUE(all(ma == floor(ma)))) 
#       stop("'ma' must only contain integer values")
#     if(length(B)!=length(ma))
#       stop("The lengths of B and ma must be equal")
#   }
#   ###########################################################
#   p <- max(ar)
#   q <- max(ma)
#   n <- length(y)
#   m <- max(p,q,na.rm=T)
#   ## initializations
#   st<- qt<- f<- ll <- dst.dqt <-nabla<-
#     dst<- df.dw <- df.dB <- 
#     dvarphi <- df.dvarphi <- rep(0,n)
#   
#   df.dA <- matrix(0,n,length(A))
#   df.dB <- matrix(0,n,length(B))
#   df.dbeta <- matrix(0,n,length(beta))
#   
#   for(i in (m+1):n)
#   {
#     f[i]  <- w + as.numeric(A%*%st[i-ar]) + 
#       as.numeric(B%*%f[i-ma]) + X[i,]%*%beta
#     qt[i]   <- linkinv(f[i])
#     st[i] <- st.func(qt[i],y[i],varphi,tau,link=link)
#     dst[i] <-st_q(qt[i],y[i],varphi,tau,link=link)
#     if(eval(qt[i])==0) dfi <- 0 else dfi <- diflink(eval(qt[i]))
#     df.dw[i] <- 1+(A*dst[i-ar])%*%df.dw[i-ar]+B%*%df.dw[i-ma]
#     df.dA[i,] <- st[i-ar]+t(A*dst[i-ar])%*%df.dA[i-ar,]+B%*%df.dA[i-ma,]
#     df.dB[i,]<- f[i-ma]+(A*dst[i-ar])%*%df.dB[i-ar,]+B%*%df.dB[i-ma,]
#     df.dbeta[i,] <- X[i,]+(A*dst[i-ar])%*%df.dbeta[i-ar,]+B%*%df.dbeta[i-ma,]
#     df.dvarphi[i]<- A%*%dvarphi[i-ar]+B%*%df.dvarphi[i-ma]
#     dvarphi[i] <- dst.varphi(qt[i],y[i],varphi,tau=tau)%*%dfi+dst[i]*df.dvarphi[i]
#   }
#   
#   Uw <- sum(nablat(qt[-c(1:m)],y[-c(1:m)],varphi,tau)*df.dw[-c(1:m)])
#   UA <- apply(
#     as.matrix(nablat(qt[-c(1:m)],y[-c(1:m)],varphi,tau)*df.dA[-c(1:m),]),
#     2,sum)
#   UB <- apply(
#     as.matrix(nablat(qt[-c(1:m)],y[-c(1:m)],varphi,tau)*df.dB[-c(1:m),]),
#     2,sum)
#   Ubeta <- apply(
#     as.matrix(nablat(qt[-c(1:m)],y[-c(1:m)],varphi,tau)*df.dbeta[-c(1:m),]),
#     2,sum)
#   Uvarphi <- sum(dll.dvarphi(qt[-c(1:m)],y[-c(1:m)],varphi,tau)+
#                    nablat(qt[-c(1:m)],y[-c(1:m)],varphi,tau)*df.dvarphi[-c(1:m)])
#   
#   if(any(beta==0)) {rval <- c(Uw,UA,UB,Uvarphi)
#   names(rval) <- c("Uw","UA","UB","Uvarphi")
#   }
#   else {rval <- rval <- c(Uw,UA,UB,Ubeta,Uvarphi)
#   names(rval) <- c("Uw",c(paste0("UA", 1:length(A))),
#                    c(paste0("UB", 1:length(B))),
#                    c(paste0("Ubeta", 1:length(beta))),
#                    "Uvarphi")
#   }
#   
#   
#   return(rval)
# }
