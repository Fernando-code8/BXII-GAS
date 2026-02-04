################################################################################
# PAPER: Quantile-based dynamic modeling of asymmetric data: a novel Burr XII
#        approach for positive continuous random variables
# SECTION: 6. Applications
# GOAL: Provides functions to fit the BXII-ARMA model
# AUTHORS: Fernando Jose Monteiro de Araujo, Renata Rojas Guerra and 
#          Fernando Arturo Pena-Ramirez
# LAST UPDATE: May 25, 2024
################################################################################

################################################################################
#######                 BXII-ARMA MODEL implementation                 #########
################################################################################

bxiiarma.fit<-function (y, ar = NA, ma = NA, X = NA,h1=6, tau = .5,diag1=0,
                        resid=3,mod=3,link = "log")
{
  if (min(y) < 0)
    stop("OUT OF SUPPORT!")
  
  if(is.ts(y)==T)  freq<-frequency(y) 
  else stop("data can be a time-series object")
  
  z<-c()
  maxit1<-10000
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  y1 <- y[(m+1):n]
  p1 <- length(ar)
  q1 <- length(ma)
  y_prev <- c(rep(NA,(n+h1)))
  error <- rep(0,n)
  eta <- rep(NA,n)
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("log", "sqrt"))){
    stats <- make.link(linktemp)
  }  else {
      stop(paste(linktemp, "link not available, available links are \"log\" and 
                 \"sqrt\""))
    }
  
  link = linktemp 
  linkfun = stats$linkfun
  linkinv = stats$linkinv 
  mu.eta = stats$mu.eta 
  diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  ynew = linkfun(y) 
  ynew_ar <- suppressWarnings(matrix(ynew,(n-1),max(p,1,na.rm=T)))
  
###########################################################3  
  if(any(is.na(ar)) == F) {
    names_phi <- c(paste("phi", ar, sep = ""))
    Z <- suppressWarnings(matrix(ynew, (n-1), p1)[m:(n-1),])} else {
      ar = p1<-0; Z <- NA } 

  if(any(is.na(ma)) == F) {
    names_theta <- c(paste("theta", ma, sep = ""))
  } else ma = q1 <- 0 
  
  if(class(X)== "numeric"){
    if(any(is.na(X)) == F){
    names_beta<-c(paste("beta", 1 : ncol(as.matrix(X)), sep = ""))
    Xm <- X[(m+1):n]
    k = 1
    c<-1 
  } else {
    k = 0 
    X <- matrix(rep(0,n), nrow = n)
    Xm <- NA
    c<-1 
  }}else{
  if(any(is.na(X)) == F){
    names_beta<-c(paste("beta", 1 : ncol(as.matrix(X)), sep = ""))
    Xm <- X[(m+1):n,]
    k = ncol(X)
    c<-1
  } else {
    k = 0 
    X <- matrix(rep(0,n), nrow = n)
    Xm <- NA
    c<-1
  }}
  
  ###FB  recorrences   #
  q_1 <- max(q1, 1)
  R <- matrix(rep(NA, (n-m)*q_1), ncol = q_1)
  k_i <- q1/q_1 
  deta.dalpha <- rep(0, n)
  deta.dbeta <- matrix(0, ncol=max(k,1), nrow=n)
  deta.dphi <- matrix(0, ncol=p1, nrow=n)
  deta.dtheta <- matrix(0, ncol=q_1, nrow=n)
  
  Xstart <- (cbind(rep(1, (n-m)), Xm, Z))
  Xstart <- matrix(apply(Xstart, 1, na.omit),nrow = (n-m),byrow = T)
  ols <- lm.fit(Xstart, ynew[(m+1) : n])$coef
  initial <- rep(0, k+p1+q1+1)
  initial[1 : (k+p1+1)] <- ols
  initial<-c(initial,c)
  
  loglik <- function(z) 
  {
    alpha <- z[1]
    if(k==0) beta = as.matrix(0) else beta = as.matrix(z[2:(k+1)])
    if(p1==0) {phi = as.matrix(0);ar=1} else phi = as.matrix(z[(k+2):(k+p1+1)]) 
    if(q1==0) theta = as.matrix(0) else theta = as.matrix(z[(k+p1+2):(k+p1+q1+1)])
    c <- z[p1+q1+2]
    Xbeta <- X%*%beta
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    
    for(i in (m+1):n)
    {
      eta[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1), ar] 
                - Xbeta_ar[(i-1), ar])%*%phi + t(theta)%*%error[i-ma]
      error[i] <- ynew[i] - eta[i] 
    }
    mu <- linkinv(eta[(m+1):n])
    s=1
    ll <- log((log(1/(1-tau))*c)/(s^(c)*log(1+(mu/s)^c)))+(c-1)*(log(y1))
          +(log(1-tau)/log(1+(mu/s)^c)-1)*((log(1+(y1/s)^c)))
    sum(ll)
  } 
  # __________________________________________________________________________ #
  escore.bxiiarma <- function(z)
    {
    alpha <- z[1]
    if(k==0)  beta = as.matrix(0) else beta = as.matrix(z[2:(k+1)])
    if(p1==0) {phi = as.matrix(0);ar=1} else phi = as.matrix(z[(k+2):(k+p1+1)]) 
    if(q1==0) theta = as.matrix(0) else  theta = as.matrix(z[(k+p1+2):(k+p1+q1+1)])
    c <- z[p1+q1+2]
    Xbeta <- X%*%beta
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    for(i in (m+1):n)
      {
      eta[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1),ar] 
              - Xbeta_ar[(i-1),ar])%*%phi + t(theta)%*%error[i-ma]
      error[i] <- ynew[i] - eta[i] 
    }
    
    mu <- linkinv(eta[(m+1):n])
    Xbeta <- X%*%beta
    for(i in 1:(n-m)){
    R[i,] <- error[i+m-ma]*k_i}

    for(i in (m+1):n)
    {
      deta.dalpha[i] <- 1 - deta.dalpha[i-ma]%*%theta
      deta.dbeta[i,] <- X[i,] - t(phi)%*%X[i-ar,] - t(theta)%*%deta.dbeta[i-ma,]
      deta.dphi[i,] <- ynew_ar[i-ar]- Xbeta[i-ar] - t(theta)%*%deta.dphi[i-ma,]
      deta.dtheta[i,] <- R[(i-m),] - t(theta)%*%deta.dtheta[i-ma,]
    }
    
    v <- deta.dalpha[(m+1):n]
    rM <- deta.dbeta[(m+1):n,]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    
    mT <- diag(mu.eta(eta[(m+1):n]))
    
    vh <- ((-c*mu^(c-1))/((1+mu^c)*log(1+mu^c)))*(1+((log(1-tau)*(log(1+y1^c)))/log(1+mu^c))) #loglik/dmu
    
    Ualpha <- t(v) %*% mT %*% vh
    Ubeta <- t(rM) %*% mT %*% vh
    Uphi <-   t(rP) %*% mT %*% vh
    Utheta <- t(rR) %*% mT %*% vh
  
    # derivada em relacao a c
    a <- as.vector(1/c+(log(y1))-((1-log(1-tau)/log(1+mu^c))*((((y1^c)*log(y1))/(1+y1^c))))-(((mu^c)*log(mu))/((1+mu^c)*log(1+mu^c)))*(1+((log(1-tau)*(log(1+y1^c)))/log(1+mu^c))))
    
    Uc <- sum(a)
    
    rval <- c(Ualpha,Ubeta,Uphi,Utheta,Uc)
    return(rval[rval!=0])
    }
  
  opt<-optim(initial, loglik, escore.bxiiarma, method = "BFGS", hessian= TRUE,
             control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
  
  z$conv <- opt$conv
  coef <- (opt$par)#[1:(1+p1+q1+k)]
  alpha <- coef[1]
  if(k==0) beta=names_beta=NULL else z$beta <- coef[2:(k+1)]
  if(p1==0) phi=names_phi=NULL else z$phi <- coef[(k+2):(k+p1+1)]
  if(q1==0) theta=names_theta=NULL else z$theta <- coef[(k+p1+2):(k+p1+q1+1)]
  c <- coef[p1+q1+2]
  
  #############################################################################
  errorhat<-rep(0,n) # E(error)=0
  etahat<-rep(NA,n)

  if( mod==2
    # any(is.na(ar)==T) && any(is.na(ma)==F)
    ){
  # Modelo MA
  for(i in (m+1):n)
  {
    etahat[i]<-alpha + (z$theta%*%errorhat[i-ma])
    errorhat[i]<- ynew[i]-etahat[i] # predictor scale
  }
  }

  if(mod==3
    # any(is.na(ar)==F) && any(is.na(ma)==F)
    ){
  # Modello ARMA
  for(i in (m+1):n)
  {
    etahat[i]<-alpha + (z$phi%*%ynew[i-ar]) + (z$theta%*%errorhat[i-ma])
    errorhat[i]<- ynew[i]-etahat[i] # predictor scale
  }
  }

  if(mod==1
    # any(is.na(ar)==F) && any(is.na(ma)==F)
    ){
  # Modelo AR
  for(i in (m+1):n)
  {
    etahat[i]<-alpha + (z$phi%*%ynew[i-ar])
    errorhat[i]<-ynew[i]-etahat[i] # predictor scale
  }
  }


  muhat <- linkinv(etahat[(m+1):n])
  y1<-y[(m+1):n]

  z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
  z$etahat <- etahat
  z$errorhat <- errorhat
  z$mustarhat <- digamma(muhat * c) - digamma((1 - muhat) * c)

  ##############################################################################
  names_par <- c("alpha",names_beta,names_phi,names_theta,"c")
  names(coef)<-names_par
  z$coeff <- coef
  J_inv <- solve(-(opt$hessian))
  z$stderror<-sqrt(diag(J_inv))
  z$zstat <- abs(z$coef/z$stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat))
  z$loglik <- opt$value
  z$counts <- as.numeric(opt$counts[1])
  z$c<-c
  z$stats <- stats
  
  
  
  if(any(is.na(X)==F))
  {
    z$k<- (p1+q1+2+k)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+2)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }
  
  # #### Forecasting
  ynew_prev <- c(ynew,rep(NA,h1))
  y_prev[1:n] <- z$fitted
  

  # X_prev<- rbind(X,X_hat) #x_hat sao as covariaveis
  if(mod==2
    # any(is.na(ar)==T) && any(is.na(ma)==F)
    ){
  # Modelo MA
  for(i in 1:h1)
  {
    ynew_prev[n+i] <- alpha +(z$theta%*%errorhat[n+i-ma])
    y_prev[n+i] <- linkinv(ynew_prev[n+i])
    errorhat[n+i] <- 0
  }
  }

  if(mod==3
    # any(is.na(ar)==F) && any(is.na(ma)==F)
    ){
  # Modelo ARMA
  for(i in 1:h1)
  {
    ynew_prev[n+i] <- alpha + (z$phi%*%ynew_prev[n+i-ar]) 
    + (z$theta%*%errorhat[n+i-ma])
    y_prev[n+i] <- linkinv(ynew_prev[n+i])
    errorhat[n+i] <- 0
  }
    }

  if(mod==1
    # any(is.na(ar)==F) && any(is.na(ma)==F)
    ){
  # Modelo AR
  for(i in 1:h1)
  {
    ynew_prev[n+i] <- alpha + (z$phi%*%ynew_prev[n+i-ar])
    y_prev[n+i] <- linkinv(ynew_prev[n+i])
  }
  }

  z$serie <- y
  # z$barma <- names_par
  z$forecast <- y_prev[(n+1):(n+h1)]
  setwd("~/GitHub/BXII-GAS/BXII_ARMA")
  source('bxii-funcs.r')

  residc <- z$resid3
  # residuals
  res1 <- y-z$fitted
  vary <- ((log(0.5)/log(1-z$fitted^c)) * beta(1+2/c, log(0.5)/log(1-z$fitted^c))
                      - (log(0.5)/log(1-z$fitted^c) 
                         * beta(1+1/c, log(0.5)/log(1-z$fitted^c)))^2 )

  z$resid1 <- (res1/sqrt(vary))[(m+1):n]

  l_tilde <- log(dbxii(y,y,z$c))
  l_hat <- log(dbxii(y,z$fitted,z$c))

  dt <- (l_tilde-l_hat)[(m+1):n]
  dt[which(dt<0)]<-0

  z$l_hat <- l_hat

  z$resid2 <- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt))

  z$resid3 <- as.vector(qnorm(pbxii(y[(m+1):n],z$fitted[(m+1):n],z$c)))

  if(resid==1) residc <- z$resid1
  if(resid==2) residc <- z$resid2
  if(resid==3) residc <- z$resid3

  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),
                              round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  z$model <- model_presentation

  # ###################################################
  # ################## GRAPHICS #######################
  # ###################################################
  if(diag1==0)
  {

    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    print(c("AIC:",round(z$aic,4)," SIC:",round(z$bic,4)," HQ:",round(z$hq,4)),
          quote=F)

    print("Residuals:",quote=F)
    print(summary(residc))

    t<-seq(-5,n+6,by=1)


    max_y<- max(c(z$fitted,y),na.rm=T)
    min_y<- min(c(z$fitted,y),na.rm=T)

    densidade<-density(residc)

    max_r<- max(residc,na.rm=T)
    min_r<- min(residc,na.rm=T)

    fim<-end(y)[1]+end(y)[2]/12


    y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))

      # Initial graphics
      w1<-4
      h11<-4
        
      postscript(file = "res_indice.eps",width = w1, height = h11,family = "Times")
      par(mar=c(5,6,4,1)+.1)
      plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
      lines(t,rep(-3,n+12),lty=2,col=1)
      lines(t,rep(3,n+12),lty=2,col=1)
      lines(t,rep(-2,n+12),lty=3,col=1)
      lines(t,rep(2,n+12),lty=3,col=1)
      dev.off()

      postscript(file = "resid_density.eps",width = w1, height = h11,family = "Times")
      par(mar=c(5,6,4,1)+.1)
      plot(densidade,ylab="Density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
      lines(densidade$x,dnorm(densidade$x),lty=2)
      legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
           pt.bg="white", lty=c(1,2), bty="n",cex=0.35)
        dev.off()

        postscript(file = "res_ACF.eps",width = w1, height = h11,family = "Times")
        par(mar=c(5,6,4,1)+.1)
        acf(residc,ylab="ACF",xlab="Lag", main="")
        dev.off()
        
        postscript(file = "res_PACF.eps",width = w1, height = h11,family = "Times")
        par(mar=c(5,6,4,1)+.1)
        pacf(residc,ylab="PACF",xlab="Lag", main="")
        dev.off()

    }else{print="You have not activated the function to generate the graphs"}
  return(z)
}


