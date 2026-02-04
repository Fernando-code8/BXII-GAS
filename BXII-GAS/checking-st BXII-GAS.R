#######################################
### calculating st for the BXII-GAS ###
#######################################
rm(list=ls())
gc()
# initial values
## static parameters
mu0=.5; c=2; tau=.25; y=3
w=-1; A= 0.1; B=0.7
## logit link
link<-VGAM::loglink
linkinv = function(x)link(x,inverse = T)
diflink<-function(x)link(x,inverse = F,deriv = 1)
diflink2<-function(x)link(x,inverse = F,deriv = 2)

ut=(1+mu0^c)
ht= log(1-tau)/log(ut)
vh <- ((-c*mu0^(c-1))/(ut*log(ut)))*(1+(ht*log(1+y^c))) ## loglik/dmu
at<-(c^2*mu0^(2*(c-1)))/((ut^2)*(log(ut)^2)) ## loglik/dmu2
########################
st0 <- (1/(at))*vh
st0

## first iteration
ft1  <- w + A*st0 #
mu1 <- linkinv(ft1) # exp(w + A*st0)

# loglik of the BXII without considering f 
s=1
loglik<-expression(log((log(1/(1-tau))*c)/(s^(c)*log(1+(mu1/s)^c)))+(c-1)*(log(y))
                   +(log(1-tau)/log(1+(mu1/s)^c)-1)*((log(1+(y/s)^c))))

# checking the st without considering f 
## checking nabla_t
ut1=(1+mu1^c)
ht1= log(1-tau)/log(ut1)
nablat<-((-c*mu1^(c-1))/(ut1*log(ut1)))*(1+(ht1*log(1+y^c)))

c(eval(D(loglik,"mu1")),nablat)


## checking S_t
ut1=(1+mu1^c)
at1<-(c^2*mu1^(2*(c-1)))/((ut1^2)*(log(ut1)^2)) ## loglik/dmu2
St<-1/at1


ht1= log(1-tau)/log(ut1)
vh1 <- ((-c*mu1^(c-1))/(ut1*log(ut1)))*(1+(ht1*log(1+y^c))) ## loglik/dmu
at1<-(c^2*mu1^(2*(c-1)))/((ut1^2)*(log(ut1)^2)) ## loglik/dmu2
########################
st1 <- (1/at1)*vh1
st1


St.p1<-function(y)(
  log(1/(1-tau))*(c*y^(c-1))/(s^(c)*log(1+(mu1/s)^c))*
    (1+(y/s)^c)^(log(1-tau)/(log(1+(mu1/s)^c))-1)#densidade
  *(((-c*mu1^(c-1))/((1+mu1^c)*log((1+mu1^c))))*(1+((log(1-tau)/log((1+mu1^c)))*log(1+y^c))))^2 #nablat^2
)

1/integrate(St.p1,0,Inf)$value
St

c(nablat*St,st1,nablat/integrate(St.p1,0,Inf)$value)

#implementação
ut=(1+mu1^c)
ht= log(1-tau)/log(ut)
at<-(c^2*mu1^(2*(c-1)))/((ut^2)*(log(ut)^2)) ## loglik/dmu2

nablat<-((-c*mu1^(c-1))/(ut*log(ut)))*(1+(ht*log(1+y^c)))
St<-1/at
nablat*St
#Paper
stpaper=-((ut*log(ut))/(c*mu1^(c-1)))*(1+ht*log(1+y^c))


##########################################

# log-likelihood with ft1
s=1
loglik_ft1<-expression(log((log(1/(1-tau))*c)/(s^(c)*log(1+(exp(ft1)/s)^c)))+(c-1)*(log(y))
           +(log(1-tau)/log(1+(exp(ft1)/s)^c)-1)*((log(1+(y/s)^c))))

# adding dmut.df in nabla_t
eval(D(loglik_ft1,"ft1"))
nablat/diflink(mu1)

###
ut1=(1+exp(ft1)^c)
ht1= log(1-tau)/log(ut1)
# adding dmut.df in St
St.p2<-function(y)(
  log(1/(1-tau))*(c*y^(c-1))/(s^(c)*log(1+(exp(ft1)/s)^c))*
    (1+(y/s)^c)^(log(1-tau)/(log(1+(exp(ft1)/s)^c))-1)*#densidade avaliada em m1= exp(ft1)
    eval(D(expression(
      (((-c*exp(ft1)^(c-1))/((1+exp(ft1)^c)*log((1+exp(ft1)^c))))*(1+((log(1-tau)/log((1+exp(ft1)^c)))*log(1+y^c)))) # nablat
      *(exp(ft1))#derivada da inversa da função de ligação
    ),"ft1"
    )
    )
)

-1/integrate(St.p2,0,Inf)$value
St*diflink(mu1)^2

# checking st of ft1
(nablat/diflink(mu1))*(diflink(mu1)^2*St)
st1*diflink(mu1)

#Paper
stpaper=-((ut*log(ut))/(c*mu1^(c-1)))*(1+ht*log(1+y^c))

## calculating the derivative with respect to varphi
dll.dc=1/c+log(y)-(((1-ht)*y^c*log(y))/(1+y^c))-((mu1^c*log(mu1))/(ut*log(ut)))*(1+ht*log(1+y^c))
eval(D(loglik,"c"))
dll.dc


