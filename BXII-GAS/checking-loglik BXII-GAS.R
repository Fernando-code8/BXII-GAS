mu=.7
c=2
tau=.5
y<-.5

s=1
loglik<-expression(log((log(1/(1-tau))*c)/(s^(c)*log(1+(mu/s)^c)))+(c-1)*(log(y))
                   +(log(1-tau)/log(1+(mu/s)^c)-1)*((log(1+(y/s)^c))))

eval(D(loglik,"mu"))

ut=(1+mu^c)
ht= log(1-tau)/log(ut)

((-c*mu^(c-1))/(ut*log(ut)))*(1+(ht*log(1+y^c)))
