#####################    Gráficos   #############################################

TC4<-paste("TC4",".pdf",sep="")
pdf(file = TC4,width = 4, height = 4,family = "Times")
{
par(mfrow=c(1,1))
par(mar=c(2.8, 2.5, 0.5, 0.5))
par(mgp=c(1.3, 0.45, 0))

# TC
omega<-TC[,1]
alpha<-TC[,2]
beta<-TC[,3]
# delta1<-TC[,4]
# delta2<-TC[,5]
psi<-TC[,4]

plot(omega,xaxs="r",type="l",
     ylab=expression("CR95%"),
     xlab=expression(plain(italic(n))),
     axes=FALSE,
     ylim=c(min(c(TC)),max(c(TC))))

lines(omega, lty=1,col=1)
lines(alpha, lty=1,col=2)
lines(beta, lty=1,col=3)
# lines(delta1, lty=1,col=4)
# lines(delta2, lty=1,col=5)
lines(psi, lty=1,col=6)

points(omega,pch=1,col=1)
points(alpha,pch=2,col=2)
points(beta,pch=3,col=3)
# points(delta1,pch=4,col=4)
# points(delta2,pch=5,col=5)
points(psi,pch=6,col=6)

legend("bottomright",c(expression(plain(hat(omega))),
# legend("topright",c(expression(plain(hat(omega))),
                    expression(plain(hat(alpha))), 
                    expression(plain(hat(beta))),
                    # expression(plain(hat(gamma)[1])),
                    # expression(plain(hat(gamma)[2])),
                    expression(plain(hat(c)))),
       lty=c(1,1,1,1,1,1
             ), pch=c(1,2,3,#4,
                      #5,
                      6),col=c(1,2,3,#4,
                               #5,
                               6),bty="n",cex=1.25)
axis(1,1:5,c(100,200,300,500,1000))
box();axis(2)
} 
dev.off()

# RB
RB4<-paste("RB4",".pdf",sep="")
pdf(file = RB4,width = 4, height = 4,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.5, 0.5, 0.5))
  par(mgp=c(1.3, 0.45, 0))
  
  # RB%
  omega<-RB[,1]
  alpha<-RB[,2]
  beta<-RB[,3]
  # delta1<-RB[,4]
  # delta2<-RB[,5]
  psi<-RB[,4]
  
  plot(omega,xaxs="r",type="l",
       ylab=expression("RB%"),
       xlab=expression(plain(italic(n))),
       axes=FALSE,
       ylim=c(min(c(RB)),max(c(RB))))
  
  lines(omega, lty=1,col=1)
  lines(alpha, lty=1,col=2)
  lines(beta, lty=1,col=3)
  # lines(delta1, lty=1,col=4)
  # lines(delta2, lty=1,col=5)
  lines(psi, lty=1,col=6)
  
  points(omega,pch=1,col=1)
  points(alpha,pch=2,col=2)
  points(beta,pch=3,col=3)
  # points(delta1,pch=4,col=4)
  # points(delta2,pch=5,col=5)
  points(psi,pch=6,col=6)
  
  legend("bottomright",c(expression(plain(hat(omega))),
# legend("topright",c(expression(plain(hat(omega))),
                         expression(plain(hat(alpha))), 
                         expression(plain(hat(beta))),
                         # expression(plain(hat(gamma)[1])),
                         # expression(plain(hat(gamma)[2])),
                         expression(plain(hat(c)))),
         lty=c(1,1,1,1,1,1
               ), pch=c(1,2,3,#4,
                       # 5,
                        6),col=c(1,2,3,#4,
                                 #5,
                                 6),bty="n",cex=1.25)
  axis(1,1:5,c(100,200,300,500,1000))
  box();axis(2)
} 
dev.off()


# RMSE
RMSE4<-paste("RMSE4",".pdf",sep="")
pdf(file = RMSE4,width = 4, height = 4,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.5, 0.5, 0.5))
  par(mgp=c(1.3, 0.45, 0))

  # # RMSE
  omega<-RMSE[,1]
  alpha<-RMSE[,2]
  beta<-RMSE[,3]
  # delta1<-RMSE[,4]
  # delta2<-RMSE[,5]
  psi<-RMSE[,4]
  
  plot(omega,xaxs="r",type="l",
       ylab=expression("RMSE"),
       xlab=expression(plain(italic(n))),
       axes=FALSE,
       ylim=c(min(c(RMSE)),max(c(RMSE))))
  
  lines(omega, lty=1,col=1)
  lines(alpha, lty=1,col=2)
  lines(beta, lty=1,col=3)
  # lines(delta1, lty=1,col=4)
  # lines(delta2, lty=1,col=5)
  lines(psi, lty=1,col=6)
  
  points(omega,pch=1,col=1)
  points(alpha,pch=2,col=2)
  points(beta,pch=3,col=3)
  # points(delta1,pch=4,col=4)
  # points(delta2,pch=5,col=5)
  points(psi,pch=6,col=6)
  
# legend("bottomright",c(expression(plain(hat(omega))),
  legend("topright",c(expression(plain(hat(omega))),
                         expression(plain(hat(alpha))), 
                         expression(plain(hat(beta))),
                         # expression(plain(hat(gamma)[1])),
                         # expression(plain(hat(gamma)[2])),
                         expression(plain(hat(c)))),
         lty=c(1,1,1,1,1,1
               ), pch=c(1,2,3,#4,
                        #5,
                        6),col=c(1,2,3,#4,
                                # 5,
                                 6),bty="n",cex=1.25)
  axis(1,1:5,c(100,200,300,500,1000))
  box();axis(2)
} 
dev.off()
