################################################################################
# PAPER: Quantile-based dynamic modeling of asymmetric data: a novel Burr XII
#        approach for positive continuous random variables
# SECTION: 2. The Burr XII ARMA model
# GOAL: Plots pdf - Burr XII distribution.
# AUTHORS: Fernando Jose Monteiro de Araujo, Renata Rojas Guerra and 
#          Fernando Arturo Pena-Ramirez
# LAST UPDATE: May 25, 2024
################################################################################


################################################################################
##################  Reparametrized Burr XII density plots     ##################
################################################################################

# pdf of BXII - mu = 0.9 (Figure 1(a))

{
rm(list = ls())
c<-5
s<-1
t<-.25
m= 0.2

f_me<-function(y){log(1/(1-t))*(c*y^(c-1))/(s^(c)*log(1+(m/s)^c))*
    (1+(y/s)^c)^(log(1-t)/(log(1+(m/s)^c))-1)}

integrate(f_me,0,Inf)



h11 <- w1 <- 4

setEPS()



postscript(file = "dens1.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx=0.25
tox= 1.75
yliminf=0
ylimsup= 5.5

c<-5
t<-.5
m= 0.9
curve(f_me,from=fromx, to=tox,
      add = FALSE, lty=1, type = "l", cex.lab=1.3, xlab = expression("y"),
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup),
      col =1, lwd = 2.0)

c<-8

curve(f_me,from=fromx, to=tox, add = TRUE, lty=2, type = "l",
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup), col = 2, lwd = 2.0)

c<-10

curve(f_me,from=fromx, to=tox, add = T, lty=3, type = "l",
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup), col = 3, lwd = 2.0)

c<-13

curve(f_me,from=fromx, to=tox, add = TRUE, lty=5, type = "l",
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup), col = 5, lwd = 2.0)
c<-15

curve(f_me,from=fromx, to=tox, add = T, lty=1, type = "l", cex.lab=1.3,
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup), col ="purple", 
      lwd = 2.0)

legend("topright", c(expression(paste(plain(c), " = 5 ")),
                     c(expression(paste(plain(c), " = 8 ")),
                       c(expression(paste(plain(c), " = 10 ")),
                         c(expression(paste(plain(c), " = 13 ")),
                           c(expression(paste(plain(c), " = 15 "))))))),
       col = c(1,2,3,5,"purple"),
       lty= c(1,2,3,5,1),
       lwd = c(3,3,3,3,3,3), bty="n", cex = .8)
dev.off()
}

# pdf of BXII - mu = 0.5 (Figure 1(b))

{
  
postscript(file = "dens2.eps",width = w1, height = h11,family = "Times")
par(mar=c(5,6,4,1)+.1)
fromx=0
tox= 1
yliminf=0
ylimsup= 3.5

c<-0.05
t<-.5
m= 0.5

curve(f_me,from=fromx, to=tox,
      add = FALSE, lty=1, type = "l", cex.lab=1.3, xlab = expression("y"),
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup),
      col =1, lwd = 2.0)


c<-0.5

curve(f_me,from=fromx, to=tox, add = TRUE, lty=2, type = "l",
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup), col = 2, lwd = 2.0)

c<-1.5

curve(f_me,from=fromx, to=tox, add = T, lty=3, type = "l",
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup), col = 3, lwd = 2.0)

c<-4

curve(f_me,from=fromx, to=tox, add = TRUE, lty=5, type = "l",
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup), col = 5, lwd = 2.0)

c<-5

curve(f_me,from=fromx, to=tox, add = T, lty=1, type = "l", cex.lab=1.3,
      ylab = expression("f(y)"),ylim =c(yliminf,ylimsup), col ="purple", 
      lwd = 2.0)

legend("topright", c(expression(paste(plain(c), " = 0.05 ")),
                     c(expression(paste(plain(c), " = 0.5 ")),
                       c(expression(paste(plain(c), " = 1.5 ")),
                         c(expression(paste(plain(c), " = 4 ")),
                           c(expression(paste(plain(c), " = 5 "))))))),
       col = c(1,2,3,5,"purple"),
       lty= c(1,2,3,5,1),
       lwd = c(3,3,3,3,3,3), bty="n", cex = .8)
dev.off()
}

