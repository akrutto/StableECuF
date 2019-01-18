#charactericticfunction (pm=1, (A)) 
ch.f1.neq.1 <- function(u) {(complex(
  modulus=exp(-(g*abs(u))^a),
  argument=d*u+((g*abs(u))^a)*b*tan(a*pi/2)*sign(u)
))}
ch.f1.eq.1 <- function(u) {exp(complex(
  real=(-(g*abs(u))),
  imaginary=d*u-g*abs(u)*b*(2/pi)*sign(u)*log(abs(u+10^(-25)))
))}
ch.f1. <- function(u) {if(a==1) {ch.f1.eq.1(u)} else {ch.f1.neq.1(u)}}
ch.f1.mod <- function(u) {if(a==1) {Mod(ch.f1.eq.1(u))} else {Mod(ch.f1.neq.1(u))}}
ch.f1.real <- function(u) {if(a==1) {Re(ch.f1.eq.1(u))} else {Re(ch.f1.neq.1(u))}}
ch.f1.im <- function(u) {if(a==1) {Im(ch.f1.eq.1(u))} else {Im(ch.f1.neq.1(u))}}

a=1
g=1
b=0
d=0
curve(exp(-(g*abs(x))^a),-100,100)
curve(Mod(ch.f1.neq.1(x)),-100,100,add=T)
curve(Re(ch.f1.neq.1(x)),-100,100,n=1001)
curve(exp(-(g*abs(x))),-100,100,n=10001,add=T)
curve(Mod(ch.f1.eq.1(x)),-100,100,n=10001,ylim=c(-1,1))
curve(d*x+g*abs(x)*b*(2/pi)*sign(x)*log(abs(x)),-100,100,n=10001)
curve(Arg(ch.f1.eq.1(x))+pi,-100,100,n=10001)
curve(Im(ch.f1.eq.1(x)),-10,10,n=10001,add=T)
curve(Re(ch.f1.eq.1(x)),-10,10,n=10001,ylim=c(-1,1))


g=10000
b=1
d=0
ch.f1.real(1)

#cumulant function (pm=1, M) 
cu.f1.neq.1 <- function(u) {(complex(
  real=-(g*abs(u))^a,
  imaginary=d*u+((g*abs(u))^a)*b*tan(a*pi/2)*sign(u)
))}
cu.f1.eq.1 <- function(u) {(complex(
  real=-(g*abs(u))^a,
  imaginary=d*u-((g*abs(u)))*b*(2/pi)*sign(u)*log(abs(u+10^(-25)))))}
cu.f1. <- function(u) {if(a==1) {cu.f1.eq.1(u)} else {cu.f1.neq.1(u)}}
cu.f1.mod <- function(u) {if(a==1) {Mod(cu.f1.eq.1(u))} else {Mod(cu.f1.neq.1(u))}}
cu.f1.real <- function(u) {if(a==1) {Re(cu.f1.eq.1(u))} else {Re(cu.f1.neq.1(u))}}
cu.f1.im <- function(u) {if(a==1) {Im(cu.f1.eq.1(u))} else {Im(cu.f1.neq.1(u))}}

cu.f1.real(1)

#Graphs: Characteristic function 
op <- par(mfcol = c(5,2),mar=c(5, 4, 1, 1) + 0.1)
b=1
g=1
d=0

a=2
u <- seq(-100, 100, by=0.1)
plot(ch.f1.(u),type='l',main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste( "Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1.8
plot(ch.f1.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1.5
#u <- seq(-100*pi, 100*pi, by=0.1)
plot(ch.f1.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1.2
plot(ch.f1.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1.1
plot(ch.f1.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1
plot(ch.f1.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )), 
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=0.9
plot(ch.f1.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))    
a=0.8
plot(ch.f1.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=0.5
plot(ch.f1.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=0.2
u <- seq(-100, 100,by=0.1)
plot(ch.f1.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
(0.5)^0.1
(2)^0.1
(1.1)^0.1
(0.9)^0.1
(0.5)^2
(2)^2
10000^0.1
#Graphs: Abs, Real, Im of characteristic function 

#main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) ))
op <- par(mfcol = c(5,2),mar=c(5, 4, 1, 1) + 0.1)
u1=-10
u2=10
b=1
g=1
d=0
a=2
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#legend(x=-7,y=1.5,legend=bquote(paste(alpha ==.(a),'   ','   ')),bty='n')
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.8
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.5
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.2
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.1
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.9
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.8
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.5
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.2
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Phi(u))),xlim=c(-6,6)) 
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')



op <- par(mfcol = c(5,2),mar=c(5, 4, 1, 1) + 0.1)
u1=-10
u2=10
b=0
g=1
d=0

a=2
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#legend(x=-7,y=1.5,legend=bquote(paste(alpha ==.(a),'   ','   ')),bty='n')
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.8
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.5
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.2
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.1
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.9
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.8
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.5
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.2
curve(cu.f1.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(Psi(u))),xlim=c(-3,3)) 
curve(cu.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=4.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')




#Cumulant function Figures
op <- par(mfcol = c(1,1))
Rpsi<-cu.f1.real
u1=-1
u2=2.2
b=1  
g=1
d=0
a=2
curve(Rpsi,u1,u2,n=10001,ylab=bquote(paste("Re(",psi,"(u))")),ylim=c(-2,0),xlim=c(-0.05,2.2))
a=1.8
curve(Rpsi,u1,u2,n=10001,add=T)
a=1.5
curve(Rpsi,u1,u2,n=10001,add=T)
a=1
curve(Rpsi,u1,u2,n=10001,add=T)
a=0.8
curve(Rpsi,u1,u2,n=10001,add=T)
a=0.5
curve(Rpsi,u1,u2,n=10001,add=T)
a=0.3
curve(Rpsi,u1,u2,n=10001,add=T)
abline(v=0.0)
abline(v=0.5)
abline(v=1)
abline(h=-1)
g=2
d=0
a=2
curve(Rpsi,u1,u2,add=T,col=2,lty=2)
a=1.8
curve(Rpsi,u1,u2,add=T,col=2,lty=2)
a=1.5
curve(Rpsi,u1,u2,add=T,col=2,lty=2)
a=1
curve(Rpsi,u1,u2,add=T,col=2,lty=2)
a=0.8
curve(Rpsi,u1,u2,add=T,col=2,lty=2)
a=0.5
curve(Rpsi,u1,u2,add=T,col=2,lty=2)
a=0.3
curve(Rpsi,u1,u2,add=T,col=2,lty=2)
abline(v=0.0)
abline(v=2)
abline(v=1)
abline(h=-1)
g=1/2
d=0
a=2
a=2
curve(Rpsi,u1,u2,add=T,col=3,lty=3)
a=1.8
curve(Rpsi,u1,u2,add=T,col=3,lty=3)
a=1.5
curve(Rpsi,u1,u2,add=T,col=3,lty=3)
a=1
curve(Rpsi,u1,u2,add=T,col=3,lty=3)
a=0.8
curve(Rpsi,u1,u2,add=T,col=3,lty=3)
a=0.5
curve(Rpsi,u1,u2,add=T,col=3,lty=3)
a=0.3
curve(Rpsi,u1,u2,add=T,col=3,lty=3)
abline(v=0.0)
abline(v=0.5)
abline(v=1)
abline(h=-1)
legend("topright",cex=0.8,legend=c('gamma=1', 'gamma=1/2', 'gamma=2'),lty=1:3,bty='n',col=1:3)




Ipsi <- cu.f.im
Ipsi(1)
u1=-2*pi
u2=2*pi
b=1
g=10
d=20
a=0.3
curve(Ipsi,u1,u2,col=7)
a=1.8
curve(Ipsi,u1,u2,add=T,col=1)
a=1.5
curve(Ipsi,u1,u2,add=T,col=1)
a=1
curve(Ipsi,u1,u2,add=T,col=1)
a=0.8
curve(Ipsi,u1,u2,add=T,col=1)
a=0.5
curve(Ipsi,u1,u2,add=T,col=1)
a=2
curve(Ipsi,u1,u2,add=T)
g=1
d=0
a=0.3
curve(Ipsi,u1,u2,col=7,add=T)
a=1.8
curve(Ipsi,u1,u2,add=T,col=2)
a=1.5
curve(Ipsi,u1,u2,add=T,col=3)
a=1
curve(Ipsi,u1,u2,add=T,col=4)
a=0.8
curve(Ipsi,u1,u2,add=T,col=5)
a=0.5
curve(Ipsi,u1,u2,add=T,col=6)
a=2
curve(Ipsi,u1,u2,add=T)


#Graphs: Abs, Real, Im of fharacteristic function 
op <- par(mfcol = c(5,2))
u1=-100
u2=100
b=1
g=1
d=0

a=2
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi(u))),xlim=c(-6,6),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.8
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi (u) )),xlim=c(-6,6),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.5
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi (u) )),xlim=c(-6,6),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.2
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi (u) )),xlim=c(-6,6),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.1
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi (u) )),xlim=c(-6,6),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1.1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi (u) )),xlim=c(-6,6),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.9
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi (u) )),xlim=c(-6,6),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.8
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi (u) )),xlim=c(-6,6),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.5
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi (u) )),xlim=c(-20,20),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.2
curve(ch.f1.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,ylab=bquote(paste(varphi (u) )),xlim=c(-20,20),
      main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )))
curve(ch.f1.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f1.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')

Y=danish
hatpsi <- function (u) {{# empirical  cumulant function real and imaginary part for resampled data
  k <- sum(cos(u*Y))	
  s <- sum(sin(u*Y))    
  r <- -log((1/N)*sqrt(k^2+s^2))#actually -r[i] as in formulas
  v <-  atan2(s,k)}#end empirical  cumulant function real and imaginary part
  return(r)}
Y=rnorm(100)
curve(r,0,1)
hatpsi(0)
x=c(0,0.5,1)
plot(x, sum(cos(x*Y)))


