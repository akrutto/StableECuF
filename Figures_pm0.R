#characterictic function (pm=0, (M)) 
ch.f0.neq.1 <- function(u) {(complex(
  modulus=exp(-(g*abs(u))^a),
  argument=-((g*abs(u))^a)*b*tan(a*pi/2)*sign(u)*((abs(g*(u+10^(-25)))^(1-a))-1)+d*u
))}
ch.f0.eq.1 <- function(u) {(complex(
  modulus=exp(-(g*abs(u))),
  argument=d*u-((g*abs(u)))*b*(2/pi)*sign(u)*log(g*abs(u+10^(-25)))))}
ch.f0. <- function(u) {if(a==1) {ch.f0.eq.1(u)} else {ch.f0.neq.1(u)}}
ch.f0.mod <- function(u) {if(a==1) {Mod(ch.f0.eq.1(u))} else {Mod(ch.f0.neq.1(u))}}
ch.f0.real <- function(u) {if(a==1) {Re(ch.f0.eq.1(u))} else {Re(ch.f0.neq.1(u))}}
ch.f0.im <- function(u) {if(a==1) {Im(ch.f0.eq.1(u))} else {Im(ch.f0.neq.1(u))}}
arctan2 <- function(u) {atan2(ch.f0.im(u),ch.f0.real(u))}
arctan2(0)

#cumulant function (pm=0, M) 
cu.f0.neq.1 <- function(u) {(complex(
  real=-(g*abs(u))^a,
  imaginary=-((g*abs(u))^a)*b*tan(a*pi/2)*sign(u)*((g*abs(u))^(1-a)-1)+d*u #abs((u+10^(-25)))^(1-a)
))}
cu.f0.eq.1 <- function(u) {(complex(
  real=-(g*abs(u)),
  imaginary=u*(-b*g*(2/pi)*log(g*abs(u))+d)
))}
cu.f0. <- function(u) {if(a==1) {cu.f0.eq.1(u)} else {cu.f0.neq.1(u)}}
cu.f0.mod <- function(u) {if(a==1) {Mod(cu.f0.eq.1(u))} else {Mod(cu.f0.neq.1(u))}}
cu.f0.real <- function(u) {if(a==1) {Re(cu.f0.eq.1(u))} else {Re(cu.f0.neq.1(u))}}
cu.f0.im <- function(u) {if(a==1) {Im(cu.f0.eq.1(u))} else {Im(cu.f0.neq.1(u))}}
Ipsi <- cu.f0.im

#check
g=0.5
d=0
d=4
b=1
a=0.2
g*cu.f0.im(1/g)
a=1
a=1.8




cu.f0.mod(1)
cu.f0.real(1)
u=0
a=2
cu.f0.im(0)
tan(pi)
#Graphs: Characteristic function 
CH=0
a=0.3
for (i in 1:length(u)){CH[i]=ch.f0.neq.1(u[i])}
plot(u,CH)

op <- par(mfcol = c(5,2),mar=c(5, 4, 1, 1) + 0.1)
b=1
g=1
d=0
a=2
u <- seq(-100, 100, by=0.1)
plot((ch.f0.(u)),type='l',main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1.8
plot(ch.f0.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1.5
plot(ch.f0.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1.2
plot(ch.f0.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1.1
plot(ch.f0.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=1
plot(ch.f0.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )), 
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=0.9
plot(ch.f0.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))    
a=0.8
plot(ch.f0.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=0.5
plot(ch.f0.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))
a=0.2
plot(ch.f0.(u),type='l',
     main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) )),
     ylab=bquote(paste("Im(",varphi,"(u))")),xlab=bquote(paste("Re(",varphi,"(u))")))

#Graphs: Abs, Real, Im of characteristic function 

#main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) ))
op <- par(mfcol = c(5,2),mar=c(5, 4, 1, 1) + 0.1)
u1=-10
u2=10
b=1
g=1
d=0
a=2
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#legend(x=-7,y=1.5,legend=bquote(paste(alpha ==.(a),'   ','   ')),bty='n')
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.8
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.5
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.2
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.1
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.9
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.8
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.5
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.2
curve(ch.f0.mod,u1, u2,n=10001,ylim=c(-0.8,1),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(varphi(u))),xlim=c(-6,6)) 
curve(ch.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(ch.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=5.0,y=0.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')


op <- par(mfcol = c(5,2),mar=c(5, 4, 1, 1) + 0.1)
u1=-10
u2=10
b=1
g=1
d=0

#Graphs: Abs, Real, Im of cumulant function 

#main=bquote(paste(alpha ==.(a),", ",beta ==.(b),", ",gamma ==.(g),", ",delta ==.(d) ))
op <- par(mfcol = c(5,2),mar=c(5, 4, 1, 1) + 0.1)
u1=-10
u2=10
b=1
g=1
d=0
a=2
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=2.7,y=1.8,bquote(paste(alpha ==.(a),'   ','   ')))
#legend(x=-7,y=1.5,legend=bquote(paste(alpha ==.(a),'   ','   ')),bty='n')
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.8
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=2.7,y=1.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.5
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=2.7,y=1.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.2
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=2.8,y=1.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1.1
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=2.8,y=1.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=1
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=2.9,y=1.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.9
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=2.9,y=1.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.8
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=2.9,y=1.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.5
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=3,y=1,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')
a=0.2
curve(cu.f0.mod,u1, u2,n=10001,ylim=c(-2,2),xlab='u',lty = 1,lwd=1,
      ylab=bquote(paste(psi(u))),xlim=c(-3,3)) 
curve(cu.f0.real,u1,u2,n=10001,add=T,lty=2,lwd=1)
curve(cu.f0.im, u1,u2,n=10001,add=T,lty=3,lwd=1.5)
text(x=2.6,y=1.8,bquote(paste(alpha ==.(a),'   ','   ')))
#abline(h=0, lwd=0.2)
#legend("topright",cex=0.8,legend=c('Absolute value', 'Real part', 'Imaginary part'),lty=1:3,bty='n')


#Cumulant function Figures KOKKU
Rpsi <- cu.f0.real
u1=-1
u2=
b=1 #not relevant for theoretical, but influence the empricial 
d=0 #not relevant
#par(mar=c(5.1,5,0,2.1))
#op <- par(mfcol = c(1,1),mar=c(3,5,1,2.1),cex.sub=1.2,mgp=c(1.5,0.5,0))
op <- par(mfcol = c(1,2),mar=c(4,4,1,1),cex.sub=1.2,mgp=c(2.5,0.7,0),oma = c(4, 1, 1, 1),cex=1.05)
a=0.2
g=1
curve(Rpsi,u1,u2,n=10001,axes=F,ylab='',ylim=c(-1.2,0),frame=T,xlim=c(0,2.5),xlab=paste("u"),lty=1,lwd=1,col=1)#col=c1,cex.axis=1.2, cex.lab=1.2
#abline(h=-1)
ticks = c(- 1.5, -1,- 0.5, 0,0.5)
axis(side = 1,at = c(-0.7,0,0.5,1,2,3.5),labels=c('',0,0.5,1.0,2.0,''))
axis(side = 2, at = ticks,labels=c('',-1, - 0.5, 0,''),las=1)
a=1.8
curve(Rpsi,u1,u2,n=10001,add=T,lty=2,col=1)
g=2
a=0.2
curve(Rpsi,u1,u2,n=10001,add=T,lty=1,lwd=1,col=1)
a=1.8
curve(Rpsi,u1,u2,n=10001,add=T,lty=2,lwd=1,col=1)
g=1/2
a=0.2
curve(Rpsi,u1,u2,n=10001,add=T,lty=1,lwd=1,col=1)
a=1.8
curve(Rpsi,u1,u2,n=10001,add=T,lty=2,lwd=1,col=1)
points((1/1), -1,col=4,pch=15,cex=1)
points((1/2), -1,col=4,pch=19,cex=1)
points((1/0.5), -1,col=4,pch=17,cex=1)

#Imagianry part. Graphics.
Ipsi <- cu.f0.im
d=0
u1=-1
u2=4
a=1.8
g=1
curve(Ipsi,u1,u2,n=10001,ylab='',axes=F,xlab=paste("u"),frame=T,xlim=c(-0.1,3.5),ylim=c(-0.2,0.2), col=1,lty=2)#,
ticks = c(-0.3,-0.2, 0,0.2,  0.3)
axis(side = 1,at = c(-0.5,0,0.5,1,2,4),labels=c('',0,0.5,1.0,2.0,''))
axis(side = 2, at = ticks,labels=c('',-0.2,0,0.2,''),las=1)
g=2
curve(Ipsi,u1,u2,n=10001,add=T,lty=2,lwd=1, col=1)
g=0.5
curve(Ipsi,u1,u2,n=10001,add=T,lwd=1,lty=2, col=1)
a=0.2
g=1
curve(Ipsi,u1,u2,n=10001, add=T,col=1,lty=1,lwd=1)#,
g=1/2
curve(Ipsi,u1,u2,n=10001,add=T,lty=1,lwd=1,col=1)
g=2
curve(Ipsi,u1,u2,n=10001,add=T,lty=1,lwd=1,col=1)
g=1
points(1/g,0,col=4,pch=15,cex=1)
g=2
points(1/g,0, col=4,pch=19,cex=1)
g=0.5
points(1/g,0,col=4,pch=17,cex=1)

#abline(h=0)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(x=-0.31,y=-0.7, legend = c(expression(paste(gamma, " = ", 2)),
                                  expression(paste(gamma, " = ", 1)),expression(paste(gamma, " = ", 0.5))), 
       xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n",  pch = c(16,15,17), col = c(4,4,4),cex = 1.1,pt.cex=c(1,1,1))
legend("bottom", legend = c(expression(paste(alpha, " = ", 0.2)),
                            expression(paste(alpha, " = ", 1.8))), 
       xpd = TRUE, horiz = T, inset = c(0, 0), bty = "n", text.width=0.15,  lwd=c(1,1),lty=c(1,2),col=c(1,1),cex = 1.1)#
##
##LÕPP
##


#Cumulant function Figures
Rpsi <- cu.f0.real
u1=-2
u2=3
b=0 #not relevant for theoretical, but influence the empricial 
d=4 #not relevant
#par(mar=c(5.1,5,0,2.1))
#op <- par(mfcol = c(1,1),mar=c(3,5,1,2.1),cex.sub=1.2,mgp=c(1.5,0.5,0))
op <- par(mfcol = c(1,2),mar=c(4,4,1,1),cex.sub=1.2,mgp=c(2.5,0.7,0),oma = c(1.5, 1, 1, 1))
g=1
c1=1
a=0.2
curve(Rpsi,u1,u2,n=10001,ylab=bquote(paste("Re(",psi,"(u))")),ylim=c(-1.2,0),xlim=c(0,2.2),xlab=paste("u"),lty=1,col=4)#col=c1,cex.axis=1.2, cex.lab=1.2
#abline(h=-1)
a=1.8
curve(Rpsi,u1,u2,n=10001,add=T,lty=1,col=51)
#a=1.5
#curve(Rpsi,u1,u2,n=10001,add=T,col=c1)
#a=1
#curve(Rpsi,u1,u2,n=10001,add=T,col=c1)
#a=0.8
#curve(Rpsi,u1,u2,n=10001,add=T,col=c1)
#a=0.5
#curve(Rpsi,u1,u2,n=10001,add=T,col=c1)
#a=0.2
#curve(Rpsi,u1,u2,n=10001,add=T,col=c1)
g=2
a=0.2
c2=1
l=1.5
lt2=2
curve(Rpsi,u1,u2,n=10001,add=T,lty=lt2,lwd=l,col=4)
a=1.8
curve(Rpsi,u1,u2,n=10001,add=T,lty=lt2,lwd=l,col=51)
#a=0.5
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt2,col=c2)
#a=1
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt2,col=c2)
#a=0.8
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt2,col=c2)
#a=0.5
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt2,col=c2)
#a=0.2
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt2,col=c2)
g=1/2
c3=1
a=0.2
l=2
lt3=3
curve(Rpsi,u1,u2,n=10001,add=T,lty=lt3,lwd=l,col=4)
a=1.8
curve(Rpsi,u1,u2,n=10001,add=T,lty=lt3,lwd=l,col=51)
#a=0.5
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt3,lwd=l,col=c3)
#a=1
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt3,lwd=l,col=c3)
#a=0.8
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt3,lwd=l,col=c3)
#a=0.5
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt3,lwd=l,col=c3)
#a=0.2
#curve(Rpsi,u1,u2,n=10001,add=T,lty=lt3,lwd=l,col=c3)
#abline(v=0.0,lwd=0.5)
#abline(v=0.5,lwd=0.5)
#abline(v=1,lwd=0.5)
#abline(v=2,lwd=0.5)
legend(y=0.05,x=1.7,legend = c(expression(paste(gamma, " = ", 0.5)),
                              expression(paste(gamma, " = ", 1)),
                              expression(paste(gamma, " = ", 2)))
                            , lty=c(lt3,1,lt2),cex=1, bty='n' ,lwd=c(2,1,1.5))
legend(x=1.8,y=-0.3,legend = c(expression(paste(alpha, " = ", 1.8)),
                               expression(paste(alpha, " = ", 0.2))),
       pch=c(15,15), bty='n',col=c(51,4) )
abline(h=-1,lwd=0.5,lty=1,col='red')



abline(v=0,lty=6)
abline(h=-0.1,lwd=0.5,col=2)
abline(h=-0.3,lwd=0.5,col=2)
legend(x=1.5,y=-0.1,legend = c(expression(paste(u[1])),
                            expression(paste(u[2]))),lty=c(1,1),cex=1.3, bty='n',lwd=c(0.5,0.5),col=c(2,2))

abline(v=0.03,lwd=0.5,col=2)
abline(v=0.09,lwd=0.5,col=2)



#Imagianry part. Graphics.

Ipsi <- cu.f0.im

#FIGURE22 ALGUS
b=1  
l05=2
l2=1.5
#op <- par(mfrow = c(1,2),mar=c(3,5,2,2.1),cex.sub=1.2,mgp=c(1.5,0.5,0))
d=0
u1=-1
u2=4
a=1.8
c18=51
g=1
curve(Ipsi,u1,u2,n=10001,ylab=bquote(paste("Im(",psi,"(u))")),xlab=paste("u"),xlim=c(-0.1,3.5),ylim=c(-0.2,0.2), col=c18)#,
points(1/g,Ipsi(1/g),col='red',pch=19)
g=2
lt2=2
Ipsi(1/g)
curve(Ipsi,u1,u2,n=10001,add=T,lty=lt2,lwd=l2, col=c18)
points(1/g,Ipsi(1/g), col='red',pch=19)
g=0.5
lt3=3
Ipsi(1/g)
curve(Ipsi,u1,u2,n=10001,add=T,lwd=l05,lty=lt3, col=c18)
points(1/g,Ipsi(1/g),col='red',pch=19)
a=0.2
g=1
c02=4
curve(Ipsi,u1,u2,n=10001, add=T,col=c02)#,
g=1/2
curve(Ipsi,u1,u2,n=10001,add=T,lwd=l05,lty=lt3,col=c02)
g=2
c2=3
lt2=2
curve(Ipsi,u1,u2,n=10001,add=T,lty=lt2,lwd=l2,col=c02)

legend(x=2.1,y=0.21, legend = c(expression(paste(gamma, " = ", 0.5)),
                                expression(paste(gamma, " = ", 1)),
                                expression(paste(gamma, " = ", 2))), lty=c(lt3,1,lt2),bty='n',cex=1,lwd=c(2,1,1.5))
legend(x=2.4,y=0.1,legend = c(expression(paste(alpha, " = ", 1.8)),
                              expression(paste(alpha, " = ", 0.2))),
       pch=c(15,15), bty='n',col=c(51,4))
mtext('(a)',outer=F,side=3,line=0.5,font=2)#expression(paste(delta^0, " = ", 0))

#TEINE JOONIS
d=1
u1=-0.8
u2=50
a=1.8
c18=51
g=1
curve(Ipsi,u1,u2,n=10001,ylab=bquote(paste("Im(",psi,"(u))")),xlab=paste("u"),xlim=c(-0,6),ylim=c(0,5), col=c18)#,
points(1/g,Ipsi(1/g),col='red',pch=19)
g=2
lt2=2
Ipsi(1/g)
curve(Ipsi,u1,u2,n=10001,add=T,lty=lt2,lwd=l2, col=c18)
points(1/g,Ipsi(1/g), col='red',pch=19)
g=0.5
lt3=3
Ipsi(1/g)
curve(Ipsi,u1,u2,n=10001,add=T,lwd=l05,lty=lt3, col=c18)
points(1/g,Ipsi(1/g),col='red',pch=19)
a=0.2
g=1
c02=4
curve(Ipsi,u1,u2,n=10001, add=T,col=c02)#,
g=1/2
curve(Ipsi,u1,u2,n=10001,add=T,lwd=l05,lty=lt3,col=c02)
g=2
c2=3
lt2=2
curve(Ipsi,u1,u2,n=10001,add=T,lty=lt2,lwd=l2,col=c02)


legend(x=0,y=5.2, legend = c(expression(paste(gamma, " = ", 0.5)),
                                expression(paste(gamma, " = ", 1)),
                                expression(paste(gamma, " = ", 2))), lty=c(lt3,1,lt2),bty='n',cex=1,lwd=c(2,1,1.5))
legend(x=0,y=3.7,legend = c(expression(paste(alpha, " = ", 1.8)),
                              expression(paste(alpha, " = ", 0.2))),
       pch=c(15,15), bty='n',col=c(51,4) )
mtext('(b)',outer=F,side=3,line=0.5
      ,font=2)#expression(paste(delta^0, " = ", 1))

#FIGURE22 LÕPP

a=0.2
d=10
g=1
curve(Ipsi,u1,u2,n=10001,ylab=bquote(paste("Im(",psi,"(u))")),lty=3,xlab=paste("u"),xlim=c(0,1),ylim=c(-0.1,1))#,
#a=1.5
#curve(Ipsi,u1,u2,n=10001,add=T,col=2)
g=1/2
curve(Ipsi,u1,u2,n=10001,add=T,lty=2)
#a=0.8
#curve(Ipsi,u1,u2,n=10001,add=T,col=4)
g=2
curve(Ipsi,u1,u2,n=10001,add=T,lwd=1.5)
#a=0.2
#curve(Ipsi,u1,u2,n=10001,add=Tcol=6)
#legend('top',legend=bquote(paste(alpha ==.(a), '   ',gamma ==.(g)  ,'   ')))
#text(y=1,x=0.5,bquote(paste(gamma ==.(g) )))
g=0.5
d=10
curve(Ipsi,u1,u2,n=10001,ylab=bquote(paste("Im(",psi,"(u))")),lty=3,xlab=paste("u"),ylim=c(-0.1,1),xlim=c(0,1))#
#a=1.5
#curve(Ipsi,u1,u2,n=10001,add=T,col=2)
d=1
curve(Ipsi,u1,u2,n=10001,add=T,lty=2)
#a=0.8
#curve(Ipsi,u1,u2,n=10001,add=T,col=4)
d=0
curve(Ipsi,u1,u2,n=10001,add=T,lwd=1.5)
#a=0.2
#curve(Ipsi,u1,u2,n=10001,add=T,col=6)
legend(y=0.4,x=0.6, legend = c(expression(paste(delta, " = ", 10)),
                              expression(paste(delta, " = ", 1)),
                              expression(paste(delta, " = ", 0))), lty=c(3,2,1),bty='n',cex=1,lwd=c(1,1,1.5))
#legend('top',legend=bquote(paste(alpha ==.(a), '   ',gamma ==.(g)  ,'   ')))
#text(y=0.9,x=0.5,bquote(paste(gamma ==.(g) )))
mtext('(b)',outer=F,side=3,line=1,font=2)

#a=0.2
op <- par(mfrow = c(1,2))
g=2
a=0.2
d=10
curve(Ipsi,u1,u2,n=10001,ylab=bquote(paste("Im(",psi,"(u))")),lty=3,xlab=paste("u"),xlim=c(0,1),ylim=c(-0.1,1))#,
#a=1.5
#curve(Ipsi,u1,u2,n=10001,add=T,col=2)
d=1
curve(Ipsi,u1,u2,n=10001,add=T,lty=2)
#a=0.8
#curve(Ipsi,u1,u2,n=10001,add=T,col=4)
d=0
curve(Ipsi,u1,u2,n=10001,add=T,lwd=1.5)
#a=0.2
#curve(Ipsi,u1,u2,n=10001,add=Tcol=6)
legend(y=0.4,x=0.6, legend = c(expression(paste(d, " = ", 0)),
                           expression(paste(d, " = ", 1)),
                           expression(paste(d, " = ", 10))), lty=c(1,2,3),bty='n',cex=1,lwd=c(1,1,1.5))
#legend('top',legend=bquote(paste(gamma ==.(g))))
#text(y=1,x=0.5,bquote(paste(gamma ==.(g) )))
mtext('(a)',outer=F,side=3,line=1,font=2)
g=0.5
d=10
curve(Ipsi,u1,u2,n=10001,ylab=bquote(paste("Im(",psi,"(u))")),lty=3,xlab=paste("u"),ylim=c(-0.1,1),xlim=c(0,1))#
#a=1.5
#curve(Ipsi,u1,u2,n=10001,add=T,col=2)
d=1
curve(Ipsi,u1,u2,n=10001,add=T,lty=2)
#a=0.8
#curve(Ipsi,u1,u2,n=10001,add=T,col=4)
d=0
curve(Ipsi,u1,u2,n=10001,add=T,lwd=1.5)
#a=0.2
#curve(Ipsi,u1,u2,n=10001,add=T,col=6)
legend('right', legend = c(expression(paste(d, " = ", 0)),
                           expression(paste(d, " = ", 1)),
                           expression(paste(d, " = ", 10))), lty=c(1,2,3),bty='n',cex=1,lwd=c(1,1,1.5))
#legend('top',legend=bquote(paste(gamma ==.(g) )))
#text(y=1,x=0.5,bquote(paste(gamma ==.(g) )))
mtext('(b)',outer=F,side=3,line=1,font=2)
