#Selection of u1 u2
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


#selection
Ipsi <- cu.f0.im
Rpsi <- cu.f0.real
Ipsi <- cu.f0.im
Rpsi(-0.5)

a=1.5
b=1 
d=0 
g=1
n=2000
# u väärtused graafikule üle mõlema
u <- c(seq(0, 1/10^8, by=1/10^10),seq(1/10^8, 1/10^6,by=1/10^8),seq(1/10^6,1/10^4,by=1/10^6),seq(1/10^4,1/10^2,by=1/10^4) ,seq(1/10^2,1.7,by=1/10^2))#,seq(1/10^2,0.1,by=1/10^3) 
u=c(-u,u)


a=1.8
r=0
v=0
Y18=rstable(n, alpha=a,beta=b, gamma=g, delta=d,param=0)
Y=Y18
#u18 <- seq(0, 2, by=1/100) #a=1.8 väärtused
for (i in 1:length(u)){# empirical  cumulant function real and imaginary part for resampled data
  k <- sum(cos(u[i]*Y))	
  s <- sum(sin(u[i]*Y))    
  r[i] <- log((1/n)*sqrt(k^2+s^2))#actually -r[i] as in formulas
  v[i] <-  atan2(s,k)}
v18=v
r18=r

a=0.2
r=0
v=0
#u02=seq(-0.0002,0.007,by=0.00001) #a=0.2 väärtused
Y02=rstable(n,alpha=a,beta=b, gamma=g, delta=d,param=0)
Y=Y02
for (i in 1:length(u)){# empirical  cumulant function real and imaginary part for resampled data
  k <- sum(cos(u[i]*Y))	
  s <- sum(sin(u[i]*Y))    
  r[i] <- log((1/n)*sqrt(k^2+s^2))#actually -r[i] as in formulas
  v[i] <-  atan2(s,k)}
v02=v
r02=r

a=1
b=1
r=0
v=0
si=0
ko=0
Y10=rstable(n, alpha=a,beta=b, gamma=g, delta=d,param=0)
Y=Y10
for (i in 1:length(u)){# empirical  cumulant function real and imaginary part for resampled data
  k <- sum(cos(u[i]*Y))/n
  s <- sum(sin(u[i]*Y))/n   
  r[i] <- log(sqrt(k^2+s^2))#actually -r[i] as in formulas
  v[i] <-  atan2(s,k)
  si[i]=s
  ko[i]=k}
v10=v
r10=r

#JOONIS 1
#REAALOSA
#par(mar=c(5.1,5,0,2.1))
g=1
d=0
u1=0
u2=2
#op <- par(mfcol = c(1,2),mar=c(3,5,1,5),cex.sub=1.2,mgp=c(2.0,0.7,0))
#op <- par(mfrow = c(1,2),mar=c(3,5,2,2.1),cex.sub=1.2,mgp=c(1.5,0.5,0))
#op <- par(mfrow = c(1,2),mar=c(3,5,2,2.1),cex.sub=1.2,mgp=c(1.5,0.5,0))
op <- par(mfcol = c(1,1
                    ),mar=c(2.9,2,0,0))
c18='blue'
c02='blue'
c10='blue'
a=0.2
p1=-0.1
p2=-0.5
curve(Rpsi,u1,u2,n=10001,ylim=c(-1.5,0),xlim=c(0,1.5), axes=F,xlab=paste("u"),frame=T,type='l',lty=2,col=c02, ylab='')#col=c1,cex.axis=1.2, cex.lab=1.2
lines(u,r10,col=1,type='l',lty=1,lwd=0.8)
lines(u,r02,col=1,type='l',lty=1,lwd=0.8)
lines(u,r18,col=1,lty=1,type='l',lwd=0.8)
abline(h=-0.5,col=2)
abline(h=-0.1,col=2)
ticks = c(-1.5,- 1, - 0.5, -0.1,0)
a=1.8
curve(Rpsi,u1,u2, add=T,xlab=paste("u"),lty=2,col=c18, ylab='')#col=c1,cex.axis=1.2, cex.lab=1.2
axis(side = 1)
axis(side = 2, at = ticks,las=1)
#abline(h=-1)
#a=0.9
a=1
curve(Rpsi,u1,u2,n=10001,add=T,lty=2,col=c18)
legend(1.1,-0.1,legend = c(expression(paste("Re(",hat(psi)[n],"(u))")),                  
                             expression(paste("Re(",psi,"(u))"))),
       lty=c(1,2), lwd=c(1,1), bty='n',col=c(1,'blue') )
text(1.33,-0.95,expression(paste(alpha, " = ", 0.2)))
text(1.31,-1.2,expression(paste(alpha, " = ", 1.0)))
text(1.33,-1.45,expression(paste(alpha, " = ", 1.8)))

p1=-0.1
p2=-0.5
u[min(which(r18<=p1))]
u[min(which(r18<=p2))]
u[min(which(r10<=p1))]
u[min(which(r10<=p2))]
u[min(which(r02<=p1))]
u[min(which(r02<=p2))]

#legend(x=0.78, y=-0.2,legend = c(expression(paste(alpha, " = ", 1.8)),
#      expression(paste(alpha, " = ", 0.2))),
#      lty=c(2,1) ,  bty='n',col=c(1,1))

#abline(h=-1,lwd=0.5,lty=1,col='red')
#text(x=0.15,y=-0.78,expression(paste(alpha, " = ", 0.2)),cex=0.8)
#text(x=0.77,y=-0.78,expression(paste(alpha, " = ", 1.8)),cex=0.8)
#abline(h=p1,lwd=0.5,lty=1,col='red')
#abline(h=p2,lwd=0.5,lty=1,col='red')
#abline(h=-1,lwd=0.5, col='grey')

#lines(u,r18,col=4,type='p',pch=20,size=1)
#points(1,-1,col='red',pch=19)
#Imaginaarosa
a=0.2
d=0
g=1
curve(Ipsi,u1,u2,add=F, axes=F,col=c02,xlab=paste("u"),ylab='', frame=T,lty=1,xlim=c(0,1.5),ylim=c(-0.2,0.3))
ticks = c(-0.2,0, 0.3)
axis(side = 1)
axis(side = 2, at = ticks,las=1)
lines(u,v02,col=1,lty=1,lwd=0.8)
lines(u,v18,col=1,lty=1,lwd=0.8)
lines(u,v10,col=1,lty=1,lwd=0.8)

a=1.8
curve(Ipsi,u1,u2,col=c18,axes=T, add=T,xlab=paste("u"),lty=1, ylab='')#col=c1,cex.axis=1.2, cex.lab=1.2

u[min(which(r18<=p1))]
u[min(which(r18<=p2))]
u[min(which(r18<=p1))]
u[min(which(r18<=p2))]
a=1
curve(Ipsi,u1,u2, col=c10,axes=T, add=T,xlab=paste("u"),lty=1, ylab='')#col=c1,cex.axis=1.2, cex.lab=1.2
points(u[min(which(r10<=p1))],v10[min(which(r10<=p1))],col=2,bg='red',pch=17,cex=1)
points(u[min(which(r10<=p2))],v10[min(which(r10<=p2))],col=2,bg='red',pch=17, cex=1)
points(u[min(which(r18<=p1))],v18[min(which(r18<=p1))],col=2,bg=2,pch=15,cex=1)
points(u[min(which(r18<=p2))],v18[min(which(r18<=p2))],col=2,bg=2,pch=15, cex=1)
points(u[min(which(r02<=p1))],v02[min(which(r02<=p1))],col=2,bg=2, pch=19,cex=1)
points(u[min(which(r02<=p2))],v02[min(which(r02<=p2))],col=2,pch=19,bg=2, cex=1)
u[min(which(r10<=p1))]
u[min(which(r10<=p2))]
#points(1,0,col='red',pch=19)
legend('topright',legend = c(expression(paste("Im(",hat(psi),"(u))")),                  
                             expression(paste("Im(",psi,"(u))"))),
       lty=c(1,1), lwd=c(2,2),bty='n',col=c(1,'grey') )



par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#legend(x=-0.2,y=-0.7, legend = c(expression(paste("Im(",psi,"(u))")),                  
#                                                             expression(paste("Im(",hat(psi),"(u))"))),
#                                    pch=c(15,15), col=c(1,4) , 
#    xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n",cex = 1)
legend("bottom", legend = c(expression(paste(alpha, " = ", 0.2)),
                            expression(paste(alpha, " = ", 1)),expression(paste(alpha, " = ", 1.8))),
       pch=c(19,17,15) ,  bty='n',col=c(2,2,2), 
       xpd = TRUE, horiz = T, inset = c(0, 0), text.width=0.15,cex = 1)#
#,
