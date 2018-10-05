
Re.cu <- function (u,X){
  uX <- X%*%u 
  N <- length(X)
  return (log((1/N)*sqrt((colSums(cos(uX)))^2+(colSums(sin(uX)))^2)))
}
Im.cu <- function (u,X){
  uX <- X%*%u 
  N <- length(X)
  return (atan2(colSums(sin(uX))/N,colSums(cos(uX))/N))
}
Re.ch <- function (u,X){
  uX <- X%*%u
  N <- length(X)
  return ((1/N)*(colSums(cos(uX))))
}
Im.ch <- function (u,X){
  uX <- X%*%u
  N <- length(X)
  return ((1/N)*(colSums(sin(uX))))
}

#Asymptotic normality
#LAMBDA
N <- length(X)
Lambda2 <- function(u1,u2)
{N <- length(X)  
L <- matrix(NA,nrow=4,ncol=4)
L[1,1] <- 1+Re.ch(2*u1,X)-2*(Re.ch(u1,X))^2 #R1R1
L[1,2] <- Re.ch(u1-u2,X)+ Re.ch(u1+u2,X) -2*Re.ch(u1,X)*Re.ch(u2,X) #R1R2
L[1,3] <- Im.ch(u1+u1,X)-0-2*Re.ch(u1,X)*Im.ch(u1,X) #R1I1
L[1,4] <- Im.ch(u1+u2,X)-Im.ch(u1-u2,X)-2*Re.ch(u1,X)*Im.ch(u2,X) #R1I2
L[2,2] <- 1+Re.ch(2*u2,X)-2*(Re.ch(u2,X))^2 #R2R2
L[2,3] <- Im.ch(u2+u1,X)-Im.ch(u2-u1,X)-2*Re.ch(u2,X)*Im.ch(u1,X) ##R2I1
L[2,4] <- Im.ch(u2+u2,X)-0-2*Re.ch(u2,X)*Im.ch(u2,X) #R1I1
L[3,3] <- 1-Re.ch(2*u1,X)-2*(Im.ch(u1,X))^2
L[3,4] <- Re.ch(u1-u2,X)- Re.ch(u1+u2,X) -2*Im.ch(u1,X)*Im.ch(u2,X) #I1I2
L[4,4] <- 1-Re.ch(2*u2,X)-2*(Re.ch(u2,X))^2 #I2I2
L[lower.tri(L)] <- t(L)[lower.tri(L)]
L <- L/(2) 
return((L))
}
#Xi
xi2 <- function(u1,u2){
xivec  <- matrix(NA,nrow=4,ncol=1)
xivec[1,1]  <- Re.ch(u1,X)/exp(2*Re.cu(u1,X))*Re.cu(u1,X)*log(u1/u2)
xivec[2,1]  <- -Re.ch(u2,X)/exp(2*Re.cu(u2,X))*Re.cu(u2,X)*log(u1/u2)
xivec[3,1]   <- Im.ch(u1,X)/exp(2*Re.cu(u1,X))*Re.cu(u1,X)*log(u1/u2)
xivec[4,1] <-  -Im.ch(u2,X)/exp(2*Re.cu(u2,X))*Re.cu(u2,X)*log(u1/u2)
return(xivec)
}

N=10000
af=1.2
bt=0.9
X <- rstable(N,alpha=af,beta=bt,param=0)
Y=X
X <- matrix(X,nrow=length(X),ncol=1)
Lambda2(0.03,0.09)
t(xi2(0.03,0.09))%*%Lambda2(0.03,0.09)%*%xi2(0.03,0.09)
AV2out(0.03,0.09)
AV2out(u[1],u[2])
AV2out <- function(u1,u2){ #for outer function
  #u
  #u1 <- u[1]
  #u2 <- u[2]
  AV<-(t(xi2(u1,u2))%*%Lambda2(u1,u2)%*%xi2(u1,u2))
  return(AV)
}

AV2out(0.03,0.09)
x=seq(1/100,1,by=1/10^2)
y=x
z <- outer(x, y, AV2out)
brk<- quantile( c(sqrt(z)),na.rm=TRUE)
cols = rev(colorRampPalette(grey(seq(0.8,0.2,  len =length(brk)-1)), space = "rgb")(length(brk)-1))
image.plot(x,y,sqrt(z),breaks=brk, col=cols, legend.shrink=1,
           lab.breaks=c(NA,NA,NA,NA,80))



AV2min <- function(u){ #for optim/genSAN function
  u1 <- u[1]
  u2 <- u[2]
  return(sqrt(abs(t(xi2(u1,u2))%*%Lambda2(u1,u2)%*%xi2(u1,u2))))
}
optim(c(0.01,0.02),AV2min,method='Nelder-Mead')
GenSA(c(0.01,0.02),AV2min,lower=c(10^(-10),10^(-10)),upper=c(1,1),control = list(max.time=6, temperature=10))
u=c(0.01380356,0.01540287)


