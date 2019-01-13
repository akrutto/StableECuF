a=1.5
b=0.5
n=1000
set.seed(4)
X <- rstable(n,alpha=a,beta=b,param=0)
Re.cu <- function (u){
  uX <- (u)%*%X  
  return (log((1/length(X))*sqrt((rowSums(cos(uX)))^2+(rowSums(sin(uX)))^2)))
}
#Re(u1)
Im.cu <- function (u){
  uX <- (u)%*%X  
  N<-lentgh(X)
  return (atan2(rowSums(sin(uX))/N,rowSums(cos(uX))/N))
}

Re.ch <- function (u){
  uX <- (u)%*%X  
  return (((1/length(X))*(rowSums(cos(uX)))))
}
Im.ch <- function (u){
  uX <- (u)%*%X  
  return (((1/length(X))*(rowSums(sin(uX)))))
}

#alpha:
f.e<-function(x,y){log(Re.cu(x)/Re.cu(y))/log(x/y)} #alpha
N <- length(X)
Lambda2 <- function(u1,u2)
{N <- length(X)  
L <- matrix(NA,nrow=4,ncol=4)
L[1,1] <- 1+Re.ch(2*u1)-2*(Re.ch(u1))^2 #R1R1
L[1,2] <- Re.ch(u1-u2)+ Re.ch(u1+u2) -2*Re.ch(u1)*Re.ch(u2) #R1R2
L[1,3] <- Im.ch(u1+u1)-0-2*Re.ch(u1)*Im.ch(u1) #R1I1
L[1,4] <- Im.ch(u1+u2)-Im.ch(u1-u2)-2*Re.ch(u1)*Im.ch(u2) #R1I2
L[2,2] <- 1+Re.ch(2*u2)-2*(Re.ch(u2))^2 #R2R2
L[2,3] <- Im.ch(u2+u1)-Im.ch(u2-u1)-2*Re.ch(u2)*Im.ch(u1) ##R2I1
L[2,4] <- Im.ch(u2+u2)-0-2*Re.ch(u2)*Im.ch(u2) #R1I1
L[3,3] <- 1-Re.ch(2*u1)-2*(Im.ch(u1))^2
L[3,4] <- Re.ch(u1-u2)- Re.ch(u1+u2) -2*Im.ch(u1)*Im.ch(u2) #I1I2
L[4,4] <- 1-Re.ch(2*u2)-2*(Im.ch(u2))^2 #I2I2 1-R22-2*(b[4]^2) R22=Re.ch(2*u2,X)Im.ch(u2,X)
L[lower.tri(L)] <- t(L)[lower.tri(L)]
L <- L/(2) 
return((L))
}


#Xi
u1=0.2
u2=0.5
X<-X05
Lambda2(0,0)

xi2 <- function(u1,u2){
  xivec  <- matrix(NA,nrow=4,ncol=1)
  xivec[1,1]  <- Re.ch(u1)/(exp(2*Re.cu(u1))*Re.cu(u1)*log(u1/u2))
  xivec[2,1]  <- -Re.ch(u2)/(exp(2*Re.cu(u2))*Re.cu(u2)*log(u1/u2))
  xivec[3,1]   <- Im.ch(u1)/(exp(2*Re.cu(u1))*Re.cu(u1)*log(u1/u2))
  xivec[4,1] <-  -Im.ch(u2)/(exp(2*Re.cu(u2))*Re.cu(u2)*log(u1/u2))
  return(xivec)
}
Lambda2(0.2,0.5) 
xi2(0.2,0.5)
u1=0.2
u2=0.5
AV2out <- function(u1,u2){ #for outer function
  #u
  #u1 <- u[1]
  #u2 <- u[2]
  AV<-(t(xi2(u1,u2))%*%Lambda2(u1,u2)%*%xi2(u1,u2))
  return(AV)
}
AV2out(0.2,0.5)
AV2out(1,2)





