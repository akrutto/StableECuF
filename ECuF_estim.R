library(stable)
memory.size(max = NA)
memory.limit(size = NA)
ugamma<-function (x){xX <- X%*%x 
return (log((1/N)*sqrt((colSums(cos(xX)))^2+(colSums(sin(xX)))^2))+1)
}
u1f<-function (x){xX <- X%*%x 
return (log((1/N)*sqrt((colSums(cos(xX)))^2+(colSums(sin(xX)))^2))+0.1)
}
u2f<-function (x){xX <- X%*%x 
return (log((1/N)*sqrt((colSums(cos(xX)))^2+(colSums(sin(xX)))^2))+0.5)
}

X <- rstable(1000,alpha=1,beta=0,gamma=1,delta=0,param=0)
X=matrix(X,nrow=length(X),ncol=1)
N=length(X)
IQR <- as.numeric(quantile(X,0.75)-quantile(X,0.25))
IQR
uniroot(u1f, c(0,1/(IQR)))
u=c(0,0)
u[1]=uniroot(u1f, c(0,1))$root
u[2]=uniroot(u1f, c(0,1))$root
u
plot(u1f,0,1)

#ESTIMATION
alfa=0
beeta1=0
beeta2=0 
gamma_std=0
gamma_data=0
delta1_std=0
delta2_std=0
delta1_data=0
delta2_data=0
hbeeta=0
hdelta=0
hatg=0
hatd=0
Y=0

#General info
Nn=c(300,500,1000,2000,2500,4000,5000,10000)#sample sizec(200,500,1000,2000,2500,4000,5000,10000)
K=100 #replications
g=1 #gamma
d=0 #delta
Bb=c(0,0.5,1)#betac(0,0.5,1)
avec=c(0.25,0.35,0.5,0.8,1,1.2,1.5,1.8,2)#alphac(0.5,1,1.5)
avec
E=matrix(NA,ncol=29,nrow=length(avec)*K*length(Nn)*length(Bb)) #tulemused
colnames(E)<-c("beta","alpha",'hatg','hatd',
               'e_halpha', 'e_hbeta', "e_hgamma",'e_hdelta',
               'N','u1','u2',
'MLe_halpha', 'MLe_hbeta', "MLe_hgamma",'MLe_hdelta',
'QBe_halpha', 'QBe_hbeta', "QBe_hgamma",'QBe_hdelta',
'Che_halpha', 'Che_hbeta', "Che_hgamma",'Che_hdelta',
'ECuF', 'ML', "QB",'EChF', 'Reu1','Reu2')
rownames(E) <- rownames(E, do.NULL = FALSE,prefix = "Obs.")#
head(E)
dim(E)
U=1 #E matrix K times counting for each avec
u=c(1,1)

#ALGUS 
for (B in 1:length(Bb)){
  b=Bb[B]
  for(n in 1:length(Nn)){
    N=Nn[n]
    for(l in 1:length(avec)){#for each alpha at fixed beta=0,1
      for (j in 1:K){ #replicates
        XX<-rstable(N,alpha=avec[l],beta=b, delta=d, gamma=g,param=0) #simulations rstable1 for b=1, rstable for b=0
        X <- matrix(XX,nrow=length(XX),ncol=1)
        TrRe <- as.numeric(quantile(XX,0.9)-quantile(XX,0.1))
                hatg<-1/uniroot(ugamma, c(0,TrRe))$root #for standard case with n=5000
        k <- sum(cos((1/hatg)*Y))	
        s <- sum(sin((1/hatg)*Y))    
        hatd <- hatg*atan2(s,k)
        ptm<-proc.time()
        
        #u=GenSA(c(0.03,0.09),AV,                                    #For Asymptotic Variance
         #       lower=c(10^(-10),10^(-10)),upper=c(1,1),
          #      control = list(max.time=60, temperature=15))$par
       u=c(0,0)
        u[1]=uniroot(u1f, c(10^-8,3*hatg))$root #for standard case with n=5000
        u[2]=uniroot(u2f, c(10^-8,3*hatg))$root #for standard case with n=5000
        r=0
        v=0
        k=0
        s=0
       # r=c(-0.1,-0.5) by assumption but we recalculate from data
        Y<-XX
        for (i in 1:2){# empirical  cumulant function real and imaginary part at u1, u2
          k <- sum(cos(u[i]*Y))	
          s <- sum(sin(u[i]*Y))    
          r[i] <- -log((1/N)*sqrt(k^2+s^2))
          v[i] <-  atan2(s,k)
          #lõpp  real and imaginary part of empirical cumulant function
        }
        ##cumulant estimates for 1.,..,K replicates
        alfa[j] <- log(r[1]/r[2])/log(u[1]/u[2])
        alfa[j] <-min(max(alfa[j],0.1),2) 
        gamma_std[j] <- exp((log(u[1])*log(r[2])- log(u[2])*log(r[1]))/log(r[1]/r[2])) #standardized
        beeta1[j] <- # alpha is not 1
          ((u[2]*v[1]-u[1]*v[2])/((gamma_std[j]^alfa[j])*(tan(pi*alfa[j]/2))*(-u[1]*(u[2]^alfa[j])+u[2]*(u[1]^alfa[j])))) 
        beeta2[j] <- #alpha is 1
          ((v[1]*u[2]-v[2]*u[1])/( (2*gamma_std[j]/pi)*u[1]*u[2]*(log(u[2])-log(u[1])) )) 
        hbeeta[j]<-beeta1[j]
        if (abs(alfa[j]-1)<0.01) {hbeeta[j]<-beeta2[j]} 
        hbeeta[j]<-min(max(hbeeta[j],-1),1) 
        c1 <- (u[1]^alfa[j])*((u[1]*gamma_std[j])^(1-alfa[j])-1) #constant for delta_std alpha NEQ 1 calculations 
        c2 <- (u[2]^alfa[j])*((u[2]*gamma_std[j])^(1-alfa[j])-1) #constant for delta_std alpha NEQ 1 calculations 
        delta1_std[j] <- # alpha is not 1
          ( v[1]*c2-v[2]*c1) / ( u[2]*(u[1]^alfa[j])-u[1]*(u[2]^alfa[j]) ) 
        c1 <- u[1]*log(gamma_std[j]*u[1])#constant for delta_std alpha=1 calculations 
        c2 <- u[2]*log(gamma_std[j]*u[2])#constant for delta_std alpha=1 calculations 
        delta2_std[j] <- #alpha is 1
          ( v[1]*c2-v[2]*c1) / ( u[1]*u[2]*(log(u[2])-log(u[1])) ) 
                hdelta[j]<-delta1_std[j] 
        if (abs(alfa[j]-1)<0.01) {hdelta[j]<-delta2_std[j]}
        E[U+j-1,24]<-(proc.time() - ptm)[1]
        E[U+j-1,28]<- r[1]
        E[U+j-1,29]<- r[2]
        E[U+j-1,1]<-b
        E[U+j-1,2]<-avec[l]
        E[U+j-1,3]<-hatg-g#errors
        E[U+j-1,4]<-hatd-d
        E[U+j-1,5]<- alfa[j]-avec[l] 
        E[U+j-1,6]<- hbeeta[j]-b
        E[U+j-1,7]<- gamma_std[j]-g
        E[U+j-1,8]<- hdelta[j]-d
        E[U+j-1,9]<- N
        E[U+j-1,10]<- u[1]
        E[U+j-1,11]<- u[2]
        ptm <- proc.time()
  ML<-stable.fit.mle(X,param=0)
        E[U+j-1,25]<-(proc.time() - ptm)[1]
        E[U+j-1,12]<-(ML[1]-avec[l])#errors
        E[U+j-1,13]<-(ML[2]-b)
        E[U+j-1,14]<-(ML[3]-g)
        E[U+j-1,15]<-(ML[4]-d)
        ptm <- proc.time()
  QB<-stable.fit.quantile(X,param=0)
        E[U+j-1,26]<-(proc.time() - ptm)[1]
        E[U+j-1,16]<-(QB[1]-avec[l])#errors
        E[U+j-1,17]<-(QB[2]-b)
        E[U+j-1,18]<-(QB[3]-g)
        E[U+j-1,19]<-(QB[4]-d)
        ptm <- proc.time()
   ChF<-stable.fit.ecf(X,gamma0=(E[U+j-1,18]+g),delta0=E[U+j-1,19],param=0)
        E[U+j-1,27]<-(proc.time() - ptm)[1]
        E[U+j-1,20]<- (ChF[1]-avec[l])#errors
        E[U+j-1,21]<- (ChF[2]-b)
        E[U+j-1,22]<- (ChF[3]-g)
        E[U+j-1,23]<- (ChF[4]-d)
     #lõpp K
      }
      U=U+K
      #lõpp avec N
    }}} 
#LÕPP

head(E) #results
dim(E)
write.csv(E,"ECuF_15.11.18.csv")
E0=E[E[,1]==0,]#beta=0
E1=E[E[,1]==1,]#beta=1
E05=E[E[,1]==0.5,]#beta=0.5
#speed N=1000, K=20
#seewave
sECuF<-tapply(E[,24],list(E[,1],E[,2],E[,9]),rms)#beta=0
sMLE<-tapply(E[,25],list(E[,1],E[,2],E[,9]),rms)
sQB<-tapply(E[,26],list(E[,1],E[,2],E[,9]),rms)
sEChF<-tapply(E[,27],list(E[,1],E[,2],E[,9]),rms)
sECuF
sMLE
sQB
sEChF
#alpha estimates
CuF2_a0<-tapply(E0[,5],list(E0[,2],E0[,9]),rms)
ChF_a0<-tapply(E0[,20],list(E0[,2],E0[,9]),rms)
ML_a0<-tapply(E0[,12],list(E0[,2],E0[,9]),rms)
QB_a0<-tapply(E0[,16],list(E0[,2],E0[,9]),rms)
CuF1_a0_u1<-tapply(E0[,10],list(E0[,2],E0[,9]),rms)
CuF1_a0_u2<-tapply(E0[,11],list(E0[,2],E0[,9]),rms)
#persp(CuF2_a0)
ChF_a0
ML_a0
QB_a0
CuF2_a1<-tapply(E1[,5],list(E1[,2],E1[,9]),rms)
ChF_a1<-tapply(E1[,20],list(E1[,2],E1[,9]),rms)
ML_a1<-tapply(E1[,12],list(E1[,2],E1[,9]),rms)
QB_a1<-tapply(E1[,16],list(E1[,2],E1[,9]),rms)
CuF1_a1_u1<-tapply(E1[,10],list(E1[,2],E0[,9]),rms)
CuF1_a1_u2<-tapply(E1[,11],list(E1[,2],E0[,9]),rms)
CuF2_a1
ChF_a1
ML_a1
QB_a1
CuF2_a05<-tapply(E05[,5],list(E05[,2],E05[,9]),rms)
ChF_a05<-tapply(E05[,20],list(E05[,2],E05[,9]),rms)
ML_a05<-tapply(E05[,12],list(E05[,2],E05[,9]),rms)
QB_a05<-tapply(E05[,16],list(E05[,2],E05[,9]),rms)
CuF2_a05
ChF_a05
ML_a05
QB_a05

CuF2_b0<-tapply(E0[,6],list(E0[,2],E0[,9]),rms)
ChF_b0<-tapply(E0[,21],list(E0[,2],E0[,9]),rms)
ML_b0<-tapply(E0[,13],list(E0[,2],E0[,9]),rms)
QB_b0<-tapply(E0[,17],list(E0[,2],E0[,9]),rms)
CuF1_b0_1<-tapply(E0[,10],list(E0[,2],E0[,9]),rms)
CuF1_b0_2<-tapply(E0[,11],list(E0[,2],E0[,9]),rms)
CuF2_b0
ChF_b0
ML_b0
QB_b0
CuF2_b05<-tapply(E05[,6],list(E05[,2],E05[,9]),rms)
ChF_b05<-tapply(E05[,21],list(E05[,2],E05[,9]),rms)
ML_b05<-tapply(E05[,13],list(E05[,2],E05[,9]),rms)
QB_b05<-tapply(E05[,17],list(E05[,2],E05[,9]),rms)
CuF2_b05
ChF_b05
ML_b05
QB_b05
CuF2_b1<-tapply(E1[,6],list(E1[,2],E1[,9]),rms)
ChF_b1<-tapply(E1[,21],list(E1[,2],E1[,9]),rms)
ML_b1<-tapply(E1[,13],list(E1[,2],E1[,9]),rms)
QB_b1<-tapply(E1[,17],list(E1[,2],E1[,9]),rms)
CuF1_b1_1<-tapply(E1[,30],list(E1[,2],E1[,9]),rms)
CuF1_b1_2<-tapply(E1[,31],list(E1[,2],E1[,9]),rms)
CuF1_b1_3<-tapply(E1[,32],list(E1[,2],E1[,9]),rms)
CuF2_b1
ChF_b1
ML_b1
QB_b1

CuF2_g05<-tapply(E05[,7],list(E05[,2],E05[,9]),rms)
ChF_g05<-tapply(E05[,22],list(E05[,2],E05[,9]),rms)
ML_g05<-tapply(E05[,14],list(E05[,2],E05[,9]),rms)
QB_g05<-tapply(E05[,18],list(E05[,2],E05[,9]),rms)
CuF2_g05
ChF_g05
ML_g05
QB_g05


CuF2_g0<-tapply(E0[,7],list(E0[,2],E0[,9]),rms)
ChF_g0<-tapply(E0[,22],list(E0[,2],E0[,9]),rms)
ML_g0<-tapply(E0[,14],list(E0[,2],E0[,9]),rms)
QB_g0<-tapply(E0[,18],list(E0[,2],E0[,9]),rms)
CuF2_g0
ChF_g0
ML_g0
QB_g0
CuF2_g1<-tapply(E1[,7],list(E1[,2],E1[,9]),rms)
ChF_g1<-tapply(E1[,22],list(E1[,2],E1[,9]),rms)
ML_g1<-tapply(E1[,14],list(E1[,2],E1[,9]),rms)
QB_g1<-tapply(E1[,18],list(E1[,2],E1[,9]),rms)
CuF2_g1
ChF_g1
ML_g1
QB_g1

CuF2_d05<-tapply(E05[,8],list(E05[,2],E05[,9]),rms)
ChF_d05<-tapply(E05[,23],list(E05[,2],E05[,9]),rms)
ML_d05<-tapply(E05[,15],list(E05[,2],E05[,9]),rms)
QB_d05<-tapply(E05[,19],list(E05[,2],E05[,9]),rms)
CuF2_d05
ChF_d05
ML_d05
QB_d05
CuF2_d0<-tapply(E0[,8],list(E0[,2],E0[,9]),rms)
ChF_d0<-tapply(E0[,23],list(E0[,2],E0[,9]),rms)
ML_d0<-tapply(E0[,15],list(E0[,2],E0[,9]),rms)
QB_d0<-tapply(E0[,19],list(E0[,2],E0[,9]),rms)
CuF2_d0
ChF_d0
ML_d0
QB_d0
CuF2_d1<-tapply(E1[,8],list(E1[,2],E1[,9]),rms)
ChF_d1<-tapply(E1[,23],list(E1[,2],E1[,9]),rms)
ML_d1<-tapply(E1[,15],list(E1[,2],E1[,9]),rms)
QB_d1<-tapply(E1[,19],list(E1[,2],E1[,9]),rms)
CuF2_d1
ChF_d1
ML_d1
QB_d1
c1=2 #Cu
c2=1 #ML
c3=1 #Ch
c4=1 #QB
l1=1 #Cu
l2=2 #ML
l3=3 #Ch
l4=4 #QB
size=1 #choice of n size in graphs
#ALgus alpha ja beeta ja gamma ja delta üle n=5000, b=0.5
op <- par(mfrow = c(2,2),mar=c(4,4,1,1),cex.sub=1.2,mgp=c(2.5,0.7,0),oma = c(2, 1, 1, 1),cex=1.1)
plot(avec,CuF2_a05[,size] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(alpha),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.25))#,xlab=bquote(paste(alpha)),xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(0.2,0.5,1,1.5,1.8,2),labels=c(0.2,0.5,1,1.5,1.8,2))
axis(2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5,1,1.1),labels=c(0,0.1,0.2,0.3,0.4,0.5,1,1.1))#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(avec,ChF_a05[,size],type='l',lty=l3,col=c3)
lines(avec,ML_a05[,size],type='l',lty=l2,col=c2)
lines(avec,QB_a05[,size],type='l',lty=l4,col=c4)

#legend('top',legend=c('n=5000','b=0.5'))
plot(avec,CuF2_b05[,size] ,type='l',col=2,axes=F,lty=1,frame=T,
     ylab=bquote(paste("RMSE(",hat(beta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.25))#,at ",beta,"=0.5" ),xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(0.2,0.5,0.8,1,1.2,1.5,1.8,2),labels=c(0.2,0.5,0.8,1,1.2,1.5,1.8,2))
axis(2,las=1,at=c(0,0.5,1,1.5),labels=c(0,0.5,1,1.5))#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
#axis(1,tick=T,las=1,at=c(0.2,0.5,0.8,1,1.2,1.5,1.8,2),labels=c(0.2,0.5,0.8,1,1.2,1.5,1.8,2))
#axis(2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5,1,1.1),labels=c(0,0.1,0.2,0.3,0.4,0.5,1,1.1))#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(avec,ChF_b05[,size],type='l',lty=l3,col=c3)
lines(avec,ML_b05[,size],type='l',lty=l2,col=c2)
lines(avec,QB_b05[,size],type='l',lty=l4,col=c4)
#legend('top',legend=c('n=5000','b=0.5'))
plot(avec,CuF2_g05[,size] ,type='l',col=2,axes=F,lty=1,frame=T,
     ylab=bquote(paste("RMSE(",hat(gamma),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.25))#,at ",beta,"=0.5" ),xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(0.2,0.5,0.8,1,1.2,1.5,1.8,2),labels=c(0.2,0.5,0.8,1,1.2,1.5,1.8,2))
axis(2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5,1,1.1),labels=c(0,0.1,0.2,0.3,0.4,0.5,1,1.1))#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(avec,ChF_g05[,size],type='l',lty=l3,col=c3)
lines(avec,ML_g05[,size],type='l',lty=l2,col=c2)
lines(avec,QB_g05[,size],type='l',lty=l4,col=c4)
#lines(avec,CuF1_g05[,1],type='l',lty=1,col=4)
#legend('top',legend=c('n=5000','b=0.5'))
plot(avec,CuF2_d05[,size] ,type='l',col=2,axes=F,lty=1,frame=T,
     ylab=bquote(paste("RMSE(",hat(delta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.25))#,at ",beta,"=0.5" ),xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(0.2,0.5,0.8,1,1.2,1.5,1.8,2),labels=c(0.2,0.5,0.8,1,1.2,1.5,1.8,2))
axis(2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5,1,1.1),labels=c(0,0.1,0.2,0.3,0.4,0.5,1,1.1))#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(avec,ChF_d05[,size],type='l',lty=l3,col=c3)
lines(avec,ML_d05[,size],type='l',lty=l2,col=c2)
lines(avec,QB_d05[,size],type='l',lty=l4,col=c4)
#lines(avec,CuF1_d05[,1],type='l',lty=1,col=4)
#legend('top',legend=c('n=5000','b=0.5'))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", legend= c('ECuF','EChF','ML','QB'),
     xpd = TRUE, horiz = T, inset = c(0, 0), bty = "n", text.width=0.15,  lty=c(l1,l3,l2,l4),col=c(c1,c3,c2,c4),cex = 1)#

#legend("bottom", legend= c('ECuF2','ECuF1','EChF','ML','QB'),
 #      xpd = TRUE, horiz = T, inset = c(0, 0), bty = "n", text.width=0.15,  lty=c(l1,1,l3,l2,l4),col=c(c1,4,c3,c2,c4),cex = 1)#




#ALGUS alpha üle Nn b=0,1
op <- par(mfcol = c(3,2),mar=c(4,4,1,1),cex=1.1,mgp=c(2.5,0.7,0),oma = c(2, 1, 1, 1))
plot(Nn,CuF2_a1[1,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(alpha),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_a1[1,],type='l',lty=l3,col=c3)
lines(Nn,ML_a1[1,],type='l',lty=l2,col=c2)
lines(Nn,QB_a1[1,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_a1_u1[1,],type='l',lty=l4,col=3)
lines(Nn,CuF1_a1_u2[1,],type='l',lty=l4,col=3)
legend('topright',legend=c(expression(paste(alpha, " =  0.3")),expression(paste(beta, " = 1"))),bty='n',horiz=T,cex=1)
plot(Nn,CuF2_a1[2,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(alpha),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_a1[2,],type='l',lty=l3,col=c3)
lines(Nn,ML_a1[2,],type='l',lty=l2,col=c2)
lines(Nn,QB_a1[2,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_a1_u1[2,],type='l',lty=l4,col=3)
lines(Nn,CuF1_a1_u2[2,],type='l',lty=l4,col=3)
legend('topright',legend=c(expression(paste(alpha, " =  1.0")),expression(paste(beta, " =  1"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_a1[3,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(alpha),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_a1[3,],type='l',lty=l3,col=c3)
lines(Nn,ML_a1[3,],type='l',lty=l2,col=c2)
lines(Nn,QB_a1[3,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_a1_u1[3,],type='l',lty=l4,col=3)
lines(Nn,CuF1_a1_u2[3,],type='l',lty=l4,col=3)

legend('topright',legend=c(expression(paste(alpha, " =  1.8")),expression(paste(beta, " =  1"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_a0[1,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(alpha),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_a0[1,],type='l',lty=l3,col=c3)
lines(Nn,ML_a0[1,],type='l',lty=l2,col=c2)
lines(Nn,QB_a0[1,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_a0_u1[1,],type='l',lty=l4,col=3)
lines(Nn,CuF1_a0_u2[1,],type='l',lty=l4,col=3)
legend('topright',legend=c(expression(paste(alpha, " =  0.3")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_a0[2,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(alpha),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_a0[2,],type='l',lty=l3,col=c3)
lines(Nn,ML_a0[2,],type='l',lty=l2,col=c2)
lines(Nn,QB_a0[2,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_a0_u1[2,],type='l',lty=l4,col=3)
lines(Nn,CuF1_a0_u2[2,],type='l',lty=l4,col=3)
legend('topright',legend=c(expression(paste(alpha, " =  1.0")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_a0[3,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(alpha),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_a0[3,],type='l',lty=l3,col=c3)
lines(Nn,ML_a0[3,],type='l',lty=l2,col=c2)
lines(Nn,QB_a0[3,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_a0_u1[3,],type='l',lty=l4,col=3)
lines(Nn,CuF1_a0_u2[3,],type='l',lty=l4,col=3)
legend('topright',legend=c(expression(paste(alpha, " = 1.8")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
#LÕPP Nn
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend= c('ECuF2','EChF','ML','QB'),
       xpd = TRUE, horiz = T, inset = c(0, 0), bty = "n", text.width=0.15,  lty=c(l1,l3,l2,l4),col=c(c1,c3,c2,c4),cex = 1)#
#,legend("bottom",  
#
#ALGUS beeta
op <- par(mfcol = c(3,2),mar=c(4,4,1,1),cex=1.1,mgp=c(2.5,0.7,0),oma = c(2, 1, 1, 1))
plot(Nn,CuF2_b1[1,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(beta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_b1[1,],type='l',lty=l3,col=c3)
lines(Nn,ML_b1[1,],type='l',lty=l2,col=c2)
lines(Nn,QB_b1[1,],type='l',lty=l4,c4)
legend('topright',legend=c(expression(paste(alpha, " =  0.3")),expression(paste(beta, " = 1"))),bty='n',horiz=T,cex=1)
plot(Nn,CuF2_b1[2,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(beta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_b1[2,],type='l',lty=l3,col=c3)
lines(Nn,ML_b1[2,],type='l',lty=l2,col=c2)
lines(Nn,QB_b1[2,],type='l',lty=l4,col=c4)
legend('topright',legend=c(expression(paste(alpha, " =  1.0")),expression(paste(beta, " =  1"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_b1[3,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(alpha),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_b1[3,],type='l',lty=l3,col=c3)
lines(Nn,ML_b1[3,],type='l',lty=l2,col=c2)
lines(Nn,QB_b1[3,],type='l',lty=l4,col=c4)
legend('topright',legend=c(expression(paste(alpha, " =  1.8")),expression(paste(beta, " =  1"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_b0[1,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(beta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_b0[1,],type='l',lty=l3,col=c3)
lines(Nn,ML_b0[1,],type='l',lty=l2,col=c2)
lines(Nn,QB_b0[1,],type='l',lty=l4,col=c4)
legend('topright',legend=c(expression(paste(alpha, " =  0.3")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_b0[2,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(beta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_b0[2,],type='l',lty=l3,col=c3)
lines(Nn,ML_b0[2,],type='l',lty=l2,col=c2)
lines(Nn,QB_b0[2,],type='l',lty=l4,col=c4)
legend('topright',legend=c(expression(paste(alpha, " =  1.0")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_b0[3,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(beta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_b0[3,],type='l',lty=l3,col=c3)
lines(Nn,ML_b0[3,],type='l',lty=l2,col=c2)
lines(Nn,QB_b0[3,],type='l',lty=l4,col=c4)
legend('topright',legend=c(expression(paste(alpha, " = 1.8")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
#LÕPP Nn
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend= c('ECuF1','EChF','ML','QB'),
       xpd = TRUE, horiz = T, inset = c(0, 0), bty = "n", text.width=0.15,  lty=c(l1,l3,l2,l4),col=c(c1,c3,c2,c4),cex = 1)#
#,legend("bottom",  
#

#
#ALGUS gamma
op <- par(mfcol = c(3,2),mar=c(4,4,1,1),cex=1.1,mgp=c(2.5,0.7,0),oma = c(2, 1, 1, 1))
plot(Nn,CuF2_g1[1,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(gamma),")")),pch=19,cex=0.7,xlab='',ylim=c(0,2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_g1[1,],type='l',lty=l3,col=c3)
lines(Nn,ML_g1[1,],type='l',lty=l2,col=c2)
lines(Nn,QB_g1[1,],type='l',lty=l4,c4)
lines(Nn,CuF1_g1[1,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  0.3")),expression(paste(beta, " = 1"))),bty='n',horiz=T,cex=1)
plot(Nn,CuF2_g1[2,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(gamma),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.5))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_g1[2,],type='l',lty=l3,col=c3)
lines(Nn,ML_g1[2,],type='l',lty=l2,col=c2)
lines(Nn,QB_g1[2,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_g1[2,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  1.0")),expression(paste(beta, " =  1"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_g1[3,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(gamma),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_g1[3,],type='l',lty=l3,col=c3)
lines(Nn,ML_g1[3,],type='l',lty=l2,col=c2)
lines(Nn,QB_g1[3,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_g1[3,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  1.8")),expression(paste(beta, " =  1"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_g0[1,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(gamma),")")),pch=19,cex=0.7,xlab='',ylim=c(0,2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_g0[1,],type='l',lty=l3,col=c3)
lines(Nn,ML_g0[1,],type='l',lty=l2,col=c2)
lines(Nn,QB_g0[1,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_g0[1,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  0.3")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_g0[2,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(gamma),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.5))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_g0[2,],type='l',lty=l3,col=c3)
lines(Nn,ML_g0[2,],type='l',lty=l2,col=c2)
lines(Nn,QB_g0[2,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_g0[2,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  1.0")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_g0[3,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(gamma),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_g0[3,],type='l',lty=l3,col=c3)
lines(Nn,ML_g0[3,],type='l',lty=l2,col=c2)
lines(Nn,QB_g0[3,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_g0[3,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " = 1.8")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
#LÕPP Nn
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend= c('ECuF2','ECuF1','EChF','ML','QB'),
       xpd = TRUE, horiz = T, inset = c(0, 0), bty = "n", text.width=0.15,  lty=c(l1,1,l3,l2,l4),col=c(c1,4,c3,c2,c4),cex = 1)#
#,legend("bottom",  
#

#
#ALGUS delta
op <- par(mfcol = c(3,2),mar=c(4,4,1,1),cex=1.1,mgp=c(2.5,0.7,0),oma = c(2, 1, 1, 1))
plot(Nn,CuF2_d1[1,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(delta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_d1[1,],type='l',lty=l3,col=c3)
lines(Nn,ML_d1[1,],type='l',lty=l2,col=c2)
lines(Nn,QB_d1[1,],type='l',lty=l4,c4)
lines(Nn,CuF1_d1[1,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  0.3")),expression(paste(beta, " = 1"))),bty='n',horiz=T,cex=1)
plot(Nn,CuF2_d1[2,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(delta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.6))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_d1[2,],type='l',lty=l3,col=c3)
lines(Nn,ML_d1[2,],type='l',lty=l2,col=c2)
lines(Nn,QB_d1[2,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_d1[2,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  1.0")),expression(paste(beta, " =  1"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_d1[3,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(delta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_d1[3,],type='l',lty=l3,col=c3)
lines(Nn,ML_d1[3,],type='l',lty=l2,col=c2)
lines(Nn,QB_d1[3,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_d1[3,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  1.8")),expression(paste(beta, " =  1"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_d0[1,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(delta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_d0[1,],type='l',lty=l3,col=c3)
lines(Nn,ML_d0[1,],type='l',lty=l2,col=c2)
lines(Nn,QB_d0[1,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_d0[1,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  0.3")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_d0[2,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(delta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.6))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_d0[2,],type='l',lty=l3,col=c3)
lines(Nn,ML_d0[2,],type='l',lty=l2,col=c2)
lines(Nn,QB_d0[2,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_d0[2,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " =  1.0")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
plot(Nn,CuF2_d0[3,] ,type='l',col=c1,axes=F,lty=l1,frame=T,
     ylab=bquote(paste("RMSE(",hat(delta),")")),pch=19,cex=0.7,xlab='',ylim=c(0,0.2))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=1,at=c(200,1000,3000,5000),labels=c(200,1000,3000,5000))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(Nn,ChF_d0[3,],type='l',lty=l3,col=c3)
lines(Nn,ML_d0[3,],type='l',lty=l2,col=c2)
lines(Nn,QB_d0[3,],type='l',lty=l4,col=c4)
lines(Nn,CuF1_d0[3,],type='l',lty=1,col=4)
legend('topright',legend=c(expression(paste(alpha, " = 1.8")),expression(paste(beta, " =  0"))),bty='n',horiz=T,cex = 1)
#LÕPP Nn
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend= c('ECuF2','ECuF1','EChF','ML','QB'),
       xpd = TRUE, horiz = T, inset = c(0, 0), bty = "n", text.width=0.15,  lty=c(l1,1,l3,l2,l4),col=c(c1,4,c3,c2,c4),cex = 1)#
#,legend("bottom",  
#


#Results for 1/gamma at -1 and delta
#dim(E)
#E[nrow(E),]
#write.csv(E,"Asym.comp12.02.csv")
#SimFile<-read.csv("simulatsioonid25.09.17.csv",header=T)

E0=E[E[,1]==0,]#beta=0
E1=E[E[,1]==1,]#beta=1
E05=E[E[,1]==0.5,]#beta=0.5
library(seewave) #for rms
ECuF1_d0<-tapply(E0[,4],list(E0[,2],E0[,9]),rms)
ECuF1_g0<-tapply(E0[,3],list(E0[,2],E0[,9]),rms)
ECuF1_d1<-tapply(E1[,4],list(E1[,2],E1[,9]),rms)
ECuF1_g1<-tapply(E1[,3],list(E1[,2],E1[,9]),rms)
ECuF1_d05<-tapply(E05[,4],list(E05[,2],E05[,9]),rms)
ECuF1_g05<-tapply(E05[,3],list(E05[,2],E05[,9]),rms)
#speed
ECuF0<-tapply(E0[,24],list(E0[,2],E0[,9]),rms)
MLE0<-tapply(E0[,25],list(E0[,2],E0[,9]),rms)
QB0<-tapply(E0[,26],list(E0[,2],E0[,9]),rms)
EChF0<-tapply(E0[,27],list(E0[,2],E0[,9]),rms)
ECuF1<-tapply(E1[,24],list(E0[,2],E0[,9]),rms)
MLE1<-tapply(E1[,25],list(E0[,2],E0[,9]),rms)
QB1<-tapply(E1[,26],list(E0[,2],E0[,9]),rms)
EChF1<-tapply(E1[,27],list(E0[,2],E0[,9]),rms)
ECuF05<-tapply(E05[,24],list(E0[,2],E0[,9]),rms)
MLE05<-tapply(E05[,25],list(E0[,2],E0[,9]),rms)
QB05<-tapply(E05[,26],list(E0[,2],E0[,9]),rms)
EChF05<-tapply(E05[,27],list(E0[,2],E0[,9]),rms)


op <- par(mfcol = c(2,2),mar=c(4,4,1,1),cex.sub=1.2,mgp=c(2.5,0.7,0),oma = c(2, 1, 1, 1))
plot(avec,ECuF1_g0[,1] ,type='l',col=1,axes=F,lty=1,
     ylab=bquote(paste("RMSE(",hat(gamma),")")),pch=19,cex=0.7,xlab=bquote(paste(alpha)),ylim=c(0,0.4))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=2,at=avec,labels=avec)
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
lines(avec,ECuF1_g0[,2],type='l')#500
lines(avec,ECuF1_g0[,3],type='l')
lines(avec,ECuF1_g0[,4],type='l')
lines(avec,ECuF1_g0[,5],type='l')
lines(avec,ECuF1_g0[,6],type='l')
lines(avec,ECuF1_g1[,1],type='l',lty=2)#200
lines(avec,ECuF1_g1[,2],type='l',lty=2)
lines(avec,ECuF1_g1[,3],type='l',lty=2)
lines(avec,ECuF1_g1[,4],type='l',lty=2)
lines(avec,ECuF1_g1[,5],type='l',lty=2)
lines(avec,ECuF1_g1[,6],type='l')
plot(avec,ECuF1_d0[,1] ,type='l',col=1,axes=F,lty=1,
     ylab=bquote(paste("RMSE(",hat(delta),")")),pch=19,cex=0.7,xlab=bquote(paste(alpha)),ylim=c(0,0.4))#,,xlim=c(0,2),,#,ylim=c(0,0.4)
axis(1,tick=T,las=2,at=c(0.1,0.2,0.3,0.4,0.5,0.8,0.9,1.0,1.1,1.2,1.5,1.8,2),labels=c(0.1,0.2,0.3,0.4,0.5,0.8,0.9,1.0,1.1,1.2,1.5,1.8,2))
axis(2,las=1)#at=c(0,0.1, 0.2, 0.3,0.4, 0.5)
#lines(avec,ECuF1_d0[,1],type='l')
lines(avec,ECuF1_d0[,2],type='l')
lines(avec,ECuF1_d0[,3],type='l')
lines(avec,ECuF1_d0[,4],type='l')
lines(avec,ECuF1_d0[,5],type='l')
lines(avec,ECuF1_d0[,6],type='l')
#lines(avec,ECuF1_d0[,7],type='l')
#lines(avec,ECuF1_d0[,8],type='l')
#lines(avec,ECuF1_d0[,9],type='l')
#lines(avec,ECuF1_d0[,10],type='l')
lines(avec,ECuF1_d1[,1],type='l',lty=2)
lines(avec,ECuF1_d1[,2],type='l',lty=2)
lines(avec,ECuF1_d1[,3],type='l',lty=2)
lines(avec,ECuF1_d1[,4],type='l',lty=2)
lines(avec,ECuF1_d1[,5],type='l',lty=2)
lines(avec,ECuF1_d1[,6],type='l',lty=2)
#lines(avec,ECuF1_d1[,7],type='l',lty=2)
#lines(avec,ECuF1_d1[,8],type='l',lty=2)
#lines(avec,ECuF1_d1[,9],type='l',lty=2)
#lines(avec,ECuF1_d1[,10],type='l',lty=2)
#lines(avec,ECuF1_d1[,11],type='l',lty=2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend = c(expression(paste(beta, " = ", 0)),
                            expression(paste(beta, " = ", 1))), 
       xpd = TRUE, horiz = T, inset = c(0, 0), bty = "n", text.width=0.15,  lty=c(1,2),col=c(1,1),cex = 1)#
#,

#Code with  u calculated and which.min
hatg=0
hatd=0
alfa=0
beeta1=0
beeta2=0 #just for the control
gamma_std=0
gamma_data=0
delta1_std=0
delta2_std=0
delta1_data=0
delta2_data=0
hbeeta=0
hdelta=0
#General info
Nn=c(5000)#sample sizec(500,1000,2000,2500,4000,5000)
K=100 #replications
g=1 #gamma
d=0 #delta
Bb=c(0.5)#beta
avec=c(0.2,0.3,0.5,0.8,1.2,1.5,1.8,2)
avec
E=matrix(NA,ncol=27,nrow=length(avec)*K*length(Nn)*length(Bb)) #tulemused
#colnames(E)<-c("beta","alpha",'hatg','hatd','e_halpha', 'e-hbeta', "e_hgamma",'e_hdelta','N','u1','u2')
#rownames(E) <- rownames(E, do.NULL = FALSE,prefix = "Obs.")#
head(E)
dim(E)
r1=-0.1
r2=-0.5
#simuleerime gamma=1 ümbruses, seega u võib olla 1 ümbruses, (r=-1), reaalsuses läbib 1/g
ustand <- c(seq(0.5, 1.5, by=1/10^2))#,seq(1/10^2,0.1,by=1/10^3) 
#uvalik r=-0.01 ja r=-0.05 jaoks vajame  0 lähedal alpha=0.2 jaoks
uvalik <- c(seq(0, 1/10^8, by=1/10^10),seq(1/10^8, 1/10^6,by=1/10^8),seq(1/10^6,1/10^4,by=1/10^6),seq(1/10^4,3/10^3,by=1/10^4) ,seq(3/10^3,10,by=1/10^2))#,seq(1/10^2,0.1,by=1/10^3) 
#otsime null lähedal, seega alates väga väikesest
U=1 #E matrix K times counting for each avec
#ALGUS 
for (B in 1:length(Bb)){
  b=Bb[B]
  for(n in 1:length(Nn)){
    N=Nn[n]
    for(l in 1:length(avec)){#for each alpha at fixed beta=0,1
      for (j in 1:K){ #replicates
        X<-rstable(N,alpha=avec[l],beta=b, delta=d, gamma=g,param=0) #simulations rstable1 for b=1, rstable for b=0
        ptm<-proc.time()
        Y<-X
        k=0
        s=0
        r=0
        v=0
        for (i in 1:length(ustand)){# empirical  cumulant function real and imaginary part for standardizing
          k <- sum(cos(ustand[i]*Y))	
          s <- sum(sin(ustand[i]*Y))    
          r[i] <- log((1/N)*sqrt(k^2+s^2))
          #v[i] <-  atan2(s,k)
          #lõpp ustand
        }
        hatg<-1/ustand[(which.min(abs(r+1)))]
        k <- sum(cos((1/hatg)*Y))	
        s <- sum(sin((1/hatg)*Y))    
        hatd <- hatg*atan2(s,k)
        #standardiseerime
        if (hatg>1){Y<-(Y-hatd)/hatg} 
        r=0
        v=0
        k=0
        s=0
        for (i in 1:length(uvalik)){# empirical cumulant function real and imaginary part for standardized data
          k <- sum(cos(uvalik[i]*Y))	
          s <- sum(sin(uvalik[i]*Y))    
          r[i] <- log((1/N)*sqrt(k^2+s^2))
          #v[i] <-  atan2(s,k)
          #lõpp u12 otsing
        }
        u=c(uvalik[which.max(r<=r1)],uvalik[which.max(r<=r2)+1])
        u
        #u=c(uvalik[min(which(r<=quantile(r,probs=0.75)))],uvalik[which.min(r[r>=(quantile(r,probs=0.75)-0.1)])])
        r=0
        v=0
        k=0
        s=0
        for (i in 1:2){# empirical  cumulant function real and imaginary part at u1, u2
          k <- sum(cos(u[i]*Y))	
          s <- sum(sin(u[i]*Y))    
          r[i] <- -log((1/N)*sqrt(k^2+s^2))#actually -r[i] as in formulas
          v[i] <-  atan2(s,k)
          #lõpp  real and imaginary part of empirical cumulant function
        }
        ##cumulant estimates for 1.,..,K replicates
        alfa[j] <- log(r[1]/r[2])/log(u[1]/u[2])
        if (alfa[j] > 2){alfa[j]<- 2} 
        gamma_std[j] <- exp((log(u[1])*log(r[2])- log(u[2])*log(r[1]))/log(r[1]/r[2])) #standardized
        gamma_data[j] <- hatg*gamma_std[j]#data
        beeta1[j] <- # alpha is not 1
          ((u[2]*v[1]-u[1]*v[2])/((gamma_std[j]^alfa[j])*(tan(pi*alfa[j]/2))*(-u[1]*(u[2]^alfa[j])+u[2]*(u[1]^alfa[j])))) 
        beeta2[j] <- #alpha is 1
          ((v[1]*u[2]-v[2]*u[1])/( (2*gamma_std[j]/pi)*u[1]*u[2]*(log(u[2])-log(u[1])) )) 
        hbeeta[j]<-beeta1[j]
        if (abs(alfa[j]-1)<0.01) {hbeeta[j]<-beeta2[j]} 
        if (hbeeta[j]> 1){hbeeta[j]<-1} 
        if (hbeeta[j]< (-1)){hbeeta[j]<- -1} 
        c1 <- (u[1]^alfa[j])*((u[1]*gamma_std[j])^(1-alfa[j])-1) #constant for delta_std alpha NEQ 1 calculations 
        c2 <- (u[2]^alfa[j])*((u[2]*gamma_std[j])^(1-alfa[j])-1) #constant for delta_std alpha NEQ 1 calculations 
        delta1_std[j] <- # alpha is not 1
          ( v[1]*c2-v[2]*c1) / ( u[2]*(u[1]^alfa[j])-u[1]*(u[2]^alfa[j]) ) 
        c1 <- u[1]*log(gamma_std[j]*u[1])#constant for delta_std alpha=1 calculations 
        c2 <- u[2]*log(gamma_std[j]*u[2])#constant for delta_std alpha=1 calculations 
        delta2_std[j] <- #alpha is 1
          ( v[1]*c2-v[2]*c1) / ( u[1]*u[2]*(log(u[2])-log(u[1])) ) 
        delta1_data[j] <- hatg*(delta1_std[j])+hatd #alpha not 1
        delta2_data[j] <- hatg*(delta2_std[j])+hatd #if alpha is in 0.99-1.01
        hdelta[j]<-delta1_data[j] 
        if (abs(alfa[j]-1)<0.01) {hdelta[j]<-delta2_data[j]}
        E[U+j-1,24]<-(proc.time() - ptm)[1]
        E[U+j-1,1]<-b
        E[U+j-1,2]<-avec[l]
        E[U+j-1,3]<-hatg-g
        E[U+j-1,4]<-hatd-d
        E[U+j-1,5]<- alfa[j]-avec[l]
        E[U+j-1,6]<- hbeeta[j]-b
        E[U+j-1,7]<- gamma_data[j]-g
        E[U+j-1,8]<- hdelta[j]-d
        E[U+j-1,9]<- N
        E[U+j-1,10]<- u[1]
        E[U+j-1,11]<- u[2]
        ptm <- proc.time()
        ML<-stable.fit.mle(X,param=0)
        E[U+j-1,25]<-(proc.time() - ptm)[1]
        E[U+j-1,12]<-(ML[1]-avec[l])
        E[U+j-1,13]<-(ML[2]-b)
        E[U+j-1,14]<-(ML[3]-g)
        E[U+j-1,15]<-(ML[4]-d)
        ptm <- proc.time()
        QB<-stable.fit.quantile(X,param=0)
        E[U+j-1,26]<-(proc.time() - ptm)[1]
        E[U+j-1,16]<-(QB[1]-avec[l])
        E[U+j-1,17]<-(QB[2]-b)
        E[U+j-1,18]<-(QB[3]-g)
        E[U+j-1,19]<-(QB[4]-d)
        E[U+j-1,26]<-(proc.time() - ptm)[1]
        ptm <- proc.time()
        ChF<-stable.fit.ecf(X,gamma0=(E[U+j-1,18]+g),delta0=E[U+j-1,19],param=0)
        E[U+j-1,27]<-(proc.time() - ptm)[1]
        E[U+j-1,20]<- (ChF[1]-avec[l])
        E[U+j-1,21]<- (ChF[2]-b)
        E[U+j-1,22]<- (ChF[3]-g)
        E[U+j-1,23]<- (ChF[4]-d)
        #lõpp K
      }
      U=U+K
      #lõpp avec N
    }}} 
#LÕPP

