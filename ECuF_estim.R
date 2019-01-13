#library(stable)
#library(stabledist)
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
        
       u=c(0,0)
        u[1]=uniroot(u1f, c(10^-8,3*hatg))$root #for standard case with n=5000
        u[2]=uniroot(u2f, c(10^-8,3*hatg))$root #for standard case with n=5000
        r=0
        v=0
        k=0
        s=0
       # r=c(-0.1,-0.5)  recalculated from data
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


