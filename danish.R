library(SMPracticals)
summary(danish)
median(danish)
N=length(danish)
q=quantile(danish,probs=0.75)-quantile(danish,probs=0.25)

u <- seq(-1, 3, by=0.001)
u <- c(seq(0, 1/10^8, by=1/10^10),seq(1/10^8, 1/10^6,by=1/10^8),seq(1/10^6, 
 1/10^4,by=1/10^6),seq(1/10^4,1/10^3,by=1/10^4),seq(1/10^3,1/10^2,by=1/10^3),seq(1/10^2,0.1,by=1/10^3),seq(0.1,1,by=1/10^2))# 
u <- c(seq(1/10^8, 1/10^6,by=1/10^8),seq(1/10^6, 
1/10^4,by=1/10^6),seq(1/10^4,1/10^3,by=1/10^4) ,seq(1/10^3,10,by=1/10^2))#,seq(1/10^2,0.1,by=1/10^3)
k=0
s=0
Y=(danish-median(danish))/q #original data
Y=danish*10^6
Y=danish
median(Y)
r=0
v=0
for (i in 1:length(u)){# empirical  cumulant function real and imaginary part for resampled data
  k <- sum(cos(u[i]*Y))	
  s <- sum(sin(u[i]*Y))    
  r[i] <- log((1/N)*sqrt(k^2+s^2))
  v[i] <- atan2(s,k)}#actually -r[i] as in formulas
#in millions DKK
plot(u,r,type='l',ylim=c(-1,0))
#Leiame esialgsed skaala ja asukoha
r[(min(which(r<=-0.1))-2):(min(which(r<=-0.1))+2)]
u[(min(which(r<=-0.1))-2):(min(which(r<=-0.1))+2)]
r[(min(which(r<=-0.5))-2):(min(which(r<=-0.5))+2)]
u[(min(which(r<=-0.5))-2):(min(which(r<=-0.5))+2)]
1/u[(min(which(r<=-1))-2):(min(which(r<=-1))+2)]
v[(min(which(r<=-1))-2):(min(which(r<=-1))+2)]
Im=mean(v[(min(which(r<=-1))-2):(min(which(r<=-1))+2)])
g1=1/u[min(which(r<=-1))]
g1
plot(u,v,type='l',ylim=c(-5,5),xlim=c(0,3/g1))
d1=g1*v[min(which(r<=-1))]
d1
d1=g1*v[(min(which(r<=-1))-2)]
d1

#Standardized data
#d1=0
#Y=(danish*-d1)/g1 #ebavajalik, sest g1<1
Y=danish
k=0
s=0
r=0
v=0
u <- c(seq(0, 1/10^8, by=1/10^10),seq(1/10^8, 1/10^6,by=1/10^8),seq(1/10^6, 
1/10^4,by=1/10^6),seq(1/10^4,1/10^3,by=1/10^4) ,seq(1/10^3,10,by=1/10^2))#,seq(1/10^2,0.1,by=1/10^3)
for (i in 1:length(u)){# empirical  cumulant function real and imaginary part for resampled data
  k <- sum(cos(u[i]*Y))	
  s <- sum(sin(u[i]*Y))    
  r[i] <- log((1/N)*sqrt(k^2+s^2))
  v[i] <- atan2(s,k)}#actually -r[i] as in formulas
plot(u,r,type='l',ylim=c(-1,0))
a=0.88
b=1
d=1.207
g=0.647
curve(Ipsi,u1,u2,n=10001,add=T,lty=1,col=c18)
curve(Rpsi,u1,u2,n=10001,add=T,lty=1,col=c18)



plot(u,v)
#u valik
u12=c(u[min(which(r<=quantile(r,probs=0.75)))],u[which.min(r[r>=(quantile(r,probs=0.75)-0.1)])]) # See sobib ka a=0.1, b=1!
u12
p1=-0.1
p2=-0.25
p1=-0.1
p2=-0.5
u121=c(u[max(which(r>=p1))],u[max(which(r>=p2))])
u121
u122=c(u[min(which(r<=p1))],u[min(which(r<=p2))])
u122
u12=(u121+u122)/2
u12
r121=c(r[max(which(r>=p1))],r[max(which(r>=p2))])
r121
r122=c(r[min(which(r<=p1))],r[min(which(r<=p2))])
r122
#reduced data
#Y=(danish)/median(danish)
#Y=(danish)/(quantile(danish,probs=0.75)-quantile(danish,probs=0.25)) 
danish=danish/10^6
Y=danish*10^6
rr1=0
median(danish)
  for (i in 1:length(u)){
  k <- sum(cos(u[i]*Y))	
  s <- sum(sin(u[i]*Y))    
  rr1[i] <- log((1/N)*sqrt(k^2+s^2))}#actually -r[i] as in formulas
#end empirical  cumulant function real and imaginary part
Y=danish/median(danish)
rr2=0
median(danish)
for (i in 1:length(u)){
  k <- sum(cos(u[i]*Y))	
  s <- sum(sin(u[i]*Y))    
  rr2[i] <- log((1/N)*sqrt(k^2+s^2))}#actually -r[i] as in formulas
#end empirical  cumulant function real and imaginary part
Y=(danish)/(quantile(danish,probs=0.75)-quantile(danish,probs=0.25)) 
rr3=0
median(Y)
for (i in 1:length(u)){
  k <- sum(cos(u[i]*Y))	
  s <- sum(sin(u[i]*Y))    
  rr3[i] <- log((1/N)*sqrt(k^2+s^2))}#actually -r[i] as in formulas
#end empirical  cumulant function real and imaginary part

#Graafik
par(mar=c(5.1,5,0,2.1))
#u <- seq(0, 1, by=0.01)
plot(u,r,type='l',ylab=bquote(paste("Re(",hat(psi)[n],"(u))")),xlim=c(0,2),ylim=c(-1,0))
lines(u,rr1)
lines(u,rr2,lty=2, col=51)
lines(u,rr3,lty=3,col=4)
#abline(h=quantile(rr,probs=0.9))#proovi danish u=0.1, ...
#abline(h=quantile(rr,probs=0.9)-0.1)
abline(h=-0.1,col=2)
#abline(h=quantile(rr,probs=0.75)-0.05)
abline(h=-0.5,col=2)
#abline(h=quantile(rr,probs=0.75)-0.15)
legend('bottomright',legend=c(expression(paste("original data in millions")),
                           expression(paste("original data reduced by data median")), expression(paste("original data reduced by data interquartile range"))),bty='n',cex=1,col=c(1,51,4),lty=c(1,2,3))
rr2[which.min(rr2<=-0.1)]
rr2[which.min(rr2<=-0.1)-1]

#Kommentaa:r u valik vaatasime rr kvantiile, 0.9, 0.8, 0.75 ja -0.1, -0.2, -0.25. Arvestades ka sim a<0.5 leidsime, et parim on 

#Histogramm kahjudele MIDAGI VIGA?
quantile(danish,probs=0.99) #VaR
mean(subset(danish, danish>=quantile(danish,probs=0.99)  )) #TVaR
plot(danish,type="h", xlab="1980-1990", ylab="Losses (millions DKK)")
par(mar=c(5,2,0,0))
hist(danish,breaks=100,freq=FALSE,main="",xlab="claims sizes (millions DKK)", ylab="")
#hist(danish,breaks=100, freq=FALSE, main="",xlab="losses (millions DKK)",xlim=c(200,300), ylim=c(0.00,5/2500)) #tail
text(
  x=140, y=0.25,cex=1.3,bquote(paste(min ==.(round(min(danish),1))," ",
                      '  1stQ' ==.(round(quantile(danish,probs=0.25),1))," ",
                           "  median" == .(round(median(danish),1))," ",
                           "  mean"==.(round(mean(danish),1))," ",
                      "  3rdQ" ==.(round(quantile(danish,probs=0.75),1)),"    ",
                      VaR[0.99]==.(round(quantile(danish,probs=0.99),1)),"    ",
                      max==.(round(max(danish),1)))))
#text(  x=140, y=0.29,cex=1.5,bquote(paste("  3rdQ" ==.(round(quantile(danish,probs=0.75),2)),"    ",
#                                    VaR[0.99]==.(round(quantile(danish,probs=0.99),2)),"    ",
#                                    max==.(round(max(danish),2)))))
#text(x=40, y=0.3, cex = 1.5, bquote(paste(S[f] == .(round(max(danish),2)),' e')))
arrows(x0=42,y0=0.1,x1=57,y1=0)
text(x=40,y=0.115,57.4)
arrows(x0=58,y0=0.1,x1=66,y1=0)
text(x=58,y=0.115,65.7)
arrows(x0=128,y0=0.1,x1=145,y1=0)
text(x=126,y=0.115,144.7)
arrows(x0=144,y0=0.1,x1=152,y1=0)
text(x=144,y=0.115,152.4)
arrows(x0=253,y0=0.1,x1=263.3,y1=0)
text(x=253,y=0.115,263.3)
legend(x=33,y=0.27,"57",box.lty=0,bg="transparent")
legend(x=47,y=0.27,"66",box.lty=0,bg="transparent")
legend(x=135,y=0.2,"145",box.lty=0)
legend(x=142,y=0.2,"152",box.lty=0)
legend(x=253,y=0.2,"263",box.lty=0)
legend("top",expression('var'[t]))
#histogrammi lıpp


#Comparison with truncated data
summary(danish)
summary(subset(danish, danish<=quantile(danish,probs=0.75)& danish>=quantile(danish,probs=0.25)))
###Cumulant Estimates
##location estimation: erinevad varinadidi valitud, j‰‰me median juurde
#d0=0
d0=median(danish)
d0
##scale estimation
#g0=as.numeric((quantile(danish,0.72)-quantile(danish,0.28))) #siin vıib mıelda ka mitte jagada ?/1.654 Seda polegi tehtud!
#g0=as.numeric((quantile(danish,0.59)-quantile(danish,0.41)))
g0=as.numeric((quantile(danish,0.75)-quantile(danish,0.25)))
#g0=as.numeric((quantile(danish,0.8)-quantile(danish,0.2)))
#g0=median(danish)
g0
#empirical cumulant function
r=0 #real part
v=0 #imaginary part
#Cumulant estimates
g0=1
d0=0
alfa=0
beeta1=0
beeta2=0 #just for the control
gamma_std=0
gamma_data=0
delta1_std=0
delta2_std=0
delta1_data=0
delta2_data=0
n=length(seq_along(danish))
H=matrix(0.,nrow=1,8) #matrix for the mean, sd of the estimates at umatrix values
H=matrix(0.,1,8)
colnames(H)<-c("u[1]","u[2]","mean(alfa)","mean(beeta1)","mean(beeta2)","mean(gamma)","mean(delta1)", "mean(delta2)")
rownames(H) <- rownames(H, do.NULL = FALSE,prefix = "Obs.")#
H
#u12=c(u[min(which(r<=quantile(r,probs=0.75)))],u[which.min(r[r>=(quantile(r,probs=0.75)-0.1)])]) # See sobib ka a=0.1, b=1!
  u=u12
  u=c(1.7e-07, 1.450e-06)
 u=c( 0.5, 1.7)
for (uvalik in 1:1){#estimates (mean) at each umatrix pair of arguments
 # uvalik=1 j=1 #prooviks
  #u <- umatrix[uvalik,,drop=FALSE]
  #u=c(0.3,0.9)
  #u on eespool leitud vt 
         for (j in 1:1)  {# bootstrap replicates k=1
          #summary(Y)
       for (i in 1:2){# empirical  cumulant function real and imaginary part for resampled data
        k <- sum(cos(u[i]*Y))	
        s <- sum(sin(u[i]*Y))    
        r[i] <- -log((1/n)*sqrt(k^2+s^2))#actually -r[i] as in formulas
        v[i] <-  atan2(s,k)
                  }#end empirical  cumulant function real and imaginary part
    
    ##cumulant estimates for 1.,..,k. replicates
    alfa[j] <- log(r[1]/r[2])/log(u[1]/u[2])
    gamma_std[j] <- exp((log(u[1])*log(r[2])- log(u[2])*log(r[1]))/log(r[1]/r[2])) #standardized
    gamma_data[j] <- abs(g0)*gamma_std[j]#data
    beeta1[j] <- # alpha is not 1
      sign(g0)*(
        (u[2]*v[1]-u[1]*v[2])/
          ((gamma_std[j]^alfa[j])*(tan(pi*alfa[j]/2))*(-u[1]*(u[2]^alfa[j])+u[2]*(u[1]^alfa[j])))) 
    beeta2[j] <- sign(g0)*((v[1]*u[2]-v[2]*u[1])/( (2*gamma_std[j]/pi)*u[1]*u[2]*(log(u[2])-log(u[1])) )) #alpha is 1
    c1 <- (u[1]^alfa[j])*((u[1]*gamma_std[j])^(1-alfa[j])-1) #constant for delta_std alpha NEQ 1 calculations 
    c2 <- (u[2]^alfa[j])*((u[2]*gamma_std[j])^(1-alfa[j])-1) #constant for delta_std alpha NEQ 1 calculations 
    delta1_std[j] <- ( v[1]*c2-v[2]*c1) / ( u[2]*(u[1]^alfa[j])-u[1]*(u[2]^alfa[j]) ) 
    c1 <- u[1]*log(gamma_std[j]*u[1])#constant for delta_std alpha=1 calculations 
    c2 <- u[2]*log(gamma_std[j]*u[2])#constant for delta_std alpha=1 calculations 
    delta2_std[j] <- ( v[1]*c2-v[2]*c1) / ( u[1]*u[2]*(log(u[2])-log(u[1])) ) 
    delta1_data[j] <- g0*(delta1_std[j])+d0 #alpha not 1
    delta2_data[j] <- g0*(delta2_std[j])+d0 # if alpha is 1 (0.8-1.2)
    }#bootstrap j
  hinnangud <- cbind(
    u[1],u[2],
    mean(alfa),
    mean(beeta1),
    mean(beeta2),
    mean(gamma_data),
    mean(delta1_data),
    mean(delta2_data))
  hinnangud
}# end estimates at u
H[1,]<-hinnangud
H
 
###TULEMUSED
> hinnangud for original data R^-1(-0.3), (-0.6)
[,1] [,2]      [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
[1,]  0.6  1.3 0.8765742 1.220174 1.359359 0.427293 1.286135 1.268588
> hinnangud for original data R^-1(-0.3), (-0.5)
[,1] [,2]      [,3]     [,4]     [,5]      [,6]     [,7]     [,8]
[1,] 0.59 1.03 0.9265803 1.270463 1.364152 0.4585426 1.272524 1.258543
> hinnangud for original data R^-1(-0.2), (-0.4)
[,1]  [,2]   [,3]     [,4]      [,5]      [,6]     [,7]     [,8]
[1,] 0.394 0.786 1.0226 1.025178 0.9961411 0.5132007 1.372826 1.378214
>hinnangud at reduced median and  u=c(0.3,0.97)
      [,1] [,2]     [,3]      [,4]    [,5]      [,6]      [,7]      [,8]
[1,]  0.3 0.97 0.936882 0.9628004 1.08025 0.4658688 1.395078 1.364352

#TULEMUSED(18.07.2017)
#Y=danish/median(danish) AFMATH poster
-sum(log(dstable(danish,alpha=0.9304694, beta=1, gamma=0.453631, delta=1.320678))) #at 0.3, 0.9
#NLL=3900.1
-sum(log(dstable(danish,alpha=0.8875012, beta=1, gamma=0.404670, delta=1.320680))) #at 0.3, 0.9
#NLL=3872.485
-sum(log(dstable(danish,alpha=0.8651273  , beta=1, gamma=0.417990, delta=1.295763))) #at 0.3, 0.9
#NLL=3872.485
Y=(danish-median(danish))/range
u=c(0.3,0.9)
> hinnangud
      [,1] [,2]      [,3]      [,4]    [,5]      [,6]     [,7]     [,8]
[1,]  0.3  0.9 0.9493122 0.9355256 1.02267 0.4752084 1.405417 1.383164
-sum(log(dstable(danish,alpha=0.9493122 ,beta=0.9355256,gamma=0.4752084 ,delta= 1.405417)))
[1] 3916.093
#u=75%+0.15  
> hinnangud
      [,1] [,2]      [,3]      [,4]     [,5]      [,6]     [,7]     [,8]
[1,] 0.25 0.73 0.9219754 0.9903556 1.158225 0.4407274 1.391911 1.343979
-sum(log(dstable(danish,alpha= 0.9219754 ,beta=0.9903556,gamma=0.4407274 ,delta= 1.343979)))
[1] 3903.535
#u=75%+0.1
> hinnangud
      [,1] [,2]      [,3]      [,4]     [,5]      [,6]     [,7]     [,8]
[1,] 0.25 0.87 0.9425835 0.9861055 1.098458 0.4665395 1.389088 1.359213
-sum(log(dstable(danish,alpha= 0.9425835 ,beta=0.9861055,gamma= 0.4665395 ,delta= 1.389088)))
[1] 3913.92
# u by 0.1
#u12=c(u[min(which(rr<=quantile(rr,probs=0.75)))],u[which.min(rr[rr>=(quantile(rr,probs=0.75)-0.1)])]) See sobib ka a=0.1, b=1!
> hinnangud #kui ei lahuta median(danish), siis komapealt sama
      [,1] [,2]      [,3]      [,4]     [,5]      [,6]     [,7]     [,8]
[1,]  0.3  0.5 0.8693076 0.9732987 1.297344 0.3829411 1.397144 1.298724
-sum(log(dstable(danish,alpha= 0.8693076 ,beta=0.9732987,gamma= 0.3829411 ,delta= 1.397144)))
[1] 3902.276
#u by 0.01
#u12=c(u[min(which(rr<=quantile(rr,probs=0.75)))],u[which.min(rr[rr>=(quantile(rr,probs=0.75)-0.1)])]) See sobib ka a=0.1, b=1!
> H
u[1] u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1 0.25 0.57  0.8931764      1.06114     1.340784   0.4052438     1.345869     1.260449
-sum(log(dstable(danish,alpha= 0.8931764      ,beta=1,gamma= 0.4052438 ,delta=  1.345869   )))
[1] 3899.836 There were 50 or more warnings (use warnings() to see the first 50)
> -sum(log(dstable(danish,alpha= 0.8931764,beta=0.99,gamma= 0.4052438 ,delta=  1.345869   )))
[1] 3864.222 #muutsime b=0.99

#TULEMUSED 29.08
#1) valime v‰ikesed u leidsime gamma on 0.333, seega ei standradiseeri, ei v‰henda. Originaalandmetele meetod
#2) teeme uue  u <- seq(0, 1, by=0.01)
u[1] u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1 0.25 0.46  0.9328928    0.9103773     1.030359   0.4414247     1.431003     1.398271
-sum(log(dstable(danish,alpha= 0.9328928,beta=(0.9103773+1.030359)/2 ,gamma= 0.4414247 ,delta=  ( 1.398271  +1.431003  )/2)))
3894.21
#EI standardiseeri, p1=0.15, p2=0.25
> H
u[1]  u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1 0.301 0.521  0.9947002    0.9359447    0.9440063   0.4984578     1.416889     1.414865
-sum(log(dstable(danish,alpha= 0.946039      ,beta=0.9551293 ,gamma=0.4742873  ,delta=   1.37096   )))
 #EI standardiseeri, p1=0.1, p2=0.3
> H
u[1] u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1 0.19 0.61   0.946039    0.9551293     1.051537   0.4742873     1.395708      1.37096
-sum(log(dstable(danish,alpha= 0.946039      ,beta=0.9551293 ,gamma=0.4742873  ,delta=   1.37096   )))
[1] 3941.456
#EI standardiseeri, p1=0.1, p2=0.4
> H
u[1] u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1 0.19 0.81  0.9573832     0.970692     1.040155   0.4880072     1.381428     1.366581

-sum(log(dstable(danish,alpha=  0.9573832     ,beta=0.970692  ,gamma=0.4880072    ,delta=    1.366581 )))

#standardiseerime (danish-median)/quartile p1=0.1, p2=0.2
> H
u[1]  u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1 0.281 0.621   0.881608    0.9867292     1.264707   0.3981428     1.389257     1.308274
-sum(log(dstable(danish,alpha=  0.881608   ,beta=0.9867292  ,gamma=0.3981428    ,delta=  1.389257)))
[1] 3871.133
#valime p1=0.05, p2=0.1
u[1]  u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1 0.151 0.291   1.075729     1.328613      1.10082   0.6297922     1.096358     1.198792
> #valime p1=0.05, p2=0.1
  > H
u[1]  u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1 0.151 0.621  0.9732871      1.17181     1.246502    0.471487     1.299539     1.275344

#Proovime ka "kırgemalt"
u12=c(u[which.max(rr<=quantile(rr,probs=0.9))],u[which.min(rr[rr>=(quantile(rr,probs=0.9)-0.1)])]) 
u12 #0.1, 0.3
> hinnangud
      [,1] [,2]     [,3]     [,4]      [,5]      [,6]     [,7]    [,8]
[1,]  0.1  0.3 1.109955 1.249976 0.9349479 0.6672947 1.162426 1.31207
-sum(log(dstable(danish,alpha=1.109955 ,beta=1,gamma= 0.6672947,delta= 1.162426)))
[1] 4476.763 #tundub et liiga 0 l‰hedal miskip‰rast ei anna head tulemust
#NB! KORDAME:
#UUS(sisemine vahemik) ja piisab u samm 0.1, mitte 0.01 ja See sobib ka a=0.1, b=1!
u12=c(u[which.max(rr<=quantile(rr,probs=0.75))],u[which.min(rr[rr>=(quantile(rr,probs=0.75)-0.1)])]) 
u12 #0.3, 0.5
[1] 3902.276
#VANA:u1 sisemine ja u2 esimene v‰ljas, aga see ei sobi kui a<0.5, b=1
u12=c(u[min(which(rr<=quantile(rr,probs=0.75)))],u[min(which(rr<=(quantile(rr,probs=0.75)-0.1)))]) 
#see on varasem l‰henemine, aga n¸¸d sean piirid ja vıtan sealt seest. enne vıtsin alumise v‰ljastpoolt.
u12 #0.3,0.6
[1] 3903.535
#Ja "kırgemalt" kui75% alustamine ei anna nii head tulemust
#Y=danish/range st mediaanita, aga jagame range mitte mediaan nagu AFMATH
#u=75%+0.1=75%+0.15=0.3,0.9 st alimine v‰ljast
> hinnangud
      [,1] [,2]      [,3]      [,4]    [,5]      [,6]     [,7]     [,8]
[1,]  0.3  0.9 0.9493122 0.9355256 1.02267 0.4752084 1.405417 1.383164
-sum(log(dstable(danish,alpha= 0.9493122 ,beta= 0.9355256,gamma= 0.4752084 ,delta=  1.405417)))
[1] 3935.388
#u=75%+0.1 alumine seest ja tuleb sama mis danish-medina(danish)/range ?
  >    u[1] u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1  0.3  0.5  0.8693076    0.9732987     1.297344   0.3829411     1.397144     1.298724

#TULEMUSED lıpp (18.07.2017)



#VANA KOODIDJUPP (afmath, Fama&Roll inspiration)
>hinnagud at reduced 50% range and u=???
> H #original data
u[1]  u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1  0.03  0.09  1.4002948     1.271524    0.3645770  1.19699926    1.2116587    1.8276779
Obs.2  0.03  0.90  1.1410112     1.319266    0.9363118  0.56205746    1.3107117    1.3840898
Obs.3  0.03  9.00  0.9517399     2.432523    2.6878696  0.24955543    0.6223326    0.6617722
Obs.4  0.03 90.00  0.7586874     4.029294    7.4681748  0.07182361    0.2813672    0.6552075
Obs.5  0.30  0.09  1.0215268     1.214224    1.1513812  0.52395934    1.2498950    1.2736727
Obs.6  0.30  0.90  1.0126710     1.048857    1.0313690  0.51554941    1.3493927    1.3524836
Obs.7  0.30  9.00  0.7821498     2.254950    2.5902248  0.29741412    0.6797065    0.7991204
Obs.8  0.30 90.00  0.5796257     3.336401    4.2361320  0.12783658    0.3704358    0.8600781
Obs.9  3.00  0.09  0.9346930     3.330654    3.6889138  0.39451937   -0.6624312   -0.6460022
Obs.10 3.00  0.90  0.7767051     7.340175    7.7158710  0.40827770   -0.4394944   -0.3954186
Obs.11 3.00  9.00  0.5575955    -6.207237   -3.6166019  0.44214616   -1.1889148   -1.0897231
Obs.12 3.00 90.00  0.3699857    -3.624883   -0.7422615  0.51024316   -1.0879598   -0.9047923
-sum(log(dstable(danish,alpha=1.0126710,beta=1  ,gamma=0.51554941 ,delta= 1.3524836 )))
> H # 59%,41% d0=59%,41%mean
u[1]  u[2] mean(alfa) cv(alfa) mean(beeta1) cv(beeta1) mean(beeta2) cv(beeta2) mean(gamma) cv(gamma) mean(delta1) cv(delta1) mean(delta2) cv(delta2)
Obs.1  0.03  0.09  1.0838610       NA    1.2573998         NA   1.00891124         NA   0.6149956        NA     1.177772         NA     1.289166         NA
Obs.2  0.03  0.90  0.9553654       NA    1.2151956         NA   1.31757062         NA   0.4047863        NA     1.322969         NA     1.319671         NA
Obs.3  0.03  9.00  0.8115847       NA    0.9887837         NA   1.28624789         NA   0.2166620        NA     1.683621         NA     1.780733         NA
Obs.4  0.03 90.00  0.6100033       NA    1.3680769         NA   3.14959559         NA   0.0549232        NA     1.702109         NA     1.909041         NA
Obs.5  0.30  0.09  0.9769780       NA    0.9571047         NA   0.99328254         NA   0.4935356        NA     1.395401         NA     1.387271         NA
Obs.6  0.30  0.90  0.8031846       NA    1.1217805         NA   1.24865908         NA   0.3951776        NA     1.320905         NA     1.315344         NA
Obs.7  0.30  9.00  0.6650906       NA    0.3417714         NA   0.30634132         NA   0.3048786        NA     1.589840         NA     1.618120         NA
Obs.8  0.30 90.00  0.4412708       NA    0.1789425         NA   0.13982952         NA   0.1418257        NA     1.653802         NA     1.679123         NA
Obs.9  3.00  0.09  0.8410448       NA    0.1960035         NA   0.22212374         NA   0.3441192        NA     1.948736         NA     1.957457         NA
Obs.10 3.00  0.90  0.7396587       NA   -2.1769261         NA  -1.84023035         NA   0.3900918        NA     1.448075         NA     1.437566         NA
Obs.11 3.00  9.00  0.4452772       NA    3.5656925         NA   0.75058251         NA   0.7759434        NA     3.348093         NA     2.553600         NA
Obs.12 3.00 90.00  0.2187451       NA    2.7993349         NA   0.02755441         NA   4.6447340        NA     6.267055         NA     2.199519         NA
-sum(log(dstable(danish,alpha=0.9769780,beta=0.99328254   ,gamma=0.4935356 ,delta=1.387271 )))
[1] 3954.37
> H #g0=28%72% d0=25%mean
       u[1]  u[2] mean(alfa) cv(alfa) mean(beeta1) cv(beeta1) mean(beeta2) cv(beeta2) mean(gamma) cv(gamma) mean(delta1) cv(delta1) mean(delta2) cv(delta2)
Obs.1  0.03  0.09  1.4567238       NA    1.2229170         NA    0.2759654         NA  1.35881961        NA     1.274093         NA     1.963700         NA
Obs.2  0.03  0.90  1.1670045       NA    1.3161614         NA    0.8515978         NA  0.58375471        NA     1.322332         NA     1.431689         NA
Obs.3  0.03  9.00  1.0107525       NA    0.8253205         NA    0.8073511         NA  0.30270537        NA     2.016327         NA     2.013025         NA
Obs.4  0.03 90.00  0.7888348       NA    1.6022316         NA    2.8064441         NA  0.07616511        NA     1.794633         NA     1.920918         NA
Obs.5  0.30  0.09  1.0646724       NA    1.2495548         NA    1.0561091         NA  0.58156021        NA     1.198748         NA     1.280913         NA
Obs.6  0.30  0.90  0.9894312       NA    0.9633131         NA    0.9791214         NA  0.50143671        NA     1.390941         NA     1.387499         NA
Obs.7  0.30  9.00  0.8476137       NA    0.1817009         NA    0.2019753         NA  0.35301986        NA     1.924245         NA     1.931972         NA
Obs.8  0.30 90.00  0.6019672       NA    0.5597144         NA    0.7553096         NA  0.12997571        NA     1.746365         NA     1.827640         NA
Obs.9  3.00  0.09  0.9504001       NA    1.2238912         NA    1.3357806         NA  0.39804829        NA     1.318908         NA     1.318163         NA
Obs.10 3.00  0.90  0.8005124       NA    1.1164654         NA    1.2278807         NA  0.39612322        NA     1.318531         NA     1.317408         NA
Obs.11 3.00  9.00  0.7574146       NA   -2.4138725         NA   -2.0313640         NA  0.39543095        NA     1.345779         NA     1.343541         NA
Obs.12 3.00 90.00  0.4065315       NA   -1.4020672         NA   -0.4320668         NA  0.38449636        NA     1.341564         NA     1.333329         NA

-sum(log(dstable(danish,alpha=0.9894312,beta=0.9791214  ,gamma=0.50143671 ,delta= 1.387499  )))
[1] 3971.957
-sum(log(dstable(danish,alpha= 1.0107525,beta=0.8073511  ,gamma=0.30270537  ,delta= 2.013025   )))
[1] 6794.179
> H #d0=0, g0=72%28%
        u[1]  u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1  0.03  0.09  1.4567238    1.2229170    0.2759654  1.35881961    1.2740930    1.9636996
Obs.2  0.03  0.90  1.1670045    1.3161614    0.8515978  0.58375471    1.3223321    1.4316890
Obs.3  0.03  9.00  1.0107525    2.4166720    2.3640547  0.30270537    0.5512813    0.5416134
Obs.4  0.03 90.00  0.7888348    4.1185542    7.2139960  0.07616511    0.2531627    0.5777811
Obs.5  0.30  0.09  1.0646724    1.2495548    1.0561091  0.58156021    1.1987477    1.2809127
Obs.6  0.30  0.90  0.9894312    0.9633131    0.9791214  0.50143671    1.3909407    1.3874985
Obs.7  0.30  9.00  0.8476137    2.1955034    2.4404801  0.35301986    0.6068004    0.7001668
Obs.8  0.30 90.00  0.6019672    3.2463286    4.3807759  0.12997571    0.3219737    0.7933649
Obs.9  3.00  0.09  0.9504001    3.8703868    4.2242215  0.39804829   -1.2650453   -1.2674020
Obs.10 3.00  0.90  0.8005124    8.8028545    9.6813172  0.39612322   -1.3058804   -1.3147358
Obs.11 3.00  9.00  0.7574146   -6.0898299   -5.1248198  0.39543095   -1.1921258   -1.1977712
Obs.12 3.00 90.00  0.4065315   -4.7367374   -1.4596925  0.38449636   -1.1902355   -1.2180583
-sum(log(dstable(danish,alpha=0.9894312,beta=0.9791214  ,gamma=0.50143671 ,delta= 1.387499  )))
[1] 3971.957
> H H #d0=0, g0=75%25%
        u[1]  u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1  0.03  0.09  1.4734096    0.9507300    0.1903415  1.46618737    1.5661181    2.1698894
Obs.2  0.03  0.90  1.1774860    1.2809001    0.7835215  0.60510204    1.3264685    1.4655224
Obs.3  0.03  9.00  1.0294026    2.5299101    2.3733314  0.32101844    0.3751211    0.3507747
Obs.4  0.03 90.00  0.8125719    3.9003504    6.4896517  0.08363647    0.3108857    0.5891143
Obs.5  0.30  0.09  1.1156650    1.2978047    0.9511464  0.67418117    1.1270760    1.2940908
Obs.6  0.30  0.90  0.9493122    0.9355256    1.0226700  0.47520840    1.4054165    1.3831637
Obs.7  0.30  9.00  0.8554495    2.3117669    2.6145907  0.36737586    0.4296177    0.5171476
Obs.8  0.30 90.00  0.6213092    2.9720569    4.1580273  0.13775188    0.3629931    0.8020223
Obs.9  3.00  0.09  0.9649712    1.2534304    1.3413518  0.40903435    1.3124259    1.3081572
Obs.10 3.00  0.90  0.8285659    1.1500995    1.2953279  0.39624372    1.3107004    1.3021246
Obs.11 3.00  9.00  0.7910485    5.4286018    4.9108462  0.39204250    1.0604634    1.0870330
Obs.12 3.00 90.00  0.4419961    4.6198162    1.9103339  0.32553329    1.0437204    1.2087646

Obs.6  0.30  0.90  0.9493122    0.9355256    1.0226700  0.47520840    1.4054165    1.3831637
-sum(log(dstable(danish,alpha=0.9493122, beta=0.9355256, gamma=0.47520840, delta=1.4054165)))
[1] 3935.389
> H #d0=0, g0=median
        u[1]  u[2] mean(alfa) cv(alfa) mean(beeta1) cv(beeta1) mean(beeta2) cv(beeta2) mean(gamma) cv(gamma) mean(delta1) cv(delta1) mean(delta2) cv(delta2)
Obs.1  0.03  0.09  1.4151740       NA    0.7779971         NA    0.1785724         NA   1.3002113        NA    1.7647334         NA    2.2563914         NA
Obs.2  0.03  0.90  1.1758103       NA    1.2659656         NA    0.7660293         NA   0.6078574        NA    1.3250897         NA    1.4759016         NA
Obs.3  0.03  9.00  1.0440666       NA    2.5113737         NA    2.2806160         NA   0.3447112        NA    0.3116926         NA    0.2747531         NA
Obs.4  0.03 90.00  0.8524601       NA    3.5377294         NA    5.0886782         NA   0.1104765        NA    0.3303973         NA    0.5893558         NA
Obs.5  0.30  0.09  1.1855222       NA    1.3604831         NA    0.8260226         NA   0.7802252        NA    1.0810394         NA    1.3548648         NA
Obs.6  0.30  0.90  0.9258032       NA    0.9672474         NA    1.1121474         NA   0.4523610        NA    1.3943584         NA    1.3551831         NA
Obs.7  0.30  9.00  0.8741230       NA    2.3242020         NA    2.5999899         NA   0.3904785        NA    0.3701640         NA    0.4488452         NA
Obs.8  0.30 90.00  0.6737712       NA    2.7048867         NA    3.4111023         NA   0.1783492        NA    0.3936155         NA    0.8279505         NA
Obs.9  3.00  0.09  0.9832502       NA    1.3220375         NA    1.3678999         NA   0.4083710        NA    1.2998943         NA    1.2966776         NA
Obs.10 3.00  0.90  0.8333979       NA    1.1788924         NA    1.3497700         NA   0.3877679        NA    1.2996926         NA    1.2858821         NA
Obs.11 3.00  9.00  0.8670737       NA    5.4924578         NA    5.2591196         NA   0.3929175        NA    0.9473172         NA    0.9695674         NA
Obs.12 3.00 90.00  0.5358575       NA    3.9931442         NA    2.0941162         NA   0.3211140        NA    1.0001707         NA    1.1729066         NA
-sum(log(dstable(danish,alpha=0.9258032, beta=1, gamma=0.4523610, delta=1.3551831)))
[1] 3914.184
-sum(log(dstable(danish,alpha=0.9258032, beta=0.9672474, gamma=0.4523610, delta=1.3943584)))
[1] 3900.117 
qstable(0.99,alpha=0.9258032, beta=0.9672474, gamma=0.4523610, delta=1.3943584)
[1] 43.87247

>H #g0=28%72%jagatud 1.65 d0=25%mean
        u[1]  u[2] mean(alfa) mean(beeta1) mean(beeta2) mean(gamma) mean(delta1) mean(delta2)
Obs.1  0.03  0.09  1.3478583    1.3335971    0.4873213  1.05778384    1.1481082    1.6583117
Obs.2  0.03  0.90  1.0814208    1.3539704    1.1266661  0.48644372    1.2905650    1.3218850
Obs.3  0.03  9.00  0.9105328    1.2079025    1.4269156  0.23266302    1.6634272    1.7065577
Obs.4  0.03 90.00  0.6812046    1.7028733    4.1129689  0.04834441    1.7443146    1.9032920
Obs.5  0.30  0.09  0.9585255    1.1453816    1.2576783  0.45921169    1.3063139    1.2703541
Obs.6  0.30  0.90  0.9496646    1.1469442    1.2129879  0.45204942    1.3050160    1.2961961
Obs.7  0.30  9.00  0.7522847    0.5882994    0.6421369  0.28934760    1.5983865    1.6358827
Obs.8  0.30 90.00  0.4942619    0.4605312    0.6138150  0.09434900    1.7097367    1.7693382
Obs.9  3.00  0.09  0.9064187    1.1583073    1.3047366  0.38894889    1.3255596    1.3455550
Obs.10 3.00  0.90  0.8148504    0.9859589    0.9778294  0.40920406    1.3193617    1.3276148 #PARIM beta2, delta2
Obs.11 3.00  9.00  0.4863387   -1.8689560   -0.7163616  0.57459100    0.8384026    0.9789747
Obs.12 3.00 90.00  0.2336798   -3.6143618   -0.1552105  1.42797113   -0.2666780    0.9523573
-sum(log(dstable(danish,alpha=0.8148504 ,   beta=0.9859589    ,  gamma=0.40920406  ,  delta=1.3193617 )))
[1] 3861.758
#PARIM, aga VaR v‰ga suur!
-sum(log(dstable(danish,alpha=0.8148504   , beta=    0.9778294 , gamma=0.40920406  ,  delta= 1.3276148 )))
[1] 3860.766
qstable(0.99,alpha=0.8148504   , beta=    0.9778294 , gamma=0.40920406  ,  delta= 1.3276148)
[1] 77.11202

##STABLE log-likelihood ja kvantiilid muud meetodid
n=length(danish)
param=4
#MLE
NLL=-sum(log(dstable(danish,alpha=0.882,beta=0.99,gamma=0.372451,delta=1.31661)))#MLE
2*NLL+2*param #MLE
2*NLL+param*log(n) #MLE
qstable(0.99,alpha=0.882,beta=0.99,gamma=0.372451,delta=1.31661)
#Quantile
-sum(log(dstable(danish,alpha=1.005,beta=0.99,gamma=0.744903,delta=1.63386)))#Quantile
qstable(0.99,alpha=1.005,beta=0.999,gamma=0.744903,delta=1.63386)
[1] 4202.418
2*4202.418+2*param 
2*4202.418+param*log(n) 
#Characteristic
-sum(log(dstable(danish,alpha=0.9491,beta=1,gamma=0.462746,delta=1.30677)))#Characteristic
qstable(0.99,alpha=0.9491,beta=1,gamma=0.462746,delta=1.30677)
[1] 4202.418
2*4202.418+2*param #MLE
2*4202.418+param*log(n) #MLE


###BOOTSTRAP Y=1:n, for (j in 1:1) asemel for (j in 1:20) ja siis:
t <- runif(n, 1/n,1) #kui uvalik=1, siis siit algus

for (r in 1:n)   {#resampling with replacing
  a <- t[r]%/%(1/n)
  Y[r] <- X[a]
}#end resampling


##MUUD JAOTUSED
# MLE for Cauchy #We formulate the log-likelihood function.
library(rmutil)
LL <- function(mu,sigma) { R = dstable(danish, mu, sigma) #
-sum(log(R))  }
# And apply MLE to estimate the two parameters 
library(stats4)
  mle(LL, start = list(mu = 0.3, sigma=1))
 
  Call:
    mle(minuslogl = LL, start = list(mu = 0.3, sigma = 1))
  
  Coefficients:
    mu     sigma 
  1.4546057 0.4926708 
  
#AIC
-2*sum(log(dcauchy(danish,location=1.4546057,scale=0.4926708)))+2*2
  qcauchy(0.99,location=1.4546057,scale=0.4926708)
#AIC  9130.979
#NLL   5039.77

 # MLE for Levy #We formulate the log-likelihood function.
  library(rmutil)
  LL <- function(mu,sigma) { R = dlevy(danish, mu, sigma) #
  -sum(log(R))  }
  # And apply MLE to estimate the two parameters 
  library(stats4)
  mle(LL, start = list(mu = 0.3, sigma=1))
  Call:
    mle(minuslogl = LL, start = list(mu = 0.3, sigma = 1))
  
  Coefficients:
    mu     sigma 
  0.2946309 1.2008186 
 #AIC
  -2*sum(log(dlevy(danish,m=0.2945,s=1.201289)))+2*2
  #AIC  10083.54
  #NLL   5039.77  
  qlevy(0.99,m=0.2945,s=1.201289)
  # MLE for Burr #We formulate the log-likelihood function.
    library(fitdistrplus)
    fitdist(as.matrix(danish)[,1], "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1))
    
#VaR
    qburr(0.99,shape1=0.087683260,14.93807786,rate= 1.08585632)  
[1] 30.9853
    qburr(0.99,shape1=0.087683260,14.93807786,1.08585632)  
    [1] 30.9853
    #AIC
-2*sum(log(dburr(danish,shape1=0.087683260,14.93807786,1.08585632)))+2*3
  [1] 7676.239
##pareto
    fitdist(as.matrix(danish)[,1], "pareto", start = NULL)  
    Fitting of the distribution ' pareto ' by maximum likelihood 
    Parameters:
      estimate Std. Error
    shape  5.164748  0.4222731
    scale 11.887356  1.1246792
    -2*sum(log(dpareto(danish,shape=5.164748,scale=11.887356)))+2*3
        qpareto(0.99,shape=5.164748,scale=11.887356) 
##Inv.Burr
        library(actuar})
        LL <- function(s1,s2,r) { R = dinvburr(danish, s1,s2,r) #
        -sum(log(R))  }
        # And apply MLE to estimate the two parameters 
        library(stats4)
        mle(LL, start = list(s1 = 1, s2=1,r=11))       
        Call:
          mle(minuslogl = LL, start = list(s1 = 1, s2 = 1, r = 11))
        
        Coefficients:
          s1         s2          r 
        111.256592   2.024997   7.096223 
  sum(log(dinvburr(danish,shape1=111.25, shape2=2.02,rate=7.09)  ))
  qinvburr(0.99,shape1=111.2565, shape2=2.024997,rate=7.096223) 
##GenPareto
  LL <- function(s1,s2,r) { R = dgenpareto(danish, s1,s2,r) #
  -sum(log(R))  }
  # And apply MLE to estimate the two parameters 
  library(stats4)
  mle(LL, start = list(s1 = 1, s2=1,r=11))
  Call:
    mle(minuslogl = LL, start = list(s1 = 1, s2 = 1, r = 11))
  
  Coefficients:
    s1         s2          r 
  2.764906 230.881676  51.382574 
 sum(log(dgenpareto(danish, 2.764906, 230.881676 , 51.382574 )))
 qgenpareto(0.99,2.764906, 230.881676 , 51.382574)
 ##LogLogistic
 library(actuar})
 LL <- function(s1,r) { R = dllogis(danish, s1,r) #
 -sum(log(R))  }
 # And apply MLE to estimate the two parameters 
 library(stats4)
 mle(LL, start = list(s1 = 1, r=11))       
 sum(log(dllogis(danish,shape=2.653,scale=1.770)  ))
 qllogis(0.99,shape=2.653,scale=1.770) 
##log-normal
  LL <- function(s1,r) { R = dlnorm(danish, s1,r) #
 -sum(log(R))  }
 # And apply MLE to estimate the two parameters 
 library(stats4)
 mle(LL, start = list(s1 = 1, r=11)) 
 Call:
   mle(minuslogl = LL, start = list(s1 = 1, r = 11))
 
 Coefficients:
   s1         r 
 0.6718793 0.7323344 
 sum(log(dlnorm(danish,0.6718793 ,0.7323344 )  ))
 qlnorm(0.99,0.6718793, 0.7323344 ) 
 ##InvPareto
 library(actuar})
 LL <- function(s1,r) { R = dinvpareto(danish, s1,r) #
 -sum(log(R))  }
 # And apply MLE to estimate the two parameters 
 library(stats4)
 mle(LL, start = list(s1 = 1, r=11))  
 Call:
   mle(minuslogl = LL, start = list(s1 = 1, r = 11))
 
 Coefficients:
   s1            r 
 114.51096654   0.01418952 
 sum(log(dinvpareto(danish,114.51096654 ,  0.01418952 )))
 qinvpareto(0.99,114.51096654  , 0.01418952 )
 ##gamma
 LL <- function(s1,r) { R = dgamma(danish, s1,r) #
 -sum(log(R))  }
 # And apply MLE to estimate the two parameters 
 library(stats4)
 mle(LL, start = list(s1 = 1, r=1)) 
 Call:
   mle(minuslogl = LL, start = list(s1 = 1, r = 1))
 
 Coefficients:
   s1         r 
 1.2579981 0.4107491 
 qgamma(0.99,1.2579981, 0.4107491 )
 

  
  library(mixtools)
  gm<-normalmixEM(danish,k=4,lambda=c(0.25,0.25,0.25,0.25),mu=c(-0.01,0.01,1,1),sigma=c(0.01,0.02,1,1))
  gm$lambda
  gm$mu
  gm$sigma
    gm$loglik
  

    
    
    
    
  