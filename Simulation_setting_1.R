#Code for simulation scenario 1
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
Rcpp::sourceCpp("MCCLC_full_categorical.cpp") #full SEM-Gibbs
#Rcpp::sourceCpp("MCCLC_hybrid_categorical.cpp") #hybrid SEM-Gibbs

L <- 4
R <- H <- 2
K <- 24
Q <- 36
ph <- c(0.625,0.375)
pr <- c(0.25,0.75)

pxwz<-matrix(c(0.4,0.3,0.2,0.1,0.15,0.05,0.35,0.45,0.05,0.55,0.15,0.25,0.3,0.05,0.5,0.15),nrow=H*R,ncol=L,byrow=T)

c1=2
c2=3
c3=4

p1<-matrix(c(0.8,0.9,0.3,0.1,0.2,0.1,0.7,0.9),nrow=c1,ncol=L,byrow=T)
p2<-matrix(c(0.7,0.8,0.2,0.3,0.3,0.2,0.8,0.7),nrow=c1,ncol=L,byrow=T)
p3<-matrix(c(0.7,0.8,0.1,0.15,0.1,0.05,0.7,0.65,0.2,0.15,0.2,0.2),nrow=c2,ncol=L,byrow=T)
p4<-matrix(c(0.8,0.1,0.7,0.2,0.05,0.7,0.1,0.7,0.15,0.2,0.2,0.1),nrow=c2,ncol=L,byrow=T)
p5<-matrix(c(0.65,0.1,0.6,0.2,0.1,0.65,0.1,0.6,0.1,0.1,0.2,0.1,0.15,0.15,0.1,0.1),nrow=c3,ncol=L,byrow=T)
p6<-matrix(c(0.75,0.2,0.7,0.1,0.05,0.6,0.1,0.7,0.05,0.05,0.1,0.15,0.15,0.15,0.1,0.05),nrow=c3,ncol=L,byrow=T)


S <- 100
for(s in 1:S){

components<-rep(0,K)
while(length(which(components==1))!=round(K*ph[1])){ #remove for simulation scheme i. (random sampling of memberships of level-2 units)
  components <- sample(1:H,prob=ph,size=K,replace=TRUE)      
}

components_out <- cbind(components_out, components)

components2<-rep(0,Q)
while(length(which(components2==1))!=round(Q*pr[1])){ #remove for simulation scheme i. (random sampling of memberships of level-2 units)
  components2 <- sample(1:R,prob=pr,size=Q,replace=TRUE)      
}

components2_out <- cbind(components2_out, components2)


nk<-13 # scenario 1.1
nq<-13
#nk<-24 # scenario 1.2
#nq<-24

n=K*nk*Q

data <- matrix(nrow=n,ncol=3) 
data[,1] <- seq(1:n)
d1<-NULL
for(k in 1:K){
  d1<-c(d1,rep(k,nk*Q))
}

data[,2]<-d1
d2<-NULL
for(k in 1:K){
  d2<-c(d2,rep(seq(1:Q),nq))
}

data[,3]<-d2

colr<-NULL
colr<-c(colr,rep(components2,nq*K))
colh<-NULL
for(i in 1:(K)){
  colh<-c(colh,rep(components[i],nk*Q))
}
datac<-cbind(data,colh,colr)


colh_out <- cbind(colh_out, colh)
colr_out <- cbind(colr_out, colr)

datacc<-datac[order(datac[,4],datac[,5]),]
data<-datacc[,1:3]

count<-c(length(which(components==1))*length(which(components2==1))*nk,length(which(components==1))*length(which(components2==2))*nk,length(which(components==2))*length(which(components2==1))*nk,length(which(components==2))*length(which(components2==2))*nk)


data2<-NULL
samples2 <- NULL
w <- 1
for(j in (1:H)){
  for(m in (1:R)){
    samples <- sample(1:L,prob=pxwz[w,],size=conteggio[w],replace=TRUE) 
    for(i in (1:count[w])){
      data2p = cbind(sample(0:(c1-1),prob=p1[,samples[i]],size=1),sample(0:(c1-1),prob=p2[,samples[i]],size=1),sample(0:(c2-1),prob=p3[,samples[i]],size=1),sample(0:(c2-1),prob=p4[,samples[i]],size=1),sample(0:(c3-1),prob=p5[,samples[i]],size=1),sample(0:(c3-1),prob=p6[,samples[i]],size=1))
      data2=rbind(data2,data2p)
    } 
    samples2<-c(samples2,samples)  
    w=w+1    
  }
}
data<-cbind(data,data2)  
datac<-cbind(data,samples2)


a<-order(data[,2],data[,3])
s<-NULL
for(i in 1:(dim(data)[1])){
  s<-rbind(s,data[a[i],])
}
s<-as.matrix(s)
data<-s


y<-list(s[1,])
for(i in (2:dim(s)[1])){
  y[[i]]<-s[i,]
}

main2(y)


}

