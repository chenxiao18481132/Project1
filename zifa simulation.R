#simulation

library(MASS)
simulation<-function(D,K,N,cluster_n,lambda){
if(cluster_n==1){
Z=matrix(rnorm(N*K),N,K)
cluster_id=rep(1,N)
}else{
  Z=matrix(0,N,K)
  cluster_id=sample(1:cluster_n,N,replace = TRUE)
  cluster_mu=matrix((runif(K*cluster_n,0,10)-5)/2,cluster_n,K)#+0.8*matrix(rep(1:cluster_n,K),cluster_n,K)
  for(i in 1:N)
    Z[i,]=mvrnorm(1,cluster_mu[cluster_id[i],],diag(rep(0.05,K)))
}
Aa=matrix(runif(D*K,-0.5,0.5),D,K)#A
mu=runif(D,2.8,3.2)
sigma=runif(D,0.2,0.4)
X<-matrix(0,N,D)
for(i in 1:N){
e=mvrnorm(n=1, rep(0, D), diag(sigma))
X[i,]=Aa%*%Z[i,]+mu+e
}

P=exp(-lambda*X^2)#p(x_ij)
h<-function(x){a=rbinom(1,1,x)
return(a)}
Y<-matrix(0,N,D)
for (i in 1:N) {
  for (j in 1:D) {
    Y[i,j]=(1-h(P[i,j]))*X[i,j]
    
  }
  
}
#sum(Y==0)
list0=list(Y=Y,Z=Z,X=X,A=Aa,sigma=sigma,mu=mu,c_id=cluster_id)
return(list0)
}
#  testing
D=20
K=2
N=200
cluster_n=3
lambda=0.1
set.seed(67)
testdata<-simulation(D,K,N,cluster_n,lambda)
Y=testdata$Y
Z=testdata$Z
X=testdata$X
A=testdata$A
sigmas=testdata$sigma
cluster_id=testdata$c_id
mu=testdata$mu
source("ZIFA_EM.R")

#TEST
a<-fitmodel(Y,2)
print(crossprod(a$lambda-lambda))
print(crossprod(a$sigmas-sigmas))
print(sum((a$A-A)^2))
print(crossprod(a$mu-mu))
km<-kmeans(a$EZ,2,nstart = 40)
zh=a$EZ
par(mfrow=c(1,2))
plot(zh,col="white",xlab="dim1",ylab="dim2",main=" E(Z)")#  use EZ to cluster the simulaiton data
points(zh[cluster_id==1,1:2],col=2,pch=0) 
points(zh[cluster_id==2,1:2],col=3,pch=6)
points(zh[cluster_id==3,1:2],col=4,pch=2)
plot(Z,col="white",,xlab="dim1",ylab="dim2",main="Z")# real cluster_id
points(Z[cluster_id==1,1:2],col=2,pch=15)
points(Z[cluster_id==2,1:2],col=3,pch=25)
points(Z[cluster_id==3,1:2],col=4,pch=17)