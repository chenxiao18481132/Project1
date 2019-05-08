
Estep<-function(A_old,sigma_old,lamda_old,mu_old,Y)
{#Input:
  #A_old,sigma_old,lamda_old,mu_old: old parameters.
  # Y: observed values.
  D=nrow(A_old)
  K=ncol(A_old)
  N=nrow(Y)
  #EX,EXZ,EX2,EZ,EZZT
  EX=matrix(0,N,D)
  EXZ=array(0,c(N,D,K))
  EX2=matrix(0,N,D)
  EZ=matrix(0,N,K)
  EZZT=array(0,c(N,K,K))
  #1.1P(Z,X)
  U=c(rep(0,K),mu_old)
  R1=cbind(diag(rep(1,K)),t(A_old))
  R2=cbind(A_old,(diag(sigma_old)+A_old%*%t(A_old)))
  R=rbind(R1,R2)
  for(i in 1:N)
  {
    Y_i=Y[i,]
    Y_is_zero=abs(Y_i)<1e-6
    #1.2P(Z,x_0,x+)
    label=(K+1):(D+K)
    Label=c(1:K,label[Y_is_zero],label[!Y_is_zero])
    sigma_ZX0=R[c(1:K,label[Y_is_zero]),c(1:K,label[Y_is_zero])]
    sigma_zx0_xplus=R[c(1:K,label[Y_is_zero]),label[!Y_is_zero]]
    sigma_xplus=R[label[!Y_is_zero],label[!Y_is_zero]]
    mu_zx0=U[c(1:K,label[Y_is_zero])]
    MU_XPLUS=U[label[!Y_is_zero]]
    #2.P(Z,X_0|X+)
    sigma_aa_ab=sigma_zx0_xplus%*%solve(sigma_xplus)
    mu_ZX0_Xplus=mu_zx0+sigma_aa_ab%*%(Y_i[!Y_is_zero]-MU_XPLUS)
    SIGMA_ZX0_Xplus=sigma_ZX0-sigma_aa_ab%*%t(sigma_zx0_xplus)
    #3.P(Z,X_0|Y)
    if(length(mu_ZX0_Xplus)>K){#X_0  exsit
      L=length(mu_ZX0_Xplus)#ZX0  LENGTH
      s_inverve=solve(SIGMA_ZX0_Xplus)#L*L
      sigma_zx0_y=solve(s_inverve+2*lamda_old*(diag(c(rep(0,K),rep(1,L-K)))))
      mu_zx0_y=sigma_zx0_y%*%s_inverve%*%mu_ZX0_Xplus
    
    #EX
      EZ[i,]=mu_zx0_y[1:K]
      EX[i,Y_is_zero]=mu_zx0_y[(K+1):L]
      EX2[i,Y_is_zero]=(mu_zx0_y[(K+1):L])^2+diag(as.matrix(sigma_zx0_y[(K+1):L,(K+1):L]))
      EZZT[i,,]=sigma_zx0_y[1:K,1:K]+mu_zx0_y[1:K]%*%t(mu_zx0_y[1:K])
      EXZ[i,Y_is_zero,]=sigma_zx0_y[(K+1):L,1:K]+mu_zx0_y[(K+1):L]%*%t(mu_zx0_y[1:K])
    }else{#latent X_0 does not exsit, 
      EZ[i,]=mu_ZX0_Xplus
      EZZT[i,,]=SIGMA_ZX0_Xplus+mu_ZX0_Xplus%*%t(mu_ZX0_Xplus)
    }
  }
  LIST=list(EZ=EZ,EX=EX,EX2=EX2,EZZT=EZZT,EXZ=EXZ)
  return(LIST)
}
#Mstep 
#Input:
#A_old,sigma_old,lamda_old,mu_old: old parameters.
# Y: observed values.
#EZ, EZZT, EX, EXZ, EX2: expectations of latent variables computed in E-step.


#Returns:
# A, mu, sigma, lamda: new values of parameters.
Mstep<-function(A_old,sigma_old,lamda_old,mu_old,Y,EZ,EX,EX2,EZZT,EXZ){
  D=nrow(A_old)
  K=ncol(A_old)
  N=nrow(Y)
  A = matrix(0,D,K)
  mu=rep(0,D)
  Y_is_zero = abs(Y) < 1e-6
  S=matrix(0,(K+1),(K+1))
  #A  mu
  #j=1:D   A:D,K   Y:N,D
  for(i in 1:N){
  S[1:K,1:K]=S[1:K,1:K]+EZZT[i,,]
  S[1:K,K+1]=S[1:K,K+1]+EZ[i,]
  
  }
  S[K+1,K+1]=N
  S[K+1,1:K]=t(S[1:K,K+1])
  #
  #j  1:D
  p1<-matrix(0,D,K)
  p2<-c(0)
  for (j in 1:D) {
   
  #EX_ij Z_i+y_ij EZ_i
    if(sum(Y_is_zero[,j])>1){
      p1[j,]=colSums(EXZ[Y_is_zero[,j],j,])+colSums(Y[!Y_is_zero[,j],j]*EZ[!Y_is_zero[,j],])
      
      p2[j]=sum(EX[Y_is_zero[,j],j])+sum(Y[!Y_is_zero[,j],j])
    }else {
      if((sum(Y_is_zero[,j])==1)){
      p1[j,]=colSums(Y[!Y_is_zero[,j],j]*EZ[!Y_is_zero[,j],])+EXZ[Y_is_zero[,j],j,]
      p2[j]=sum(EX[Y_is_zero[,j],j])+sum(Y[!Y_is_zero[,j],j])
      
    }else{
      p1[j,]=colSums(Y[,j]*EZ[,])
      p2[j]=sum(Y[,j])
    }}
  p=c(p1[j,],p2[j])
  theta=solve(S)%*%p
  A[j,]=theta[1:K]
  mu[j]=theta[K+1]
  }
  ##update lamda 
 
  
  fr<-function(lambda){
    
    #Input:
    # lambda
    
    # return  L_lam :L(lambda) 
    y_squared=Y^2
    Y_is_zero = abs(Y) < 1e-6
    e = exp(-lambda *y_squared)
    
    #e[1-e<0]=1
    log_exp_lambda = log(1 -e)#log(1-exp(-lambda y^2_ij))
    log_exp_lambda[is.infinite(log_exp_lambda)]=-10000000
   
    L_lam = -sum(Y_is_zero * (-lambda*EX2 ) + (1 - Y_is_zero) * log_exp_lambda)
    
    return(L_lam)}
  grr<-function(lambda){
    y_squared=Y^2
    Y_is_zero = abs(Y) < 1e-6
    e = exp(-lambda *y_squared)
   
    exp_ratio = e / (1 -e)
    exp_ratio[is.infinite(exp_ratio)]=10000000
    grad = -sum(Y_is_zero * (-EX2) + (1 - Y_is_zero) * y_squared * exp_ratio)
    return(grad)
  }#  grad: gradient
  
  lambda=optim(c(lamda_old), fn=fr, gr=grr, method = "BFGS")$par#minimize
  ##updata sigma
  sigma=c(numeric(0))
  for(j in 1:D){
    
    E_mu_z_i=0#E(Aj z_i+u_j)^2  sum
    for(i in 1:N){
    E_mu_z_i=E_mu_z_i+A[j,]%*%EZZT[i,,]%*%A[j,]+2*mu[j]*A[j,]%*%EZ[i,]+(mu[j])^2
    }
    if(sum(Y_is_zero[,j])>0){
   # p1=colSums(EXZ[Y_is_zero[,j],j,])+colSums(Y[!Y_is_zero[,j],j]*EZ[!Y_is_zero[,j],])
    p3=sum(EX2[Y_is_zero[,j],j])+sum((Y[!Y_is_zero[,j],j])^2)
    #p2=sum(EX[Y_is_zero[,j],j])+sum(Y[!Y_is_zero[,j],j])
    }else{
      #p1=colSums(Y[,j]*EZ[,])
      #p2=sum(Y[,j])
      p3=sum((Y[,j])^2)
    }
    #EX_ij +y_ij 
    
    #EX^2_ij +y^2_ij 
    
    sigma[j]=(-2*A[j,]%*%p1[j,]+p3-2*mu[j]*p2[j]+E_mu_z_i)/N
  }
  list2=list(sigma=sigma,lambda=lambda,mu=mu,A=A)
}
###initialize
initialization<-function(Y,K){
N=nrow(Y)
# dimension of X   N*D
# dimension of Z   N*K
Y_copy=Y
mu_old=colMeans(Y)
Y_copy<-scale(Y_copy,scale=FALSE)
fa1 <- factanal(Y_copy, factors = K, rotation="varimax")
A_old=as.matrix(fa1$loadings)
sigma_old=diag(cov(Y_copy)-A_old%*%t(A_old))
#ini  lambda

lambda_old=runif(1,20,30)
lists<-list(lambda=lambda_old,A=A_old,sigma=sigma_old,mu=mu_old)
return(lists)}
fitmodel<-function(Y,k,iter=100){# k: size of latent Z
  set.seed(356)
  converge=FALSE
  ini<-initialization(Y,k)
  A_o=ini$A
  sigma_o=ini$sigma
  mu_o=ini$mu
  lambda_o=ini$lambda
  list1=Estep(A_o,sigma_o,lambda_o,mu_o,Y)
  
  for (e in 1:iter) {
    lambda1=lambda_o
    sigmas1=sigma_o
    mu1=mu_o
    A1=A_o
    list2=Mstep(A_o,sigma_o,lambda_o,mu_o,Y,list1$EZ,list1$EX,list1$EX2,list1$EZZT,list1$EXZ)
    
    A_o=list2$A
    sigma_o=list2$sigma
    lambda_o=list2$lambda
    mu_o=list2$mu
    gap_new=c((crossprod(lambda_o-lambda1))/crossprod(lambda_o), (crossprod(sigma_o-sigmas1))/crossprod(sigma_o), (sum((A_o-A1)^2))/sum((A_o)^2),(crossprod(mu_o-mu1))/crossprod(mu_o))
    if(all((gap_new)<0.00001))
    {converge=TRUE
    break
    
    }
    
    list1=Estep(A_o,sigma_o,lambda_o,mu_o,Y)
  }
  listt<-list(EZ=list1$EZ, A=A_o,sigmas=sigma_o,lambda=lambda_o,mu=mu_o,converge=converge,gap_new=gap_new)
  return(listt)
}
options(warn=-1)